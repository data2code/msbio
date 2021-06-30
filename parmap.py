#!/usr/bin/env python
from __future__ import absolute_import
from __future__ import print_function
import multiprocessing
import multiprocessing.pool
import traceback
import util
import gc
import time
import random
import sys
import os
from six.moves import range
import socket
import zmq
import gridmap
from gridmap.data import zloads, zdumps
from collections import defaultdict
from pprint import pprint
import logging
import tqdm

GRID_MAP_DEBUG_LEVEL=logging.ERROR
#GRID_MAP_DEBUG_LEVEL=logging.INFO
#GRID_MAP_DEBUG_LEVEL=logging.DEBUG

def get_ip():
    host_name = socket.gethostname()
    ip_address = socket.gethostbyname(host_name)
    for _, _, _, _, (ip, _) in socket.getaddrinfo(socket.getfqdn(), 0):
        if ip != '127.0.0.1':
            ip_address = ip
            break
    else:
        util.warn_msg('IP address for JobMonitor server is '
                      '127.0.0.1.  Runners on other machines will be'
                      ' unable to connect.')
        ip_address = '127.0.0.1'
    return ip_address

def get_tcp():
    context = zmq.Context()
    socket= context.socket(zmq.REP)
    ip=get_ip()
    port = socket.bind_to_random_port(f'tcp://{ip}', min_port=6001, max_port=6999, max_tries=100)
    context.destroy()
    return f'tcp://{ip}:{port}'

class Server:

    def __init__(self, url=None):
        self.url=url
        # this is the master, it binds
        self.context = zmq.Context()
        #https://learning-0mq-with-pyzmq.readthedocs.io/en/latest/pyzmq/patterns/client_server.html
        # REP block on receive
        self.zsocket= self.context.socket(zmq.REP)
        self.zsocket.bind(self.url)
        self.hostname=socket.gethostname()
        self.ip=get_ip()

    def listen(self):
        msg=zloads(self.zsocket.recv())
        self.zsocket.send(zdumps("ok"))
        #print("Received message", msg)
        return msg

    def __del__(self):
        if self.context is not None:
            self.context.destroy()

class Client:

    def __init__(self, url=None):
        self.url=url
        # this is the master, it binds
        self.context = zmq.Context()
        #https://learning-0mq-with-pyzmq.readthedocs.io/en/latest/pyzmq/patterns/client_server.html
        # REP block on receive
        self.zsocket= self.context.socket(zmq.REQ)
        self.zsocket.connect(self.url)
        self.hostname=socket.gethostname()
        self.ip=get_ip()

    def talk(self, data):
        msg={"hostname": self.hostname, "ip":self.ip, "data": data}
        self.zsocket.send(zdumps(msg))
        msg = zloads(self.zsocket.recv())
        #print("Sent message", msg)
        return msg

    def __del__(self):
        if self.context is not None:
            self.context.destroy()

def master_wrapper(f_master, cache):
    def f():
        server=Server(cache['url'])
        cache['current_total']=max(cache['total'], 1)
        # if total is 0, we will change the current_total on the fly
        pg=tqdm.tqdm(total=cache['current_total'], position=0)
        cache['pg']=pg
        while True:
            msg=server.listen()
            f_master(msg["data"], cache)
        del server
    return f

def worker_wrapper(params):
    opt, others=params
    out=others[0](others[1])
    phone=Client(opt["url"])
    id=opt["id"]
    phone.talk({"id":id, "cnt":1})
    del phone
    return out

def default_master(msg, cache):
    cache['count']+=msg['cnt']
    cache['counts'][msg['id']]+=1
    #print(msg)
    if cache['total']==0 and cache['count']>=cache['current_total']:
        cache['current_total']=max(int(1.5*cache['current_total']), cache['current_total']+1)
        cache['pg'].total=cache['current_total']
        cache['pg'].refresh()
    cache['pg'].update(msg['cnt'])
    #if (time.time()-cache['last_report'])>=cache['interval'] or cache['count']>=cache['total']:
    #    cache['pg'].check(cache['count'])
    #    cache['last_report']=time.time()
    #    #print(cache['counts'])

def run_cmd(X):
    """wrapper that runs a shell command on a node"""
    if type(X) in (list, tuple):
        opt, s_cmd = X
        s_cmd=f"""PARMAP_URL={opt['url']} PARMAP_TASK_ID={opt['id']} {s_cmd}"""
    else:
        s_cmd=X
        opt_progress=None
    out=util.unix(s_cmd, l_error=True, l_print=True).strip()
    return out

def run_func(X):
    """wrapper that runs a py_file.func_name on a node, pass a task as args"""
    if len(X)==2: # there is opt_progress
        opt_progress, (py_file, func_name, task) = X
    else:
        py_file, func_name, task = X
        opt_progress=None
    import importlib
    s=os.path.abspath(py_file)
    if os.path.dirname(s) not in sys.path:
        sys.path=[os.path.dirname(s)]+sys.path
    pck, ext=os.path.splitext(os.path.basename(s))
    myapp=importlib.import_module(pck)
    func = getattr(myapp, func_name)
    if opt_progress is None:
        out=func(task)
    else:
        out=func((opt_progress, task))
    return out

def cluster_map(f_worker, tasks, opt=None):
    """opt can be used to overwrite defaults, such as:
    num_slots: number of cores per job
    mem_free: can also used to limit the number of jobs per server
    """
    logging.captureWarnings(True)
    logging.basicConfig(format=('%(asctime)s - %(name)s - %(levelname)s - ' +
                                '%(message)s'), level=GRID_MAP_DEBUG_LEVEL)

    hostname=socket.getfqdn()
    if hostname.endswith('.gnf.org'):
        queue="gnf"
        qsub_kw=os.environ.get("QSUB_KW",'-p -5')
    elif hostname.endswith('.novartis.net'):
        queue='default.q'
        qsub_kw=os.environ.get("QSUB_KW",'-l h_rt=14399')
    else: # aws
        queue='all.q'
        qsub_kw=''

    max_processes = len(tasks)
    myopt={'queue':queue, 'qsub_kw':qsub_kw, 'temp_dir':'/tmp', 'max_processes':max_processes,
        'mem_free':'4G', 'name':'gridmap_job', 'num_slots': 1, 'white_list':None, 'quiet':False}
    if opt is not None:
        myopt.update(opt)
    out=gridmap.grid_map(f_worker, tasks, **myopt)
    return out

def map(f_worker, tasks, f_master=None, n_total=0, n_CPU=0, local=True, cluster_opt=None):
    """
        f_worker(*tasks, runner_opts={'url', 'id'}): do the actual work and use runner_opts to report progress
        f_master() receive messages from worker and report progress

        for cluster use, if you have trouble picking the function, you may create shell command and then
        use run_cmd as f_worker, argument will be a shell command string
        This is how we run CellProfiler
    """
    if len(tasks)==0: return []
    url=get_tcp() # we cannot bind here, b/c once forked we have two servers listening to the same port

    if f_master is None:
        if not local:
            # this does not work for cluster, as f_worker was not pickled
            # for cluster, gridmap's default time estimation works, so you don't need this in principle
            n_total=len(tasks)
        else:
            # we wrap with worker_wrapper and use default_master for counting
            n_total=len(tasks)
            tasks=[ (f_worker, task) for task in tasks ]
            #print(tasks)
            f_master=default_master
            f_worker=worker_wrapper

    if f_master is not None:
        #pg=tqdm.tqdm(total=n_total)
        cache={'total':n_total, 'pg':None, 'count':0, 'counts':defaultdict(int), 'url':url }

        srv=multiprocessing.Process(target=master_wrapper(f_master, cache))
        srv.start()

    if n_CPU==0:
        n_CPU=max(multiprocessing.cpu_count()-2,1) if local else len(tasks)
    n_CPU=min(n_CPU, len(tasks))
    out=[]

    def wrap_task(tasks):
        if f_master is not None:
            tasks2=[]
            for i,task in enumerate(tasks):
                opt={"url":url, "id":i}
                tasks2.append( (opt, task) )
            return tasks2
        return tasks

    if local:
        if n_CPU==1:
            opt={"url":url, "id":0}
            # without pooling, easier for debugging
            out=[ f_worker((opt, task)) for task in tasks ]
        else:
            tasks2=wrap_task(tasks)
            pl=multiprocessing.Pool(n_CPU)
            out=pl.map(f_worker, tasks2)
            pl.close()
            pl.join()
    else:
        tasks2=wrap_task(tasks)
        pg_batch=None
        #if len(tasks2)>n_CPU and f_master is None:
        #    pg_batch=tqdm.tqdm(total=n_total, desc="All Tasks", position=2)
        #for i_batch,task_batch in enumerate(util.split(tasks2, chunk_size=n_CPU)):
        opt={'max_processes':n_CPU}
        if cluster_opt is not None: opt.update(opt)
        out=cluster_map(f_worker, tasks2, opt=opt)
        #    out.extend(cluster_map(f_worker, task_batch))
        #    if pg_batch is not None:
        #        pg_batch.update(len(task_batch))

    if f_master is not None:
        srv.terminate()
    return out

def _calc(params):
    """parameter should be just one tuple/list
    The first thing is to flatten the tuple/list to individual variables
    This is an example function calculate the sum from a to b"""
    opt_progress, (a,b)=params # the parameter is (start, end) passed in as tasks
    #print(opt, a, b)
    phone=Client(opt_progress["url"])
    id=opt_progress["id"]
    sum=0
    for i in range(a, b):
        sum+=i
        time.sleep(0.1)
        #print("To send message", {"id":id, "i": i, "sum":sum})
        phone.talk({"id":id, "cnt": 1})
    del phone
    return sum

def _calc2(params):
    (a,b)=params
    sum=0
    for i in range(a, b):
        sum+=i
        time.sleep(0.1)
    return sum

if __name__=="__main__":

    tasks=[(i,i+30) for i in range(10)]
    #out=map(_calc, tasks, default_master, 0, n_CPU=1)
    out=map(_calc, tasks, default_master, 30*10, n_CPU=4)
    #out=map(_calc2, tasks, n_CPU=5)
    #out=map(_calc, tasks, default_master, 30*10, n_CPU=0, local=False)
    #out=map(_calc2, tasks, n_CPU=3, local=False)
    #out=map(run_cmd, [ "hostname" for i in range(10) ], local=False)
    #out=map(run_cmd, [ """bash -c 'echo "$PARMAP_URL <>  $PARMAP_TASK_ID"'""" for i in range(3) ], f_master=default_master, local=False)
    #out=map(run_func, [ ("/depts/ChemInfo/p/python/lib/parmap.py", "_calc2", (i, i+30)) for i in range(10)], local=False)
    #out=map(run_func, [ ("/depts/ChemInfo/p/python/lib/parmap.py", "_calc", (i, i+30)) for i in range(10)], default_master, 30*10, local=False)
    print(out)
