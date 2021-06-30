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
#import sys
import os
from six.moves import range

#sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)

## no longer needed. http://bytes.com/topic/python/answers/552476-why-cant-you-pickle-instancemethods
## This code enables static or instance methods to be pickled, therefore can be used in multiprocessing calls

class MP:
    """Mutliprocessing class that starts workers and processes input. Two advantages of this class compared to other solutions: (1) it can parallelize any method, be it a function, an instance method, a static method, an anonymous function, a method within a method, etc. (2) it has the mechanism to start worker processes when memory footprints is small, and have workers continuously process taskes, unless we call stop()."""
    def __init__(self, DEBUG=False, QUIT=True):
        """DEBUG: boolean, default False. If True print out lots of debugging messages
        QUIT: boolean, default True. Once map() finishes, the workers quit by default. If False, workers stay for future tasks."""
        self.n_CPU=self.n_use=0
        self.q_in = []
        # each process has its own queue
        self.q_out = None
        self.q_flag = None
        self.c_proc= {}
        self.registry={}
        self.DEBUG=DEBUG
        self.QUIT=QUIT
        self.mgr=None
        self.lock=None
        self.work_status=None
        self.n_running=None # use a list, so that fetch routine can modify it (if we do n_running=0, fetch cannot m

    def start(self, f=None, n_CPU=0, as_daemon=True):
        """Launch workers.
        f: method, task for worker, default None, if None, worker runs a wrapper task, which expect the first argument to be a registered method name
        n_CPU: number of workers to launch.
        as_daemon: defaults to True, so the parent will terminate children,
                but this prevent the children from launching new processes
            If you launch processes within f(), as_daemon should be False
        """
        if self.has_started():
            util.warn_msg('Already started, stop first and start new!')
            self.stop()
        self.n_CPU=n_CPU if n_CPU>0 else multiprocessing.cpu_count()-1
        self.n_use=self.n_CPU # use all CPUs by default, but user could ask to use fewer
        self.q_in = [multiprocessing.Queue(1) for i in range(self.n_CPU)]
        self.q_out =[multiprocessing.Queue() for i in range(self.n_CPU)]
        self.c_proc={} #bug, if c_proc is not inside manager, each process will see {}, {1proc}, {2procs} ..., as they are created
        self.f=f
        if self.f is None:
            self.f=self.wrapper
        # put work_status, n_running into manager, so visible by all processes, see the same copy
        self.mgr=multiprocessing.Manager()
        #self.c_proc=self.mgr.dict({})
        self.work_status=self.mgr.list([False for i in range(self.n_CPU)])
        # remember which process is working on my job
        #self.has_my_job=[False for i in range(self.n_CPU)]
        self.has_my_job=None
        # whether job is done, results ready for reading
        self.job_is_done=self.mgr.list([False for i in range(self.n_CPU)])
        self.n_running=self.mgr.list([0, 0])
        # n_running [0] is # of running process, [1] # of map() calls running
        # use a list, so that fetch routine can modify it (if we do n_running=0, fetch cannot m
        self.lock=self.mgr.Lock()
        self.mail=multiprocessing.Condition() # notify all processes when a worker finishes
        procs=[multiprocessing.Process(target=self.spawn(self.f),args=(self.q_in[i], self.q_out[i], self.mail, self.job_is_done)) for i in range(self.n_CPU)]
        for p in procs:
            p.daemon = as_daemon
            p.start()
            self.c_proc[p.pid]=p
            self.n_running[0]+=1
            if self.DEBUG: print("DEBUG>>> %s  %d %d" % (str(p), p.pid, self.n_running[0]))
        if self.DEBUG: print("DEBUG> self.c_proc: ", self.c_proc)

    def CPU_used(self, n_use=0):
        """Specify number of CPUs to use by default, has to be no more than self.n_CPU"""
        self.n_use=self.n_CPU
        if n_use>0:
            self.n_use=min(n_use, self.n_CPU)

    def register(self, name, func):
        """Register a task for workers.
        name: registered name
        func: method pointer
        no return: registry updated"""
        self.registry[name]=func

    def register_dict(self, c):
        """Register tasks for workers.
        c: dict[str]=method, keys are registered names
        no return: registry updated"""
        for k,v in c.items():
            self.register(k, v)

    def get_func(self, name):
        """return method: previously registered under the name"""
        return self.registry[name]

    def wrapper(self, args):
        """Internal default wrapper function to be called if start(f=None)"""
        name=args[0]
        args=args[1:]
        if callable(name):
            return name(*args)
        else:
            if name not in self.registry:
                util.error_msg('Function %s were not registered!' % name)
            return self.registry[name](*args)

    def spawn(self, f):
        """internal helper function"""
        def func(q_in,q_out,mail,job_is_done):
            while True:
                (i, i_worker, x) =q_in.get()
                try:
                    if i is None:
                        # garbage collect to release memory
                        gc.collect()
                        q_out.put((None, multiprocessing.current_process().pid))
                    else:
                        q_out.put((i, f(x)))
                except:
                    #myfile=open("ERROR.parallel.txt", "a")
                    #myfile.write("\nException in process %d\n" % multiprocessing.current_process().pid)
                    #myfile.write(traceback.format_exc())
                    #myfile.close()
                    s_msg=str(multiprocessing.current_process().pid)+"\n"+traceback.format_exc()
                    q_out.put((None, s_msg))
                finally:
                    job_is_done[i_worker]=True
                    with mail:
                        mail.notify_all()
                    if i is None:
                        break
        return func

    def map(self, tasks, n_CPU=0, l_quit=None, on_task_end=None):
        """X: list[input tuple], list of input parameters. If workers were started with f=None, each element in X in passed to the wrapper task. In that case, we expect X to be a tuple (or list) and the first element of X must be either the method pointer or its registered name. However, many methods, such as instance method or func within a func, cannot be pickled, therefore the method cannot be send over the pipe. We should pre-register such methods and call them by name. Need example later.
        l_quit: boolean, default None, if specified, controls whether workers quit or not after tasks are processed. If not, workers wait for future tasks.
        return list"""
        # very similar to the idea in https://stackoverflow.com/questions/3288595/multiprocessing-how-to-use-pool-map-on-a-function-defined-in-a-class, author klaus se
        l_quit=l_quit if l_quit is not None else self.QUIT
        #if self.is_busy():
        #    util.error_msg('Works are still busy!')
        if n_CPU==0: n_CPU=self.n_use # defaults to n_use
        n_CPU=min(n_CPU, self.n_CPU) # one could start 8 CPUs, but only use 4 for mapping, if the work takes lots of memory
        res=[]
        if len(tasks)==0 and not l_quit: return res
        #print '=============', self.c_proc
        if not self.has_started() and len(tasks)>0:
            util.warn_msg('Please start processes first, no worker is running!')
            util.warn_msg('However, we will process the task with ONE cpu!!!')
            return [self.wrapper(x) for x in X]

        if len(tasks)>0 and n_CPU==0:
            return [self.wrapper(x) for x in X]

        s_pid=str(multiprocessing.current_process().pid)
        has_my_job=[False for i in range(self.n_CPU)]

        def engine():
            print("=================================================================")
            print("PID: ", str(multiprocessing.current_process().pid))
            print("WORK STATUS: ", self.work_status)
            print("HAS MY JOB: ", has_my_job)
            print("JOB IS DONE: ", self.job_is_done)
            print("N_RUNNING: (%d, %d) " % (self.n_running[0], self.n_running[1]))
            print("=================================================================")

        def is_busy():
            return sum(has_my_job)>0

        def process_out(out):
            i, x = out
            if i is None:
                self.n_running[0]-=1
                # I modify the original code, so that we can join the process and release it as soon as possible
                if type(x) is str:
                    print("Exception> "+x)
                    exit()
                else:
                    if self.DEBUG: print("Progress: %d processes remaining. Stopping %d" % (self.n_running[0], x))
                    #print self.c_proc.keys()
                    self.c_proc[x].join()
                    del self.c_proc[x]
                    if self.DEBUG: print("Progress: process %d stopped." % x)
            else:
                res.append(out)
                if self.DEBUG: print("Progress: %d of %d item calculated." % (len(res), len(tasks)))

        def fetch(l_lock=False):
            while is_busy():
                l_fetch_something=False
                for i_worker in range(self.n_CPU):
                    if has_my_job[i_worker] and self.job_is_done[i_worker]:
                        try:
                            (i,x)=self.q_out[i_worker].get()
                            process_out((i, x))
                            self.n_running[1]-=1
                            self.work_status[i_worker]=False
                            has_my_job[i_worker]=False
                            self.job_is_done[i_worker]=False
                            l_fetch_something=True
                            if self.DEBUG:
                                engine()
                            if i is not None and on_task_end is not None: # callback function, which may add new tasks that becomes feasible
                                on_task_end(tasks, i, x)
                            with self.mail:
                                self.mail.notify_all()
                        except Exception as e:
                            print("ERROR> Fail to fetch results from worker: %d" % i_worker)
                            print(traceback.format_exc())
                            return
                if not l_fetch_something:
                    if l_lock:
                        with self.mail:
                            self.mail.wait(timeout=8.5+random.random()*3)
                    else:
                        return

        if self.DEBUG: print("DEBUG> self.n_running: ", self.n_running)
        if self.DEBUG: print("DEBUG> self.c_proc.keys(): ", list(self.c_proc.keys()))
        ###ZHOU FEB16,2016
        #self.n_running[1]+=1
        ###
        task_id=0
        while True:
            if task_id==len(tasks): break
            x=tasks[task_id]
            if self.DEBUG: print("fetch job entry %d " % task_id)
            #print self.work_status
            self.lock.acquire()
            j=util.index(False, self.work_status) # find an idle worker
            l_put_something=False
            if j>=0 and sum(has_my_job)<n_CPU: #j>=0 and j<n_CPU: # we only use up to n_CPU, even if there are more workers
                #print "assing job to %d" % j
                self.work_status[j]=True # flag it as busy
                has_my_job[j]=True
                if self.DEBUG: print("DEBUG> self.c_proc.keys(): ", list(self.c_proc.keys()))
                ###ZHOU FEB16,2016
                self.n_running[1]+=1
                ###
                self.lock.release()
                self.q_in[j].put((task_id,j,x)) # assign task
                task_id+=1
                if self.DEBUG: print("Progress: send input %d of %d items." % (task_id, len(X)))
                l_put_something=True
                if self.DEBUG:
                    print(">>>A2")
                    engine()
            else:
                self.lock.release()
            # we constantly removing items from the output queue, so that the process can release some memory
            fetch()
            if not l_put_something:
                with self.mail:
                    self.mail.wait(timeout=8.5+random.random()*3)
                    fetch()

        while (is_busy()):
            fetch(True)
        if self.DEBUG:
            print(">>>A3")
            engine()
        self.lock.acquire()
        ###ZHOU FEB16,2016
        #self.n_running[1]-=1
        ###
        if self.DEBUG:
            print(">>>QUIT="+("True" if l_quit else "False"))
            print(">>>n_running[1]=%d" % self.n_running[1])
        if l_quit and self.n_running[1]==0: # I am the last one running map()
            if self.DEBUG:
                print(">>>A4")
                engine()
            for i in range(self.n_CPU):
                self.q_in[i].put((None,i,None))
                self.work_status[i]=True
                self.n_running[1]+=1
                has_my_job[i]=True
            self.lock.release()
            while (is_busy()):
                fetch(True)
                if self.DEBUG:
                    engine()
        else:
            self.lock.release()
        if self.DEBUG:
            for i,x in enumerate(res):
                print('>>>A4 ', i, type(x), type(x[0]), type(x[1]))
        res.sort(key=lambda x: x[0])
        return [x for i,x in res]

    def has_started(self):
        return len(self.c_proc)>0

    def stop(self):
        if not self.has_started(): return
        self.map([], l_quit=True)

# http://stackoverflow.com/questions/3288595/multiprocessing-using-pool-map-on-a-function-defined-in-a-class
def parmap(f, L, n_CPU=0, DEBUG=False):
    """Starts n_CPU workers that knows how to do task f on input list L
    f: method, task
    L: list[input arguments], input list
    n_CPU: int default 0. number of workers, 0 means number of physical CPUs on the machine
    DEBUG: boolean, default False. If Ture, print debugging messages
    return list"""
    if n_CPU==1:
        return [ f(x) for x in L ]
    mp=MP(DEBUG=DEBUG)
    mp.start(f=f, n_CPU = n_CPU)
    return mp.map(L)

def parprep(f=None, n_CPU=0, DEBUG=False):
    """Starts n_CPU workers that knows how to do task f on input list L
    f: method, task
    L: list[input arguments], input list
    n_CPU: int default 0. number of workers, 0 means number of physical CPUs on the machine
    DEBUG: boolean, default False. If Ture, print debugging messages
    return list"""
    mp=MP(DEBUG=DEBUG)
    mp.start(f=f, n_CPU = n_CPU)
    return mp

def proxy(f):
    return f

def map(f, tasks, n_CPU=0, progress=False, min_interval=30):
    """min_interval, do not report progress more frequent that 30 sec"""
    if len(tasks)==0: return []
    if n_CPU==0:
        n_CPU=multiprocessing.cpu_count()
    if n_CPU<=1:
        if not progress:
            return [ f(args) for args in tasks ]
        else:
            out=[]
            pg=util.Progress(len(tasks))
            i_start=time.time()
            for i, args in enumerate(tasks):
                out.append(f(args))
                if (time.time()-i_start)>=min_interval or i==len(tasks)-1:
                    pg.check(i+1)
                    i_start=time.time()
            return out
    n_CPU=min(n_CPU, len(tasks))
    pl=multiprocessing.Pool(n_CPU)
    if progress:
        # https://stackoverflow.com/questions/5666576/show-the-progress-of-a-python-multiprocessing-pool-map-call
        out=[]
        pg=util.Progress(len(tasks))
        i_start=time.time()
        #for i, _ in enumerate(pl.imap(f, tasks), 1):
        for i, _ in enumerate(pl.imap(f, tasks)):
            out.append(_)
            #print(">>>", i, time.time(), i_start, time.time()-i_start)
            if (time.time()-i_start)>=min_interval or i==len(tasks)-1:
                pg.check(i+1)
                i_start=time.time()
    else:
        out=pl.map(f, tasks)
    pl.close()
    pl.join()
    return out

def test_f(a):
    time.sleep(2)
    return a*a

def test_grand(a):
    import parmap
    return parmap.map(test_f, range(4), n_CPU=4)
    #map(f_worker, tasks, f_master=None, n_total=0, n_CPU=0, local=True):

def on_task_end(tasks, task_id, task_out):
    """task_id just finished with result task_out, we may add a new job to tasks"""
    if task_id in (3,4,12,16):
        #print("TaskID: ", task_id)
        tasks.append( (test_f, task_id*10) )
        #print("New task queue: ", [x[1] for x in tasks ])

if __name__=="__main__":

    # this is an example to show how to use on_task_end callback to add new tasks on the fly
    # when tasks have topological dependency, new tasks can be added only after some previous tasks are completed
    # on_task_end will track what has completed and add new tasks accordingly
    mp=MP(DEBUG=False, QUIT=True)
    # as_daemon=False allows us to launch multiple processes within the function
    mp.start(n_CPU=4)
    X=[ (test_f, x) for x in range(20)]
    out=mp.map(X, on_task_end=on_task_end)
    print(out)
    out=mp.map(X, on_task_end=None)
    print(out)

    mp.start(n_CPU=4, as_daemon=False)
    X=[ (test_grand, x) for x in range(20)]
    out=mp.map(X, on_task_end=on_task_end)
    print(out)

    exit()

    def calc(X):
        """parameter should be just one tuple/list
        The first thing is to flatten the tuple/list to individual variables
        This is an example function calculate the sum from a to b"""
        a,b=X # the parameter is (start, end) passed in as tasks
        sum=0
        for i in range(a, b+1):
            sum+=i
            time.sleep(0.5)
        return sum

    tasks=[(i,i+30) for i in range(10)]
    out=map(calc, tasks, n_CPU=2, progress=True, min_interval=5)
    print(out)
