#!/usr/bin/env python
from __future__ import absolute_import
from __future__ import print_function
from py2cytoscape.data.cyrest_client import CyRestClient
from py2cytoscape.data.cynetwork import CyNetwork
from py2cytoscape.data.style_client import StyleClient
from py2cytoscape.data.style import Style
from py2cytoscape import cyrest
import xgmml
import requests
import json
import pandas as pd
import util
import traceback
import re

def load_json(s_file):
    return json.loads(util.read_string(s_file))

def save_json(s_file, data):
    util.save_string(s_file, json.dumps(data))

def save_network_js(s_file, data):
    if type(data) is dict: # only one network
        name=data["data"]["name"]
        data={name: data}
    else:
        # a list of networks
        data={x["data"]["name"]:x for x in data }
    s='var networks={};'.format(json.dumps(data))
    util.save_string(s_file, s)

def save_style_js(s_file, data, default_style=None):
    if type(data) is dict: # only one style
        data=[data]
    if default_style is None:
        default_style=data[0]["title"]
    s='var styles={};'.format(json.dumps(data))
    util.save_string(s_file, s)

def save_cyjs_app(s_file, network_name, style_name=None, cyjs_path=None, s_template='CyJS.template.html'):
    s=util.read_string(s_template)
    s=re.sub('@NETWORK_NAME@', network_name, s)
    style_name=network_name+".style" if style_name is None else style_name
    s=re.sub('@STYLE_NAME@', style_name, s)
    if cyjs_path is None:
        cyjs_path="http://metascape.org/gp/Content/CyJS"
    s=re.sub('@CYJS_PATH@', cyjs_path, s)
    util.save_string(s_file, s)

def cynetwork_get_name(self):
    return self.to_json()["data"]["name"]
CyNetwork.get_name=cynetwork_get_name

class CytoscapeClient:
    """This is a wrapper to the py2cytoscape wrapper. The local installation of Py2cytoscape has been modified for it to work"""
    def __init__(self, host="localhost", port="1234", version="v1"):
        self.host=host
        self.port=port
        self.version=version
        self.cy = CyRestClient(ip=self.host, port=self.port, version=self.version)
        self.__url = 'http://' + host + ':' + str(port) + '/' + version + '/'

    def gc(self):
        """
        Garbage collection
        """
        try:
            response = requests.get(self.__url+"gc/")
        except Exception as e:
            print('Could not garbadge collect: ' + str(e))
        else:
            return response

    def status(self):
        """
        return free memory left, -1 if not successful
        """
        try:
            X = requests.get(self.__url).json()
        except Exception as e:
            print('Could not get status from cyREST: ' + str(e))
            return -1
        return max(X['memoryStatus']['maxMemory']-X['memoryStatus']['usedMemory'], X['memoryStatus']['freeMemory'])

    def random_id(self, prefix="MY_"):
        import random
        s="ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789"
        return prefix+''.join([s[random.randint(0, len(s)-1)] for i in range(10)])

    def cynet_from_xgmml(self, s_file, name=None, s_layout="force-directed", l_bundle=True, t_xy=None, l_digraph=False):
        """
        Create a Cynetwork object based on an xgmml file
        return CyNetwork
        """
        import os
        if not os.path.exists(s_file):
            util.warn_msg('Missing file: %s' % s_file)
            return None
        if not l_digraph:
            net=xgmml.Network(s_file, name=name)
        else:
            import digraph
            net=digraph.Digraph(s_file, name=name)
        if t_xy is not None:
            net.set_XY(t_xy)
            s_layout=None
        S=net.T_node.header()
        #print(str(net))
        #print(S, name, s_layout, l_bundle)
        return self.cynet_from_network(net, name=name, s_layout=s_layout, l_bundle=l_bundle)

    def cynet_from_network(self, net, name=None, s_layout="force-directed", l_bundle=True):
        """
        Create a Cynetwork object based on a xgmmal.Network object
        """
        data=net.to_json()
        return self.cynet_from_data(data, name=name, s_layout=s_layout, l_bundle=l_bundle)

    def cynet_from_data(self, data, name=None, s_layout="force-directed", l_bundle=True):
        """
        Create a Cynetwork object based on json data
        return CyNetwork
        """
        if name is None:
            name=data['data']['name']
        cynet=None
        for trial in range(3): # something Cytoscape errors, try up to three times
            #print("Trial: ", trial)
            try:
                if cynet is None:
                    cynet=self.cy.network.create(name=name, data=data)
                else:
                    break
            except:
                cynet=None
                print("Fail at trial %d" % (trial+1))
                print("JSON data>", data)
                print(traceback.format_exc())
        if cynet is None:
            import sys
            util.warn_msg('Error during plot_go_network(): %s' % sys.exc_info()[0])
            print(traceback.format_exc())

        id=cynet.get_id()

        if s_layout is not None:
            #print "apply layout %s" % s_layout
            self.cy.layout.apply(name=s_layout, network=cynet)
        if l_bundle:
            self.cy.edgebundling.apply(cynet)
            #print "apply bundle %s" % s_layout
        self.cy.layout.fit(cynet)
        return cynet

    def delete_cynet(self, cynet=None):
        """None means delete all networks"""
        if cynet is None:
            self.cy.network.delete_all()
        else:
            self.cy.network.delete(cynet)

    def get_network(self, suid):
        """Get xgmml.Network with X,Y coordinates"""
        cynet=self.cy.network.create(suid=suid)
        net=xgmml.Network.from_json(cynet.to_json())
        net.T_node.drop([x for x in ['graphics_x','graphics_y','SUID','id_original'] if x in net.T_node.header() ], axis=1, inplace=True)
        t_xy=self.cynet_get_XY(cynet)
        if t_xy is not None:
            t_xy.drop('Gene', axis=1, inplace=True)
            t_xy.rename2({'x':'graphics_x', 'y':'graphics_y'})
            net.T_node=net.T_node.merge(t_xy, left_on='id', right_on='id', how='left')
        return net

    def cynet_get_XY(self, cynet):
        cynet=self.get_cynet(cynet)
        data=cynet.get_first_view()
        if data is None: return None # view does not exist
        nodes=data['elements']['nodes']
        X=[]
        for node in nodes:
            # id is an internal id, Gene is our key
            X.append({'id': str(node['data']['id']), 'x':node['position']['x'], 'y':node['position']['y'], 'Gene':node['data']['Gene']})
        t_xy=pd.DataFrame(X)
        return t_xy

    def get_cynet(self, cynet):
        """make sure we return a cynet object"""
        if type(cynet) in (str, int):
            return self.cy.network.create(suid=int(cynet))
        return cynet

    def cynet_save(self, cynet, s_file="x.png"):
        def save_image(s_file, data):
            with open(s_file, 'wb') as f:
                f.write(data)
                f.flush()

        cynet=self.get_cynet(cynet)
        s_ext=s_file.upper()
        if s_ext.endswith('.PNG'):
            save_image(s_file, cynet.get_png())
        elif s_ext.endswith('.PDF'):
            save_image(s_file, cynet.get_pdf())

    def get_style(self, name):
        if type(name) is str:
            return Style(name)
        return name

    def delete_style(self, name):
        self.cy.style.delete(self.get_style(name))

if __name__=="__main__":
    name="My_TEST"
    host="localhost"
    port="1234"
    cc=CytoscapeClient(host=host, port=port)
    cc.gc()
    cc.cy.network.delete_all()
    cc.cy.style.delete_all()
    #cynet=cc.cynet_from_data(load_json("GONetwork.json"), name=name)
    cynet=cc.cynet_from_xgmml("GONetwork.xgmml", name=name)
    style_json=load_json('ColorByCluster.json')
    style=cc.cy.style.create(name="ColorByCluster", original_style=style_json)
    style_name=style.get_name()
    print(style)
    print(style_name)
    cc.cy.style.apply(style, network=cynet)
    cc.cy.edgebundling.apply(cynet)
    cc.cy.layout.fit(cynet)
    #cc.cy.style.apply(style, network=cynet)
    s_cys="C:\\Temp\\MyTest.cys"
    cc.cy.session.save(s_cys)
    S_id=cc.cy.network.get_all()
    if len(S_id):
        cynet=cc.get_cynet(S_id[0])
        data=cynet.get_first_view()
        data["data"]["shared_name"]=data["data"]["name"]=name
        save_json(name+".json", data)
        style_json=cc.cy.style.get("ColorByCluster", data_format='cytoscapejs')
        save_json(name+".style.json", [style_json])
        save_network_js(name+".js", data)
        save_style_js(name+".style.js", style_json, "ColorByCluster")
        save_cyjs_app("GONetwork.html", "GONetwork", "H:/tmp/CyJS")

    cc.cy.session.open(s_cys)
    S_id=cc.cy.network.get_all()
    print(S_id)
    for id in S_id:
        cynet=cc.get_cynet(id)
        print(cynet, id)
        s_name=cynet.get_name()
        s_file=s_name+'.png'
        cc.cynet_save(cynet, s_file)


