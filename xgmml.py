#!/usr/bin/env python
from __future__ import absolute_import
from __future__ import print_function
import numpy as np
import pandas as pd
import re
import util
import os
import xml.etree.ElementTree as ET
import datetime as dt
from scipy.sparse import dok_matrix
import hashlib
import six
from six.moves import range
from six.moves import zip
from xml.sax.saxutils import escape

class XGMML(object):

    def __init__(self):
        self.T_node=None
        self.T_edge=None
        self.name="untitled"

    def parse(self, s_file):
        tree=ET.parse(s_file)
        root=tree.getroot()
        self.name=root.attrib['label']
        c_name={}
        nodes=[]
        for node in root:
            if not node.tag.endswith('node'): continue
            id=node.attrib['id']
            c_name[id]=node.attrib['label']
            c={}
            #c['_id']=id
            for att in node:
                if att.tag.endswith('graphics'):
                    for k,v in att.attrib.items():
                        c['graphics_'+k]=v
                    continue
                elif not att.tag.endswith('att'):
                    continue
                v=att.attrib.get('value', None)
                ty=att.attrib['type']
                if ty=='integer':
                    v=int(v) if v is not None else 0
                elif ty=='real':
                    v=float(v) if v is not None else 0.0
                c[att.attrib['name']]='' if pd.isnull(v) else v
            nodes.append(c)
        self.T_node=pd.DataFrame(nodes)
        if 'Gene' in self.T_node.header():
            self.T_node['Gene']=self.T_node['Gene'].astype(str)

        edges=[]
        for edge in root:
            if not edge.tag.endswith('edge'): continue
            id_A=edge.attrib['source']
            id_B=edge.attrib['target']
            gene_A=id_A
            gene_B=id_B
            ty='pp'

            if 'label' in edge.attrib:
                m=re.search(r'^(\S+)\s+\((\S+)\)\s+(\S+)', edge.attrib['label'])
                if m: gene_A, ty, gene_B=m.groups()
            c={}
            c['Gene_A']=gene_A
            c['Name_A']=c_name[id_A]
            c['Gene_B']=gene_B
            c['Name_B']=c_name[id_B]
            c['TYPE']='Direct' if ty=='pp' else 'Indirect'
            for att in edge:
                if not att.tag.endswith('att'): continue
                name=att.attrib['name']
                if name in ('canonicalName', 'interaction'): continue
                v=att.attrib.get('value', None)
                ty=att.attrib['type']
                if ty=='integer':
                    v=int(v) if v is not None else 0
                elif ty=='real':
                    v=float(v) if v is not None else 0.0
                c[name]='' if pd.isnull(v) else v
            edges.append(c)
        self.T_edge=pd.DataFrame(edges)
        # xgmml may contains many visualization columns
        self.T_node.dropna(axis=1, how='all', inplace=True)
        self.T_edge.dropna(axis=1, how='all', inplace=True)

    def set_name(self, s_name):
        self.name=s_name

    def set_tables(self, T_node, T_edge):
        self.T_node=T_node
        self.T_edge=T_edge
        if len(T_node):
            self.T_node['Gene']=util.sarray2sarray(self.T_node['Gene'])
        if len(T_edge):
            self.T_edge['Gene_A']=util.sarray2sarray(self.T_edge['Gene_A'])
            self.T_edge['Gene_B']=util.sarray2sarray(self.T_edge['Gene_B'])

    def print_header(self, s_name=None):
        now=dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        s_name=s_name or self.name or "Untitled"
        s='''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<graph label="%s" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:xlink="http://www.w3.org/1999/xlink" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:cy="http://www.cytoscape.org" xmlns="http://www.cs.rpi.edu/XGMML" >
  <att name="documentVersion" value="1.1"/>
  <att name="networkMetadata">
    <rdf:RDF>
      <rdf:Description rdf:about="http://www.cytoscape.org/">
        <dc:type>Protein-Protein Interaction</dc:type>
        <dc:description>N/A</dc:description>
        <dc:identifier>N/A</dc:identifier>
        <dc:date>%s</dc:date>
        <dc:title>Cytoscape Network</dc:title>
        <dc:source>http://www.cytoscape.org/</dc:source>
        <dc:format>Cytoscape-XGMML</dc:format>
      </rdf:Description>
    </rdf:RDF>
  </att>
  <att type="string" name="backgroundColor" value="#ffffff"/>
  <att type="real" name="GRAPH_VIEW_ZOOM" value="1.0"/>
  <att type="real" name="GRAPH_VIEW_CENTER_X" value="300.0"/>
  <att type="real" name="GRAPH_VIEW_CENTER_Y" value="300.0"/>
''' % (escape(s_name), now)
        return s

    def print_xgmml(self, t_noa, t_sif, s_idCol, s_nameCol, S_nodeCol, S_edgeCol, l_remove_loop=False):
        S_name=t_noa.header()
        n=len(S_name)
        if len(t_sif) and t_noa[s_idCol].dtype is not np.dtype(object):
            t_noa[s_idCol]=t_noa[s_idCol].astype(str)
        if len(t_sif) and t_sif['Gene_A'].dtype is not np.dtype(object):
            t_sif['Gene_A']=t_sif['Gene_A'].astype(str)
        if len(t_sif) and t_sif['Gene_B'].dtype is not np.dtype(object):
            t_sif['Gene_B']=t_sif['Gene_B'].astype(str)
        S_all=set(list(t_sif['Gene_A'])+list(t_sif['Gene_B'])) if len(t_sif) else []

        c_type={}
        for i in range(n):
            ty=t_noa.iloc[:,i].dtype
            if np.issubdtype(ty, np.integer):
                ty='integer'
            elif np.issubdtype(ty, np.floating):
                ty='real'
            else:
                ty='string'
            c_type[S_name[i]]=ty

        s=""
        id=0
        c_id={}
        if s_nameCol not in S_name: s_nameCol=s_idCol
        #t_noa.to_csv('t.csv')
        #print t_noa[:3]
        #print t_noa.col_types()
        #print c_type

        for i in range(len(t_noa)):
            s_gene=t_noa[s_idCol].iloc[i]
            if t_noa[s_idCol].iloc[i] not in S_all: continue
            if s_gene not in c_id:
                id+=1
                c_id[s_gene]=str(id)
            s_name=t_noa[s_nameCol].iloc[i]
            #print s_name, type(s_name),s_gene, type(s_gene),c_id[s_gene],type(c_id[s_gene])
            s+="  <node label=\""+escape(s_name)+"\" id=\""+c_id[s_gene]+"\">\n"
            s+="    <att type=\"string\" name=\"canonicalName\" value=\""+escape(s_name)+"\"/>\n"
            for s_col in S_name:
                if s_col in S_nodeCol:
                    #in S_listCol or c_type[s_col]!='string':
                    v=t_noa[s_col].iloc[i]
                    v='NA' if pd.isnull(v) else str(v)
                    s+="    <att type=\""+c_type[s_col]+"\" name=\""+escape(s_col)+"\" value=\""+escape(str(t_noa[s_col].iloc[i]), entities={'"':'&quot;'})+"\"/>\n"
                #else:
                #    print t_noa[s_col].iloc[i]
                #    S2=re.split(r',\s*', t_noa[s_col].iloc[i])
                #    s+="  <att type=\"list\" name=\""+s_col+"\">\n"
                #    for s2 in S2:
                #        s+="    <att type=\"string\" name=\""+s_col+"\" value=\""+s2+"\"/>\n"
                #    s+="  </att>\n"
            S_graph=[x for x in S_name if x.startswith('graphics_')]
            if S_graph:
                s+="    <graphics "+' '.join([x[9:]+"=\""+str(t_noa[x].iloc[i])+"\"" for x in S_graph])+"/>\n"
            else:
                x=np.random.rand()*300.0+150
                y=np.random.rand()*300.0+150
                # do not specify h and w, Cytoscape 3.x will consider node size locked and need to Remove Byp. in order to resize
                #s+="    <graphics type=\"ELLIPSE\" h=\"10.0\" w=\"10.0\" x=\""+str(x)+"\" y=\""+str(y)+"\"/>\n"
                s+="    <graphics type=\"CIRCLE\" x=\""+str(x)+"\" y=\""+str(y)+"\"/>\n"
            s+="  </node>\n"

        S_name=t_sif.header()
        n=len(S_name)
        c_type={}
        for i in range(n):
            ty=t_sif.iloc[:,i].dtype
            if ty is np.dtype(int):
                ty='int'
            elif ty is np.dtype(float):
                ty='real'
            else:
                ty='string'
            c_type[S_name[i]]=ty
        for i in range(len(t_sif)):
            s_A=t_sif['Gene_A'].iloc[i]
            s_type='pp' if t_sif['TYPE'].iloc[i]=='Direct' else 'ppp'
            s_B=t_sif['Gene_B'].iloc[i]
            idA=c_id[s_A]
            idB=c_id[s_B]
            if l_remove_loop and idA==idB: continue  # remove loop edge
            s+="  <edge label=\""+s_A+" ("+s_type+") "+s_B+"\" source=\""+idA+"\" target=\""+idB+"\">\n"
            s+="    <att type=\"string\" name=\"canonicalName\" value=\""+s_A+" ("+s_type+") "+s_B+"\"/>\n"
            s+="    <att type=\"string\" name=\"interaction\" value=\""+s_type+"\"/>\n"
            for s_col in S_name:
                if s_col in S_edgeCol:
                    v=t_sif[s_col].iloc[i]
                    v='NA' if pd.isnull(v) else str(v)
                    s+="    <att type=\""+c_type[s_col]+"\" name=\""+s_col+"\" value=\""+escape(str(t_sif[s_col].iloc[i]), entities={'"':'&quot;'})+"\"/>\n"
            s+="  </edge>\n"
        s+="</graph>"
        return s

    def save(self, s_file=None, s_graph=None, l_remove_loop=False):
        if s_file is not None:
            fname, ext = os.path.splitext(s_file)
            if ext != '.xgmml': fname=s_file
            s_name=os.path.split(fname)[1]
        else:
            fname=s_graph or self.name or 'Untitled'
            s_name=fname
        ext='.xgmml'
        s_out=self.print_header(s_name)
        S_node=[]
        S_edge=[]
        S_REMOVE_NODE=set(["label","id","canonicalName"])
        S=self.T_node.header()
        S_node.extend([s for s in S if s not in S_REMOVE_NODE and not s.startswith('graphics')])
        S_REMOVE_EDGE=set(["label","Gene_A","Name_A","InteractionType","Gene_B","Name_B","canonicalName","InteractionType","TYPE"])
        S=self.T_edge.header()
        S_edge.extend([s for s in S if s not in S_REMOVE_EDGE])
        s_out+=self.print_xgmml(self.T_node, self.T_edge, "Gene", "Symbol", S_node, S_edge, l_remove_loop=l_remove_loop)
        util.save_list(fname+ext, s_out, s_end='\n')

    def add_node_attr(self, t_attr, s_key="Gene"):
        #S=self.T_node.header()
        self.T_node=pd.merge(self.T_node, t_attr, left_on="Gene", right_on=s_key, how="left")
        #S2=self.T_node.header()
        #print S2
        #for s in S2:
        #    if s not in S:
        #        self.T_node[s].fillna('NA', inplace=True)

    def add_edge_attr(self, t_attr, S_key=None):
        S_key=S_key or ["Gene_A", "InteractionType", "Gene_B"]
        #S=self.T_edge.header()
        self.T_edge=pd.merge(self.T_edge, t_attr, left_on=["Gene_A", "InteractionType", "Gene_B"], right_on=S_key, how="left")
        #S2=self.T_node.header()
        #for s in S2:
        #    if s not in S:
        #        self.T_node[s].fillna('NA', inplace=True)

    def to_network(self, allow_indirect=False):
        return Network(self.T_edge, allow_indirect=allow_indirect, name=self.name, T_node=self.T_node, s_noa=None)

#  Reimplementation based on MCODE Java Code version 1.2
#    v1.2 is written in the README file, although the release version for Cytoscape was 1.3.1
#    http://chianti.ucsd.edu/svn/csplugins/trunk/mskcc/gbader/mcode/src/csplugins/mcode/
#    svn co http://chianti.ucsd.edu/svn/csplugins/trunk/mskcc/gbader/mcode/
#  Key Improvements:
#    (1) add Decompose after Haircut, although it does not seem to change results too much
#        why MOCDE in Cytoscape gives a cluster containing no seed?
#    (2) add an optimization routine to get better clusters than original implementation
#        when there are multiple starting nodes with the same score, the order of how these nodes are explored
#        makes a difference. This also cause the algorithm to be non-deterministic. Different MCODE runs may ends up with
#        different results.  What I did is to loop through all the nodes with the same score, then choose the node
#        that gives the best-scoring network. Then repeat and identify the next best node to pick
#        yes, slower, but it finds better-scoring clusters.
#

class NodeInfo(object):

    def __init__(self):
      self.density = 0.0        #neighborhood density
      self.numNodeNeighbors = 0 #number of node nieghbors
      self.coreLevel = 0        #e.g. 2 = a 2-core
      self.coreDensity = 0.0    #density of the core neighborhood
      self.nodeNeighbors=[]     #stores node indices of all neighbors
      self.score=0.0            #node score

    def score_node(self, degreeCutoff):
        if (self.numNodeNeighbors > degreeCutoff):
            self.score = self.coreDensity*self.coreLevel
        else:
            self.score = 0.0

    def __str__(self):
        s=">density: %.2f; numNeighbors: %d; coreLevel: %.2f; coreDensity:%.2f; score: %.2f" % \
            (self.density, self.numNodeNeighbors, self.coreLevel, self.coreDensity, self.score)
        return s

    def clone(self):
        x=NodeInfo()
        x.density=self.density
        x.numNodeNeighbors=self.numNodeNeighbors
        x.coreLevel=self.coreLevel
        x.coreDensity=self.coreDensity
        x.nodeNeighbors=self.nodeNeighbors[:]
        x.score=self.score
        return x

class Network(object):

    @staticmethod
    def canonicalize_table(t):
        """Make sure gid_A<=gid_B"""
        mask=t.Gene_A>t.Gene_B
        tmp=t[mask].copy()
        if len(tmp)==0: return t
        t.loc[tmp.index, 'Gene_A']=tmp.Gene_B
        t.loc[tmp.index, 'Gene_B']=tmp.Gene_A
        return t

    def from_table(self, t_edge):
        self.data={}
        if "TYPE" not in t_edge.header():
            t_edge['TYPE']=['Direct']*len(t_edge)
        t=t_edge
        if not self.allow_indirect:
            t=t_edge[~ t_edge.TYPE.isin(["Indirect","ppp"])]
        nodes=set(t_edge.Gene_A)|set(t_edge.Gene_B)
        idx=util.index('SCORE', [x.upper() for x in t.header()])
        score=np.ones(len(t)) if idx<0 else t.iloc[:, idx].values
        data={k:{} for k in nodes}
        [ (data[k].__setitem__(v,c) or data[v].__setitem__(k,c)) for k,v,c in zip(t.Gene_A, t.Gene_B, score) ]
        self.data=data
        return data
        #for i in range(len(t_edge)):
        #    if not self.allow_indirect and t_edge['TYPE'].iloc[i] in ["Indirect","ppp"]: continue
        #    s1=t_edge['Gene_A'].iloc[i]
        #    s2=t_edge['Gene_B'].iloc[i]
        #    if s1 not in self.data: self.data[s1]={}
        #    if s2 not in self.data: self.data[s2]={}
        #    score=1 if idx<0 else t_edge.iat[i, idx]
        #    self.data[s1][s2]=score
        #    self.data[s2][s1]=score

    def create_node(self):
        #S_nodes=[]
        #for k,v in self.data.items():
        #    S_nodes.append(k)
        #    S_nodes.extend(list(v.keys()))
        #S_nodes=list(set(S_nodes))
        # data is symmetrical
        S_nodes=list(self.data.keys())
        self.T_node=pd.DataFrame({'Gene':S_nodes})

    def create_edge(self):
        rows=[{'Gene_A':a, 'Gene_B':b, 'TYPE':'Direct', 'SCORE':c } \
            for a,v in self.data.items() for b,c in v.items() if a<b ]
        self.T_edge=pd.DataFrame(rows)

    def __init__(self, x=None, allow_indirect=False, name='Untitled', T_node=None, s_noa=None, premade_T_edge=None, skip_copy=False):
        """premade_T_edge: if not None and x is dict, it will be used as self.T_edge for speed-up
            for human ppi, it takes 10 secs to clone data dict, so skip_copy can speed it up
        """
        self.allow_indirect=allow_indirect
        self.name=name or 'Untitled'
        self.data=None
        self.T_node=None
        self.T_edge=None

        if x is None:
            self.data={}
            return

        if isinstance(x, Network):
            import copy
            self.data=copy.deepcopy(x.data)
            if x.T_node is not None: self.T_node=x.T_node.copy()
            if x.T_edge is not None: self.T_edge=x.T_edge.copy()
        elif type(x) is dict:
            if not skip_copy:
                import copy
                self.data=copy.deepcopy(x)
            else:
                self.data=x
            # create empty T_node and T_edge
            if T_node is None:
                self.create_node()
            else:
                self.T_node=T_node if skip_copy else T_node.copy()
            if premade_T_edge is not None:
                self.T_edge=premade_T_edge if skip_copy else premade_T_edge.copy()
            else:
                self.create_edge()
        elif isinstance(x, pd.DataFrame):
            self.T_edge=Network.canonicalize_table(x.copy())
            self.from_table(self.T_edge)
            if T_node is not None:
                self.T_node=T_node.copy()
            else:
                self.create_node()
        elif type(x) is str: # filename
            fname, ext = os.path.splitext(x)
            if (ext==".sif"):
                t=pd.read_csv(x, sep="\t")
                self.from_table(t)
                self.T_edge=t.copy()
                if t_noa:
                    self.T_node=pd.read_csv(s_noa, sep="\t")
                    self.T_node['Gene']=self.T_node['Gene'].astype(str)
            elif (ext==".xgmml"):
                xg=XGMML()
                xg.parse(x)
                self.from_table(xg.T_edge)
                self.T_node=xg.T_node.copy()
                self.T_edge=xg.T_edge.copy()
                self.name=xg.name
        if len(self.T_edge)>0 and 'Indirect' in self.T_edge['TYPE'].unique():
            self.allow_indirect=True

    def is_node(self, node):
        return node in self.data

    def is_empty(self):
        if not self.data: return True
        return len(self.data)==0

    def nodes(self):
        return list(self.data.keys())

    def node_MD5(self):
        """Return a MD5 representing sort node ID"""
        X=sorted(self.nodes())
        if util.is_python3():
            return hashlib.md5((" ".join(X)).encode("utf-8")).digest()
        else:
            return hashlib.md5(" ".join(X)).digest()

    def __contains__(self, node):
        return node in self.data

    def nof_nodes(self):
      return len(self.nodes())

    def nof_edges(self):
        n=0
        for k,neighbors in self.data.items():
            for x in neighbors:
                if k>x: n+=1
        return n

    def are_neighbors(self, node1, node2):
        return node2 in self.data.get(node1, {})

    def neighbors(self, node):
        if self.is_node(node):
            return list(self.data[node].keys())
        return []

    # need test
    def remote_neighbors(self, S_node, i_hops=1, include_start=False):
        if type(S_node) is not list:
            S_node=[S_node]
        S_node=[x for x in S_node if x in self] # in case node is not in network
        if not S_node: return {}
        c_dist={ node:0 for node in S_node }
        if i_hops is None: i_hops=self.nof_nodes() # infinite hops
        if i_hops==0:
            if include_start:
                return c_dist
            return {}
        S_frontier=S_node
        for i in range(i_hops):
            S_frontier=[y for x in S_frontier for y in self.neighbors(x) if y not in c_dist]
            if not S_frontier: break
            for x in S_frontier:
                c_dist[x]=i+1
        if not include_start:
            for node in S_node:
                del c_dist[node]
        return c_dist

    def nodes_reachable_by(self, S_start, S_end=None, i_hops=1):
        '''identify the nodes in S_start and S_end that can be connected within i_hops
        return (SubsetOfStart, IntermediateNodes, SubsetofEnd)'''
        S_end = S_end or S_start
        c_start=set(S_start)
        c_end=set(S_end)
        S_share=list(c_start.intersection(c_end))
        if i_hops==0:
            return (S_share, [], S_share)
        C_neighbors1=self.remote_neighbors(S_start, i_hops=1, include_start=True)
        if i_hops==1:
            S_e=[ x for x in S_end if x in C_neighbors1 ]
            C_neighbors2=self.remote_neighbors(S_end, i_hops=1, include_start=True)
            S_b=[ x for x in S_start if x in C_neighbors2 ]
            return (S_b, [], S_e)
        # i_hops > 1
        (S_b, S_m, S_e) = self.nodes_reachable_by(list(C_neighbors1.keys()), S_end, i_hops-1)
        c_m=set(S_m)
        S_m2=[ x for x in S_b if x not in c_start and x not in c_m ]
        S_b=[ x for x in S_b if x in c_start ]
        S_m.extend(S_m2)
        return (S_b, S_m, S_e)

    def paths_between(self, S_start, S_end=None, i_hops=1, l_indirect=False):
        (S_b, S_m, S_e) = self.nodes_reachable_by(S_start, S_end, i_hops)
        if not l_indirect:
            return self.network(S_b+S_m+S_e)
        else: # set all indirect path to ppp
            net=self.subnetwork(S_b+S_e) # get all direct edges
            S_node=net.nodes()
            t_edge=net.T_edge
            rows=[]
            S_new_node=[]
            for s_a in S_node:
                for s_b in S_node:
                    if s_a>=s_b: continue
                    if self.are_neighbors(s_a, s_b): continue
                    if (self.shortest_path(s_a, s_b, i_hops)):
                        rows.append({'Gene_A':s_a, 'Gene_B':s_b, 'TYPE':'Indirect'})
                        S_new_node.extend([s_a, s_b])
            if rows:
                S_node=list(set(S_node+S_new_node))
                t=pd.DataFrame(data={'Gene':S_node})
                t_node=pd.merge(t, self.T_node, left_on='Gene', right_on='Gene', how='left')
                t_edge=pd.concat([t_edge, pd.DataFrame(rows)], axis=0, ignore_index=True)
                tmp_edge=self.T_edge.copy()
                tmp_edge=tmp_edge.drop('TYPE', axis=1)
                t_edge=pd.merge(t_edge, tmp_edge, left_on=['Gene_A','Gene_B'], right_on=['Gene_A','Gene_B'], how='left')
                net=Network(t_edge, T_node=t_node, allow_indirect=True)
            return net

    def shortest_path(self, s_start, s_end, max_hops=None):
        if s_start==s_end:
            return [(s_start)]
        if max_hops is None: max_hops=self.nof_nodes() # infinite hops
        if max_hops==0: return []
        if self.are_neighbors(s_start, s_end):
            return [(s_start, s_end)]
        elif max_hops==1:
            return []
        S_neighbor=self.neighbors(s_start)
        S_path=[]
        for s in S_neighbor:
            path=self.shortest_path(s, s_end, max_hops-1)
            if path: S_path.extend(path)
        if not S_path: return []
        hops=min([len(x) for x in S_path])
        S_path=[(s_start,)+path for path in S_path if len(path)==hops]
        return S_path

    def shortest_path_to_targets(self, s_start, S_end, max_hops=None):
        S_path=[]
        for s_end in S_end:
            path=self.shortest_path(s_start, s_end, max_hops)
            if path: S_path.extend(path)
        if not S_path: return []
        S_path=sorted(S_path, key=lambda x: len(x))
        hops=len(S_path[0])
        for i in range(1, len(S_path)):
            if len(S_path[i])!=hops: break
        S_path=S_path[:i]
        return S_path

    def all_neighbors(self, node):
        """Pull out all nodes connectable to node, used to decompose a network into components"""
        n=self.nof_nodes()
        c=self.remote_neighbors(node, n)
        return list(c.keys())

    def delete_node(self, node):
        S=self.neighbors(node)
        del self.data[node]
        [ self.data[s].__delitem__(node) for s in S ]
        DEL=[node]
        for s in S:
            if s in self.data and len(self.data[s])==0:
                del self.data[s]
                DEL.append(s)
        self.T_node=self.T_node[ ~ self.T_node.Gene.isin(DEL)].copy()
        self.T_edge=self.T_edge[~ ((self.T_edge.Gene_A.isin(DEL))|(self.T_edge.Gene_B.isin(DEL))) ].copy()
        if self.allow_indirect and sum(self.T_edge['TYPE'].isin('Indirect'))==0:
            self.allow_indirect=False

    def degree(self, node):
        return len(self.neighbors(node))

    def subnetwork(self, S_node, l_keep_as_singleton=False):
        """if l_keep_as_singleton is true, the node, if was in the original network, will be retained as a singleton node, even if it has no more neighbors"""
        S_node=[x for x in S_node if x in self.data]
        c_keep=set(S_node)
        tmp=pd.DataFrame({"Gene_A":list(c_keep)})
        T_edge=self.T_edge.merge(tmp, left_on='Gene_A', right_on='Gene_A')
        tmp.rename2({"Gene_A":"Gene_B"})
        T_edge=T_edge.merge(tmp, left_on='Gene_B', right_on='Gene_B')
        if l_keep_as_singleton:
            S_keep=set(T_edge.Gene_A)|set(T_edge.Gene_B)
            S_add=[x for x in c_keep if x not in S_keep ]
            if len(S_add):
                tmp=pd.DataFrame({"Gene_A":S_add, "Gene_B":S_add})
                tmp['SCORE']=1
                tmp['TYPE']='Direct'
                T_edge=pd.concat([T_edge, tmp], ignore_index=True)
        tmp=pd.DataFrame({"Gene": list(set(T_edge.Gene_A)|set(T_edge.Gene_B))})
        T_node=self.T_node.merge(tmp, left_on="Gene", right_on="Gene")
        return Network(T_edge, allow_indirect=self.allow_indirect, T_node=T_node)

    def combine_network(self, others=None):
        '''Merge this network with another network.  All nodes and edges are kept.'''
        if type(others) != list:
            others=[others]
        data1=[self.T_node]
        data2=[self.T_edge]
        for another_net in others:
            data1.append(another_net.T_node)
            data2.append(another_net.T_edge)
            if another_net.allow_indirect:
                self.allow_indirect=True
        self.T_node=pd.concat(data1, ignore_index=True)
        self.T_node.drop_duplicates(['Gene'], inplace=True)
        self.T_edge=pd.concat(data2, ignore_index=True)
        # remove indirect edge if direct already exists
        self.T_edge.sort_values(['Gene_A','Gene_B','TYPE','SCORE'], ascending=[True,True,True,False], inplace=True)
        self.T_edge.drop_duplicates(['Gene_A','Gene_B'], inplace=True)
        self.from_table(self.T_edge)
        # since TYPE = direct will appear first in the sorted table

    # a network may consists of multiple disconnected components (or after haircut)
    def decompose(self):
        subnets=[]
        all_nodes=self.nodes()
        c_seen={}
        for node in all_nodes:
            if node in c_seen: continue
            S=[node]+self.all_neighbors(node)
            for s in S:
                c_seen[s]=True
            subnets.append(self.subnetwork(S))
        L=[(x, x.nof_nodes(), x.nof_edges()) for x in subnets]
        # return the biggest network first
        L=sorted(L, key=lambda x: (-x[1], -x[2]))
        L=[ x[0] for x in L ]
        return L

    def __str__(self):
        s="Nodes: "+str(self.nof_nodes())+"\n"
        s+="Edges: "+str(self.nof_edges())+"\n"
        if util.is_python3():
            import io
            output = io.StringIO()
        else:
            import cStringIO
            output=cStringIO.StringIO()
        self.T_edge.to_csv(output, index=False, sep="\t")
        return s+output.getvalue()+"\n"

    def to_xgmml(self, s_file, l_remove_loop=False):
        xg=XGMML()
        xg.set_tables(self.T_node, self.T_edge)
        xg.save(s_file, self.name, l_remove_loop=l_remove_loop)

    def to_sif(self, s_file):
        s_file, s_ext=os.path.splitext(s_file)
        self.T_node.to_csv(s_file+'.noa', sep="\t", index=False)
        self.T_edge.to_csv(s_file+'.sif', sep="\t", index=False)

    @staticmethod
    def set_XY_(net, t_xy):
        # table t_xy has columns "id", "x", "y", "Gene", typically produced by cytoscape.py cynet_get_XY
        net.T_node.drop([x for x in ['graphics_x','graphics_y'] if x in net.T_node.header() ], axis=1, inplace=True)
        if "id" in t_xy.header():
            t_xy.drop("id", axis=1, inplace=True)
        t_xy.rename2({'x':'graphics_x', 'y':'graphics_y'})
        net.T_node=net.T_node.merge(t_xy, left_on='Gene', right_on='Gene', how='left')
        net.T_node['graphics_x'].fillna(0)
        net.T_node['graphics_y'].fillna(0)

    def set_XY(self, t_xy):
        return Network.set_XY_(self, t_xy)

    @staticmethod
    def to_json_(net):
        """To Cytoscape json format"""
        data={'data': {}, 'elements': { 'nodes': [], 'edges': [] }}
        data['data']['name']=net.name
        data['data']['selected']=True
        data["format_version"]="1.0"
        data["generated_by"]="xgmml.py"
        data["target_cytoscapejs_version"]="~2.1"

        id_type='string'
        if len(net.T_node)==sum(net.T_node['Gene'].astype(str).apply(lambda x:re.match(r'-?\d+$', x) is not None)):
            id_type='integer' # ID is actually integer

        id=0
        c_id={}
        for i,r in net.T_node.iterrows():
            id+=1
            c_id[r['Gene']]=str(id) if id_type=='string' else str(r['Gene'])
        #else:
        #    for i,r in self.T_node.iterrows():
        #        c_id[r['Gene']]=r['id']

        nodes=data['elements']['nodes']
        S=[x.upper() for x in net.T_node.header()]
        S_label=[x for x in ['Symbol','SYMBOL','canonicalName','Name','NAME','Label','LABEL'] if x in S]
        for i,r in net.T_node.iterrows():
            c=r.to_dict()
            c['id']=c['name']=c_id[r['Gene']]
            #if len(S_label):
            #    c['label']=r[S_label[0]]
            #c={k:v for k,v in c.items() if k not in S_REMOVE_NODE}
            c_node={'data': c}
            if 'graphics_x' in c and 'graphics_y' in c:
                c_node['position']={'x':float(c['graphics_x']), 'y':float(c['graphics_y']) }
                c.pop('graphics_x', None)
                c.pop('graphics_y', None)
                c.pop('graphics_type', None)
            nodes.append(c_node)
        edges=data['elements']['edges']
        #cnt=0
        for i,r in net.T_edge.iterrows():
            c=r.to_dict()
            #cnt+=1
            #c['id']=c['label']=str(cnt)
            ##c['source']=c_id[r['Gene_A']]
            ##c['target']=c_id[r['Gene_B']]
            ia=r['Gene_A'] if r['Gene_A'] in c_id.keys() else r['Name_A']
            ib=r['Gene_B'] if r['Gene_B'] in c_id.keys() else r['Name_B']
            c['source']=c_id[ia]
            c['target']=c_id[ib]
            c['interaction']=r['TYPE']
            #c['label']=c['Name_A']+" ("+c['TYPE']+") "+c['Name_B']
            #c={k:v for k,v in c.items() if k not in S_REMOVE_EDGE}
            edges.append({'data': c})
        return data

    def to_json(self):
        return Network.to_json_(self)

    @staticmethod
    def from_json(data):
        nodes=data['elements']['nodes']
        edges=data['elements']['edges']
        s_name=data['data']['name']
        S_n=[node['data'] for node in nodes]
        S_e=[edge['data'] for edge in edges]
        T_node=pd.DataFrame(S_n)
        T_edge=pd.DataFrame(S_e)
        S=T_node.header()
        if 'Gene' not in S:
            T_node['Gene']=T_node['id']
        if 'Symbol' not in T_node.header():
            T_node.rename2({'canonicalName':'Symbol'})
        S=set(T_node.header())
        T_node.drop([x for x in ['name','selected','shared_name'] if x in S], axis=1, inplace=True)
        S=T_edge.header()
        if 'Gene_A' not in S:
            T_edge.rename2({'source':'Gene_A'})
        if 'Gene_B' not in S:
            T_edge.rename2({'target':'Gene_B'})
        if 'TYPE' not in S:
            T_edge.rename2({'interaction':'TYPE'})
        S=T_edge.header()
        T_edge.drop([x for x in ['source','target','interaction','SUID','canonicalName','id','label','selected','shared_interaction','shared_name'] if x in S], axis=1, inplace=True)
        return Network(T_edge, name=s_name, T_node=T_node)

    def add_node_attr(self, t_attr, s_key="Gene"):
        t_attr[s_key]=t_attr[s_key].astype(str)
        self.T_node=pd.merge(self.T_node, t_attr, left_on="Gene", right_on=s_key, how="left")

    @staticmethod
    def add_node_degree(netwk):
        """Add a DEGREE node attribute to the network"""
        netwk.T_node["DEGREE"]=[ netwk.degree(x) for x in netwk.T_node["Gene"] ]

    def add_a_node_attr(self, s_attr, c_attr, s_NULL=""):
        """c_attr: dict key is Gene, value is type
        s_attr: attribute name"""
        self.T_node[s_attr]=self.T_node['Gene'].apply(lambda x: c_attr.get(x, s_NULL))

    def del_node_attr(self, S_attr=None):
        if S_attr is None: return
        S_attr=[x for x in S_attr if x in self.T_node.header()]
        if len(S_attr):
            self.T_node.drop(S_attr, axis=1, inplace=True)

    def del_edge_attr(self, S_attr=None):
        if S_attr is None: return
        S_attr=[x for x in S_attr if x in self.T_edge.header()]
        if len(S_attr):
            self.T_edge.drop(S_attr, axis=1, inplace=True)

    @staticmethod
    def overconnected(input_network, S_node, min_links=2, p_cutoff=0.01, min_enrichment=0):
        """Find nodes that overconnect to S_node, with at least min_links, hyper p<=0.01
        This implements GeneGo R overconnected method"""
        # original net
        S_node=[x for x in S_node if x in input_network ]
        N=input_network.nof_nodes()
        n2=len(S_node)
        S_new=list(input_network.remote_neighbors(S_node, i_hops=1, include_start=True).keys())
        import stats
        S_old=set(S_node)
        data=[]
        for x in S_new: # also score existing old nodes
            n1=input_network.degree(x) # total nof links in input_network
            n=len(set(input_network.neighbors(x)) & S_old)
            if n<min_links: continue # too few
            logP=np.log10(max(stats.hyper(n, N, n1, n2), 1e-100))
            if logP>np.log10(p_cutoff): continue
            ef=n*N*1.0/n2/n1
            if ef<=min_enrichment: continue;
            z=stats.ZScore_GeneGo(n, N, n1, n2)
            data.append({"Node":x, "#PPI":N, "#Hits":n2, "#LinksInPPI":n1, "#LinksToHits":n, "Enrichment":ef, "Z-score":z, "LogP":logP})
        if len(data):
            t=pd.DataFrame(data)
            t['Log(q-value)']=np.log10(np.clip(
                stats.adjust_p(np.power(10, t.LogP.values), N=len(S_new), method="BH"), 1e-100, 1.0))
            t.sort_values(['LogP','Enrichment'], ascending=[True,False], inplace=True)
            t['InHits']=t.Node.apply(lambda x: x in S_old)
            t=t.reindex(columns=["Node","InHits","#PPI","#Hits","#LinksInPPI","#LinksToHits","Enrichment","Z-score","LogP","Log(q-value)"])
            return t
        return None

    @staticmethod
    def propagation(input_network, prior = [], alpha = 0.8, l1norm_cutoff = 1E-6, remove_source = False, smoothing = False, col_name = None, use_edge_weight=False, l_Otsu_split=False):
        """ Network propagation algorithm [1].
        prior - list of seeds (subset of values in T_node['Gene'] column)
        alpha - (o <= alpha <= 1) - weights the importance of 2 constraints.
                Values closer to 0 conservative to a more conservative behaviour.
        l1norm_cutoff - convergence criteria
        remove_source - if True, source will be removed after the first iterations.
        In this case, the final score doe not depend on alpha, end equivalent to the case of alpha = 0.
        smoothing - if True, a smoothing step will be applied. Not recommended.
        col_name - if a string is supplied, a new attribute will be added to the nodes table,
                   the value of the attribute is the score.


        Returns a tuple:
            (1) dictionary of scores, where keys are node ids, and values are scores.
            (2) network with a new node attribute, id col_name is not None,
                None otherwise.

        [1] Vanunu et al, Associating genes and protein complexes with disease via network propagation, 2010.

        4/26/2015: Yingyao: when l_Otsu_split is True, we first limit the size of node to 3*#original nodes,
        then apply Otsu's threshold method to split nodes into two groups based on the scores
        This seems to be a reasonable way to separate the dense nodes from the rest
        """
        network = Network(input_network)
        # remove orphan nodes
        for node in network.T_node.Gene.values.tolist():
            if len(network.neighbors(node)) == 0:
                network.delete_node(node)
        network.T_edge.reset_index (drop = True, inplace = True)
        network.T_node.reset_index (drop = True, inplace = True)
        n_nodes = len(network.T_node)

        # map node names to index and back
        idx_to_name = dict(zip(network.T_node.index.values.tolist(), network.T_node.Gene.values.tolist()))
        name_to_idx = dict(zip(network.T_node.Gene.values.tolist(), network.T_node.index.values.tolist()))

        # initialize graph matrix and other values
        matrix = dok_matrix((n_nodes,n_nodes), dtype = float)
        if use_edge_weight:
            col_idx=network.T_edge.col_index('SCORE')
            if col_idx<0: # SCORE is nto found
                use_edge_weight=False
        for _, row in network.T_edge.iterrows():
            i = name_to_idx[row['Gene_A']]
            j = name_to_idx[row['Gene_B']]
            if use_edge_weight:
                matrix[i, j]=matrix[j, i]=row['SCORE']
            else:
                # no weights assumed so far. 1 should be replaced by the weight if available.
                matrix[i, j] = 1
                matrix[j, i] = 1
        diag = matrix.sum(1)
        norm = np.sqrt(diag * diag.transpose())
        w_prime = matrix.todense() / norm

        # propagation routine
        # (it's useful to keep it as a separate routine in case we
        # decide to add a permutation test later -
        # we would just need to shuffle prior_values and call the function)

        def propagate(prior_values):
            ft1 = prior_values.copy()
            y = ft1 * (1.0 - alpha)
            ft = np.zeros(n_nodes)
            continue_propagation = True
            while continue_propagation:
                ft = alpha * ft1 * w_prime + y
                if (abs(ft - ft1)).sum() < l1norm_cutoff:
                    continue_propagation = False
                else:
                    ft1 = ft
                if remove_source:
                    y = ft1 * (1.0 - alpha)
            if smoothing:
                ft = 1. * ft1 * w_prime
            return np.array(ft)[0]

        # running propagation
        prior_values = np.zeros(n_nodes)
        for node in prior:
            if node in name_to_idx: # only if node in the network
                prior_values[name_to_idx[node]] = 1
        score = propagate(prior_values)
        final_score = dict(zip(input_network.T_node.Gene.values, np.zeros(len(input_network.T_node))))

        for idx in range(len(score)):
            final_score[idx_to_name[idx]] = score[idx]

        # identify densely connected nodes, 4/26/2015, YZ
        if l_Otsu_split:
            import stats
            R=np.array(list(final_score.values()))
            R.sort()
            R=R[::-1]
            R=R[: len(prior)*3]
            cutoff_score=stats.Otsu_threshold(R)
            final_score={k:v for k,v in final_score.items() if v>=cutoff_score}
            input_network=input_network.subnetwork(list(final_score.keys()))

        # adding a new attribute to the input network
        if not col_name is None:
            df = pd.DataFrame({col_name : list(final_score.values()), "Gene" : list(final_score.keys())})
            input_network.add_node_attr(df, s_key = "Gene")
            return final_score, input_network

        return final_score, None

if __name__=="__main__":
    x=XGMML()
    x.parse('~/Cytoscape/RM_.xgmml')
    #print x.T_node
    #print x.T_edge
    ## pretend to add expression
    T_attr=x.T_node.reindex(columns=["Gene"])
    T_attr['Activity']=np.random.randn(len(T_attr))
    x.add_node_attr(T_attr)
    x.save("test.xgmml", "mygraph")

    net=Network("RM_.xgmml")
    print(net)
    print(net.neighbors('2'))
    print(net.subnetwork(['2','6275', '10525', '348', '7276']))
    print(net.remote_neighbors('2',3))
    print(net.all_neighbors('2'))
    S=net.decompose()
    print("Subnetworks:\n")
    for s in S:
        print(s)

