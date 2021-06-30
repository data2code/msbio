#!/usr/bin/env python
from __future__ import absolute_import
from __future__ import print_function
import numpy as np
import pandas as pd
import re
import util
import os
import entrez as ez
from xgmml import *
import parallel
import db
from six.moves import range
import setting

class MCODECluster(Network):

    def __init__(self, network, seedNode=None, score=0.0):
        self.seedNode=seedNode
        self.score=score or 0.0
        super(MCODECluster, self).__init__(network)

    def __str__(self):
        s="SeedNode: "+self.seedNode+"\n"
        s+="Score: "+str(self.score)+"\n"
        s+=super(MCODECluster, self).__str__()
        return s

class MCODE(Network):

    #the parameters used for this instance of the algorithm
    params = {
        'includeLoops': False,
        'degreeCutoff': 2,
        'kCore': 2, # kCore must be greater than 1
        'maxDepthFromStart': 100,
        'nodeScoreCutoff': 0.2,
        'fluff': False,
        'haircut': True,
        'fluffNodeDensityCutoff': 0.1
    }

    def get_max_score(self):
        return self.C_nodeByScore[0][0]

    @staticmethod
    def calc_density(gpInputGraph, includeLoops=False):
        if (gpInputGraph.is_empty()): return -1.0
        loopCount=0
        if includeLoops:
            S_node=gpInputGraph.nodes()
            for node in S_node:
                if (gpInputGraph.are_neighbors(node, node)): loopCount+=1
        n=gpInputGraph.nof_nodes()
        possibleEdgeNum = n**2
        actualEdgeNum = gpInputGraph.nof_edges() - loopCount
        return actualEdgeNum*1.0/possibleEdgeNum

    def score_network(self, network):
        numNodes = network.nof_nodes()
        density = MCODE.calc_density(network, self.params['includeLoops'])
        score = density * numNodes
        return score

    @staticmethod
    def get_KCore(gpInputGraph, k):
        if (gpInputGraph is None or gpInputGraph.is_empty()):
            util.error_msg("GetKCore(): no input network!")
        #filter all nodes with degree less than k until convergence
        firstLoop = True
        gpOutputGraph = None
        while True:
            numDeleted = 0
            S_node = gpInputGraph.nodes()
            alCoreNodes=[node for node in S_node if len(gpInputGraph.data[node]) >= k]
            # ZHOU 3/16/2018
            if len(alCoreNodes)<k: return None
            #
            if (len(S_node)>len(alCoreNodes) or firstLoop):
                gpOutputGraph = gpInputGraph.subnetwork(alCoreNodes)
                if (gpOutputGraph.is_empty()): return None
                #iterate again, but with a new k-core input graph
                gpInputGraph = gpOutputGraph
                firstLoop = False
            else:
                #stop the loop
                break
        return gpOutputGraph

    def get_highest_KCore(self, gpInputGraph):
        S_md5=[]
        if self.l_cache:
            s_md5=gpInputGraph.node_MD5()
            S_md5.append(s_md5)
            if s_md5 in self.cache_kcore:
                k=self.cache_kcore[s_md5]
                gpPrevCore=MCODE.get_KCore(gpInputGraph, k)
                return {'k':self.cache_kcore[s_md5], 'network':gpPrevCore}

        gpCurCore=gpPrevCore=None
        ### ZHOU 3/16/2018 tries to speed things up, find the max possible k
        R_degree=np.array([gpInputGraph.degree(x) for x in gpInputGraph.nodes()])
        if len(R_degree)==0: return {'k':0, 'network':gpPrevCore}
        # https://stackoverflow.com/questions/26984414/efficiently-sorting-a-numpy-array-in-descending-order/26984520
        R_degree[::-1].sort()
        # if k is a max, then there must be at least k nodes with degrees >= k
        lb=ub=R_degree[-1]
        if lb==0: exit()
        # all degrees in tmp are possible k, so we use max(tmp)
        tmp=R_degree[R_degree<=np.arange(1,len(R_degree)+1)]
        if len(tmp)>0:
            ub=max(tmp)
            # Be aware if candidates [5 5 5 2 2 2], max k can be 3, not in the candidate list
            # the above tmp will give ub=2, so we need to find the closest failure
            tmp2=R_degree[R_degree>ub]
            if len(tmp2): ub=min(tmp2)-1
        gpPrevCore=MCODE.get_KCore(gpInputGraph, lb)
        while (lb<ub): # lb is already a solution, ub has not been explored yet
            k=max(lb+1, int(lb*0.3+ub*0.7)) # empirically seems a bit better to bias towards ub
            #print("Try: ", k, "[", lb, ub, "]")
            gpCurCore = MCODE.get_KCore(gpInputGraph, k)
            if gpCurCore is None or gpCurCore.is_empty():
                #print("fail")
                ub=k-1
            else:
                gpPrevCore = gpCurCore
                gpInputGraph=gpCurCore # let's shrink the search space
                # use cache to avoid recomputing.
                if self.l_cache:
                    s_md5=gpInputGraph.node_MD5()
                    S_md5.append(s_md5)
                    if s_md5 in self.cache_kcore:
                        k=self.cache_kcore[s_md5]
                        gpPrevCore=MCODE.get_KCore(gpInputGraph, k)
                        lb=ub=k
                lb=k
        #while True:
        #    gpCurCore = MCODE.get_KCore(gpInputGraph, i)
        #    if gpCurCore is None or gpCurCore.is_empty(): break
        #    gpPrevCore = gpCurCore
        #    gpInputGraph = gpCurCore
        #    i+=1
        #    print("try ",i)
        #k = i-1
        #print("Answer ", lb)
        if self.l_cache:
            self.cache_kcore.update({x:lb for x in S_md5})
        return {'k':lb, 'network':gpPrevCore}
        #in the last iteration, gpCurCore is null (loop termination condition)

    #@staticmethod
    def calc_node_info(self, node, degreeCutoff=None):
        k = self.degree(node)
        neighbors = self.neighbors(node)
        s_md5=""
        #print("::::",node, ":::", k, "::::")
        #sw=util.StopWatch()
        if (k < 2):
            nodeInfo = NodeInfo()
            if (k == 1):
                nodeInfo.coreLevel = 1
                nodeInfo.coreDensity = 1.0
                nodeInfo.density = 1.0
                nodeInfo.numNodeNeighbors = len(neighbors); #########
                nodeInfo.nodeNeighbors = neighbors; ####
                # why ignore neighbor when k==1 in the original code???
        else:
            gpNodeNeighborhood = self.subnetwork(neighbors+[node])
            #sw.check('subnetwork')
            if (gpNodeNeighborhood.is_empty()):
                util.error_msg("In calc_node_info(): gpNodeNeighborhood was None.")
            #calculate the node information for each node
            if self.l_cache:
                s_md5=gpNodeNeighborhood.node_MD5()
                if s_md5 in self.cache_info:
                    #self.hit+=1
                    nodeInfo=self.cache_info[s_md5].clone()
                    nodeInfo.nodeNeighbors=neighbors
                    return nodeInfo
            nodeInfo = NodeInfo()
            #density
            nodeInfo.density = MCODE.calc_density(gpNodeNeighborhood, self.params['includeLoops'])
            #w.check('density')
            nodeInfo.numNodeNeighbors = len(neighbors)
            #calculate the highest k-core
            c = self.get_highest_KCore(gpNodeNeighborhood)
            #w.check('kcore')
            k = c['k']
            gpCore = c['network']
            nodeInfo.coreLevel = k
            if (gpCore is not None and not gpCore.is_empty()):
                nodeInfo.coreDensity = MCODE.calc_density(gpCore, self.params['includeLoops'])
            #w.check('cacl_density')
            #record neighbor array for later use in cluster detection step
            nodeInfo.nodeNeighbors = neighbors
        if degreeCutoff: nodeInfo.score_node(degreeCutoff)
        if self.l_cache:
            self.cache_info[s_md5]=nodeInfo
        return nodeInfo

    def score_graph(self):
        self.C_nodeInfo = {}
        self.C_nodeByScore = []
        S_node = self.nodes()
        rows=[]
        #sw=util.StopWatch()
        if self.CPU<=1:
            for node in S_node:
                nodeInfo = self.calc_node_info(node, self.params['degreeCutoff'])
                self.C_nodeInfo[node]=nodeInfo
                rows.append({'Node':node, 'Score':nodeInfo.score, 'Density':nodeInfo.density, 'numNodeNeighbors':nodeInfo.numNodeNeighbors})
        else:
            def f(X):
                return self.calc_node_info(X[0], X[1])
            #mp=parallel.MP()
            #mp.start(f, n_CPU=self.CPU)
            L=[ (x, self.params['degreeCutoff']) for x in S_node]
            out=parallel.parmap(f, L, n_CPU=self.CPU)
            #out=mp.map(L)
            for i,node in enumerate(S_node):
                self.C_nodeInfo[node]=out[i]
                rows.append({'Node':node, 'Score':out[i].score, 'Density':out[i].density, 'numNodeNeighbors':out[i].numNodeNeighbors})
        #sw.check('Done scoring')
        t=pd.DataFrame(rows)
        t=t.sort_values(['Score', 'Density', 'numNodeNeighbors', 'Node'], ascending=[False, False, False, True])
        grps=t.groupby(by='Score')
        self.C_nodeByScore=[ (score,list(grp['Node'])) for score,grp in grps]
        self.C_nodeByScore.sort(key=lambda x: x[0])
        self.C_nodeByScore.reverse()

    def get_cluster_core_internal(self, startNode, c_nodeSeen, startNodeScore, currentDepth, myCluster, nodeScoreCutoff, maxDepthFromStart):
        #base cases for recursion
        if (startNode in c_nodeSeen): return
        c_nodeSeen[startNode]=True
        if (currentDepth > maxDepthFromStart): return
        #don't exceed given depth from start node

        #Initialization
        neighbors = self.C_nodeInfo[startNode].nodeNeighbors
        #neighbors.sort()
        #print "A:"+startNode
        #print c_nodeSeen
        for node in neighbors:
            #go through all currentNode neighbors to check their core density for cluster inclusion
            #print "Neigh:"+node
            if (node in c_nodeSeen): continue
            if (self.C_nodeInfo[node].score >= (startNodeScore - startNodeScore * nodeScoreCutoff)):
                myCluster.append(node)
                #try to extend cluster at this node
                self.get_cluster_core_internal(node, c_nodeSeen, startNodeScore, currentDepth + 1, myCluster, nodeScoreCutoff, maxDepthFromStart)

    def get_cluster_core(self, startNode, c_nodeSeen, nodeScoreCutoff, maxDepthFromStart):
        myCluster = []
        self.get_cluster_core_internal(startNode, c_nodeSeen, self.C_nodeInfo[startNode].score, 1, myCluster, nodeScoreCutoff, maxDepthFromStart)
        return self.subnetwork(myCluster+[startNode])

    def fluff_cluster_boundary(self, myCluster, c_nodeSeen):
        #create a temp list of nodes to add to avoid concurrently modifying 'cluster'
        nodesToAdd = []
        #Keep a separate internal nodeSeenHashMap because nodes seen during a fluffing should not be marked as permanently seen,
        #they can be included in another cluster's fluffing step.
        c_nodeSeenInternal = {}
        #add all current neighbour's neighbours into cluster (if they have high enough clustering coefficients) and mark them all as seen
        S_node=myCluster.nodes()
        for node in S_node:
            neighbors = self.C_nodeInfo[node].nodeNeighbors
            for nb in neighbors:
                if (nb in c_nodeSeen): continue
                if (nb in c_nodeSeenInternal): continue
                if (self.C_nodeInfo[nb].density > self.params['fluffNodeDensityCutoff']):
                    nodesToAdd.append(nb)
                    c_nodeSeenInternal[nb]=True

        #Add fluffed nodes to cluster
        if (len(nodesToAdd)>0):
            return self.subnetwork(S_node+nodesToAdd)
        return myCluster

    def filter_cluster(self, gpClusterGraph):
        if (gpClusterGraph.is_empty()): return True
        #filter if the cluster does not satisfy the user specified k-core
        gpCore = MCODE.get_KCore(gpClusterGraph, self.params['kCore'])
        if (gpCore is None or gpCore.is_empty()): return True
        return False

    @staticmethod
    def haircut_cluster(myCluster):
        #get 2-core
        gpCore = MCODE.get_KCore(myCluster, 2)
        if (gpCore is not None or not gpCore.is_empty()):
            #clear the cluster and add all 2-core nodes back into it
            #must add back the nodes in a way that preserves gpInputGraph node indices
            # we cannot do myCluster=gpCore, which will not change myCluster outside
            return gpCore
        return myCluster

    def find_clusters(self, l_decompose=True, l_optimized=True):
        if (self.is_empty()):
            util.error_msg("In find_Clusters(): input network is empty!")
        if (not len(self.C_nodeInfo.keys()) or not len(self.C_nodeByScore)):
            util.error_msg("In find_Clusters(): C_nodeInfo or C_nodeByScore is None.")
        C_results=[]
        cnt=0
        #initialization
        c_nodeSeen= {} #key is nodeIndex, value is true/false
        c_nodeSeenSnapshot={}
        findingTotal = len(self.C_nodeInfo.keys())
        rows=[]
        for score,alNodesWithSameScore in self.C_nodeByScore:
            if not l_optimized or len(alNodesWithSameScore)<=1:
                for currentNode in alNodesWithSameScore:
                    if currentNode in c_nodeSeen: continue
                    alCluster = self.get_cluster_core(currentNode, c_nodeSeen, self.params['nodeScoreCutoff'], self.params['maxDepthFromStart'])
                    if (alCluster is not None and not alCluster.is_empty()):
                        #make sure seed node is part of cluster, if not already in there
                        if (not self.filter_cluster(alCluster)):
                            if (self.params['haircut']): alCluster=MCODE.haircut_cluster(alCluster)
                            if (self.params['fluff']): alCluster=self.fluff_Cluster_boundary(alCluster, c_nodeSeen)
                            if l_decompose:
                                c_components=alCluster.decompose()
                            else:
                                c_components=[alCluster]
                            for comp in c_components:
                                cnt+=1
                                score=self.score_network(comp)
                                C_results.append(MCODECluster(comp, currentNode, score))
                                rows.append({'ID':cnt, 'Score':score, 'NofNode':comp.nof_nodes(), 'SeedScore':self.C_nodeInfo[currentNode].score})
            else:
                def f(X):
                    tmp_rows=[]
                    c_stack={}
                    currentNode=X[0]
                    c_nodeSeenCopy=X[1].copy()
                    #if currentNode in c_nodeSeen: continue
                    alCluster = self.get_cluster_core(currentNode, c_nodeSeenCopy, self.params['nodeScoreCutoff'], self.params['maxDepthFromStart'])
                    if (alCluster is not None and not alCluster.is_empty()):
                        #make sure seed node is part of cluster, if not already in there
                        if (not self.filter_cluster(alCluster)):
                            if (self.params['haircut']): alCluster=MCODE.haircut_cluster(alCluster)
                            if (self.params['fluff']): alCluster=self.fluff_Cluster_boundary(alCluster, c_nodeSeenCopy)
                            if l_decompose:
                                c_components=alCluster.decompose()
                            else:
                                c_components=[alCluster]
                            for k,comp in enumerate(c_components):
                                score=self.score_network(comp)
                                tmp_rows.append({'ID':currentNode, 'Score':score, 'NofNode':comp.nof_nodes(), 'SeedScore':self.C_nodeInfo[currentNode].score, 'ComponentIndex':k})
                                c_stack[currentNode]={'nodeSeen':c_nodeSeenCopy, 'components':c_components}
                    return (tmp_rows, c_stack)

                while (len(alNodesWithSameScore)):
                    tmp_rows=[]
                    c_stack={}
                    L=[ (x, c_nodeSeen) for x in alNodesWithSameScore if x not in c_nodeSeen ]
                    #if self.CPU<=1:
                    #    out=[f(x) for x in L]
                    #else:
                    #    mp=parallel.MP()
                    #    mp.start(f, n_CPU=self.CPU)
                    #    out=mp.map(L)
                    out=parallel.parmap(f, L, n_CPU=self.CPU)
                    for X in out:
                        tmp_rows.extend(X[0])
                        c_stack.update(X[1])
                    tmp=pd.DataFrame(tmp_rows)
                    if len(tmp):
                        tmp=tmp.sort_values(['Score','NofNode','SeedScore', 'ID'], ascending=[False, False, False, True])
                        bestNode=tmp['ID'].iloc[0]
                        c_nodeSeen=c_stack[bestNode]['nodeSeen']
                        for comp in tmp_rows:
                            if comp['ID']!=bestNode: continue
                            compIdx=comp['ComponentIndex']
                            cnt+=1
                            C_results.append(MCODECluster(c_stack[bestNode]['components'][compIdx], bestNode, comp['Score']))
                            rows.append({'ID':cnt, 'Score':comp['Score'], 'NofNode':comp['NofNode'], 'SeedScore':self.C_nodeInfo[bestNode].score})

                        alNodesWithSameScore=[ x for x in alNodesWithSameScore if x !=bestNode]
                    else:
                        for x in c_stack:
                            for s in x.nodeSeen.keys():
                                c_nodeSeen[s]=True
                        alNodesWithSameScore=[]
        C_sorted=[]
        t=pd.DataFrame(rows)
        if len(t):
            t=t.sort_values(['Score','NofNode','SeedScore'], ascending=[False, False, False])
            for i in range(len(t)):
                C_sorted.append(C_results[t['ID'].iloc[i]-1])
        return C_sorted

    def to_MCODE_table(self, S_mcode_clusters):
        rows=[]
        for i,c in enumerate(S_mcode_clusters):
            S_nodes=c.nodes()
            for node in S_nodes:
                ty='Seed' if node==c.seedNode else 'Clustered'
                rows.append({'Cluster':i+1, 'Score':c.score, 'Type':ty, 'Gene':node})
        if len(rows)==0:
            return None
        t=pd.DataFrame(rows)
        if 'Symbol' in self.T_node.header():
            c_name={self.T_node['Gene'].iloc[i] : self.T_node['Symbol'].iloc[i] for i in range(len(self.T_node))}
            t['Symbol']=t['Gene'].map(c_name)
        return t

    @staticmethod
    def MCODE_label(network, s_col_name='MCODE_LABEL'):
        """Label nodes in the network by their MCODE cluster IDs, great for coloring nodes"""
        network=Network(network) # makes a copy
        L=network.decompose()
        c_attr={}
        for j,net in enumerate(L):
            mc=MCODE(net)
            mc.params['hariCut']=True
            components=mc.find_clusters(True, True)
            for i,c in enumerate(components):
                S_nodes=c.nodes()
                for x in S_nodes:
                    if x not in c_attr:
                        c_attr[x]="N%dC%d" % (j+1, i+1)
                    else:
                        c_attr[x]+=" N%dC%d" % (j+1, i+1)
        network.add_a_node_attr(s_col_name, c_attr)
        return network

    def __init__(self, network, n_CPU=0, l_cache=True):
        self.C_nodeInfo = None
        #key is the node name, value is a NodeInfo instance
        C_nodeByScore = None
        #a collection of array, {{Na, Nb}, {Nc}, {Nd,Ne} ...}, where nodes are sorted by descending score
        # nodes with the same NodeInfo.score are group in one array
        super(MCODE, self).__init__(network)
        #"Scoring all nodes in the network ..."
        self.CPU=n_CPU
        self.l_cache=l_cache
        #self.hit=0
        self.cache_info={}
        self.cache_kcore={}
        self.score_graph()
        #for c,v in self.C_nodeInfo.items():
        #    print c, v


class Cache(object):
    DATA_DIR=setting.ppi['DATA_DIR']
    ppi_data={'LOCAL':{}, 'GPDB':{}, 'HISTORY':{}}
    ppi_node={'LOCAL':{}, 'GPDB':{}, 'HISTORY':{}}
    ppi_edge={'LOCAL':{}, 'GPDB':{}, 'HISTORY':{}}
    CUTOFF_PHYS=132
    CUTOFF_COMB=187
    VERSION=setting.ppi.get('VERSION',2) # 1: without STRING DB, 2: with STRING DB, database scheme changed

    @staticmethod
    def gene2node(S_gene, con=None):
        if con is None: con=db.DB('METASCAPE')
        t_node=con.sql_in("SELECT gid Gene,source_id Symbol from gid2source_id t where gid in (", ") and t.id_type_id=1", util.rarray2iarray(S_gene))
        t_node['Gene']=t_node.Gene.astype(str)
        if len(S_gene)!=len(t_node):
            util.warn_msg("Strange, gene ID has no symbol?")
            t=pd.DataFrame({'Gene':list(S_gene)})
            t_node=t.merge(t_node, left_on='Gene', right_on='Gene', how='left')
            X=t_node.Symbol.isnull()
            #print(t_node.loc[X][:10])
            if X.any():
                t_node.loc[X,'Symbol']=t_node.loc[X,'Gene']
        return t_node

    @staticmethod
    def df2data(t, con=None):
        nodes=set(t.Gene_A)|set(t.Gene_B)
        data={ k:{} for k in nodes }
        [ (data[k].__setitem__(v,c) or data[v].__setitem__(k,c)) for k,v,c in zip(t.Gene_A, t.Gene_B, t.SCORE) ]
        return (data, Cache.gene2node(nodes, con=con))

    @staticmethod
    def get(l_use_GPDB=True, S_DB=None, tax_id=9606):
        """In VERSION=2, S_DB is a string, one of "PHYSICAL_CORE","PHYSICAL_ALL","COMBINED_CORE","COMBINED_ALL"
        getting a phyiscal db will populate both PHYSICAL_CORE and PHYSICAL_ALL
        getting a combined db will populate all four databases
        """
        S_DB=S_DB or Cache.get_DB(l_use_GPDB)
        if Cache.VERSION==1: # in version one we merge all db data in S_DB
            S_DB.sort()
            s_db=":".join(S_DB)
            if not (tax_id in Cache.ppi_data['HISTORY'] and s_db in Cache.ppi_data['HISTORY'][tax_id]):
                s_key=Cache.key(l_use_GPDB)
                Cache.load(tax_id=tax_id, l_use_GPDB=l_use_GPDB, S_DB=S_DB)
                data=None
                out_node=[]
                for x in S_DB:
                    #print ">>>>>>>>>", S_DB, x, Cache.ppi_data[s_key][tax_id].keys()
                    c=Cache.ppi_data[s_key][tax_id].get(x, {})
                    if data is None:
                        data=c
                    else:
                        for k in c.keys():
                            for v,score in c[k].items():
                                if k not in data:
                                    data[k]=c[k].copy()
                                else:
                                    data[k][v]=max(score, data[k].get(v,0))
                    out_node.append(Cache.ppi_node[s_key][tax_id].get(x, pd.DataFrame()))
                t_node=pd.concat(out_node, ignore_index=True)
                t_node.drop_duplicates('Gene', inplace=True)
                if tax_id not in Cache.ppi_data['HISTORY']:
                    Cache.ppi_data['HISTORY'][tax_id]={}
                    Cache.ppi_node['HISTORY'][tax_id]={}
                    Cache.ppi_edge['HISTORY'][tax_id]={}
                Cache.ppi_data['HISTORY'][tax_id][s_db]=data
                Cache.ppi_node['HISTORY'][tax_id][s_db]=t_node
        else: # In VERSION 2, each entry in S_DB is its own collection
            s_db=S_DB
            #print(tax_id, list(Cache.ppi_data['HISTORY'].keys()), list(Cache.ppi_data['HISTORY'][tax_id].keys()))
            if not (tax_id in Cache.ppi_data['HISTORY'] and s_db in Cache.ppi_data['HISTORY'][tax_id]):
                Cache.load(tax_id=tax_id, l_use_GPDB=True, S_DB=s_db)
        return (Cache.ppi_data['HISTORY'][tax_id][s_db], Cache.ppi_node['HISTORY'][tax_id][s_db], \
                Cache.ppi_edge['HISTORY'][tax_id].get(s_db, None))

    @staticmethod
    def info():
        for s_key in ('LOCAL','GPDB','HISTORY'):
            print(">Databases: %s" % s_key)
            for tax_id in Cache.ppi_data[s_key].keys():
                print("TAX_ID=%d (%s)" % (tax_id, ez.Cache.C_TAX_NAME.get(tax_id, "UNKNOWN")))
                for s_db in Cache.ppi_data[s_key][tax_id].keys():
                    print("Source: %s" % s_db)
                    print("PPI_DATA=%d" % len(Cache.ppi_data[s_key][tax_id][s_db]))
                    print("PPI_NODE=%d" % len(Cache.ppi_node[s_key][tax_id][s_db]))
                    print("PPI_EDGE=%d" % len(Cache.ppi_edge[s_key][tax_id][s_db]))
            print("")

    @staticmethod
    def unload(tax_id, l_use_GPDB):
        s_key=Cache.key(l_use_GPDB)
        if tax_id in Cache.ppi_data[s_key]:
            del Cache.ppi_data[s_key][tax_id]
            del Cache.ppi_node[s_key][tax_id]

    @staticmethod
    def key(l_use_GPDB):
        return 'GPDB' if l_use_GPDB else 'LOCAL'

    @staticmethod
    def get_DB(l_use_GPDB=True):
        if Cache.VERSION==1:
            DEFAULT_DB=["BioGrid","InWeb_IM","OmniPath"] if l_use_GPDB else ["BHMRRS","CORUM","Prolexys","Chanda"] # String
        else:
            DEFAULT_DB=setting.ppi.get('DEFAULT_DB', ["PHYSICAL_CORE","PHYSICAL_ALL","COMBINED_CORE","COMBINED_ALL"][2])
        return DEFAULT_DB

    @staticmethod
    def load(tax_id=9606, l_use_GPDB=True, S_DB=None, entrez=None):
        """tax_id is None, defaults to 9606, if 0, means load all supported species,
        entrez is only used in local mode to accelerate Symbol retrieval"""
        sw=util.StopWatch()
        if Cache.VERSION==2:
            if S_DB is None: S_DB="PHYSICAL_CORE"
            if type(S_DB)!=str: util.error_msg("S_DB must be a string in VERSION 2")
            s_db=S_DB
            fn=setting.ppi.get('STRING_PATH', os.path.join(os.path.dirname(__file__),"STRING/Interaction.csv.gz"))
            mydb=db.DB('METASCAPE')
            if tax_id==0:
                S_tax_id=ez.Cache.C_TAX_ID.values()
            else:
                S_tax_id=[tax_id]
            data=[]
            for i_tax_id in S_tax_id:
                fn=setting.ppi.get('STRING_PATH', os.path.join(os.path.dirname(__file__), f"STRING/Interaction.{i_tax_id}.csv.gz"))
                if os.path.exists(fn):
                    t=util.read_csv(fn, dtype={'gid_A':str, 'gid_B':str})
                    if "PHYSICAL" in s_db:
                        t=t[t.interaction_type_id==11].copy()
                    t.rename2({'gid_A':'Gene_A', 'gid_B':'Gene_B', 'tax_id_A':'tax_id'})
                    sw.check(f"data loaded from {fn}")
                else:
                    if i_tax_id>0:
                        if "PHYSICAL" in s_db:
                            t=mydb.from_sql("SELECT gid_A Gene_A,gid_B Gene_B,interaction_type_id,score_physical,score_combined,tax_id_A tax_id,support from interaction where tax_id_A=? and interaction_type_id=11", params=[i_tax_id])
                        else:
                            t=mydb.from_sql("SELECT gid_A Gene_A,gid_B Gene_B,interaction_type_id,score_physical,score_combined,tax_id_A tax_id,support from interaction where tax_id_A=?", params=[i_tax_id])
                    else:
                        if "PHYSICAL" in s_db:
                            t=mydb.from_sql("SELECT gid_A Gene_A,gid_B Gene_B,interaction_type_id,score_physical,score_combined,tax_id_A tax_id,support from interaction where interaction_type_id=11")
                        else:
                            t=mydb.from_sql("SELECT gid_A Gene_A,gid_B Gene_B,interaction_type_id,score_physical,score_combined,tax_id_A tax_id,support from interaction")
                    t['Gene_A']=t.Gene_A.astype(str)
                    t['Gene_B']=t.Gene_B.astype(str)
                if sum(t.Gene_A>t.Gene_B):
                    util.info_msg("Genes not order by str, canonicalize required!")
                    t=Network.canonicalize_table(t) # since we change type to str, we need to reorder it
                data.append(t)
            if len(data)==1:
                t=data[0]
            else:
                t=pd.concat(data, ignore_index=True)
            #sw.check("Canonicalized")
            t['TYPE']='Direct'
            sw.check("Start processing each tax_id")
            S_tax_id=t.tax_id.unique()
            for tax_id in S_tax_id:
                #for tax_id,t_v in t.groupby('tax_id'):
                #sw.check("ENTER GROUPBY")
                if tax_id not in Cache.ppi_data['HISTORY']:
                    Cache.ppi_data['HISTORY'][tax_id]={}
                    Cache.ppi_node['HISTORY'][tax_id]={}
                    Cache.ppi_edge['HISTORY'][tax_id]={}
                if "COMBINED" in s_db:
                    tmp=t.loc[t.tax_id==tax_id, ['Gene_A','Gene_B','TYPE','score_combined','support']].copy()
                    #sw.check("COPY")
                    tmp.rename2({'score_combined':'SCORE'})
                    data,t_node=Cache.df2data(tmp, con=mydb)
                    #sw.check("DICT")
                    Cache.ppi_data['HISTORY'][tax_id]["COMBINED_ALL"]=data
                    Cache.ppi_node['HISTORY'][tax_id]["COMBINED_ALL"]=t_node
                    Cache.ppi_edge['HISTORY'][tax_id]["COMBINED_ALL"]=tmp
                    #sw.check("Combined all")
                    tmp=tmp[tmp.SCORE>=Cache.CUTOFF_COMB].copy()
                    #sw.check("FILTER")
                    data,t_node=Cache.df2data(tmp, con=mydb)
                    #sw.check("DICT2")
                    Cache.ppi_data['HISTORY'][tax_id]["COMBINED_CORE"]=data
                    Cache.ppi_node['HISTORY'][tax_id]["COMBINED_CORE"]=t_node
                    Cache.ppi_edge['HISTORY'][tax_id]["COMBINED_CORE"]=tmp
                    #tmp=t_v[t_v.interaction_type_id==11]
                tmp=t.loc[(t.tax_id==tax_id) & (t.interaction_type_id==11)]
                tmp=tmp[['Gene_A','Gene_B','TYPE','score_physical','support']].copy()
                tmp.rename2({'score_physical':'SCORE'})
                #sw.check("Combined core")
                data,t_node=Cache.df2data(tmp, con=mydb)
                Cache.ppi_data['HISTORY'][tax_id]["PHYSICAL_ALL"]=data
                Cache.ppi_node['HISTORY'][tax_id]["PHYSICAL_ALL"]=t_node
                Cache.ppi_edge['HISTORY'][tax_id]["PHYSICAL_ALL"]=tmp
                #sw.check("Physical all")
                tmp=tmp[tmp.SCORE>=Cache.CUTOFF_COMB].copy()
                data,t_node=Cache.df2data(tmp, con=mydb)
                Cache.ppi_data['HISTORY'][tax_id]["PHYSICAL_CORE"]=data
                Cache.ppi_node['HISTORY'][tax_id]["PHYSICAL_CORE"]=t_node
                Cache.ppi_edge['HISTORY'][tax_id]["PHYSICAL_CORE"]=tmp
                #sw.check("Physical core")
                sw.check(f"processed :{tax_id}")
                t=t.loc[t.tax_id!=tax_id]
            return

        S_DB=S_DB or Cache.get_DB(l_use_GPDB)
        if tax_id is None:
            util.error_msg('tax_id must be an int, or 0 means all supported species')
        tax_id=abs(tax_id)
        s_key=Cache.key(l_use_GPDB)
        S_tax_id=[]
        if not l_use_GPDB:
            if tax_id not in (0,9606):
                util.error_msg('Local database only supports human!')
            tax_id=9606
            if tax_id in Cache.ppi_data[s_key]:
                S_DB=[x for x in S_DB if x not in Cache.ppi_data[s_key][tax_id]]
                if len(S_DB)==0: return
            S_tax_id=[tax_id]
            T=[]
            for filename in S_DB:
                print("loading PPI database: "+filename+" ...")
                if os.path.isfile(filename):
                    t=pd.read_csv(filename)
                    t['ds']=filename
                    T.append(t)
                elif os.path.isfile(Cache.DATA_DIR+filename+".csv"):
                    t=pd.read_csv(Cache.DATA_DIR+filename+".csv")
                    t['ds']=filename
                    T.append(t)
                else:
                    util.warn_msg('PPI database ' + filename + ' not found.')
            if len(T)>1:
                t=pd.concat(T, axis=0, ignore_index=True)
            else:
                t=T[0]
            t=t[(t.Gene_A!=t.Gene_B) & (t.Score>=0.5)].copy()
            eg=entrez
            if eg is None:
                eg=ez.EntrezGene(tax_id=tax_id)
            else:
                eg.load_organism(tax_id=tax_id)
            c_seen={}
            t.index=list(range(len(t)))
            t['Gene_A']=t.Gene_A.astype(str)
            t['Gene_B']=t.Gene_B.astype(str)
            S_gene_A=t.Gene_A.tolist()
            S_gene_B=t.Gene_B.tolist()
            for i in range(len(t)):
                gene_A=S_gene_A[i]
                gene_B=S_gene_B[i]
                if gene_A not in c_seen:
                    c_seen[gene_A]=eg.fix_gene_id(gene_A)
                S_gene_A[i]=c_seen[gene_A]
                if S_gene_A[i] is None: continue
                if gene_B not in c_seen:
                    c_seen[gene_B]=eg.fix_gene_id(gene_B)
                S_gene_B[i]=c_seen[gene_B]
            t['Gene_A']=S_gene_A
            t['Gene_B']=S_gene_B
            t=t[~(t.Gene_A.isnull() | t.Gene_B.isnull())].copy()
            t.index=list(range(len(t)))
            t['tax_id']=tax_id
        else:
            mydb=db.DB('METASCAPE')
            if tax_id>0 and tax_id in Cache.ppi_data[s_key]:
                S_DB=[x for x in S_DB if x not in Cache.ppi_data[s_key][tax_id]]
                if len(S_DB)==0: return
            if tax_id>0:
                print("loading PPI database from database for tax_id: %d ..." % tax_id)
                t=mydb.sql_in("SELECT gid_A Gene_A,gid_B Gene_B,0 Score,tax_id_A tax_id,ds from interaction where interaction_category!='genetic' and gid_A!=gid_B and tax_id_A=tax_id_B and tax_id_A=? and ds in (", ")", S_DB, params_before=[tax_id])
                S_tax_id=[tax_id]
            else:
                #ZZZ modify in the future, to obtain the list of all supported tax_id
                t=mydb.from_sql('SELECT DISTINCT tax_id FROM gid2source_id')
                S_tax_id=[x for x in t.tax_id.astype(int).tolist() if x not in Cache.ppi_data[s_key]]
                if len(S_tax_id):
                    s_tax_id=",".join(util.iarray2sarray(S_tax_id))
                    print("loading PPI database for tax_id: %s ..." % s_tax_id)
                    t=mydb.sql_in("SELECT gid_A Gene_A,gid_B Gene_B,0 Score,tax_id_A tax_id,ds from interaction where interaction_category!='genetic' and gid_A!=gid_B and tax_id_A=tax_id_B and ds in (", ")", S_DB)
                    #t=mydb.sql_in("SELECT gid_A Gene_A,gid_B Gene_B,0 Score,tax_id_A tax_id,ds from interaction where interaction_category!='genetic' and gid_A!=gid_B and tax_id_A=tax_id_B and tax_id_A in ("+s_tax_id+") and ds in (", ")", S_DB)
                else:
                    t=pd.DataFrame()
            if len(t):
                t['Gene_A']=t.Gene_A.astype(str)
                t['Gene_B']=t.Gene_B.astype(str)
                if sum(t.Gene_A>t.Gene_B):
                    t=Network.canonicalize_table(t) # since we change type to str, we need to reorder it

        for x in S_tax_id:
            #print ">>>>>>>>>>>>>>>>>>>>>>>", x
            if x not in Cache.ppi_data[s_key]:
                Cache.ppi_data[s_key][x]={}
                Cache.ppi_node[s_key][x]={}
            for y in S_DB:
                Cache.ppi_data[s_key][x][y]={}
                Cache.ppi_node[s_key][x][y]=pd.DataFrame()

        if len(t)==0: return
        for k,t_v in t.groupby(['tax_id','ds']):
            #print ">>>", k, len(t_v)
            #t_v=t_v.copy()
            if k[0] not in S_tax_id: continue
            data={}
            t_node=None
            #t_v=t_v.copy()
            #t_v.index=list(range(len(t_v)))
            #for i in t_v.index:
                #if i%1000==0: print i
            for row in t_v.itertuples():
                gene_A=row.Gene_A #t_v.ix[i,'Gene_A']
                gene_B=row.Gene_B #t_v.ix[i,'Gene_B']
                score=row.Score #t_v.ix[i,'Score']
                if gene_A not in data:
                    data[gene_A]={gene_B:score}
                else:
                    data[gene_A][gene_B]=max(score, data[gene_A].get(gene_B,0))
                if gene_B not in data:
                    data[gene_B]={gene_A:score}
                else:
                    data[gene_B][gene_A]=max(score, data[gene_B].get(gene_A,0))
            Cache.ppi_data[s_key][k[0]][k[1]]=data
            S_gene=list(data.keys())
            if l_use_GPDB:
                t_node=Cache.gene2node(S_gene, con=mydb)
            else:
                t_node=eg.gene_sarray_to_table(S_gene, l_description=False)
            Cache.ppi_node[s_key][k[0]][k[1]]=t_node

# YZHOU: for InWeb_IM, their web GUI uses a threshold for score
#From: Rasmus Borup Hansen [mailto:rbh@intomics.com]
#Sent: Friday, February 03, 2017 4:22 AM
#Subject: Re: Interaction not shown in InBio Map
#
#To make a long story short: We've tried a number of different strategies for choosing a cutoff, and right now the web interface uses 0.156.
#
#Best,
#
#Rasmus

class PPI(Network):

    def __init__(self, tax_id=9606, l_use_GPDB=False, S_DB=None):
        """tax_id is None, defaults to 9606, if 0, means load all species
        Warning: S_DB is set in Cache.load(), so preload Cache if you want to use different database"""
        self.tax_id=tax_id
        data, t_node, t_edge=Cache.get(tax_id=tax_id, l_use_GPDB=l_use_GPDB, S_DB=S_DB)
        print("PPI databases loaded")
        super(PPI, self).__init__(data, T_node=t_node, name='proteome', premade_T_edge=t_edge, skip_copy=True)


if __name__=="__main__":

    #Cache.load(tax_id=9606, S_DB='PHYSICAL_CORE')
    #Cache.load(tax_id=9606, S_DB='COMBINED_CORE')
    sw=util.StopWatch()
    Cache.load(tax_id=9606, S_DB='COMBINED_CORE')
    Cache.info()
    sw.check('Loaded')
    #Cache.load(tax_id=0, l_use_GPDB=True)
    #Cache.load(tax_id=0, S_DB=['BioGrid','GeneGO'], l_use_GPDB=True)
    #Cache.info()
    #exit()
    ppi=PPI(l_use_GPDB=True, tax_id=9606)
    sw.check('Ready')
    exit()
    ppi=PPI(l_use_GPDB=True, tax_id=9606)
    print(list(Cache.ppi_data['GPDB'].keys()))
    #ppi.T_node.to_csv('t1.csv')
    #ppi.T_edge.to_csv('t2.csv')
    print(ppi.data['132884'])
    S_node=['132884','191','537']
    test=ppi.subnetwork(S_node)
    print(test.nof_nodes())
    exit()

    ## example
    S_node=util.read_list('~/RM_Hits.txt')
    test=ppi.subnetwork(S_node)
    test.to_xgmml('RM_.xgmml')
    exit()
    S_node=util.read_list('~/CM_Hits.txt')
    test=ppi.subnetwork(S_node)
    test.to_xgmml('CM_.xgmml')

    exit()
    #print ppi.T_node[:5]
    #print ppi.T_edge[:5]
    test=ppi.subnetwork(S_node)
    #print test
    exit()

    mc=MCODE(net)
    #print mc.C_nodeByScore
    mc.params['hairCut']=True
    c=mc.find_clusters(True, True)
    print(mc.to_MCODE_table(c))
    for i,cp in enumerate(c):
        print(">>> Rank "+str(i)+" <<<")
        cp.to_xgmml('out/test'+str(i))
        S=cp.nodes()
        for node in S:
            nodeInfo=mc.C_nodeInfo[node]
            print("Node=> "+node)
            print(nodeInfo)
