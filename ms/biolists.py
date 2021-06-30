#!/usr/bin/env python
from __future__ import absolute_import
from __future__ import print_function
import pandas as pd
import stats
import numpy as np
import util
import go
import re
import brewer2mpl
import json
import os
import time
from . import msobject
import glob
import tempfile
import cluster
import traceback
import parallel
import copy
import math
import ppi
from six.moves import range
import cytoscape as cyto

#import mahotas

class NamedList(msobject.MSObject):
    '''List where elements can be access either by name or by index'''

    def __init__(self, values=None, names=None, data=None):
        #print ">>>> ", items, ids, '<<<<<<'
        self.values=None
        self.names=None
        self.c_map=None
        self.current=None
        values = values or []
        names = names or []
        if data is not None and isinstance(data, NamedList):
            self.values=data.values[:]
            self.names=data.names[:]
            self.c_map=data.c_map.copy()
        else:
            self.values=values[:]
            if names:
                self.names=names[:]
            else:
                self.names=[
                        getattr(v, 'id', getattr(v, 'name', str(i)))
                        if v is not None else str(i)
                        for i,v in enumerate(values)]
        self.c_map={s:i for i,s in enumerate(self.names)}
        if len(self.c_map)<len(self.names):
            c=util.unique_count(self.names)
            s=", ".join([k for k,v in c.items() if v>1])
            util.error_msg("Duplicate names: %s" % s)

    def __getitem__(self, index):
        if type(index) is int:
            return self.values[index]
        elif index in self.c_map:
            return self.values[self.c_map[index]]
        else:
            util.error_msg('Name: '+index+' is not found!')

    def add(self, name, value, index=-1):
        if index<0: index=len(self.names)+index+1
        self.values[index:index]=[value]
        self.names[index:index]=[name]
        for i in range(index, len(self.names)):
            self.c_map[self.names[i]]=i

    def __contains__(self, name):
        return name in self.names

    def remove(self, name):
        if name not in self.c_map:
            util.error_msg('Name: '+name+' not in the named list!')
        idx=self.c_map[name]
        del self.values[idx]
        del self.names[idx]
        del self.c_map[name]
        for i in range(idx, len(self.names)):
           self.c_map[self.names[i]]=i
        return idx

    def clear(self):
        self.values=[]
        self.names=[]
        self.c_map={}

    def __len__(self):
        return len(self.names)

    def __iter__(self):
        self.current=0
        return self

    # python2
    def next(self):
        if self.current==len(self.names):
            raise StopIteration
        else:
            self.current+=1
            return (self.names[self.current-1], self.values[self.current-1])

    # python3
    def __next__(self):
        if self.current==len(self.names):
            raise StopIteration
        else:
            self.current+=1
            return (self.names[self.current-1], self.values[self.current-1])

    def __str__(self):
        S=[]
        for i,s in enumerate(self.names):
            if i>=5: break # don't print too many
            S.append(s+':'+str(self.values[i]))
        return '\n'.join(S)

class GeneList(msobject.MSObject):

    def __init__(self, name, S_gene=None):
        self.name=name
        if S_gene is None:
            self.data=set()
        else:
            self.data=set(S_gene)
        # make sure Gene IDs are in str type
        for x in self.data:
            if type(x) is not str:
                self.data={str(int(x)) for x in self.data}
            break
        self.data={x for x in self.data if re.match('^\d+$', x) }

    def is_empty(self):
        if self.data is None: return True
        return len(self.data)==0

    def name(self):
        return self.name

    def rename(self, name):
        self.name=name

    def __contains__(self, s_gene):
        return s_gene in self.data

    def genes(self):
        return list(self.data)

    def __len__(self):
        return len(self.data)

    def __str__(self):
        s='GeneList: '+self.name+' #:'+str(len(self.data))
        return s

    def to_str(self):
        return self.name+"\n"+",".join(sorted(list(self.data)))

    def go_analysis(self, go, S_go=None, SRC_GENE=None, min_overlap=3, min_enrichment=0, p_cutoff=0.01, n_CPU=0, S_go_category=None, l_background_by_ontology=False, max_list_size=0):
    # require go object to be passed in, which mostly had databases pre-cached
        if max_list_size>0 and len(self.data)>max_list_size:
            return GOList(None, self.name, genelist=self, go_opts={'SRC_GENE':SRC_GENE, 'min_overlap':min_overlap, 'min_enrichment':min_enrichment, 'p_cutoff':p_cutoff, 'n_CPU':n_CPU, 'S_go_category':S_go_category, 'l_background_by_ontology':l_background_by_ontology})
        return GOList(go.analysis(self.genes(), S_go=S_go, SRC_GENE=SRC_GENE, min_overlap=min_overlap, min_enrichment=min_enrichment, p_cutoff=p_cutoff, n_CPU=n_CPU, S_go_category=S_go_category, l_background_by_ontology=l_background_by_ontology), self.name, genelist=self, go_opts={'SRC_GENE':SRC_GENE, 'min_overlap':min_overlap, 'min_enrichment':min_enrichment, 'p_cutoff':p_cutoff, 'n_CPU':n_CPU, 'S_go_category':S_go_category, 'l_background_by_ontology':l_background_by_ontology})

    @staticmethod
    def mcode_analysis(net, l_MCODE_optimize=False, min_size=3, max_size=500, n_CPU=0, l_connect_in_merge=True):
        n_CPU=1 # not much benefit for parallelize
        print("Find MCODE components ...")
        data=[]
        id_mcode=0
        #mcode_nets=[]
        S_sort_mcode=[]
        c_cluster={}
        c_score={}
        c_type={}
        out_mcode=[]
        t_mcode=None

        s_name=net.name
        netcomps=net.decompose()
        nets=[]
        for i,x in enumerate(netcomps):
            x.name=s_name+"_SUB"+str(i+1)
            if GeneList.filter_net(x, min_size, max_size) is not None:
                nets.append(x)
        for i,x in enumerate(nets):
            mc=ppi.MCODE(x, n_CPU=n_CPU, l_cache=True)
            mc.params['hariCut']=True
            mc.params['fluff']=False # fluff must be False, so that each gene uniquely belongs to one MCODE cluster
            components=mc.find_clusters(True, l_MCODE_optimize)
            if len(components):
                print("Found %d MCODE components" % len(components))
            if len(components)==0: continue
            out_mcode.extend(components)
            t=mc.to_MCODE_table(components)
            t['Cluster']+=id_mcode
            t['Network']=x.name
            data.append(t)
            for j in t.index:
                c_cluster[t.loc[j,'Gene']]=t.loc[j,'Cluster']
                c_score[t.loc[j,'Gene']]=t.loc[j,'Score']
                c_type[t.loc[j,'Gene']]=t.loc[j,'Type']
            #mcode_nets.extend(components)
            for j,cp in enumerate(components):
                S_sort_mcode.append([ j+1+id_mcode, cp.nof_nodes() ])
            id_mcode+=len(components)

        mcode_net=None
        if len(data):
            t_mcode=pd.concat(data, ignore_index=True) if len(data) else None
            # rename cluster ID, so that larger cluster has smaller cluster ID
            S_sort_mcode=sorted(S_sort_mcode, key=lambda x: -x[1])
            #print S_sort_mcode
            c_map_cluster={x[0]:(i+1) for i,x in enumerate(S_sort_mcode)}
            #print c_map_cluster
            c_cluster={k:c_map_cluster.get(v, 0) for k,v in c_cluster.items() }
            #print c_cluster
            t_mcode['Cluster']=t_mcode['Cluster'].apply(lambda x: c_map_cluster.get(x,0))
            if l_connect_in_merge:
                mm=net.subnetwork(list(t_mcode.Gene))
            else:
                mm=out_mcode[0]
                if len(out_mcode)>1:
                    mm.combine_network(out_mcode[1:])
            mm.name=s_name+"_MCODE_ALL"
            mm.add_a_node_attr('MCODE_CLUSTER_ID', c_cluster, s_NULL=0)
            mm.add_a_node_attr('MCODE_SCORE', c_score, s_NULL=0)
            mm.add_a_node_attr('MCODE_TYPE', c_type)
            mcode_net=mm
            #S_gene=set(t_mcode.Gene)
            #t_evi['EvidenceMCODE']=t_evi.Gene.apply(lambda x: 1 if x in S_gene else 0)
        if len(c_cluster):
            net.add_a_node_attr('MCODE_CLUSTER_ID', c_cluster, s_NULL=0)
            net.add_a_node_attr('MCODE_SCORE', c_score, s_NULL=0)
            net.add_a_node_attr('MCODE_TYPE', c_type)
        return(mcode_net, t_mcode)

    @staticmethod
    def filter_net(net, min_size=3, max_size=500):
        n=net.nof_nodes()
        if min_size>0 and n<min_size: return None
        if max_size>0 and n>max_size: return None
        return net

    def ppi_analysis(self, myppi, s_name=None, l_MCODE=True, min_size=3, max_size=500, l_overconnect=False, l_propagation=False, opt_over=None, n_CPU=0, l_MCODE_optimize=False, l_connect_in_merge=True, l_indirect_PPI=False, indirect_size=10, max_list_size=0):
        """l_connect_in_merge, if true, we may add ppi that connect MCODE components, if false, each MCODE component remains isolated."""
        #n_CPU=1 # not much benefit for parallelize
        sw=util.StopWatch("GeneList::ppi_analysis")
        if len(self.data)==0 or len(self.genes())<min_size: return (None, None, None, None, None)
        # add new logic on Jan 29, 2018, YZ
        if max_list_size>0 and len(self.genes())>max_list_size: return (None, None, None, None, None)
        #print self.genes()
        #print type(self.genes()[0])
        #myppi.T_edge.to_csv('t.csv', index=False)
        net=myppi.subnetwork(self.genes())
        if s_name is None:
            s_name=self.name
        #print("NAME>>>>>>>>>>##############", s_name, net.nof_nodes())
        if net.nof_nodes()>0:
            # add new logic on Jan 29, 2018, YZ
            print("Size of total network: %d" % net.nof_nodes())
            # comment out on Mar 11, 2019
            #if net.nof_nodes()>max_size: return (None, None, None, None, None)
            net.name=s_name
            nets=[net]
        if l_indirect_PPI and net.nof_nodes()<=indirect_size: # network too thin, allow indirect interactions
            S_genes=self.genes()
            net_indirect=net.paths_between(S_genes, S_genes, i_hops=2, l_indirect=True)
            #print ">>>>>>>>>>>>>>>>>>>>>>>>", len(S_genes), net_indirect.nof_nodes()
            if net_indirect.nof_nodes()>0:
                net_indirect.name=s_name+"_OneHop"
                nets.append(net_indirect)
            else:
                return (None, None, None, None, None)
        elif net.nof_nodes()==0:
            return (None, None, None, None, None)

        t_evi=pd.DataFrame({'Gene':net.nodes()})
        t_evi['EvidencePPI:%s' % s_name]=1
        #netcomps=net.decompose()
        #sw.check('Overall Network')
        #print "Network Pieces Found: %d" % len(netcomps)
        ##l_onePiece=len(netcomps)==1
        ##if l_onePiece:
        ##    nets=[net]
        ##else:
        ##    nets=[net]+netcomps
        #for i,x in enumerate(netcomps):
        #    #if i==0:
        #    #    net.name=s_name # use the original network name
        #    #else:
        #    x.name=s_name+"_SUB"+str(i+1)
        #    if GeneList.filter_net(x, min_size, max_size) is not None:
        #        nets.append(x)
        out=nets[:]
        #print "Remaining %d sizable network components." % (len(nets)-1)
        t_mcode=mcode_net=None
        #print ">>>>>>>>>", l_MCODE
        if l_MCODE:
            (mcode_net, t_mcode)=GeneList.mcode_analysis(net, l_MCODE_optimize=l_MCODE_optimize, min_size=min_size, max_size=max_size, n_CPU=n_CPU, l_connect_in_merge=l_connect_in_merge)
            if mcode_net is not None:
                out.append(mcode_net)
                S_gene=set(t_mcode.Gene)
                t_evi['EvidenceMCODE:%s' % s_name]=t_evi.Gene.apply(lambda x: 1 if x in S_gene else 0)
        t_over=None
        DEFAULT_OVER={'min_links':2, 'p_cutoff':0.01, 'min_enrichment':3.0}
        opt=DEFAULT_OVER.copy()
        if opt_over is not None:
            opt.update(opt_over)
        if l_overconnect:
            print("Find overconnected nodes ...")
            data=[]
            for i,x in enumerate(nets):
                if re.search(r'_SUB\d+$', x.name):
                    continue # only calculate once for the overall network
                t=myppi.overconnected(myppi, net.nodes(), **opt)
                if t is not None and len(t)>0:
                    t['Network']=x.name
                data.append(t)
            t_over=pd.concat(data, ignore_index=True) if len(data) else None
            sw.check('Overconnect done')
        if l_propagation:
            print("Find new nodes by propagation ...")
            for i,x in enumerate(nets):
                if re.search(r'_SUB\d+$', x.name):
                    continue # only calculate once for the overall network
                c, tmp=myppi.propagation(myppi, prior=x.nodes(), l_Otsu_split=True)
                if filter_net(tmp, min_size, max_size) is not None:
                    x.name=x.name+"_PROP"
                    out.append(x)
            sw.check('Propagation done')
        sw.check('Network DONE')
        for x in out:
            ppi.Network.add_node_degree(x)
        if t_evi is not None:
            t_evi['Gene']=t_evi.Gene.astype(str)
        return (out, t_mcode, t_over, net, t_evi)

    @staticmethod
    def save_network_files(nets, t_mcode, t_over, s_output_dir="", s_output_xgmml=None):
        if s_output_xgmml is None:
            s_output_xgmml=s_output_dir
        if nets is not None and len(nets)>0:
            for x in nets:
                x.to_xgmml(os.path.join(s_output_xgmml, x.name))
        if t_mcode is not None and len(t_mcode):
            util.df2sdf(t_mcode).to_csv(os.path.join(s_output_dir, "MCODE.csv"), index=False)
        if t_over is not None and len(t_over):
            util.df2sdf(t_over).to_csv(os.path.join(s_output_dir, "PPI_OverConnect.csv"), index=False)

class GeneLists(NamedList):

    def __init__(self, genelists=None):
        # genelists, a list of GeneList objects
        genelists=[] if genelists is None else genelists
        if isinstance(genelists, GeneList):
            genelists=[genelists]
        names=[]
        values=[]
        for x in genelists:
            if len(x)==0:
                util.warn_msg('Empty gene list ignored: %s' % x.name)
                continue
            names.append(x.name)
            values.append(x)
        self.t_mem=None
        super(GeneLists, self).__init__(values, names)

    def rename(self, name, new_name):
        i=util.index(name, self.names)
        if i>=0:
            self.values[i].rename(new_name)
            self.names[i]=new_name
            self.c_map[new_name]=self.c_map[name]
            del self.c_map[name]
            return
        util.warn_msg('gene list: '+name+' is not found!')

    def is_empty(self):
        if len(self.values)==0: return True
        l_non_empty=False
        for x in self.values:
            if len(x)>0:
                l_non_empty=True
                break
        return not l_non_empty

    def merge(self):
        S_gene=set()
        for x in self.values:
            S_gene |= x.data
        return GeneList(name='Merged', S_gene=S_gene)

    def membership(self, s_gene):
        return [ 1 if s_gene in x else 0 for x in self.values ]

    def membership_table(self, l_gene_as_index=True):
        if self.t_mem is not None:
            df=self.t_mem.copy()
            if l_gene_as_index:
                df.index=df.Gene
            else:
                df.index=list(range(len(df)))
            return df
        S_gene=self.merge().genes()
        if len(S_gene)==0:
            df=pd.DataFrame([], columns=["_MEMBER_"+x.name for x in self.values])
            df['_PATTERN_']=[]
            df['_WEIGHTED_RANK_']=[]
            df['_RANK_']=[]
            df['Gene']=[]
            util.warn_msg('All gene lists are empty')
            return df
        M=np.transpose(np.array([ np.array([1 if g in x.data else 0 for g in S_gene], dtype=np.int32) for x in self.values ]))
        S_pat=[ "M"+"".join([str(x) for x in M[i,:]]) for i in range(M.shape[0]) ]
        R_w=np.array([1.0/len(x) for x in self.values])
        R_weighted_rank=(M*R_w).sum(axis=1)
        R_rank=M.sum(axis=1)
        df=pd.DataFrame(M, index=S_gene, columns=["_MEMBER_"+x.name for x in self.values])
        df['Gene']=S_gene
        df['_PATTERN_']=S_pat
        df['_WEIGHTED_RANK_']=R_weighted_rank
        df['_RANK_']=R_rank
        df.sort_values(by='_RANK_', ascending=False, inplace=True)
        if not l_gene_as_index:
            df.index=list(range(len(df)))
        self.t_mem=df
        #print self.t_mem
        return df

    def overlap(self, n_total=None):
        """n_total: total number of genes in background"""
        if n_total is None:
            n_total=go.GO.TOTAL_GENE_COUNT
        N=len(self.names)
        M_cnt=np.ones([N,N], dtype=int)
        M_p=np.zeros([N,N])
        for i in range(N):
            S1=self.values[i].data
            n1=len(S1)
            M_cnt[i,i]=n1
            for j in range(i+1, N):
                S2=self.values[j].data
                n2=len(S2)
                n=len(S1 & S2)
                M_cnt[i,j]=M_cnt[j,i]=n
                M_p[i,j]=M_p[j,i]=stats.hyper(n, n_total, n1, n2)
        t_cnt=pd.DataFrame(M_cnt, columns=self.names)
        t_cnt.index=self.names
        t_p=pd.DataFrame(M_p, columns=self.names)
        t_p.index=self.names
        return (t_cnt, t_p)

    def _go_analysis(self, go, S_go=None, SRC_GENE=None, min_overlap=3, min_enrichment=0, p_cutoff=0.01, n_CPU=0, S_go_category=None, l_background_by_ontology=False, max_list_size=0):
        # require go object to be passed in, which mostly had databases pre-cached
        return GOLists([x.go_analysis(go, S_go=S_go, SRC_GENE=SRC_GENE, min_overlap=min_overlap, min_enrichment=min_enrichment, p_cutoff=p_cutoff, n_CPU=n_CPU, S_go_category=S_go_category, l_background_by_ontology=l_background_by_ontology, max_list_size=max_list_size) for x in self.values])

    def go_analysis(self, go, S_go=None, SRC_GENE=None, min_overlap=3, min_enrichment=0, p_cutoff=0.01, n_CPU=0, S_go_category=None, l_background_by_ontology=False, max_list_size=0):
        """New version, better parallelism"""
        c_hitlist={}
        c_go={}
        sw=util.StopWatch("GeneLists::go_analysis")
        for x in self.values:
            if max_list_size>0 and len(x.data)>max_list_size: continue
            if len(x.data)<min_overlap: continue
            c_hitlist[x.name]=x.genes()
            c_go[x.name]=S_go
        out=[]
        if len(c_hitlist):
            t_go=go.analysis(c_hitlist, S_go=c_go, SRC_GENE=SRC_GENE, min_overlap=min_overlap, min_enrichment=min_enrichment, p_cutoff=p_cutoff, n_CPU=n_CPU, S_go_category=S_go_category, l_background_by_ontology=l_background_by_ontology)
            if t_go is not None:
                for k,t_v in t_go.groupby('Name'):
                    t_v=t_v.copy()
                    t_v.drop('Name', axis=1, inplace=True)
                    out.append(GOList(t_v, k, genelist=self.values[self.c_map[k]], go_opts={'SRC_GENE':SRC_GENE, 'min_overlap':min_overlap, 'min_enrichment':min_enrichment, 'p_cutoff':p_cutoff, 'n_CPU':n_CPU, 'S_go_category':S_go_category, 'l_background_by_ontology':l_background_by_ontology}))
        return GOLists(out)

    def __str__(self):
        s="GeneLists: %d" % len(self)
        for i in range(len(self)):
            s+="\n"+str(self.values[i])
        return s

    def ppi_analysis(self, myppi, s_name=None, l_MCODE=True, min_size=3, max_size=500, l_overconnect=False, l_propagation=False, opt_over=None, n_CPU=0, l_MCODE_optimize=False, l_individual=True, l_merge=True, l_indirect_PPI=False, indirect_size=10, max_list_size=0):
        if len(self.values)==1:
            (out, t_mcode, t_over, nets, t_evi)=self.values[0].ppi_analysis(myppi, s_name=None, l_MCODE=l_MCODE, min_size=min_size, max_size=max_size, l_overconnect=l_overconnect, l_propagation=l_propagation, opt_over=opt_over, n_CPU=n_CPU, l_MCODE_optimize=l_MCODE_optimize, l_connect_in_merge=False, l_indirect_PPI=l_indirect_PPI, indirect_size=indirect_size, max_list_size=max_list_size)
            # GeneList returns one network obj, instead of a dictionary
            nets={} if nets is None else {nets.name : nets}
            if out is None:
                out=[]
            return (out, t_mcode, t_over, nets, t_evi)
        import ppi
        nets={}
        t_mcode=t_over=net=None
        S_mcode=[]
        S_over=[]
        out=[]
        t_evi=None
        S_evi=[]

        if l_merge:
            gl=self.merge()
            (out, t_mcode, t_over, net, tmp_evi)=gl.ppi_analysis(myppi, s_name=s_name, l_MCODE=l_MCODE, min_size=min_size, max_size=max_size, l_overconnect=l_overconnect, l_propagation=l_propagation, opt_over=opt_over, n_CPU=n_CPU, l_MCODE_optimize=l_MCODE_optimize, l_connect_in_merge=False, l_indirect_PPI=l_indirect_PPI, indirect_size=indirect_size)
            if net is not None:
                nets[s_name]=net
            #if out is not None:
            #    for x in out:
            #        print ">>>", x.name, x.nof_nodes()
            if out is None: out=[]
            if out is not None and len(out):
                t_mem=self.membership_table()
                S=[x for x in t_mem.header() if x.startswith('_MEMBER_') ]
                for i,x in enumerate(S):
                    c_dict={}
                    # _MEMBER_* columns might have been modified by Circos plot into 0.3333
                    for idx in t_mem.index:
                        c_dict[str(t_mem.loc[idx,'Gene'])]=1 if t_mem.loc[idx, x]>=0.99 else 0
                    for y in out:
                        y.add_a_node_attr('#COUNT_%03d' % (i+1), c_dict)

            if t_mcode is not None and len(t_mcode)>0:
                S_mcode.append(t_mcode)
            if l_overconnect and t_over is not None and len(t_over)>0:
                S_over.append(t_over)

        def f(gl):
            # one gene list only uses ONE cpu
            return gl.ppi_analysis(myppi, s_name=gl.name, l_MCODE=l_MCODE, min_size=min_size, max_size=max_size, l_overconnect=l_overconnect, l_propagation=l_propagation, opt_over=opt_over, n_CPU=1, l_MCODE_optimize=l_MCODE_optimize, l_connect_in_merge=False, l_indirect_PPI=l_indirect_PPI, indirect_size=indirect_size, max_list_size=max_list_size)

        if l_individual:
            #if n_CPU<=1:
            #    S_out=[f(x) for x in self.values]
            #else:
            #    L=[x for x in self.values]
            #    mp=parallel.MP()
            #    mp.start(f, n_CPU=n_CPU)
            #    S_out=mp.map(L)
            S_out=parallel.parmap(f, self.values, n_CPU)

            for i,x in enumerate(self.values):
                (out2, t_mcode2, t_over2, net, tmp_evi)=S_out[i]
                #if out2 is not None:
                #    for x in out2:
                #        print ">>>", x.name, x.nof_nodes()
                if out2 is not None and len(out2):
                    out.extend(out2)
                if t_mcode2 is not None and len(t_mcode2)>0:
                    S_mcode.append(t_mcode2)
                if l_overconnect and t_over2 is not None and len(t_over2)>0:
                    S_over.append(t_over2)
                if net is not None:
                    nets[x.name]=net
                if tmp_evi is not None:
                    tmp_evi.rename2({'EvidencePPI': 'EvidencePPI:%s' % x.name, 'EvidenceMCODE': 'EvidenceMCODE:%s' % x.name})
                    if t_evi is None:
                        t_evi=tmp_evi
                    else:
                        t_evi=t_evi.merge(tmp_evi, left_on='Gene', right_on='Gene', how='outer')
                        #print t_evi.header()
            #nets[s_name]=net

        if len(S_mcode):
            t_mcode=pd.concat(S_mcode, ignore_index=True)
        if len(S_over):
            t_over=pd.concat(S_over, ignore_index=True)
        if t_evi is not None:
            t_evi.fillna(0, inplace=True)
            t_evi=t_evi.astype(int)
            t_evi['Gene']=t_evi.Gene.astype(str)
        return (out, t_mcode, t_over, nets, t_evi)

class GOList(msobject.MSObject):

    def __init__(self, t_go=None, name='', genelist=None, go_opts=None):
        """genelist if provided must share the same name as the GOList. This will memorize the original gene list where the GOList was derived from"""
        self.S_total=set()
        if t_go is None or len(t_go)==0:
            self.data=None
            self.S_total=set()
            self.c_go=set()
        else:
            self.data=t_go
            if 'GeneID' in t_go.header(): # parent GO list may not care about GeneID
                self.S_total=set("|".join(t_go.GeneID).split("|"))
            self.c_go=set(t_go.GO)
        self.name=name
        self.genelist=genelist
        self.go_opts={} if go_opts is None else go_opts.copy()
        self.net=None
        self.go_cluster=None
        self.t_mem=None

    def clone(self):
        return GOList(t_go=self.data, name=self.name, genelist=self.genelist, go_opts=self.go_opts)

    def is_empty(self):
        if self.data is None: return True
        return len(self.data)==0

    def GOs(self, S_col=None):
        if self.data is None:
            return None
        if not S_col:
            return list(self.data.GO)
        else:
            if 'GO' not in S_col:
                S_col=['GO']+S_col
            return self.data.reindex(columns=S_col)

    def __contains__(self, s_go):
        return s_go in self.c_go

    def extract_parent_go(self):
        t_go=go.GO.extract_parent_go(self.data)
        return GOList(t_go=t_go, name=self.name, genelist=self.genelist, go_opts=self.go_opts)

    def combine(self, golist2):
        """Combine two golist into one, they should have the same name"""
        if golist2.name!=self.name:
            util.error_msg('Cannot combine lists of different names: {} {}!'.format(self.name, golist2.name))
        if len(golist2)==0:
            return
        if self.__len__()==0:
            self.data=golist2.data.copy()
            self.S_total=golist2.S_total.copy()
            self.c_go=golist2.c_go.copy()
            self.go_opts=golist2.go_opts.copy()
        else:
            self.data=pd.concat([self.data, golist2.data], ignore_index=True)
            self.data.sort_values('LogP', inplace=True)
            self.c_go|=golist2.c_go
            self.S_total|=golist2.S_total
            self.go_opts.update(golist2.go_opts)

    def is_clustered(self):
        # has been analyzed by go.GO_Cluster
        if self.data is None: return False
        return 'FirstInGroupByLogP' in self.data.header()

    def cluster(self, similarity=0.3, n_CPU=0, max_terms=0, l_go_selective=False):
        """If max_terms are specified, we only cluster to top max_terms, sorted by LogP,
        this is to make sure it returns within reasonable time"""
        self.net=None
        if self.data is None:
            util.warn_msg('Cannot cluster None!')
            return None
        if max_terms>0 and len(self.data)>max_terms:
            if l_go_selective:
                t_go=self.data.sort_values(['GiniIndex','LogP'], ascending=[False, True])
            else:
                t_go=self.data.sort_values('LogP')
            if len(t_go)>max_terms:
                t_go=t_go[:max_terms].copy()
        else:
            t_go=self.data
        #self.go_cluster=go.GO_Cluster(list(self.S_total), c_go=self.data, n_CPU=n_CPU)
        # if there are too many terms, only cluster the top max_terms
        self.go_cluster=go.GO_Cluster(list(self.S_total), c_go=t_go, n_CPU=n_CPU)
        self.t_mem=self.go_cluster.data
        out=self.go_cluster.cluster(similarity=similarity, l_go_selective=l_go_selective)
        # we send a copy of t_go to cluster, so we need to add clustering columns back
        S_old=self.data.header()
        S_new=[x for x in out.header() if x not in S_old]
        self.data=self.data.merge(out[['GO']+S_new], left_on='GO', right_on='GO', how='left')
        self.data['GROUP_ID'].fillna(0, inplace=True) # incase we truncate and only cluster a subset of terms
        return out

    def network(self, max_clusters=20, max_members=10, max_nodes=250, l_go_selective=False):
        if not self.is_clustered():
            util.error_msg('Please cluster first!')
        return self.go_cluster.network(max_clusters=max_clusters, max_members=max_members, max_nodes=max_nodes, l_go_selective=l_go_selective)

    def membership(self, n_CPU=0, l_cluster=True, max_cluster=20, S_exclude=None):
        """Membership table of gene by term
        The format for column name is: GRP3|19_GO:0006826|-3.573|iron ion transport
        cluster id 'GRP3|' is omitted, if terms has not been clustered
        _SCORE_GRP3_ column is the frequency of the gene appears in terms under the given cluster GRP3
        _SCORE_ column name is used, if terms has not been clustered.
        if l_cluster is True, genes ordered based on hierachical clustering
            if terms are clustered, _SCORE_GRP*_ matrix is used
            otherwise, term membership is used
        S_exclude: remove terms fall into this list/set
        """
        sw=util.StopWatch("GOList::membership")
        if self.data is None:
            util.warn_msg('No gene in GO list!')
            return None
        l_keep_all=S_exclude is None or len(S_exclude)==0
        if l_keep_all:
            S_exclude=set()
        else:
            S_exclude=set([re.sub(r'^\d+_', '', x) for x in S_exclude])
        data=self.data
        if not l_keep_all:
            data=self.data[self.data.GO.apply(lambda x: x not in S_exclude)]
        if self.t_mem is None:
            self.go_cluster=go.GO_Cluster(list(self.S_total), c_go=data, n_CPU=n_CPU)
            self.t_mem=self.go_cluster.data
            # t_mem is a binary membership matrix Gene X GO
        l_clustered=self.is_clustered()
        t=self.t_mem.copy()
        if not l_keep_all:
            t.drop(S_exclude, axis=1, inplace=True)
        c_grp={}
        if not l_clustered:
            t['_SCORE_']=self.t_mem.mean(axis=1)
        else:
            # parallelize this will be slower, 0.1 sec versus 0.9 sec
            data=data[data['GROUP_ID']>0] # in case there are too many terms, some are not clusters (as we only cluster the top 1000 terms)
            for k,t_v in data.groupby('GROUP_ID'):
                if max_cluster>0 and k>max_cluster: continue
                # SCORE is the percentage terms a gene appears within a group
                t['_SCORE_GRP%d_' % k]=t[list(t_v.GO)].mean(axis=1)
                c_grp[k]=list(t_v.GO)
        sw.check('Score calculation')

        t.reset_index(inplace=True)
        t.rename2({'index':'Gene'})
        sw.check('Get membership matrix')

        if l_cluster:
            import fastcluster
            import matplotlib as mpl
            mpl.rcParams['pdf.fonttype'] = 42
            mpl.use('Agg')
            #import scipy.cluster.hierarchy as clst

            def f(X):
                t=X[0]
                S_col=X[1]
                s_filter=X[2]
                k=X[3]
                if s_filter is not None:
                    t=t[t[s_filter]>0]
                else:
                    # filter out GO not in any cluster (we only keep top 20)
                    t=t[t.sum(axis=1)>0]
                M=t[S_col].values
                #import ms.msobject as mo
                #mo.MSObject.dump_object(M, 'M')
                #import sys
                #sys.setrecursionlimit(10000)
                if len(S_col)>1:
                    Zc=fastcluster.linkage(M.T, method='average', metric='euclidean', preserve_input=True)
                    den_r=cluster.FastCluster.linkage2order(Zc, M.T) #clst.dendrogram(Zc)
                    #print s_filter, k, len(S_col), S_col, "<<<<<<<<<<>>>>>>>>>>>>>", den_r['leaves']
                    S_col=[ S_col[x] for x in den_r]
                if len(t)>1:
                    Zr=fastcluster.linkage(M, method='average', metric='euclidean', preserve_input=True)
                    #den_r=clst.dendrogram(Zr)
                    den_r=cluster.FastCluster.linkage2order(Zr, M)
                    S_row=[ t.index[x] for x in den_r]
                else:
                    S_row=t.index
                return (k, S_row, S_col)

            L=[]
            if not l_clustered:
                L.append([self.t_mem, self.t_mem.header(), None, None])
            else:
                S_col=[x for x in t.header() if x.startswith('_SCORE_GRP') ]
                L.append([t, S_col, None, None])
                for k,v in c_grp.items():
                    if max_cluster>0 and k>max_cluster: continue
                    L.append([t, v, '_SCORE_GRP%d_' % k, k])

            #if n_CPU<=1:
            #    out=[f(x) for x in L]
            #else:
            #    mp=parallel.MP()
            #    mp.start(f, n_CPU=n_CPU)
            #    out=mp.map(L)

            # for each cluster, we cluster the membership matrix gene-wise (row) and term-wise (col)
            out=parallel.parmap(f, L, n_CPU)
            S_term=[]
            S_sort=[]
            for k,S_row,S_col in out:
                S_term+=S_col
                if k is None:
                    s_sort="_SORT_"
                else:
                    s_sort="_SORT_GRP%d_" % k
                t[s_sort]=0
                t.loc[S_row, s_sort]=list(range(1, len(S_row)+1))
                S_sort.append(s_sort)
            S_header=[ x for x in t.header() if x not in set(S_term+S_sort) ]+S_term+S_sort
            t=t[S_header].copy()
            sw.check('Cluster genes and columns for each cluster.')

        c_rename={}
        S_del=[]
        for i,r in data.iterrows():
            s_go=r['GO']
            logP=r['LogP']
            s_cat=r['CategoryID']
            s_des=r['Description']
            if l_clustered and max_cluster>0 and r['GROUP_ID']>max_cluster:
                S_del.append(s_go)
            else:
                s_new="GRP%d|" % r['GROUP_ID'] if l_clustered else ""
                s_new+=str(s_cat)+"_"+s_go+"|"+("%.3f" % logP)+"|"+s_des
                c_rename[s_go]=s_new
        if len(S_del):
            t.drop(S_del, axis=1, inplace=True)
        t.rename2(c_rename)
        return t

    def gene_with_evidence(self, min_size=3, max_size=100, logp=-2.0, S_GO=None):
        S_gene=set()
        if self.data is None:
            util.warn_msg('No gene in GO list!')
            return S_gene
        tmp=self.data[(self.data['#GeneInGO']<=max_size) & (self.data['#GeneInGO']>=min_size) & (self.data['LogP']<=logp)]
        if S_GO is not None:
            S_GO=set(S_GO)
            tmp=tmp[tmp.GO.apply(lambda x: x in S_GO)]
        s="|".join(tmp.GeneID.tolist())
        s=re.sub(r'\s+', '', s)
        S_gene=set([x for x in s.split("|") if x!=''])
        return S_gene

    def heatmap(self, s_out, t_membership, max_clusters=20, max_members=10, max_nodes=50, S_label=None):
        """t_membership is the membership table generated by golists"""
        if not self.is_clustered():
            util.error_msg('golist_clsutered is not clustered yet!')
        s_out, s_ext=os.path.splitext(s_out)
        S_col=[x for x in t_membership.header() if x.startswith('_MEMBER_')]
        if len(S_col)<=1:
            util.warn_msg('Need at least two gene lists to plot heatmap, only %d found!' % len(S_col))
            # delete s_out files
            S_suffix=[".png",".cdt",".gtr",".atr",".kdt",".jtv",".pdf",".input"]
            S_file = glob.glob(s_out+".*")
            for s_file in S_file:
                s_name, s_ext = os.path.splitext(s_file)
                if s_ext in S_suffix:
                    os.remove(s_file)
            return
        c_go_subset=go.GO_Cluster.sample_rows(self.data, max_clusters=max_clusters, max_members=max_members, max_nodes=max_nodes)
        S_go=list(c_go_subset)
        c_map={re.sub(r'^\d+_', '', x):x for x in t_membership.index}
        S_idx=[c_map[x] for x in S_go]
        data=t_membership.loc[S_idx, S_col].values
        if S_label is None:
            S_label=[x.replace('_MEMBER_', '') for x in S_col]
        S_des=[ t_membership.loc[i,'GO']+':'+t_membership.loc[i,'Description'] for i in S_idx ]
        import fastcluster
        Zr=fastcluster.linkage(data, method='average', metric='euclidean', preserve_input=True)
        Zc=fastcluster.linkage(data.T, method='average', metric='euclidean', preserve_input=True)
        import cluster
        fc=cluster.FastCluster(data, S_col=S_label, S_row=S_go, S_description=S_des, Zr=Zr, Zc=Zc)
        fc.plot(s_out+".png", colormap=None)
        fc.save(s_out, l_norm_row=False)

    def summary_plot(self, s_out, t_membership=None, max_clusters=20, S_label=None, l_go_id=True, cluster_column=True):
        """t_membership is the membership table generated by golists containing p-value columns _LogP_..."""
        #if not self.is_clustered():
        #    util.warn_msg('golist_clustered is not clustered yet!')
        s_out, s_ext=os.path.splitext(s_out)
        S_suffix=[".png",".cdt",".gtr",".atr",".kdt",".jtv",".pdf",".input"]
        S_file = glob.glob(s_out+".*")
        for s_file in S_file:
            s_name, s_ext = os.path.splitext(s_file)
            if s_ext in S_suffix:
                os.remove(s_file)
        if self.data is None or len(self.data)==0:
            return
        l_bargraph=True
        S_col=[x for x in t_membership.header() if x.startswith('_LogP_')]
        if len(S_col)==0: return
        if len(S_col)>1:
            l_bargraph=False
        if self.is_clustered():
            c_go_subset=go.GO_Cluster.representative_rows(self.data, max_clusters=max_clusters)
        else:
            if max_clusters>0:
                c_go_subset=set(self.data['GO'][:max_clusters])
            else:
                c_go_subset=set(self.data['GO'])
        #c_go_subset=go.GO_Cluster.sample_rows(self.data, max_clusters=max_clusters, max_members=max_members, max_nodes=max_nodes)
        S_go=list(c_go_subset)
        c_map={re.sub(r'^\d+_', '', x):x for x in t_membership.index}
        S_idx=[c_map[x] for x in S_go]

        data_out=t_membership.loc[S_idx, ['GO','Description']+S_col].copy()
        data=t_membership.loc[S_idx, S_col].values
        if S_label is None:
            S_label=[x.replace('_LogP_', '') for x in S_col]
        if l_go_id:
            S_des=[ "%s: %s" % (t_membership.loc[i,'GO'], t_membership.loc[i,'Description']) for i in S_idx ]
        else:
            S_des=list(t_membership.loc[S_idx,'Description'])
        if l_bargraph:
            s_col=S_col[0]
            R_logP=t_membership.loc[S_idx, s_col].values
            t=pd.DataFrame(data={'LogP':R_logP, 'Description':S_des})
            t.sort_values('LogP', ascending=False, inplace=True)
            n=len(t)
            t.index=list(range(n))
            n=len(t)
            P=[-2, -3, -4, -6, -10, -20]
            C=brewer2mpl.get_map("YlOrBr", 'sequential', len(P), reverse=False).hex_colors[:]
            t['COLOR']='#FFFFFF'
            for i in t.index:
                p=t.loc[i,'LogP']
                for j,x in enumerate(P):
                    if p>x:
                        clr=C[max(j-1,0)]
                        break
                if p<=-20: clr=C[-1]
                t.loc[i,'COLOR']=clr

            import matplotlib as mpl
            mpl.rcParams['pdf.fonttype'] = 42
            mpl.use('Agg')
            import matplotlib.pyplot as plt
            # we will compute the height of the figure, so that the width of the bar are always about the same
            # it is not pretty when the bar is too fat or too thin, when there are too few or too many bars
            fig_width=6
            bar_sz=0.2
            margin_x=0.05
            fig_height=(bar_sz*n)+(margin_x*fig_width*2)
            margin_y=margin_x*fig_width/fig_height
            plt.figure(figsize=(fig_width, fig_height))
            ###
            plt.clf()
            plt.barh(t.index, -t.LogP.values, align='center', color=t.COLOR, edgecolor='#636363')
            plt.xlabel('-log10(P)') # TeX is not installed r'$-log_{10}P$')
            ax=plt.gca()
            ax.yaxis.tick_right()
            plt.yticks(t.index, t.Description)
            plt.xlim(left=0)
            plt.ylim(bottom=-1, top=n)
            xmin,xmax=plt.xlim()
            for x in [2,4,6,10,20]:
                plt.axvline(x=x, zorder=0, color='#bdbdbd')
            plt.xlim(right=xmax)
            ax.tick_params(axis='y', which='both',length=0)
            plt.margins(x=margin_x, y=margin_y)
            plt.savefig(s_out+".png", bbox_inches='tight')
            plt.savefig(s_out+".pdf", bbox_inches='tight')
            plt.close()
        else:
            def scale_color(x):
                return min(-x/20.0,1.0)

            def scale_distance(x):
                if x>-2: return 0.0
                if x>-3: return 2.0
                if x>-4: return 2.5
                if x>-6: return 3.0
                if x>-10: return 3.5
                if x>-20: return 4.0
                return 4.5

            import pydendroheatmap as pyhm
            cm=pyhm.DendroHeatMap.color_by_pvalue()
            vecscale=np.vectorize(scale_distance)
            data_dist=vecscale(data)

            import fastcluster
            Zr=fastcluster.linkage(data_dist, method='average', metric='euclidean', preserve_input=True)
            if cluster_column:
                Zc=fastcluster.linkage(data_dist.T, method='average', metric='euclidean', preserve_input=True)
            else:
                Zc=None

            vecscale=np.vectorize(scale_color)
            data_color=vecscale(data)

            import cluster
            den_r=cluster.FastCluster.linkage2order(Zr)
            data_out.index=range(len(data_out))
            data_out=data_out.loc[den_r[::-1]]
            data_out.to_csv(s_out+".csv", index=False)
            fc=cluster.FastCluster(data_color, S_col=S_label, S_row=S_go, S_description=S_des, Zr=Zr, Zc=Zc)
            fc.plot(s_out+".png", colormap=cm, l_normalize_for_color=False, l_legend_pvalue=True, l_pdf=True)
            fc.save(s_out, l_norm_row=False)

    def __len__(self):
        if self.data is None: return 0
        return len(self.data)

    def __str__(self):
        return 'GOList: '+self.name+' #:'+str(self.__len__())

    def save(self, s_out):
        if self.data is not None and len(self.data)>0:
            util.df2sdf(self.data).to_csv(s_out, index=False)

class GOLists(NamedList):

    def __init__(self, golists=None):
        # golists, a list of GOList objects
        self.has_genelist=True
        golists=[] if golists is None else golists
        names=[]
        values=[]
        self.go_opts={}
        for x in golists:
            if len(x)==0: continue
            names.append(x.name)
            values.append(x)
            if x.genelist is None: self.has_genelist=False
            self.go_opts.update(x.go_opts)
        super(GOLists, self).__init__(values, names)

    def merge(self, go=None, max_terms=0):
        sw=util.StopWatch("GOLists::merge")
        if len(self.values)==0:
            return GOList(t_go=None, name='Merged', go_opts=self.go_opts)

        #if len(self.values)==1:
        #    return self.values[0].clone()
        t_mem=self.membership_table(l_go_as_index=False, S_matrix=['LogP'])
        data=[x.data for x in self.values if x.data is not None]
        if len(data)==0:
            return GOList(t_go=None, name='Merged', go_opts=self.go_opts)
        t_go=pd.concat(data, ignore_index=False)
        t_go['_ORIGINAL_']=True

        if max_terms>0 and len(t_mem)>max_terms:
            # only the top max_terms
            S_logP=[x for x in t_mem.header() if x.startswith('_LogP_')]
            t_mem['_MIN_LOGP_']=t_mem[S_logP].min(axis=1)
            t_mem.sort_values('_MIN_LOGP_', inplace=True)
            t_mem=t_mem[:max_terms].copy()
            S_go=set(t_mem.GO)
            t_go=t_go[t_go.GO.apply(lambda x: x in S_go)]
            t_mem.drop('_MIN_LOGP_', axis=1, inplace=True)
        #S_matrix=['LogP']
        #sw.check("t_go: %d, t_mem: %d" % (len(t_go), len(t_mem)))
        if self.has_genelist:
            # for GO found multiple times, we recalculate p-value by merging those gene lists
            S_mem=[x for x in t_mem.header() if x.startswith('_MEMBER_')]
            if go is None:
                util.error_msg('method argument "go" must be provided!')
            c_hitlists={} # mem_pattern: values are corresponding merged hit lists
            c_hit_go={}  # mem_pattern: values are a list of GOs need to be recalculated
            for i in t_mem.index:
                s_go=t_mem.loc[i,'GO']
                if s_go in go.GO_CATEGORY:
                    s_go=str(go.GO_CATEGORY[s_go])+"_"+s_go # add prefix
                if t_mem.loc[i, '_RANK_']<=1: continue # not found by multiple lists
                s_mem=t_mem.loc[i,'_PATTERN_']
                if s_mem not in c_hit_go:
                    c_hit_go[s_mem]=[]
                c_hit_go[s_mem].append(s_go)
                if s_mem not in c_hitlists:
                    S_hit=set()
                    for j,x in enumerate(s_mem):
                        if x=='1':
                            S_hit|=self.values[j-1].genelist.data
                    c_hitlists[s_mem]=list(S_hit)
                    #print ">>>", s_mem, len(S_hit)
            l_batch_mode=False
            sw.check('c_hitlists: %d' % len(c_hitlists))
            if l_batch_mode: # 4.0 second for example
                data=[]
                for s_mem,S_go in c_hit_go.items():
                    tmp=go.analysis(c_hitlists[s_mem], S_go=S_go, **self.go_opts)
                    if tmp is not None:
                        data.append(tmp)
                data.append(t_go)
                t_go=pd.concat(data, ignore_index=True)
            else: # this mode seems faster, 1.4 second
                self.go_opts["p_cutoff"]=0.5
                tmp=go.analysis(c_hitlists, S_go=c_hit_go, **self.go_opts)
                sw.check('Finish P-value Update')
                if tmp is not None:
                    tmp.drop('Name', axis=1, inplace=True)
                    t_go=pd.concat([t_go, tmp], ignore_index=True)
            t_go['_ORIGINAL_']=t_go['_ORIGINAL_'].apply(lambda x: False if pd.isnull(x) else True)
        t_go.sort_values(['CategoryID','GO','_ORIGINAL_','LogP'],inplace=True) # need to get unique GO entries
        #t_go.to_csv('output/merge.csv', index=False)
        iB=iE=0
        I=[]
        n=len(t_go)
        t_go.index=list(range(n))
        for i in range(1, n+1):
            if i>=n or (t_go.loc[i,'CategoryID']!=t_go.loc[i-1,'CategoryID']) or (t_go.loc[i, 'GO']!=t_go.loc[i-1, 'GO']):
                iE=i-1
                if iB==iE:
                    pass # GO only occur once
                else:
                    if not t_go.loc[iB, '_ORIGINAL_'] and t_go.loc[iB, 'LogP']<=t_go.loc[iB+1, 'LogP']:
                    # First entry is merged enty (not necessary, as the merged entry might have a p-value so bad that it does not return
                        pass
                    else:
                        if not t_go.loc[iB, '_ORIGINAL_']:
                            iB=iB+1 # do not use the merged entry
                        c_Hits={}
                        S_hits=("|".join(t_go.loc[iB:iE,'Hits'])).split('|')
                        S_id=("|".join(t_go.loc[iB:iE,'GeneID'])).split('|')
                        c_Hits={x:S_id[j] for j,x in enumerate(S_hits)}
                        S_hits=list(c_Hits.keys())
                        S_hits.sort()
                        GeneID=[c_Hits[x] for x in S_hits]
                        s_hits="|".join(S_hits)
                        s_id="|".join(GeneID)
                        t_go.loc[iB,'LogP']=t_go.loc[iB:iE, 'LogP'].min()
                        t_go.loc[iB,'Enrichment']=t_go.loc[iB:iE, 'Enrichment'].max()
                        t_go.loc[iB,'Z-score']=t_go.loc[iB:iE, 'Z-score'].max()
                        t_go.loc[iB,'Hits']=s_hits
                        t_go.loc[iB,'GeneID']=s_id
                        t_go.loc[iB,'#GeneInGOAndHitList']=len(S_hits)
                        t_go.loc[iB,'#GeneInHitList']=t_go.loc[iB,'#GeneInHitList'].max()
                        t_go.loc[iB,'#TotalGeneInLibrary']=t_go.loc[iB,'#TotalGeneInLibrary'].max()
                        t_go.loc[iB,'#GeneInGO']=t_go.loc[iB,'#GeneInGO'].max()
                I.append(iB)
                iB=i
        t_go=t_go.loc[I]
        t_mem['GO']=t_mem['GO'].apply(lambda x: re.sub(r'^\d*_','',x))
        t_mem.drop('Description', inplace=True, axis=1)
        t_go=t_mem.merge(t_go, left_on='GO', right_on='GO')
        t_go.sort_values(['LogP','Enrichment','#GeneInGOAndHitList'], ascending=[True,False,False], inplace=True)
        return GOList(t_go=t_go, name='Merged', go_opts=self.go_opts)

    def extract_parent_go(self):
        return GOLists([ x.extract_parent_go() for x in self.values ])

    def membership_table(self, l_go_as_index=True, S_matrix=None):
        """additional S_matrix type can be "LogP" and/or "Enrichment", e.g., S_matrix=["LogP","Enrichment"]
        """
        if S_matrix is None:
            S_matrix=[]
        if type(S_matrix) is str:
            S_matrix=[S_matrix]
        S_matrix=[x for x in S_matrix if x in ('LogP','Enrichment')]
        S_go=[]
        c_map={}
        c_des={}
        data=[]
        for x in self.values:
            if len(x)==0: continue
            S_go.extend(list(x.data.GO))
            for i in x.data.index:
                s_go=x.data.loc[i,'GO']
                c_map[s_go]=str(x.data.loc[i,'CategoryID'])+'_'+s_go
                c_des[s_go]=x.data.loc[i,'Description']
            tmp=x.data[['GO','LogP','Enrichment']].copy()
            tmp['NAME']=x.name
            tmp['ONE']=1
            tmp['GO']=tmp['GO'].apply(lambda x:c_map[x])
            data.append(tmp)
        if len(data)==0:
            t=pd.DataFrame(columns=['NAME','ONE','GO','LogP','Enrichment'])
        else:
            t=pd.concat(data, ignore_index=True)

        if len(S_go)==0:
            df=pd.DataFrame(data=[], columns=["_MEMBER_"+x.name for x in self.values])
            df['_PATTERN_']=[]
            df['_RANK_']=[]
            df['GO']=[]
            df['Description']=[]
            return df
        S_go=util.unique(S_go)
        #t_mem=t.pivot_table('ONE', index='GO', columns='NAME')
        #t_mem.fillna(0, inplace=True)
        #M=t_mem[self.names].values
        M=np.transpose(np.array([ np.array([1 if g in x.c_go else 0 for g in S_go], dtype=np.int32) for x in self.values ]))
        S_pat=[ "M"+"".join([str(int(x)) for x in M[i,:]]) for i in range(M.shape[0]) ]
        R_rank=M.sum(axis=1)
        df=pd.DataFrame(M, index=[c_map[x] for x in S_go], columns=["_MEMBER_"+x.name for x in self.values])
        for x in S_matrix:
            t_mem=t.pivot_table(x, index='GO', columns='NAME')
            t_mem.fillna(0, inplace=True)
            c_rename={ y:'_'+x+'_'+y for y in self.names }
            t_mem.rename2(c_rename)
            df=df.join(t_mem, how='left')

        df['GO']=S_go
        df['Description']=df['GO'].apply(lambda x: c_des.get(x,x))
        df['_PATTERN_']=S_pat
        df['_RANK_']=R_rank
        df.sort_values('_RANK_', ascending=False, inplace=True)
        if not l_go_as_index:
            df.index=list(range(len(df)))
        return df

    def gene_evidence_table(self, min_size=3, max_size=100, logp=-2, S_GO=None):
        """If S_GO is provided, we only consider s_go in S_GO"""
        t_evi=None
        out=[]
        for x in self.values:
            S_gene=x.gene_with_evidence(min_size=min_size, max_size=max_size, logp=logp, S_GO=S_GO)
            if len(S_gene)==0: continue
            tmp=pd.DataFrame({'Gene':list(S_gene)})
            tmp['Value']=1
            tmp['Name']='EvidenceGO:%s' % x.name
            out.append(tmp)
        if len(out):
            tmp=pd.concat(out, ignore_index=True)
            t_evi=tmp.pivot_table('Value', index='Gene', columns='Name')
            t_evi.fillna(0, inplace=True)
            t_evi=t_evi.astype(int)
            t_evi.reset_index(inplace=True)
            t_evi['Gene']=t_evi.Gene.astype(str)
        return t_evi

    def save(self, s_out):
        data=[]
        for x in self.values:
            if x.data is not None and len(x.data)>0:
                tmp=x.data.copy()
                tmp['GeneList']=x.name
                data.append(tmp)
        if len(data)>0:
            t=pd.concat(data, ignore_index=True)
            util.df2sdf(t).to_csv(s_out, index=False)

    def combine(self, golists2):
        """Combine two golists into one, if a go list name appear in both golists, concatenate their data table
        This is to combine the results of go analysis with non-L1000 with L1000"""
        for x in golists2.values:
            if x.data is None: continue
            if x.name not in self.names:
                self.add(x.name, x)
            else:
                idx=self.c_map[x.name]
                gl1=self.values[idx]
                gl1.combine(x)

def get_color_label(t_go_mem, S_name):
    S_header=[x[8:] for x in t_go_mem.header() if x.startswith('_MEMBER_') ]
    c_col={k:i for i,k in enumerate(S_name)}
    S_idx=[ c_col.get(x, -1) for x in S_header]
    return S_idx

def net_annotation(net, genelists=None):
    """genelists can be None: for single-list, or GeneLists object, or a membership table"""
    t=net.T_node
    if len(t)==0: return
    t.index=list(range(len(t)))
    if genelists is not None:
        if type(genelists) is GeneLists:
            t_mem=genelists.membership_table(l_gene_as_index=True)
        else:
            t_mem=genelists
        if 'Gene' in t_mem.header():
            if t_mem.col_type('Gene')!='s':
                t_mem['Gene']=t_mem['Gene'].apply(lambda x: util.r2i2s(x))
            t_mem=t_mem.set_index('Gene')
        #S_header=[x for x in t_mem.header() if x.startswith('_MEMBER_')]
        S_header=[x for x in t.header() if x.startswith('_MEMBER_')]
        # T_node may not have the same # of _MEMBER_* cols as t_mem, b/c they are sorted and
        # not all gene list have enriched terms, or enriched terms survived the picking
        if '_RANK_' not in t.header():
            t['_RANK_']=1
        if '_PATTERN_' not in t.header():
            t['_PATTERN_']="M1"
        s_pat=t['_PATTERN_'].values[0]
        n=len(s_pat)-1
        for i in range(n):
            t['#COUNT_%03d' % (i+1)]=0

        for i in t.index:
            # For GPEC results, RANK=1, but is complied from hits coming from multiple gene lists
            #if t.ix[i, '_RANK_']==1:
            #    t.ix[i, '#COUNT_%03d' % (t.ix[i, '_PATTERN_'].find("1"))]=t.ix[i, '#GeneInGOAndHitList']
            #else:
                S_gene=t.loc[i,'GeneID'].split('|')
                for j,s in enumerate(t.loc[i, '_PATTERN_']):
                    if s=="1":
                        t.loc[i, '#COUNT_%03d' % j]=t_mem.loc[S_gene, S_header[j-1]].sum()
        # just in case S_gene not found in t_mem, then it introduces 'NaN', which crash Cytoscape
        for i in range(n):
            t['#COUNT_%03d' % (i+1)].fillna(0, inplace=True)

    # pick one per cluster to show labels
    n=len(t)
    iB=iE=0
    t['CLUSTER_LABEL']=' '
    for i in t.index:
        if t.loc[i, 'FirstInGroupByLogP']>0:
            t.loc[i, 'CLUSTER_LABEL']=t.loc[i, 'Description']
#    # order for picking labels
#    c_category={
#        31:1, # GeneGo process
#        27:2, # GeneGo Pathway
#        19:3, # GO Biological Process
#        24:4, # KEGG BP
#        11:5, # Canonical pathway
#        23:6, # Hallmark
#        15:7, # BioCarta
#        20:8, # GO CC
#        25:9, # KEGG CC
#        21:10, # GO MF
#        26:11, # KEGG MF
#        6:12, # Reactome
#        28:13, # Disease assoc
#        13:14, # Chemical perturbation
#        3:15, # Immuno
#        4:16, # Oncogenic
#        29:17 # Drug target
#    }
#    for i in range(1, n+1):
#        if i>=n or t.ix[i,'GROUP_ID']!=t.ix[i-1,'GROUP_ID']:
#            iE=i-1
#            if iB==iE:
#                t.ix[iB, 'CLUSTER_LABEL']=t.ix[iB, 'Description']
#            else:
#                S_rank=np.array([ c_category.get(int(t.ix[j,'CategoryID']), 100) for j in range(iB, iE+1)], dtype=int)
#                idx=S_rank.argmin()
#                t.ix[iB+idx, 'CLUSTER_LABEL']=t.ix[iB+idx, 'Description']
#            iB=i
    t.sort_values(['CLUSTER_LABEL','GROUP_ID','LogP'], inplace=True)
    # sort by label, so that nodes with label to displayed will be shown on top, labels will not be covered
    # by other non-label nodes
    S_order=["canonicalName","Description","Category","LogP","Log(q-value)","GROUP_ID","Hits","Enrichment","Z-score","#GeneInGOAndHitList","#GeneInGO","#GeneInHitList","#TotalGeneInLibrary","%InGO","STDV %InGO","BestLogPInGroup","BestEnrichmentInGroup","FirstInGroupByEnrichment","FirstInGroupByLogP","GeneID","_RANK_","_PATTERN_","CLUSTER_LABEL"]
    S_header=t.header()
    S=set(S_header)
    S0=[x for x in S_order if x in S]
    S0+=[x for x in S_header if x.startswith('#COUNT_')]
    S=set(S0)
    S0+=[x for x in S_header if x not in S]
    net.T_node=t[S0].copy()

class CytoPlot:

    def __init__(self, s_json_dir=None, host="localhost", port="1234", on_windows=False, add_logo=True):
        self.cc=None
        self.on_windows=on_windows
        self.add_logo=add_logo
        try:
            self.cc=cyto.CytoscapeClient(host=host, port=port)
        except:
            util.warn_msg('Cannot connect to Cytoscape server!')
            print(traceback.format_exc())

        if s_json_dir is None:
            s_json_dir=os.path.join(os.path.dirname(__file__), "cytoscape")
        self.s_json_dir=s_json_dir
        import ms.mssetting
        self.cyjs_src=ms.mssetting.cytoscapejs['SOURCE']

    def change_edge_color(self, X, edge_color=None):
        if edge_color is None:
            return
        for x in X["defaults"]:
            if x.get("visualProperty", "")=="EDGE_STROKE_UNSELECTED_PAINT":
                x["value"]=edge_color


    def get_ColorByPValue(self, brewer_name="YlOrBr", P=None, edge_color=None, s_json_name="ColorByPValue.json", n_node=100):
        X=json.loads(util.read_string(os.path.join(self.s_json_dir, s_json_name)))
        for x in X["mappings"]:
            if x.get("visualProperty", "")=="NODE_FILL_COLOR":
                if P is None: P=[-20.0,  -10.0, -6.0, -4.0, -3.0, -2.0]
                C=brewer2mpl.get_map(brewer_name, 'sequential', len(P), reverse=True).hex_colors[:]
                #C.insert(0, C[0])
                #print len(P), len(C)
                vp=[]
                for i in range(len(C)):
                    vp.append({'lesser':C[i], 'equal':C[i], 'value':P[i], 'greater':C[i]})
                x["points"]=vp
        self.change_edge_color(X, edge_color)
        self.adjust_label_size(X, n_node=n_node)
        return X

    def get_ColorByCounts(self, n=3, brewer_name=None, edge_color=None, s_json_name="ColorByCounts.json", S_color=None, n_node=100):
        if n>12: util.warn_msg("n exceeds 12, recycle color!")
        X=json.loads(util.read_string(os.path.join(self.s_json_dir, s_json_name)))
        for x in X["defaults"]:
            if x.get("visualProperty", "")=="NODE_CUSTOMGRAPHICS_1":
                if S_color is None:
                    if brewer_name is not None:
                        C=brewer2mpl.get_map(brewer_name, 'qualitative', min(max(n,3),12)).hex_colors[:]
                    else:
                        C=CytoPlot.get_qualitative_colors(n)
                    C=C[:n]
                else:
                    C=S_color
                if len(C)<n:
                    C=[C[i%len(C)] for i in range(n)]
                col=["#COUNT_%03d" % (i+1) for i in range(n)]
                x["value"]="org.cytoscape.PieChart:"+json.dumps({"cy_range":[1, n], "cy_colors": C, "cy_dataColumns": col, "cy_borderWidth":0.0})
        self.change_edge_color(X, edge_color)
        self.adjust_label_size(X, n_node=n_node)
        return X

    @staticmethod
    def get_qualitative_colors(n):
        if n>21: util.warn_msg("n exceeds 21, recycle color!")
        if n<=9:
            C=brewer2mpl.get_map("Set1", 'qualitative', max(n,3)).hex_colors[:]
            C=C[:n]
        elif n<=12:
            C=brewer2mpl.get_map("Paired", 'qualitative', n).hex_colors[:]
        else:
            C=brewer2mpl.get_map("Set1", 'qualitative', 9).hex_colors[:]
            C+=brewer2mpl.get_map("Set3", 'qualitative', min(n-9, 12)).hex_colors[:]
            if len(C)<n:
                C=[C[i%len(C)] for i in range(n)]
        return C

    def adjust_label_size(self, X, n_node=100):
        for x in X["defaults"]:
            if x.get("visualProperty", "")=="NODE_LABEL_FONT_SIZE":
                x["value"]=max(math.ceil(20*math.sqrt(n_node/100.0)-0.5),4.0)
                #print ">>>>>>>>>>>", n_node, x["value"]

    def get_ColorByCluster(self, n=10, brewer_name=None, edge_color=None, max_cluster=0, s_json_name="ColorByCluster.json", n_node=100):
        """max_cluster, maximum number of colors, if n>max_cluster, the extra cluster are colored by grey"""
        X=json.loads(util.read_string(os.path.join(self.s_json_dir, s_json_name)))
        m=n
        if max_cluster>0:
            m=min(n, max_cluster)
        for x in X["mappings"]:
            if x.get("visualProperty", "")=="NODE_FILL_COLOR":
                if brewer_name is None:
                    C=CytoPlot.get_qualitative_colors(m)
                else:
                    C=brewer2mpl.get_map(brewer_name, 'qualitative', m).hex_colors[:]
                vp=[{"key":"0.0", "value":"#BCBDDC"}]
                for i,c in enumerate(C):
                    vp.append({"key":str(i+1.0), "value": c})
                for i in range(m, n):
                    vp.append({"key":str(i+1.0), "value": "#eeeeee"})
                x["map"]=vp
        self.change_edge_color(X, edge_color)
        self.adjust_label_size(X, n_node=n_node)
        return X

    def style_remove_node_label(self, X):
        X2=copy.copy(X)
        if "style" in X2:
            X2["style"][0]["css"]["content"]=""
        return X2

    def get_default_style(self):
        #return self.cc.cy.style.get("default", "cytoscapejs")
        X=json.loads(util.read_string(os.path.join(self.s_json_dir, 'default.json')))
        return X[0]

    def get_DefaultStyle(self, edge_color=None, l_label=False):
        if l_label:
            X=json.loads(util.read_string(os.path.join(self.s_json_dir, 'DefaultStyleLabel.json')))
        else:
            X=json.loads(util.read_string(os.path.join(self.s_json_dir, 'DefaultStyle.json')))
        self.change_edge_color(X, edge_color)
        return X

    def apply_style(self, name, style_json, s_xgmml, s_out_dir=".", r_delay=0.0, S_fmt=None, prefix="", l_overwrite=True, l_bundling=True):
        # Lesson: we should creat a new network for each style, otherwise, result is random
        if os.path.exists(os.path.join(s_out_dir, name+'.png')):
            os.remove(os.path.join(s_out_dir, name+'.png'))
        if type(s_xgmml) is str:
            n_try=3
            i_try=1
            while True:
                cynet=self.cc.cynet_from_xgmml(s_xgmml, name=prefix+re.sub('NoLabel','', name), l_bundle=False)
                if cynet is not None or i_try>=n_try: break
                i_try+=1
                time.sleep(0.5)
            if cynet is None:
                util.error_msg('Cannot load network from: %s' % s_xgmml)
            net_obj=self.cc.get_network(cynet.get_id())
            #print type(net_obj)
            # json file rename columns starting with '#' into '_'
            c_known=set(['_GeneInGOAndHitList','_GeneInHitList','_TotalGeneInLibrary','_GeneInGO'])
            c_rename={ x: re.sub(r'^_', '#', x) for x in net_obj.T_node.header() if x in c_known or re.search(r'^_COUNT_\d+$', x) }
            if len(c_rename):
                net_obj.T_node.rename2(c_rename)
        else:
            # must be a network objects with coordinates
            net_obj=s_xgmml
            #print "Net obj passed in"
            cynet=self.cc.cynet_from_network(s_xgmml, name=prefix+re.sub('NoLabel','', name), s_layout=None, l_bundle=False)
        #print ">>>>>>>>", cynet.get_id(), "<<<<<<<<"
        #print ">>ALL STYLES>>>>>>>>>>>>>>", self.cc.cy.style.get_all()
        #print ">>Create Style>>>>>>>>>>>>", name

        if name in self.cc.cy.style.get_all() and not l_overwrite:
            style_name=name
            #print "style %s already exists" % style_name
        else:
            self.cc.delete_style(name)
            style_name=self.cc.cy.style.create(name=name, original_style=style_json).get_name()
            #print ">> >> ", name, style_name
        #print ">>After Creation>>>>>>>>>>", self.cc.cy.style.get_all(), style_name
        mystyle=self.cc.cy.style.get(style_name, "cytoscapejs")
        #if len(mystyle)==0:
        #    print(style_name, mystyle)
        #    exit()
        # move fit to the beginning, as if I put fit() to the end, it seems to screw up style for cynet_save()
        #self.cc.cy.layout.fit(cynet)
        style=self.cc.get_style(style_name)
        self.cc.cy.style.apply(style, network=cynet)
        if l_bundling: self.cc.cy.edgebundling.apply(cynet)
        self.cc.cy.layout.fit(cynet)
        # repeat to recover
        #print "applying style %s" % style_name
        self.cc.cy.style.apply(style, network=cynet)
        #self.cc.cy.layout.fit(cynet)
        return (cynet, mystyle, net_obj)
        #time.sleep(r_delay)
        #if "png" in S_fmt:
        #    self.cc.cynet_save(cynet, os.path.join(s_out_dir, prefix+name+'.png'))
        #if "pdf" in S_fmt:
        #    self.cc.cynet_save(cynet, os.path.join(s_out_dir, prefix+name+'.pdf'))

    @staticmethod
    def fix_network_for_js(data):
        for net in data:
            for node in net["elements"]["nodes"]:
                c={}
                for k,v in node["data"].items():
                    if k.startswith("_COUNT_"):
                        c[k]=v
                if len(c):
                    total=sum(c.values())
                    for k,v in c.items():
                        node["data"][k+"_PCT"]=v*100/total
        return data

    @staticmethod
    def fix_style_for_js(data, S_color, l_PPI=False):
        for style in data:
            for x in style["style"]:
                if x["selector"]=="node" and "css" in x:
                    if l_PPI and 'NoLabel' not in style['title']:
                        x["css"]["content"]="data(Symbol)";
                    x["css"]["border-width"]=4.0
                    if "ColorByCounts" in style["title"]:
                        x["css"]["pie-size"] = "100%"
                        x["css"]["border-width"]=0.0
                        for i,c in enumerate(S_color):
                            x["css"]["pie-{:d}-background-color".format(i+1)]=c
                            x["css"]["pie-{:d}-background-size".format(i+1)]="data(_COUNT_{:03d}_PCT)".format(i+1)

            if l_PPI:
                style["style"].extend([
                    {
                        "selector" : "node[DEGREE<=5]",
                        "css" : {
                          "width": 20.0,
                          "height": 20.0
                        }
                    },
                    {
                        "selector" : "node[DEGREE>5][DEGREE<20]",
                        "css" : {
                          "width": "mapData(DEGREE,5,20,35.0,50.0)",
                          "height": "mapData(DEGREE,5,20,35.0,50.0)",
                        }
                    },
                    {
                        "selector" : "node[DEGREE>=20]",
                        "css" : {
                          "width": 50.0,
                          "height": 50.0
                        }
                    }])

            else:
                style["style"].extend([
                    {
                        "selector" : "node[_GeneInGOAndHitList<=5]",
                        "css" : {
                          "width": 20.0,
                          "height": 20.0
                        }
                    },
                    {
                        "selector" : "node[_GeneInGOAndHitList>5][_GeneInGOAndHitList<20]",
                        "css" : {
                          "width": "mapData(_GeneInGOAndHitList,5,20,35.0,50.0)",
                          "height": "mapData(_GeneInGOAndHitList,5,20,35.0,50.0)",
                        }
                    },
                    {
                        "selector" : "node[_GeneInGOAndHitList>=20]",
                        "css" : {
                          "width": 50.0,
                          "height": 50.0
                        }
                    }])

        return data


    def plot_go_network(self, s_xgmml, n_list, n_cluster, s_out_dir=".", S_fmt=None, s_session=None, S_label=None, S_color=None, n_node=100):
        try:
            self.cc.gc()
            self._plot_go_network(s_xgmml, n_list, n_cluster, s_out_dir=s_out_dir, S_fmt=S_fmt, s_session=s_session, S_label=S_label, S_color=S_color, n_node=n_node)
        except Exception as e:
            if 'Restart Please' in str(e):
                raise e
            else:
                import sys
                util.warn_msg('error during plot_go_network(): %s\n%s' % (s_xgmml, sys.exc_info()[0]))
                print(traceback.format_exc())

    def plot_ppi_network(self, S_xgmml, n_list, n_cluster, s_out_dir=".", S_fmt=None, s_session=None, S_label=None, S_color=None, S_by_count=None, max_cluster=20, S_node=None, s_legend_prefix=None):
        try:
            self.cc.gc()
            self._plot_ppi_network(S_xgmml, n_list, n_cluster, s_out_dir=s_out_dir, S_fmt=S_fmt, s_session=s_session, S_label=S_label, S_color=S_color, S_by_count=S_by_count, max_cluster=max_cluster, S_node=S_node, s_legend_prefix=s_legend_prefix)
        except Exception as e:
            if 'Restart Please' in str(e):
                raise e
            else:
                import sys
                util.warn_msg('Error during plot_ppi_network(): %s\n%s' % ("\n".join(S_xgmml), sys.exc_info()[0]))
                print(traceback.format_exc())

    def _plot_ppi_network(self, S_xgmml, n_list, n_cluster, s_out_dir=".", S_fmt=None, s_session=None, S_label=None, S_color=None, S_by_count=None, max_cluster=20, S_node=None, s_legend_prefix=None):
        """Make sure Cytoscape is already running, and Cytoscape window should not be minimized, otherwise, plots may not be updated correctly."""
        #if s_session is None:
        #    util.error_msg('Session file must be specified, in order to plot, as we now generate cys first, then load in to make plots!')
        self.cc.gc()
        self.cc.cy.network.delete_all()
        self.cc.cy.style.delete_all()
        l_bundling=False # no need for this in PPI
        c_style={}
        style_clr=self.get_ColorByCluster(n_cluster, edge_color='#54278F', max_cluster=max_cluster, s_json_name="PPIColorByCluster.json")
        #self.cc.cy.style.create(name='PPIColorByCluster', data=style_clr).get_name()
        #x=self.cc.cy.style.get('PPIColorByCluster', "cytoscapejs")
        #c_style['PPIColorByCluster']=x
        style_clr2=self.get_ColorByCluster(n_cluster, edge_color='#54278F', max_cluster=max_cluster, s_json_name="PPIColorByClusterNoLabel.json")
        #self.cc.cy.style.create(name='PPIColorByClusterNoLabel', data=style_clr2).get_name()
        #x=self.cc.cy.style.get('PPIColorByClusterNoLabel', "cytoscapejs")
        #c_style['PPIColorByClusterNoLabel']=x

        r_delay=0.0
        S_json=[]
        if S_node is None:
            S_node=[100]*len(S_xgmml)
        #print "INNNNNNNNNNN", S_xgmml, S_node
        for i,s_xgmml in enumerate(S_xgmml):
            #print(">>>>>>>>>>", i, s_xgmml)
            prefix, s_ext=os.path.splitext(os.path.basename(s_xgmml))
            l_colorByCluster=True
            l_colorByCounts=False
            net_obj=None
            if S_by_count is not None and s_xgmml in S_by_count:
                #print "^^^^^^^^^^", i, s_xgmml
                s_style="PPIColorByCounts" if prefix.endswith('_MCODE_ALL') else 'PPIColorByCountsNoLabel'
                #print n_list, s_style, S_color, S_node, i
                style=self.get_ColorByCounts(n=n_list, edge_color='#54278F', s_json_name=s_style+".json", S_color=S_color, n_node=S_node[i])
                #print "AAAAAAAAAAAAAAAAAAAA",  self.cc.cy.style.get_all(), s_style
                (cynet, x, net_obj)=self.apply_style(s_style, style, s_xgmml, s_out_dir=s_out_dir, S_fmt=S_fmt, r_delay=r_delay, prefix=prefix+"_", l_overwrite=True, l_bundling=l_bundling)
                if s_style not in c_style: c_style[s_style]=x
                l_colorByCounts=True
                #if '_MCODE_ALL' not in prefix:
                #    l_colorByCluster=False # no need to color by MCODE for the full cluster at this point
            if l_colorByCluster:
                s_style="PPIColorByCluster" if prefix.endswith('_MCODE_ALL') else 'PPIColorByClusterNoLabel'
                style=style_clr if prefix.endswith('_MCODE_ALL') else style_clr2
                #print "$$$$$$$$$$$", s_style, l_colorByCluster, l_colorByCounts, type(net_obj), type(style), s_xgmml, prefix
                if 'NoLabel' not in s_style:
                    self.adjust_label_size(style, n_node=S_node[i])
                l_overwrite=s_style not in c_style
                #print ">>c_style>>>>>>>", c_style.keys()
                y=s_xgmml
                if l_colorByCounts:
                    y=net_obj
                #print type(y)
                #if type(y) is str: print y, s_style, prefix
                #print "BBBBBBBBBBBBBBBBBBBBB",  self.cc.cy.style.get_all(), s_style
                (cynet, x, net_obj)=self.apply_style(s_style, style, y, s_out_dir=s_out_dir, S_fmt=S_fmt, r_delay=r_delay, prefix=prefix+"_", l_overwrite=l_overwrite, l_bundling=l_bundling)
                if s_style not in c_style: c_style[s_style]=x
                #print "CCCCCCCCCCCCCC",  self.cc.cy.style.get_all(), s_style
            json_data=cynet.get_first_view() #to_json())
            s_name=os.path.join(s_out_dir, prefix+".json")
            util.save_string(s_name, json.dumps(json_data))
            S_json.append(json_data)
        S_style=[ v for k,v in c_style.items() ]
        S_style.append(self.get_default_style())
        S_style=[self.style_remove_node_label(x) for x in S_style]
        cyto.save_style_js(os.path.join(s_out_dir, "PPINetwork.style.js"), CytoPlot.fix_style_for_js(S_style, S_color=S_color, l_PPI=True))
        cyto.save_network_js(os.path.join(s_out_dir, "PPINetwork.js"), CytoPlot.fix_network_for_js(S_json))
        s_template=os.path.join(self.cyjs_src, 'CyJS.template.html')
        cyto.save_cyjs_app(os.path.join(s_out_dir, "PPINetwork.html"), "PPINetwork", None, "../CyJS", s_template)

        with open(os.path.join(s_out_dir, "PPI.style.json"), "w") as f:
            json.dump(S_style, f)

        if s_session is not None:
            #import ms.report as rpt
            s_cys=s_session
            if self.on_windows: #rpt.is_cytoscape_on_WIN():
                s_cys=util.format_path(s_session, 'win')
            self.cc.cy.session.save(s_cys)
        if S_fmt is not None:
            self.plot_session(s_session=s_session, s_out_dir=s_out_dir, S_fmt=S_fmt, S_label=S_label, S_color=S_color, n_cluster=n_cluster, max_cluster=max_cluster, l_try_memory_first=True, s_legend_prefix=s_legend_prefix)

    def _plot_go_network(self, s_xgmml, n_list, n_cluster, s_out_dir=".", S_fmt=None, s_session=None, S_label=None, S_color=None, max_cluster=0, n_node=100, l_try_memory_first=False):
        """Make sure Cytoscape is already running, and Cytoscape window should not be minimized, otherwise, plots may not be updated correctly."""
        if S_fmt is None: S_fmt=['png']
        #if s_session is None:
        #    util.error_msg('Session file must be specified, in order to plot, as we now generate cys first, then load in to make plots!')
        self.cc.gc()
        self.cc.cy.network.delete_all()
        self.cc.cy.style.delete_all()
        l_bundling=True
        S_style=[]
        r_delay=0.0
        S_cynet=[]
        S_name=[]
        style=self.get_ColorByCluster(n=n_cluster, edge_color='#54278F', max_cluster=0, n_node=n_node)
        (cynet, x, net_obj)=self.apply_style("ColorByCluster", style, s_xgmml, s_out_dir=s_out_dir, S_fmt=S_fmt, r_delay=r_delay, l_bundling=l_bundling)
        # get cluster labels
        S_legend_label=[""]*n_cluster
        for i in net_obj.T_node.index:
            if len(net_obj.T_node.loc[i,'CLUSTER_LABEL'])>1:
                S_legend_label[int(round(net_obj.T_node.loc[i,'GROUP_ID']-1))]=net_obj.T_node.loc[i,'CLUSTER_LABEL']
        S_cynet.append(cynet)
        S_style.append(x)
        S_name.append("ColorByCluster")
        style=self.get_ColorByPValue(edge_color='#54278F', n_node=n_node)
        (cynet, x, net_obj)=self.apply_style("ColorByPValue", style, net_obj, s_out_dir=s_out_dir, S_fmt=S_fmt, r_delay=r_delay, l_bundling=l_bundling)
        S_cynet.append(cynet)
        S_style.append(x)
        S_name.append("ColorByPValue")
        if n_list>1:
            style=self.get_ColorByCounts(n=n_list, edge_color='#54278F', S_color=S_color, n_node=n_node)
            #S_style.append(style)
            (cynet, x, net_obj)=self.apply_style("ColorByCounts", style, net_obj, s_out_dir=s_out_dir, S_fmt=S_fmt, r_delay=r_delay, l_bundling=l_bundling)
            S_cynet.append(cynet)
            S_style.append(x)
            S_name.append("ColorByCounts")
        s_name_network="GONetwork"
        if s_session is not None:
            #import ms.report as rpt
            s_cys=s_session
            if self.on_windows: #rpt.is_cytoscape_on_WIN():
                s_cys=util.format_path(s_session, 'win')
            self.cc.cy.session.save(s_cys)
            S_id=self.cc.cy.network.get_all()
            if len(S_id):
                cynet=self.cc.get_cynet(S_id[0])
                #s_name, s_ext=os.path.splitext(util.format_path(s_session))
                s_name, s_ext=os.path.splitext(util.format_path(s_xgmml))
                s_name_network=os.path.basename(s_name)
                data=cynet.get_first_view()
                data["data"]["shared_name"]=data["data"]["name"]=s_name_network
                s_json=json.dumps(data)
                util.save_string(s_name+".json", s_json)
                cyto.save_network_js(s_name+".js", CytoPlot.fix_network_for_js([data]))
        if len(S_style)>0:
            S_style.append(self.get_default_style())
            S_style=[self.style_remove_node_label(x) for x in S_style]
            s_file=s_xgmml.replace('.xgmml', '')
            with open(s_file+".style.json", "w") as f:
                json.dump(S_style, f)
            cyto.save_style_js(s_file+".style.js", CytoPlot.fix_style_for_js(S_style, S_color=S_color, l_PPI=False))
            s_template=os.path.join(self.cyjs_src, 'CyJS.template.html')
            cyto.save_cyjs_app(s_file+".html", s_name_network, None, "../CyJS", s_template)

        self.plot_session(s_session=s_session, s_out_dir=s_out_dir, S_fmt=S_fmt, S_label=S_label, S_color=S_color, n_cluster=n_cluster, max_cluster=max_cluster, l_try_memory_first=True, s_legend_prefix=S_legend_label)

    @staticmethod
    def is_good_network_image(s_png, min_fraction=0.1, min_regions=0):
        # bad images have most nodes as ovals (not circles) and all nodes a colored in certain red/blue
        # so we first check how many pixles are in red/blue, compared to non-white
        # if the fraction is too low, it's likely a good image
        # then we check how many red/blue regions, it has to exceed min_regions
        # (if we have 100 node, we probably expect 50 of them are vislable, so min_regions can be set to 50)
        # Then we check the aspect ratio, the width and heigth of the max object should be close to 1.0, a circle
        DEBUG=True
        #import scipy.misc
        #M=mahotas.imread(s_png)
        import imageio
        M=imageio.imread(s_png)
        #M*=255
        h,w,z=M.shape
        r=M[:,:,0]
        g=M[:,:,1]
        b=M[:,:,2]
        # assume a bad image has red ovals (# sometimes are blue circle)
        # bad image has most of the foreground (ovals) colored as red or blue
        fg=h*w-np.sum(np.logical_and(r==255, g==255, b==255))
        m=np.logical_and(r==200, g==0, b==0)
        if np.sum(m)<fg*min_fraction:
            # try blue color
            if DEBUG: util.info_msg('Not enough red')
            m=np.logical_and(r==0, g==153, b==204)
            if np.sum(m)<fg*min_fraction:
                if DEBUG: util.info_msg('Not enough blue')
                return True
        #mahotas.imsave('mask.png', m)
        import scipy.ndimage
        labeled,nr_regions = scipy.ndimage.label(m)
        # expecting certain # of circles/ovals
        if (nr_regions==0) or (min_regions>0 and nr_regions<min_regions):
            if DEBUG: util.info_msg('No enough nodes, nr_regions=%d, min_region=%d' % (nr_regions, min_regions))
            return True
        sz_x=np.zeros(nr_regions)
        sz_y=np.zeros(nr_regions)
        # find the oval with max sizes
        for i in range(nr_regions):
            X,Y=np.where(labeled==i+1)
            sz_x[i]=X.max()-X.min()
            sz_y[i]=Y.max()-Y.min()
        # if y is larger than x, it's an oval, not a circle, therefore image is bad
        max_y=sz_y.max()
        max_x=sz_x.max()
        ratio=max_y/max_x
        if DEBUG: util.info_msg('Max: X=%d y=%d ratio=%.2f' % (max_x, max_y, ratio))
        if max_x>5 and max_y>5 and ratio > 1.2 or ratio < 0.8:
            return False
        return True

    def plot_session(self, s_session=None, s_out_dir=".", S_fmt=None, S_network=None, S_label=None, S_color=None, n_cluster=20, max_cluster=0, l_try_memory_first=False, s_legend_prefix="Cluster"):
        import shutil
        import ms.report as rpt
        sw=util.StopWatch("Cytoplot::plot_session")
        if S_color is None:
            S_color=CytoPlot.get_qualitative_colors(len(S_label))
        for trial in range(3):
            if S_network is not None and len(S_network)==0: break
            #print("Trial %d: " % trial, s_session)
            if s_session is None or (trial==0 and l_try_memory_first):
                print("using existing session")
                pass
            elif s_session is not None:
                s_cys=s_session
                if self.on_windows: #rpt.is_cytoscape_on_WIN():
                    s_cys=util.format_path(s_session, 'win')
                #if "\\" in s_session: # windows path
                #    s_cys=s_session.replace("\\", "\\\\")
                print("loading session file: ", s_cys)
                self.cc.cy.session.open(s_cys)
            S_id=self.cc.cy.network.get_all() # every time we get a different ID
            S_bad_net=[]
            for id in S_id:
                cynet=self.cc.get_cynet(id)
                s_name=cynet.get_name()
                if (S_network is not None) and (s_name not in S_network): continue
                l_good=True
                S_png=[]
                if "png" in S_fmt:
                    s_file=os.path.join(s_out_dir, s_name+'.png')
                    self.cc.cynet_save(cynet, s_file)
                    # maybe there is no View created, somehow, thereis no expored image
                    if not os.path.exists(s_file):
                        util.warn_msg('File cannot be created: %s' % s_file)
                        continue
                    if os.path.getsize(s_file)==0:
                        util.warn_msg('File cannot be created: %s' % s_file)
                        os.remove(s_file)
                        continue
                    # check if file is in the right format
                    try:
                        rpt.crop(s_file, right=100)
                        l_good=CytoPlot.is_good_network_image(s_file)
                    except:
                        l_good=False
                        print("Host=", self.cc.host, "Port=", self.cc.port, "Memory:", self.cc.cy.status())
                        print(s_file, os.stat(s_file).st_size)
                        if os.stat(s_file).st_size<=500:
                            s=""
                            with open(s_file, 'rb') as f:
                                s=f.read().decode(errors='ignore')
                            print("png content: "+s)
                            #{"data":{},"errors":[{"status":500,"type":"urn:cytoscape:ci:cyrest-core:v1:error-handling:errors:0","message":"Uncaught exception while processing resource [qtp270404096-96]: No rendering engine for {\"network\":2018878, \"view\":2022703}","link":"file:/home/cytoscape/CytoscapeConfiguration/3/framework-cytoscape.log"}]}
                            if 'Uncaught exception while processing resource' in s and 'No rendering engine for' in s:
                                raise NameError('Cytoscape Internal Problem, Restart Please!')
                    #sw.check('Save png: %s, check good/bad: %s' % (s_file, 'GOOD' if l_good else 'BAD'))
                    if not l_good:
                        print(">>Found Bad Image: %s, trial=%d" % (s_name, trial))
                        S_bad_net.append(s_name)
                        shutil.copyfile(s_file, s_file.replace('.png', '.trial'+str(trial+1)+'.png'))
                    if (l_good or trial==2): # add legend watermark
                        if self.add_logo: # add_metascape_logo():
                            S_png.append(os.path.join(os.path.dirname(__file__), "watermark", "Watermark.png"))
                        if 'ColorByCounts' in s_name and S_label is not None:
                            s_png=os.path.join(s_out_dir, s_name+'.legend.png')
                            rpt.make_legend(S_color, S_label, s_png)
                            S_png.append(s_png)
                        elif 'ColorByPValue' in s_name:
                            S_png.append(os.path.join(os.path.dirname(__file__), "watermark", "legend_P.png"))
                        elif 'ColorByCluster' in s_name:
                            s_png=os.path.join(s_out_dir, s_name+'.legend.png')
                            if type(s_legend_prefix) is dict:
                                s_legend=s_legend_prefix.get(re.sub(r'_(PPI)?ColorByCluster','', s_name), None)
                            else:
                                s_legend=s_legend_prefix
                            #print "LEGEND>>>>>", s_legend
                            if s_legend is not None:
                                rpt.make_legend_by_cluster(n_cluster, s_png, max_cluster=max_cluster, s_label=s_legend, max_width=300)
                                S_png.append(s_png)
                        if len(S_png):
                            rpt.make_legend_png(s_file, S_png, l_clean_up=False)
                            #sw.check('Make legend: %s' % " ".join(S_png))
                if (l_good or trial==2) and "pdf" in S_fmt:
                    #print s_name+'.pdf'
                    s_pdf=os.path.join(s_out_dir, s_name+'.pdf')
                    s_legend=os.path.join(s_out_dir, s_name+'.legend.pdf')
                    self.cc.cynet_save(cynet, s_pdf)
                    #sw.check('Save pdf: %s, check good/bad: %s' % (s_pdf, 'GOOD' if l_good else 'BAD'))
                    #print s_legend, S_png, s_pdf
                    # comment out the following, b/c it takes 5-8 seconds to add legend and logo to PDF
                    if len(S_png):
                        rpt.make_legend_pdf(s_legend, S_png, l_clean_up=False)
                        #sw.check('Make legend: %s' % s_legend)
                        rpt.add_watermark_pdf([s_pdf], s_out=s_pdf, s_wm=s_legend, l_clean_up=True)
                        for x in S_png:
                            if x.startswith(s_out_dir):
                                os.remove(x)
                        #sw.check('Make new pdf: %s' % s_pdf)

            S_network=S_bad_net
            #if s_session is not None:
            #    self.cc.cy.network.delete_all()
            print("Bad networks, trial=", S_network, trial)
        #    print S_network, trial
        #self.cc.cy.network.delete_all()

    #def plot_network(self, s_xgmml, s_out_dir=".", S_fmt=None, s_session=None):
    #    """No longer in use.
    #    Make sure Cytoscape is already running, and Cytoscape window should not be minimized, otherwise, plots may not be updated correctly."""
    #    if S_fmt is None: S_fmt=['png']
    #    self.cc.cy.network.delete_all()
    #    style=self.get_DefaultStyle(edge_color='#54278F')
    #    r_delay=0.0
    #    self.apply_style("DefaultStyle", style, s_xgmml, s_out_dir=s_out_dir, S_fmt=S_fmt, r_delay=r_delay)
    #    style=self.get_DefaultStyle(edge_color='#54278F', l_label=True)
    #    self.apply_style("DefaultStyleLabel", style, s_xgmml, s_out_dir=s_out_dir, S_fmt=S_fmt, r_delay=r_delay)
    #    if s_session is not None:
    #        self.cc.cy.session.save(s_session)
    #    self.cc.cy.network.delete_all()


class DendroHeatmapPlot:

    @staticmethod
    def plot_membership(t_membership, s_out, s_symbol=None, S_label=None):
        S_col=[x for x in t_membership.header() if x.startswith('_MEMBER_')]
        S_exp=[x[8:] for x in S_col]
        if S_label is not None:
            S_exp=S_label
        if len(S_exp)<2:
            util.warn_msg("Too few columns to cluster!")
            return
        if len(t_membership)<2: util.error_msg("Too few rows to cluster!")
        if s_symbol is not None and (s_symbol in t_membership.header):
            S_gene=list(t_membership.loc[:, s_symbol])
        else:
            S_gene=list(t_membership.index)
        X=t_membership.loc[:, S_col].values
        import cluster
        cluster.FastCluster.quick_plot(X, s_out, S_gene, S_exp)

class PiePlot:

    @staticmethod
    def brighter(S_hex, i_increase=50, i_max=250):
        #code copied from http://stackoverflow.com/questions/4296249/how-do-i-convert-a-hex-triplet-to-an-rgb-tuple-and-back
        """Make input list of hex colors brighter by 50/255. Use to generate background colors based on foreground colors."""
        _NUMERALS = '0123456789abcdefABCDEF'
        _HEXDEC = {v: int(v, 16) for v in (x+y for x in _NUMERALS for y in _NUMERALS)}
        LOWERCASE, UPPERCASE = 'x', 'X'

        def rgb(triplet):
            return _HEXDEC[triplet[0:2]], _HEXDEC[triplet[2:4]], _HEXDEC[triplet[4:6]]

        def triplet(rgb, lettercase=LOWERCASE):
            return format(rgb[0]<<16 | rgb[1]<<8 | rgb[2], '06'+lettercase)

        def hex2hex(x, i_increase):
            if x.startswith('#'): x=x[1:]
            r,g,b=rgb(x)
            r,g,b=(min(r+i_increase, i_max), min(g+i_increase, i_max), min(b+i_increase, i_max))
            return '#'+triplet((r,g,b))

        return [hex2hex(x, i_increase) for x in S_hex]

    @staticmethod
    def plot(N_total, N_in_GO, I_n_hit, I_n_in_GO, s_out, S_color=None, S_fmt=['png', 'pdf']):
        import matplotlib as mpl
        mpl.rcParams['pdf.fonttype'] = 42
        mpl.use('Agg')
        import matplotlib.pyplot as plt
        if not type(I_n_hit) is list:
            I_n_hit=[I_n_hit]
            I_n_in_GO=[I_n_in_GO]
        n=len(I_n_hit) # number of hit lists
        fig=plt.figure()
        C_bg=brewer2mpl.get_map("Greys", 'sequential', 5).hex_colors[:]
        if S_color is None:
            C_yes=CytoPlot.get_qualitative_colors(n)
            C_no=PiePlot.brighter(C_yes, 80, 220)
        else:
            C_yes=S_color
            C_no=PiePlot.brighter(C_yes, 80, 220)
        n_row=int(round(math.sqrt(n*3.0/4.0)))
        n_col=int(math.ceil(n*1.0/n_row))
        for i in range(n):
            n_in_GO=I_n_in_GO[i]
            n_hit=I_n_hit[i]
            ax=fig.add_subplot(n_row, n_col, i+1)
            ax.set_aspect('equal')
            wc=ax.pie([N_in_GO, N_total-N_in_GO], colors=[C_bg[4], C_bg[1]], startangle=90, radius=1.0)
            plt.text(0.05, 0.75, '{0:.2%}, {1:d}'.format(N_in_GO*1.0/N_total, N_in_GO), fontsize=8)
            for w in wc[0]:
                w.set_edgecolor('#FFFFFF')
            wc=ax.pie([n_in_GO, n_hit-n_in_GO], colors=[C_yes[i % n], C_no[i % n]], startangle=90, radius=0.6)
            plt.text(0.05, 0.3, '{0:.2%}, {1:d}'.format(n_in_GO*1.0/n_hit, n_in_GO), fontsize=8)
            for w in wc[0]:
                w.set_edgecolor('#FFFFFF')
            p=stats.hyper(n_in_GO, N_total, n_hit, N_in_GO)
            ax.set_title('p=%.2g' % p, fontsize=10)
        if s_out.endswith('.png'): s_out=s_out[:-4]
        for x in S_fmt:
            plt.savefig(s_out+'.'+x, bbox_inches='tight')
        plt.close()

if __name__=="__main__":
    import ms.biolists as bl

    # test PPI network analysis
    sw=util.StopWatch()
    if not os.path.exists('output'):
        os.mkdir('output')
    bl.PiePlot.plot(10000, 200, [300], [20], 'output/membership.png')
    bl.PiePlot.plot(10000, 200, [300, 200, 250, 400], [20,5, 8, 10], 'output/membership.png')
    # build gene lists
    ### We upload three hit lists
    S=["ChandaConfirmedHits.csv", "Elledge 284.csv", "MercksiRNA.csv"]
    input=[]
    s_dir=os.path.dirname(bl.__file__)
    for x in S:
        t=pd.read_csv(s_dir+'/datasets/'+x)
        input.append(bl.GeneList(x, list(t['Gene ID'])))
    # create a GeneLists object
    lists=bl.GeneLists(input)

    ### Single list PPI analysis
    import entrez as ez
    myez=ez.EntrezGene(l_use_GPDB=True)
    mygo=go.GO(l_use_GPDB=True, entrez=myez)
    import ppi
    myppi=ppi.HumanPPI(S_DB=["GeneGo"], minScore=0.5, entrez_gene=myez)

    cp=bl.CytoPlot()

    onelist=input[0]
    netwks=onelist.ppi_analysis(myppi, "output", "ListOne")
    cp.plot_network("output/ListOne.xgmml", s_out_dir="output", S_fmt=["png","pdf"], s_session="C:\\Temp\\default.cys")
    ###

    sw.check("Gene List Loaded")
    # overlap analysis
    t_membership=lists.membership_table(l_gene_as_index=True)
    sw.check("Membership")
    t_membership.to_csv('output/x.csv')
    # plot dendrogram and Circos overlap plot
    #bl.DendroHeatmapPlot.plot_membership(t_membership, 'output/hit.png', S_label=['Chanda', 'Elledge', 'Merck'])
    #sw.check("Dendrogram")

    from . import circos
    cir=circos.Circos(BIN=None)
    cir.plot_membership(t_membership, outputdir="output", outputfile="CircosPlot", S_label=['Chanda', 'Elledge', 'Merck'])
    #util.unix("/bin/mv -f circos/CircosPlot.png circos/CircosPlot.svg output", l_print=False)
    #sw.check("Circos")

    # perform GO-enrichment analysis
    #mygo=go.GO(l_use_GPDB=True)
    print("Enrichment ...")
    sw.check("Load GO")
    golists=lists.go_analysis(mygo, n_CPU=12)
    sw.check("Enrichment")

    # cluster and combined three GO lists into one
    print("Merge ...")
    golist=golists.merge(go=mygo)

    # plot new Circos and heatmap by considering indirect-GO overlap
    t_go=golist.data
    if t_go is None or len(t_go)==0:
        util.warn_msg('No enrichmed GO terms for dataset')
    else:
        t_go=t_go[(t_go.LogP<=-3) & (t_go['#GeneInGO']<=100) & (t_go.Enrichment>=2.0)]
        t_go=t_go[['GeneID']]
        sw.check("Merge GO")
        cir.plot_membership(t_membership, t_go=t_go, outputdir="output", outputfile="CircosPlot_GO", S_label=['Chanda', 'Elledge', 'Merck'])

        #util.unix("/bin/mv -f circos/CircosPlot_GO.png circos/CircosPlot_GO.svg output", l_print=False)
        #sw.check("Circos")
        # side effect, t_membership now contains 0.5 elements
        bl.DendroHeatmapPlot.plot_membership(t_membership, 'output/hit_by_go.png', S_label=['Chanda', 'Elledge', 'Merck'])
        sw.check("Dendrogram")

        # Cluster enriched GO terms into GO network
        t_go=golist.cluster(n_CPU=12, similarity=0.3)
        t_go_mem=golists.membership_table()
        golist.heatmap('output/go_heatmap', t_go_mem, S_label=['Chanda', 'Elledge', 'Merck'])
        net=golist.network(max_clusters=20, max_members=10, max_nodes=250)
        # need to add gene counts
        print("Add node count ...")
        sw.check("GO Clustering")
        bl.net_annotation(net, t_membership)
        t_go.to_csv('output/x.csv', index=False)
        s_xgmml='output/GONetwork'
        net.to_xgmml(s_xgmml)
        n_list=len(lists)
        n_cluster=20
        cp=bl.CytoPlot()
        cp.plot_go_network(s_xgmml+".xgmml", n_list, n_cluster, s_out_dir="output", S_fmt=["png","pdf"], s_session="C:\\Temp\\test.cys")
        sw.check("CytoPlot")
