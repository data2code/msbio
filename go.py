#!/usr/bin/env python
from __future__ import absolute_import
from __future__ import print_function
from __future__ import division
import numpy as np
import pandas as pd
import re
import util
import os
import entrez as ez
import stats
import parallel
import xgmml
import db
import random
from six.moves import range
import setting

class Cache(object):
    DATA_DIR=setting.go['DATA_DIR']
    # share by all tax_id
    CATEGORY={'LOCAL':None, 'GPDB':None, 'L1k':None}
    GO_DESCRIPTION={'LOCAL':None, 'GPDB':None, 'L1k':None}
    GO_CATEGORY={'LOCAL':None, 'GPDB':None, 'L1k':None}
    PARENT={'LOCAL':None, 'GPDB':None, 'L1k':None}
    # per tax_id
    TOTAL_GENE_COUNT={'LOCAL':{}, 'GPDB':{}, 'L1k':{}}
    ALL_GENE={'LOCAL':{}, 'GPDB':{}, 'L1k':{}}
    CATEGORY_COUNT={'LOCAL':{}, 'GPDB':{}, 'L1k':{}}
    GO_GENE_ENRICH={'LOCAL':{}, 'GPDB':{}, 'L1k':{}}
    GO_GENE={'LOCAL':{}, 'GPDB':{}, 'L1k':{}}
    GO_PARENT={'LOCAL':{}, 'GPDB':{}, 'L1k':{}}
    GENE_GO={'LOCAL':{}, 'GPDB':{}, 'L1k':{}}
    N_TRIVIAL=800

    @staticmethod
    def dump():
        def swap(c):
            if 'LOCAL' in c:
                c['LOCAL'], c['GPDB']=c['GPDB'], c['LOCAL']
            return c

        util.dump_object([swap(Cache.CATEGORY), swap(Cache.GO_DESCRIPTION), swap(Cache.GO_CATEGORY), \
            swap(Cache.PARENT), swap(Cache.TOTAL_GENE_COUNT), swap(Cache.ALL_GENE), swap(Cache.CATEGORY_COUNT), \
            swap(Cache.GO_GENE_ENRICH), swap(Cache.GO_GENE), swap(Cache.GO_PARENT), swap(Cache.GENE_GO)],
            os.path.join(Cache.DATA_DIR, "GO.pickle.gz"))

    @staticmethod
    def get(l_use_GPDB=True, tax_id=9606, l_L1k=False, S_go_category=None):
        if tax_id==-9606: tax_id=9606
        s_key=Cache.key(l_use_GPDB, l_L1k)
        if not l_use_GPDB:
            X=util.load_object(os.path.join(Cache.DATA_DIR, "GO.pickle.gz"))
            (Cache.CATEGORY, Cache.GO_DESCRIPTION, Cache.GO_CATEGORY, Cache.PARENT, \
            Cache.TOTAL_GENE_COUNT, Cache.ALL_GENE, Cache.CATEGORY_COUNT, \
            Cache.GO_GENE_ENRICH, Cache.GO_GENE, Cache.GO_PARENT, Cache.GENE_GO) = X

        if l_L1k and tax_id!=9606:
            util.error_msg('L1k is only for tax_id 9606!')
        if tax_id not in Cache.TOTAL_GENE_COUNT[s_key]:
            Cache.load(tax_id=tax_id, l_use_GPDB=l_use_GPDB, l_L1k=l_L1k, S_go_category=S_go_category)
            #if l_L1k:
            #    Cache.loadL1k()
            #else:
            #    Cache.load(tax_id=tax_id, l_use_GPDB=l_use_GPDB)
        return (Cache.CATEGORY[s_key],
            Cache.GO_DESCRIPTION[s_key],
            Cache.GO_CATEGORY[s_key],
            Cache.PARENT[s_key],
            # per tax_id, above are shared across tax_id
            Cache.TOTAL_GENE_COUNT[s_key][tax_id],
            Cache.ALL_GENE[s_key][tax_id],
            Cache.CATEGORY_COUNT[s_key][tax_id],
            Cache.GO_GENE_ENRICH[s_key][tax_id],
            Cache.GO_GENE[s_key][tax_id],
            Cache.GO_PARENT[s_key][tax_id],
            Cache.GENE_GO[s_key][tax_id]
            )

    @staticmethod
    def info():
        for s_key in ('LOCAL','GPDB', 'L1k'):
            print(">Databases: %s" % s_key)
            print("CATEGORY=%d" % (0 if Cache.CATEGORY[s_key] is None else len(Cache.CATEGORY[s_key])))
            print("GO_DESCRIPTION=%d" % (0 if Cache.GO_DESCRIPTION[s_key] is None else len(Cache.GO_DESCRIPTION[s_key])))
            print("GO_CATEGORY=%d" % (0 if Cache.GO_CATEGORY[s_key] is None else len(Cache.GO_CATEGORY[s_key])))
            print("PARENT=%d" % (0 if Cache.PARENT[s_key] is None else len(Cache.PARENT[s_key])))
            for tax_id in Cache.TOTAL_GENE_COUNT[s_key].keys():
                print("TAX_ID=%d (%s)" % (tax_id, ez.Cache.C_TAX_NAME.get(tax_id, "UNKNOWN")))
                print("TOTAL_GENE_COUNT=%d" % Cache.TOTAL_GENE_COUNT[s_key][tax_id])
                print("ALL_GENE=%d" % len(Cache.ALL_GENE[s_key][tax_id]))
                print("CATEGORY_COUNT=%d" % len(Cache.CATEGORY_COUNT[s_key][tax_id]))
                print("GO_GENE_ENRICH=%d" % len(Cache.GO_GENE_ENRICH[s_key][tax_id]))
                print("GO_GENE=%d" % len(Cache.GO_GENE[s_key][tax_id]))
                print("GO_PARENT=%d" % len(Cache.GO_PARENT[s_key][tax_id]))
                print("GENE_GO=%d" % len(Cache.GENE_GO[s_key][tax_id]))
            print("")

    @staticmethod
    def unload(tax_id, l_use_GPDB, l_L1k):
        if tax_id==-9606: tax_id=9606
        s_key=Cache.key(l_use_GPDB, l_L1k)
        if tax_id in Cache.TOTAL_GENE_COUNT[s_key]:
            del Cache.CATEGORY_COUNT[s_key][tax_id]
            del Cache.TOTAL_GENE_COUNT[s_key][tax_id]
            del Cache.ALL_GENE[s_key][tax_id]
            del Cache.GO_GENE_ENRICH[s_key][tax_id]
            del Cache.GO_GENE[s_key][tax_id]
            del Cache.GO_PARENT[s_key][tax_id]
            del Cache.GENE_GO[s_key][tax_id]

    @staticmethod
    def key(l_use_GPDB, l_L1k):
        if l_L1k: return "L1k"
        return 'GPDB' if l_use_GPDB else 'LOCAL'

    @staticmethod
    def load(tax_id=9606, l_use_GPDB=True, user_go=None, l_L1k=False, S_go_category=None):
        """tax_id is None, defaults to 9606, if 0, means load all supported species,
        entrez_gene is only used in local mode to accelerate Symbol retrieval
        If you only need very few GO categories, you can specify it in S_go_category, l_L1k will be ingored.
        """
        if tax_id is None:
            util.error_msg('tax_id must be an int, or 0 mans all supported species')
        tax_id=abs(tax_id)
        s_key=Cache.key(l_use_GPDB, l_L1k=l_L1k)
        if tax_id!=0 and tax_id in Cache.TOTAL_GENE_COUNT[s_key]: return
        S_tax_id=[]

        # performance optimization
        if l_L1k: return Cache.loadL1k()

        if not l_use_GPDB:
            if tax_id not in (0,9606):
                util.error_msg('Local database only supports human!')
            tax_id=9606
            if tax_id in Cache.TOTAL_GENE_COUNT[s_key]: return
            S_tax_id=[tax_id]
        else:
            mydb=db.DB('METASCAPE')
            if tax_id>0:
                S_tax_id=[tax_id]
            else:
                t=mydb.from_sql('SELECT DISTINCT tax_id FROM gid2source_id')
                S_tax_id=[x for x in t.tax_id.astype(int).tolist() if x not in Cache.TOTAL_GENE_COUNT[s_key]]
        if len(S_tax_id)==0: return
        s_tax_id=",".join(util.iarray2sarray(S_tax_id))
        print("Load %s GO database for tax_id: %s ..." % (s_key, s_tax_id))

        if l_use_GPDB:
            s_where_L1k="(term_category_id>=91 AND term_category_id<=94)" if l_L1k else "(term_category_id<91 OR term_category_id>94)"
            if S_go_category is not None:
                if type(S_go_category) is not list: S_go_category=[S_go_category]
                s_where_L1k="(term_category_id in ("+",".join([str(int(x)) for x in S_go_category])+"))"
            if Cache.CATEGORY[s_key] is None:
                t=mydb.from_sql("SELECT term_category_id,category_name FROM term_category where "+s_where_L1k)
                Cache.CATEGORY[s_key] = { k:v for k,v in zip(t.term_category_id, t.category_name) } #{t.loc[i,'term_category_id']:t.ix[i,'category_name'] for i in t.index}
                t=mydb.from_sql("SELECT t.term_id GO,term_name AS DESCRIPTION,term_category_id CATEGORY_ID FROM term t where "+s_where_L1k)
                X=t.DESCRIPTION.isnull()
                if sum(X):
                    t.loc[X, 'DESCRIPTION']=t.loc[X, 'GO']
                #if not util.is_python3():
                #    t['DESCRIPTION']=t['DESCRIPTION'].apply(lambda x: unicode(x, encoding="ISO-8859-1", errors='ignore')) # L1000 has micro Mol
                Cache.GO_DESCRIPTION[s_key]=dict(zip(t.GO, t.DESCRIPTION))
                t['CATEGORY_ID']=t['CATEGORY_ID'].astype(int)
                Cache.GO_CATEGORY[s_key]={re.sub(r'^\d+_', '', row.GO):int(row.CATEGORY_ID) for row in t.itertuples() }
            if tax_id==0:
                t=mydb.from_sql("SELECT COUNT(*) as N,tax_id FROM annotation a where a.annotation_type_id=3 AND content='protein-coding' group by tax_id")
            else:
                t=mydb.sql_in("SELECT COUNT(*) as N,tax_id FROM annotation a where a.annotation_type_id=3 AND content='protein-coding' and tax_id in (", ") group by tax_id", S_tax_id)
            Cache.TOTAL_GENE_COUNT[s_key]=dict(zip(t.tax_id, t.N))
            if tax_id==0:
                t=mydb.from_sql("SELECT term_id GO,gids GENES,tax_id FROM term2gids where "+s_where_L1k)
            else:
                t=mydb.sql_in("SELECT term_id GO,gids GENES,tax_id FROM term2gids WHERE "+s_where_L1k+" and tax_id in (", ")", S_tax_id)
            #tmp=t[t.GO.apply(lambda x: x.startswith('6'))]
            #print tmp[:4]

        else:
            DATA_FILE=setting.go['DATA_FILE']
            #TAX_ID,GeneID
            t_gene=pd.read_csv(DATA_FILE)
            t_gene=t_gene[t_gene.TAX_ID==tax_id]
            C_GENE=set(t_gene['GeneID'].astype(str).tolist())
            Cache.TOTAL_GENE_COUNT[s_key][tax_id]=len(C_GENE)

            if user_go is not None:
                if os.path.isfile(user_go):
                    if user_go.upper().endswith(".CSV"):
                        t=pd.read_csv(user_go)
                    else:
                        t=pd.read_table(user_go)
            elif os.path.isfile(Cache.DATA_DIR+"AllAnnotations.tsv"):
                t=pd.read_csv(Cache.DATA_DIR+"AllAnnotations.tsv", sep='\t')
            if t is None:
                util.error_msg('No GO Annotations available.')
            #GO  TYPE    GENES   DESCRIPTION
            S=util.unique(t.TYPE)
            Cache.CATEGORY[s_key] = dict(zip(S, S))
            Cache.GO_CATEGORY[s_key]=dict(zip(t.GO, t.TYPE))
            Cache.GO_DESCRIPTION[s_key]=dict(zip(t.GO, t.DESCRIPTION))
            t['tax_id']=tax_id

        for x in S_tax_id:
            Cache.ALL_GENE[s_key][x]=set()
            Cache.GENE_GO[s_key][x]={}
            Cache.GO_GENE[s_key][x]={}
            Cache.GO_PARENT[s_key][x]={}
            Cache.CATEGORY_COUNT[s_key][x]={}
            Cache.GO_GENE_ENRICH[s_key][x]=set()

        if l_use_GPDB:
            PARENT_DISTANCE=1
            if Cache.PARENT[s_key] is None:
                t_parent=mydb.from_sql(f"SELECT term_id FROM term2term tt WHERE tt.parent_term_id='GO:0008150' and distance={PARENT_DISTANCE}")
                S_go=["19_"+x for x in t_parent.term_id]
                t_parent=mydb.sql_in("SELECT term_id, term_name FROM term WHERE term_id in (", ")", S_go)
                C_PARENT={ k:v for k,v in zip(t_parent.term_id.tolist(), t_parent.term_name.tolist()) }
                Cache.PARENT[s_key]=C_PARENT
            else:
                S_go=list(Cache.PARENT[s_key].keys())
            tt=mydb.sql_in("SELECT term_id parent_term_id,id_count,tax_id FROM term2gids t WHERE t.term_id in (", ") and tax_id in ("+s_tax_id+") ORDER BY tax_id,id_count", S_go)
            S_go=[ x[3:] for x in S_go ]
            t2=mydb.sql_in("SELECT term_id,parent_term_id FROM term2term WHERE parent_term_id in (", ")", S_go)
            t2['term_id']=t2.term_id.apply(lambda x: "19_"+x)
            t2['parent_term_id']=t2.parent_term_id.apply(lambda x: "19_"+x)
            t2=t2.merge(tt, left_on='parent_term_id', right_on='parent_term_id')
            for tax_id2,t_v in t2.groupby('tax_id'):
                t_v.sort_values(["term_id","id_count","parent_term_id"], inplace=True)
                t_v.drop_duplicates(["term_id"], inplace=True)
                Cache.GO_PARENT[s_key][tax_id2]={ k:v for k,v in zip(t_v.term_id.tolist(), t_v.parent_term_id.tolist()) }
        #sw=util.StopWatch("AAAAAAA")
        for tax_id2,t_v in t.groupby('tax_id'):
            #t_v=t_v.copy()
            GENE_GO={}
            GO_GENE={}
            GO_GENE_ENRICH=set()
            ALL_GENE=set()
            CATEGORY_COUNT={}
            s_cat=0
            S_genes=[ (row.GO, row.GENES.split(",")) for row in t_v.itertuples() ]
            if not l_use_GPDB:
                S_genes=[ (x, [y for y in Y if (y in C_GENE)]) for x,Y in S_genes ]
            GO_GENE={x: set(Y) for x,Y in S_genes if (len(Y)>0 and len(Y)<=Cache.N_TRIVIAL) }
            GO_GENE_ENRICH=set(GO_GENE.keys())
            if l_use_GPDB:
                for x in GO_GENE_ENRICH:
                    if re.sub(r'^\d+_', '', x) not in Cache.GO_CATEGORY[s_key]:
                        print(">>>>>>>>>>>>>>>>>>>", x, s_key, re.sub(r'^\d+_', '', x))
                        exit()
                S_cat=[ Cache.GO_CATEGORY[s_key][re.sub(r'^\d+_','', x)] for x in GO_GENE_ENRICH ]
            else:
                S_cat=[ Cache.GO_CATEGORY[s_key][x] for x in GO_GENE_ENRICH ]
            CATEGORY_COUNT=util.unique_count(S_cat)
            # reduce is slower
            #ALL_GENE=reduce(lambda a,b : a|b, GO_GENE.values())
            ALL_GENE=set([x for Y in GO_GENE.values() for x in Y])

            #for row in t_v.itertuples():
            ##for i in t_v.index:
            #    s_go=row.GO #t_v.ix[i, 'GO']
            #    S_genes=row.GENES.split(",") #t_v.ix[i, 'GENES'].split(",")
            #    if not l_use_GPDB:
            #        ### warning, gene ids not recognized are treated as tax ID 0!!!
            #        S_genes=[s for s in S_genes if s in C_GENE]
            #    if len(S_genes)==0: continue
            #    if len(S_genes)<=Cache.N_TRIVIAL:
            #        GO_GENE_ENRICH.add(s_go)
            #        if l_use_GPDB:
            #            s_cat=Cache.GO_CATEGORY[s_key].get(re.sub(r'^\d+_','',s_go), 0)
            #        CATEGORY_COUNT[s_cat]=CATEGORY_COUNT.get(s_cat, 0)+1
            #    GO_GENE[s_go]=set(S_genes)
            #    ALL_GENE.update(GO_GENE[s_go])
            #sw.check("TTTTTTTTT "+str(tax_id))
            for k,v in GO_GENE.items():
                for s_gene in v:
                    if s_gene not in GENE_GO:
                        GENE_GO[s_gene]={k}
                    else:
                        GENE_GO[s_gene].add(k)
            Cache.ALL_GENE[s_key][tax_id2]=ALL_GENE
            Cache.GENE_GO[s_key][tax_id2]=GENE_GO
            Cache.TOTAL_GENE_COUNT[s_key][tax_id2]=max(Cache.TOTAL_GENE_COUNT[s_key][tax_id2], len(GENE_GO))
            Cache.CATEGORY_COUNT[s_key][tax_id2]=CATEGORY_COUNT
            Cache.GO_GENE[s_key][tax_id2]=GO_GENE
            Cache.GO_GENE_ENRICH[s_key][tax_id2]=GO_GENE_ENRICH

        if l_L1k:
            s_path=setting.go['L1000_PATH']
            S_gene=util.read_list(s_path+'/L1kAllGenes.txt')
            Cache.ALL_GENE[s_key][tax_id]=set(S_gene)
            Cache.TOTAL_GENE_COUNT[s_key][tax_id]=len(S_gene)

    @staticmethod
    def loadL1k():
        """Load L1000 terms"""
        sw=util.StopWatch()
        print("Loading L1k terms ...")
        tax_id=9606
        s_key="L1k"
        s_path=setting.go['L1000_PATH']
        S_gene=util.read_list(s_path+"/L1kAllGenes.txt")
        Cache.TOTAL_GENE_COUNT[s_key][tax_id]=len(S_gene)
        t1=pd.read_csv(s_path+"/Term2Gid_L1000_PhaseI.csv")
        t2=pd.read_csv(s_path+"/Term2Gid_L1000_PhaseII.csv")
        t=pd.concat([t1, t2], ignore_index=True)
        t['category_id']=t['category_id'].astype(int)
        #t=t1[:600000]
        Cache.ALL_GENE[s_key][tax_id]=set(S_gene)
        Cache.GENE_GO[s_key][tax_id]={}
        Cache.GO_GENE[s_key][tax_id]={}
        Cache.CATEGORY_COUNT[s_key][tax_id]={}
        Cache.GO_GENE_ENRICH[s_key][tax_id]=set()
        sw.check('Loaded CSV')
        s_type='BP'
        Cache.CATEGORY[s_key]={91:'L1000 shRNA',92:'L1000 Cpd',93:'L1000 cDNA',94:'L1000 Lgnd'}
        Cache.CATEGORY_COUNT[s_key][tax_id]={}
        Cache.GO_CATEGORY[s_key]={}
        Cache.GO_DESCRIPTION[s_key]={}
        CATEGORY_COUNT={}
        t['GO_ID']=[ str(row.category_id)+"_"+row.term_id for row in t.itertuples() ]
        Cache.GO_CATEGORY[s_key]=dict(zip(t.GO_ID, t.category_id))
        Cache.GO_DESCRIPTION[s_key]=dict(zip(t.GO_ID, t.term_name))
        CATEGORY_COUNT=util.unique_count(t.category_id)
        Cache.GO_GENE[s_key][tax_id]={ row.GO_ID:set(row.gids.split(",")) for row in t.itertuples() }
        for s_go_id,S_gene in Cache.GO_GENE[s_key].items():
            for s_gene in S_gene:
                if s_gene not in Cache.GENE_GO[s_key][tax_id]:
                    Cache.GENE_GO[s_key][tax_id][s_gene]={s_go_id}
                else:
                    Cache.GENE_GO[s_key][tax_id][s_gene].add(s_go_id)

        ##for i in t.index:
        #for row in t.itertuples():
        #    i_cat=row.category_id #t.ix[i,'category_id']
        #    s_go=row.term_id #t.ix[i, 'term_id']
        #    s_des=row.term_name #t.ix[i,'term_name']
        #    S_gene=row.gids.split(",") #t.ix[i,'gids'].split(',')
        #    s_go_id=str(i_cat)+"_"+s_go
        #    Cache.GO_CATEGORY[s_key][s_go_id]=i_cat
        #    Cache.GO_DESCRIPTION[s_key][s_go_id]=s_des
        #    CATEGORY_COUNT[i_cat]=CATEGORY_COUNT.get(i_cat, 0)+1
        #    #if len(S_genes)<=Cache.N_TRIVIAL:
        #    Cache.GO_GENE[s_key][tax_id][s_go_id]=set(S_gene)
        #    for s_gene in S_gene:
        #        if s_gene not in Cache.GENE_GO[s_key][tax_id]:
        #            Cache.GENE_GO[s_key][tax_id][s_gene]={s_go_id}
        #        else:
        #            Cache.GENE_GO[s_key][tax_id][s_gene].add(s_go_id)
        Cache.CATEGORY_COUNT[s_key][tax_id]=CATEGORY_COUNT
        Cache.GO_GENE_ENRICH[s_key][tax_id]=set(Cache.GO_GENE[s_key][tax_id].keys())
        Cache.GO_PARENT[s_key][tax_id]={}
        sw.check("Done L1k")
        print('Done L1k loading')

class GO(object):

    GO_ROOT={"BP":"GO:0008150","MF":"GO:0003674","CC":"GO:0005575"}
    #GO_DESCRIPTION={"BP":"Biological Process","MF":"Molecular Function","CC":"Cellular Component"}
    # used for l_use_GPDB GO Categories
    # GO: BP:19, MF:21, CC:20
    # GeneGo: Pathway Map:27, Go Processes:31, Drug Target:29, Disease:28
    # Custom: custom gene sets
    # KEGG: pathway:24, MF: 26, CC:25
    # MSigDB: Pathway:11, BioCarta:15, Hallmark:23, Reactome:6, Onc Set:4, Immu Set:3, Chem/Genetics 13

    def get_source_id(self,term_id):
        return re.sub(r'^\d*_','',term_id)

    def get_category_id(self,term_id):
        if re.search(r'^\d+_', term_id):
            t = term_id.split('_')
            return int(t[0])
        else:
            return self.GO_CATEGORY.get(term_id, None)
        #the line below breaks on hsa_M00055
        #return int(t[0]) if len(t)>1 else None

    def get_categories(self):
        if self.GPDB:
            mydb=db.DB('METASCAPE')
            t=mydb.from_sql("SELECT term_category_id,category_name,ds data_source FROM term_category")
            return t
        return None

    @staticmethod
    def qc_categories(S_cat):
        S=set([34,35,36,91,92,93,94])
        S_qc=[x for x in S_cat if x in S]
        S_go=[x for x in S_cat if x not in S]
        return (S_go, S_qc)

    @staticmethod
    def get_go_categories():
        mydb=db.DB('METASCAPE')
        t=mydb.from_sql("select c.term_category_id,c.category_name,c.category_group,c.category_group_name_membership from term_category c where c.used_in_enrichment='Y' order by category_group,term_category_id")
        return t

    def __init__(self, tax_id=None, l_use_GPDB=True, entrez=None, l_L1k=False, r_random=0, S_go_category=None):
        """If S_go_category is provided, only load those"""
        self.eg=entrez
        self.GPDB=l_use_GPDB
        self.tax_id=tax_id
        self.L1k=l_L1k
        (self.CATEGORY, self.GO_DESCRIPTION, self.GO_CATEGORY, self.PARENT, self.TOTAL_GENE_COUNT, self.ALL_GENE, self.CATEGORY_COUNT, self.GO_GENE_ENRICH, self.GO_GENE, self.GO_PARENT, self.GENE_GO)=Cache.get(tax_id=tax_id, l_use_GPDB=l_use_GPDB, l_L1k=l_L1k, S_go_category=S_go_category)
        if r_random>0:
            self.GO_GENE_ENRICH=random.sample(self.GO_GENE_ENRICH, int(len(self.GO_GENE_ENRICH)*r_random))

    def is_L1000(self):
        return self.L1k

    def go_description(self, s_go):
        if s_go in self.GO_DESCRIPTION:
            return self.GO_DESCRIPTION[s_go]
        return s_go

    def filter_genes_by_go(self, s_go, S_genes):
        return [ x for x in S_genes if s_go in self.GENE_GO.get(x, []) ]
        #if s_go in self.GO_GENE:
        #    return list(set(S_genes).intersection(self.GO_GENE[s_go]))
        #else:
        #    return []

    def go_size(self, s_go):
        #return len([True for v in self.GENE_GO.values() if s_go in v ])
        return len(self.GO_GENE.get(s_go, []))

    def analysis_go(self, s_go, S_hit, N_total=0, SRC_GENE=None, min_overlap=3, p_cutoff=0.01):
        c={'GO':s_go, '#TotalGeneInLibrary':N_total, '#GeneInGO':0, '#GeneInHitList':0, '#GeneInGOAndHitList':0, 'LogP':0.0, 'Enrichment':0}
        #if SRC_GENE is not None:
        #    print "SRC_GENE: "+str(len(SRC_GENE))
        S_gene=self.GO_GENE[s_go]
        if len(S_gene)>=Cache.N_TRIVIAL:
            return None
        if not N_total:
            N_total=len(self.GENE_GO.keys()) #len(self.ALL_GENE), only count genes that has GO annotation
        if SRC_GENE is not None:
            S_gene=S_gene.intersection(SRC_GENE)
            S_hit=set(S_hit).intersection(SRC_GENE)
        else:
            S_hit=set(S_hit)
        S_hit=self.ALL_GENE.intersection(S_hit)
        c['#GeneInGO']=len(S_gene)
        c['#GeneInHitList']=len(S_hit)
        if c['#GeneInGO']<min_overlap or c['#GeneInHitList']<min_overlap:
            return None
        S_both=S_gene.intersection(S_hit)
        c['#GeneInGOAndHitList']=len(S_both)
        if c['#GeneInGOAndHitList']<min_overlap:
            return None
        c['%InGO']=c['#GeneInGOAndHitList']*100.0/c['#GeneInHitList']
        q=min(max(c['%InGO']/100, 1.0/c['#GeneInHitList']), 1-1.0/c['#GeneInHitList'])
        c['STDV %InGO']=np.sqrt(q*(1-q)/c['#GeneInHitList'])*100
        c['Enrichment']=c['%InGO']/100.0*N_total/c['#GeneInGO']
        S=[int(x) for x in S_both]
        S.sort()
        c['GeneID']='|'.join([str(x) for x in S])
        c['LogP']=np.log10(max(
            stats.hyper(c['#GeneInGOAndHitList'], N_total, c['#GeneInGO'], c['#GeneInHitList']), 1e-100))
        # GeneGo Z-score definition
        c['Z-score']=stats.ZScore_GeneGo(c['#GeneInGOAndHitList'], N_total, c['#GeneInGO'], c['#GeneInHitList'])
        if c['LogP']>np.log10(p_cutoff): return None
        return c

    def analysis_go_RSA(self, s_go, S_hit, S_score, N_total=0, SRC_GENE=None, min_overlap=3, p_cutoff=0.01, l_keep_most=True):
        """Input is a list of hits with scores, the smaller the score, the better!
        We then iteratively try different score cutoffs and use the cutoff that produce the best P-value (Bonferroni corrected)"""
        c={'GO':s_go, '#TotalGeneInLibrary':N_total, '#GeneInGO':0, '#GeneInHitList':0, '#GeneInGOAndHitList':0, 'Cutoff':None, '#HitRemain':0, '#HitInGORemain':0, 'LogP':0.0, 'Enrichment':0}
        #if SRC_GENE is not None:
        #    print "SRC_GENE: "+str(len(SRC_GENE))
        S_gene=self.GO_GENE[s_go]
        if len(S_gene)>=Cache.N_TRIVIAL:
            return None
        if not N_total:
            N_total=len(self.GENE_GO) #len(self.ALL_GENE), only count genes that has GO annotation
        t_hit=pd.DataFrame(data={'Hit':S_hit, 'Score':S_score})
        if SRC_GENE is not None:
            S_gene=S_gene.intersection(SRC_GENE)
            t_hit=t_hit[ t_hit.Hit.apply(lambda x: x in SRC_GENE) ]
            S_hit=set(t_hit.Hit)
        else:
            S_hit=set(S_hit)
        t_hit.sort_values('Score', inplace=True)
        c['#GeneInGO']=len(S_gene)
        c['#GeneInHitList']=len(S_hit)
        if c['#GeneInGO']<min_overlap or c['#GeneInHitList']<min_overlap:
            return None
        S_both=S_gene.intersection(S_hit)
        c['#GeneInGOAndHitList']=len(S_both)
        if c['#GeneInGOAndHitList']<min_overlap:
            return None

        I_index=np.arange(len(t_hit))[t_hit.Hit.apply(lambda x: x in S_gene).values]
        I_rank=stats.RSA_rank(t_hit.Score.values, I_index)
        rslt=stats.RSA_score(I_rank, N_total, l_BonferroniCorrection=True, l_keep_most=l_keep_most, p_cutoff=p_cutoff)
        c['#HitInGORemain']=rslt["cutoff"]+1
        if c['#HitInGORemain']<min_overlap: return None
        c['#HitRemain']=I_rank[rslt["cutoff"]]
        c['Cutoff']=t_hit.Score.values[rslt["cutoff"]]
        c['%InGO']=c['#HitInGORemain']*100.0/c['#HitRemain']
        q=min(max(c['%InGO']/100, 1.0/c['#HitRemain']), 1-1.0/c['#HitRemain'])
        c['STDV %InGO']=np.sqrt(q*(1-q)/c['#HitRemain'])*100
        c['Enrichment']=c['%InGO']/100.0*N_total/c['#GeneInGO']
        S=[int(x) for x in S_both]
        S.sort()
        c['GeneID_All']='|'.join([str(x) for x in S])
        S=[int(x) for x in list(t_hit.Hit[: rslt["cutoff"]+1])]
        S.sort()
        c['GeneID']='|'.join([str(x) for x in S])
        c['LogP']=rslt['logP']
        return c

    def go_count(self, S_hit, S_go=None):
        """return a dict of GO and the number of genes appear in each GO
            if S_go is provided, only counts for those go terms"""
        c={}
        if S_go is not None: S_go=set(S_go)
        for x in S_hit:
            Y=self.GENE_GO.get(x, [])
            if S_go is not None: Y = set(Y).intersection(S_go)
            for y in Y:
                c[y]=c.get(y,0)+1
        return c

    def gene_count(self, S_go, S_gene=None):
        """return a dict of Gene and the number of GOs appear for each gene
           if S_gene is provided, only counts for those genes
        """
        c={}
        if S_gene is not None: S_gene=set(S_gene)
        for x in S_go:
            Y=self.GO_GENE.get(x, [])
            if S_gene is not None: Y = set(Y).intersection(S_gene)
            for y in self.GO_GENE.get(x, []):
                c[y]=c.get(y,0)+1
        return c

    def membership_count(self, S_go, S_gene):
        """return a dict of GO and the set of genes fall into those GO"""
        return self.go_count(S_gene, S_go)
        #c=self.go_count(S_gene)
        #if type(S_go)!=set:
        #    S_go=set(S_go)
        #c={ k:v for k,v in c.items() if k in S_go }
        #return c

    def membership_go_genes(self, S_go, S_gene):
        c={}
        S_gene=set(S_gene)
        for x in S_go:
            S=set(self.GO_GENE.get(x, [])).intersection(S_gene)
            if len(S): c[x]=list(S)
        return c
        #gene_go = { x:self.GENE_GO.get(x, []) for x in S_gene}
        #c={}
        #for k,v in gene_go.items():
        #    for g in v:
        #        if g not in S_go:
        #            continue
        #        if g not in c:
        #            c[g] = []
        #        c[g].append(k)
        #return c

    def analysis(self, S_hit, S_score=None, S_go=None, SRC_GENE=None, min_overlap=3, min_enrichment=0, p_cutoff=0.01, n_CPU=0, l_rsa_keep_most=True, S_go_category=None, l_background_by_ontology=False):
        """If Score is None, just run GO enrichment test, if Score is provided, the smaller score represents hits are more reliable.
        An iterative enrichment test, i.e., RSA routine is applied, see analysis_go_RSA is applied.
            S_go_category: a set of categories to use, useful for gene prioritization project.
                By default, we recommend [31,19,11,15,27,24]
        Special usage: both S_hit and S_go are dict, {'a':S_hit1, 'b':S_hit2}, {'a':S_go1, 'b':S_go2}
        Which is a short cut to run analysis(S_hit1, S_go1) and analysis(S_hit2, S_go2)
        i.e., we analysis multiple hit lists, each has its own S_go list.
        We pool them together, so that we can use CPUs more effectively, e.g, len(S_go1)==2 and len(S_go2)==10, we can run them in 12 CPUs once, instead of in two batches.

        """

        def go_filtered(S_go, S_go_category):
            return [x for x in S_go if self.get_category_id(x) in S_go_category]

        S_all_go_filtered=[]
        def all_go_filtered(S_go_category):
            if len(S_all_go_filtered)==0:
                S_go=self.GO_GENE_ENRICH
                if S_go_category is None:
                    S_all_go_filtered.append(list(S_go))
                else:
                    S_all_go_filtered.append(go_filtered(S_go, S_go_category))
            return S_all_go_filtered[0]

        N_go=0
        if S_go_category is not None and len(S_go_category)>0:
            # hard code for now, to be fixed later
            if type(S_go_category) in (int, str):
                S_go_category=[S_go_category]
            S_go_category={int(x) for x in S_go_category if self.CATEGORY_COUNT.get(x,0)>0 }
            for x in S_go_category:
                N_go+=self.CATEGORY_COUNT[x]
        else:
            N_go=sum(self.CATEGORY_COUNT.values())

        l_multi_list=type(S_hit) is dict
        if S_go is None:
            if l_multi_list:
                S_go={}
                for k in S_hit.keys():
                    S_go[k]=all_go_filtered(S_go_category)
            else:
                S_go=all_go_filtered(S_go_category)
        else:
            if l_multi_list:
                for k in S_hit.keys():
                    if S_go.get(k, None) is None:
                        S_go[k]=all_go_filtered(S_go_category)
                    else:
                        S_go[k]=go_filtered(S_go[k], S_go_category)
            else:
                S_go=go_filtered(S_go, S_go_category)

        if SRC_GENE is not None:
            if type(SRC_GENE) is list:
                SRC_GENE=set(SRC_GENE)
            SRC_GENE=self.ALL_GENE.intersection(SRC_GENE) # remove genes from background, if it is not in self.ALL_GENE
            N_total=len(SRC_GENE) #self.ALL_GENE.intersection(SRC_GENE))
        elif l_background_by_ontology:
            # GeneGo uses this
            if l_multi_list:
                X=set()
                for x in S_go.values():
                    X.add(set(x))
                src_genes=self.gene_count(list(X))
            else:
                src_genes=self.gene_count(S_go)
            N_total=len(src_genes)
            SRC_GENE=set(src_genes.keys())
        else:
            if self.is_L1000():
                N_total=len(self.ALL_GENE)
            else:
                N_total=len(self.GENE_GO) #len(self.ALL_GENE), only count genes that has GO annotation
            #N_total=len(self.ALL_GENE)
        # prefiltering uninteresting GO terms
        # already converted to multiple hit list situation
        #sw=util.StopWatch()
        L=[] # list of (S_hit, s_go)

        def spread_input(S_hit, S_go, key):
            #S_hit, S_go, key=(X[0], X[1], X[2])
            # may not worth it
            #c_cnt=self.go_count(S_hit, S_go)
            #S_go=[s_go for s_go in S_go if c_cnt.get(s_go,0)>=min_overlap ]
            # minimum size
            MIN_BATCH=2000
            S_go2=util.split(S_go, chunk_size=MIN_BATCH)
            return [(key, S_hit, x) for x in S_go2]

        #sw.check('To spreadout')
        if l_multi_list:
            #mp=parallel.MP()
            #m=1 if len(S_hit)<=3 else n_CPU
            #mp.start(f=spread_input, n_CPU=m)
            #L=[(X, S_go[k], k) for k,X in S_hit.items() if len(X)>=min_overlap]
            #out=mp.map(L)
            #L=[y for X in out for y in X]
            L=[]
            for k,X in S_hit.items():
                if len(X)<min_overlap: continue
                L.extend(spread_input(X, S_go[k], k))
                random.shuffle(L)
        else:
            if len(S_hit)>=min_overlap:
                L=spread_input(S_hit, S_go, 'Default')

        if self.eg is None:
            self.eg=ez.EntrezGene(tax_id=self.tax_id, l_use_GPDB=self.GPDB)
        if n_CPU==0: n_CPU=1
        #print ">>>>>>>>>>>>>>", len(L)
        S_chunk=util.split(L, n_chunk=n_CPU)
        #sw.check('Spreadout tasks: %d' % len(L))

        def analyze(L):
            """L is a list of [[s_name, S_hit, s_go]], s_go can also be a list"""
            rslt=[]
            #p=util.Progress(len(L))
            i=0
            import multiprocessing
            s_pid=str(multiprocessing.current_process().pid)
            for s_name, S_hit, S_go in L:
                i+=1
                #if (i % 50000): p.check(i, s_pid)
                if type(S_go) is str: S_go=[S_go]
                for s_go in S_go:
                    if s_go not in self.GO_GENE: continue
                    if S_score is None:
                        c=self.analysis_go(s_go, S_hit, N_total, SRC_GENE=SRC_GENE, min_overlap=min_overlap, p_cutoff=p_cutoff)
                    else:
                        c=self.analysis_go_RSA(s_go, S_hit, S_score, N_total, SRC_GENE=SRC_GENE, min_overlap=min_overlap, p_cutoff=p_cutoff, l_keep_most=l_rsa_keep_most)
                    if c is None:
                        continue
                    c['Name']=s_name
                    if min_enrichment>0 and c['Enrichment']<min_enrichment: continue
                    if p_cutoff<1 and 10**c['LogP']>p_cutoff: continue
                    c['Description']= self.go_description(s_go)
                    S_gene=c['GeneID'].split('|')
                    S_symbol=[self.eg.C_GENENEW[x] if x in self.eg.C_GENENEW else x for x in S_gene]
                    c['Hits']='|'.join(S_symbol)
                    if 'GeneID_All' in c:
                        S_gene=c['GeneID_All'].split('|')
                        S_symbol=[self.eg.C_GENENEW[x] if x in self.eg.C_GENENEW else x for x in S_gene]
                        c['Hits_All']='|'.join(S_symbol)
                    if self.GPDB:
                        c['CategoryID'] = self.get_category_id(c['GO'])
                        c['Category'] = self.CATEGORY.get(self.get_category_id(c['GO']))
                        c['GO'] = self.get_source_id(c['GO'])
                    rslt.append(c)
            return rslt
        out=parallel.parmap(analyze, S_chunk, n_CPU=n_CPU)
        #if n_CPU>1:
        #    mp=parallel.MP()
        #    mp.start(f=analyze, n_CPU=n_CPU)
        #    out=mp.map(S_chunk)
        #else:
        #    out=[analyze(x) for x in S_chunk]

        #mp.start(n_CPU=n_CPU)
        #sw.check('P-value Calculation')
        #sw.check('P-value Calculation Done')
        rslt=[]
        for x in out:
            if len(x): rslt.extend(x)

        if len(rslt):
            #sw.check('Length: %d' % len(rslt))
            t=pd.DataFrame(rslt)
            #sw.check('Table DONE')
            S=[str(a)+"_"+b for a,b in zip(t.CategoryID.tolist(), t.GO.tolist()) ]
            t['PARENT_GO']=[ '' if (x not in self.GO_PARENT) else self.GO_PARENT[x]+' '+self.PARENT[self.GO_PARENT[x]] for x in S ]
            if S_score is None:
                t=t.sort_values(['LogP','Enrichment','#GeneInGOAndHitList'], ascending=[True,False,False])
                cols = ['Name','GO','Description','PARENT_GO','LogP','Enrichment','Z-score','#TotalGeneInLibrary',
                        '#GeneInGO','#GeneInHitList','#GeneInGOAndHitList','%InGO','STDV %InGO','GeneID','Hits']
            else:
                t=t.sort_values(['LogP','Enrichment','#HitInGORemain','#GeneInGOAndHitList'], ascending=[True,False,False,False])
                cols = ['Name','GO','Description','PARENT_GO','LogP','Enrichment','Z-score','#TotalGeneInLibrary',
                    '#GeneInGO','#HitRemain','#HitInGORemain','Cutoff','#GeneInHitList','#GeneInGOAndHitList','%InGO','STDV %InGO','GeneID','Hits','GeneID_All','Hits_All']
            if self.GPDB:
                #cols.insert(1,'field1')
                cols.insert(1,'CategoryID')
                cols.insert(1,'Category')
            #sw.check('sorted DONE')
            t=t.reindex(columns=cols)
            # FDR
            #print ">>> N_go: ", N_go
            #sw.check('reindex DONE')
            t['Log(q-value)']=np.log10(np.clip(stats.adjust_p(np.power(10, t.LogP.values), N=N_go), 1e-100, 1.0))
            #sw.check('q-value DONE')
            if not l_multi_list:
                t.drop('Name', axis=1, inplace=True)
            return t
        else:
            return None

    def key_terms(self, t_old, t_new, t_union, S_old=None, t_over=None):
        """Look for terms presented in both list and the p-valeu is even better in the combined list
        S_old is the set of genes in Old set, if None, set to all genes in t_old table
        This method is to be used by analyze_key_terms."""

        if t_old is None or t_new is None or t_union is None:
            return None
        print("Old: %d, New: %d, Union: %d" % (len(t_old), len(t_new), len(t_union)))
        if S_old is None:
            S_old=set([y for x in [ t_old.loc[i, 'GeneID'].split("|") for i in t_old.index ] for y in x])
        elif type(S_old) is list:
            S_old=set(S_old)
        t_old=t_old[["GO","LogP"]]
        t_old.rename2({"LogP":"LogP_Hit"})
        t_new=t_new[["GO","LogP"]]
        t_new.rename2({"LogP":"LogP_OverConnect"})
        #t_old.rename2({"LogP":"LogP_Union"})
        t=pd.merge(t_union, t_old, on="GO")
        if len(t)==0: return None
        t=pd.merge(t, t_new, on="GO")
        if len(t)==0: return None
        t.sort_values(["LogP", "Enrichment"], ascending=[True, False], inplace=True)
        S_id1=[]
        S_id2=[]
        S_name1=[]
        S_name2=[]
        S_count=[]
        for i in t.index:
            S_id=t.loc[i, 'GeneID'].split("|")
            S_name=t.loc[i, "Hits"].split("|")
            Idx={ j for j,id in enumerate(S_id) if id in S_old }
            S_id1.append("|".join([ x for j,x in enumerate(S_id) if j in Idx ]))
            S_name1.append("|".join([ x for j,x in enumerate(S_name) if j in Idx ]))
            S_id2.append("|".join([ x for j,x in enumerate(S_id) if j not in Idx ]))
            S_name2.append("|".join([ x for j,x in enumerate(S_name) if j not in Idx ]))
            S_count.append("%d|%d" % (len(Idx), t.loc[i, "#GeneInGOAndHitList"]-len(Idx)))
        # separate gene lists into those in Hits and those not in Hits
        t["GeneID_InHit"]=S_id1
        t["Gene_InHit"]=S_name1
        t["GeneID_InOverConnect"]=S_id2
        t["Gene_InOverConnect"]=S_name2
        t["GeneCount_InHit_OverConnect"]=S_count
        S_over=list(t_over.Node) if t_over is not None else None
        other={'OldGO':t_old, 'NewGO':t_new, 'OverConnect':t_over, 'OldGeneList':list(S_old), 'OverGeneList':S_over}
        return (t, other)

    def analyze_key_terms(self, S_hit, ppi=None, p_cutoff=0.01, min_enrich_expand=0, min_links=2, n_CPU=0, S_intgo_category=None, l_background_by_ontology=True):
        """min_enrich_expand is the minimum enrichment factor for construct overconnected list."""
        if ppi is None:
            util.error_msg('Network must be provided: ppi!')
        t_old=self.analysis(S_hit, p_cutoff=p_cutoff, n_CPU=n_CPU, S_go_category=S_go_category, l_background_by_ontology=l_background_by_ontology)
        t_over=ppi.overconnected(ppi, S_hit, min_links=min_links, p_cutoff=p_cutoff, min_enrichment=min_enrich_expand)
        S_new=list(t_over.Node)
        t_new=self.analysis(S_new, p_cutoff=p_cutoff, n_CPU=n_CPU, S_go_category=S_go_category, l_background_by_ontology=l_background_by_ontology)
        S_union=list(set(S_hit+S_new))
        t_union=self.analysis(S_union, p_cutoff=p_cutoff, n_CPU=n_CPU, S_go_category=S_go_category, l_background_by_ontology=l_background_by_ontology)
        return self.key_terms(t_old, t_new, t_union, t_over)

    @staticmethod
    def extract_parent_go(t_go):
        if t_go is None: return None
        if 'PARENT_GO' not in t_go.header():
            util.error_msg('PARENT_GO column missing!')
        t=t_go[t_go.PARENT_GO!='']
        S_sort=['Name','PARENT_GO','LogP'] if 'Name' in t_go.header() else ['PARENT_GO','LogP']
        t=t.sort_values(S_sort)
        # for each parent_go, each source name, only keep the best LogP
        t=t.drop_duplicates(S_sort[:-1]).copy()
        t.rename2({'GO':'CHILD_GO', 'Description':'CHILD_Description'})
        t['GO']=[ x[3:13] for x in t.PARENT_GO ]
        t['Description']=[ x[14:] for x in t.PARENT_GO ]
        t.drop('PARENT_GO', axis=1, inplace=True)
        return t

class GO_Cluster:

    def __init__(self, S_all_gene, c_go=None, n_CPU=20):
        """S_all_gene, list[str] list of all genes went through GO analysis
        c_go, dict, keys are GO IDs, and values are list of genes in that GO
            c_go can be the table as output of the GO analysis, contains GO, GeneID or GeneID_All"""
        self.data=pd.DataFrame(data={"Gene":util.unique(S_all_gene)})
        self.data.set_index("Gene", drop=True, inplace=True)
        self.t_go=None
        self.n_CPU=n_CPU if n_CPU>0 else 1
        sw=util.StopWatch("GO_Cluster")
        if type(c_go) is not dict: # dataframe
            self.t_go=c_go
            s_col='GeneID_All' if 'GeneID_All' in self.t_go.header() else 'GeneID'
            #c_go={ k:v for k,v in zip(self.t_go.GO, self.t_go[s_col].split('|')) } #{ self.t_go.GO[i]:self.t_go[s_col][i].split('|') for i in self.t_go.index }
            c_go={ k:v for k,v in zip(self.t_go.GO, self.t_go[s_col].str.split('|')) } # YaZ

        L=util.split(list(c_go.keys()), self.n_CPU)
        def f(x):
            t, S_go = x[0], x[1]
            for k in S_go:
                for g in c_go[k]:
                    t.loc[g, k]=1
            return t

        #mp=parallel.MP()
        #mp.start(f, n_CPU=len(L))
        L=[( self.data.copy(), x) for x in L]
        #out=mp.map(L)
        out=parallel.parmap(f, L, n_CPU=len(L))
        self.data=pd.concat(out, axis=1)
        #for k,v in c_go.items():
        sw.check('Done membership...')
        #    for g in v:
        #        self.data.ix[g, k]=1
        self.data.fillna(value=0, inplace=True)
        self.S_GO=self.data.header()
        self.DM=None  # distance matrix
        self.similarity=None

    def cluster(self, similarity=0.3, l_go_selective=False):
        """If self.t_go is None, it cluster the self.data - membership matrix and return a dict
        return dict: key: GO names, value: group id
        If self.t_go is a table, it returns a table with additional columns:
            Group_ID,FirstInGroupByEnrichment,FirstInGroupByLogP,URL
        """
        sw=util.StopWatch("GO_Cluster::cluster")
        #K=stats.kappa_stat(self.data.values)
        #T_edge=pd.DataFrame({'Gene_A':[],'Gene_B':[],'TYPE':[],'SCORE':[]})
        M=self.data.values
        n,m=M.shape
        print("Matrix size: %d genes x %d GOs" % (n, m))
        S_go=self.data.header()
        out=self.t_go
        if m==0:
            self.DM=np.zeros(0)
            if out is None:
                return {}
            return None
        if m==1:
            # only 1, no need to cluster
            if out is None:
                return {S_go[0]:1}
            out['GROUP_ID']=1
            out['FirstInGroupByEnrichment']=1
            out['FirstInGroupByLogP']=1
            if l_go_selective: out['FirstInGroupByGini']=1
            self.DM=np.zeros(0)
            return out
            #T_edge=pd.DataFrame({'Gene_A':[S_go[0]],'Gene_B':[S_go[0]],'TYPE':['Direct'],'SCORE':[1.0]})
        if self.DM is None:
            self.DM=stats.kappa_stat(self.data.values, n_CPU=self.n_CPU)
        sw.check("Kappa done ...")
        #import ms.msobject as mo
        #mo.MSObject.dump_object(self.DM, s_name='untitled', s_cache_dir=".")
        import scipy.cluster.hierarchy as clst
        import fastcluster
        Z=fastcluster.linkage(1.0-self.DM, method='average')
        S=clst.fcluster(Z, 1-similarity, criterion='distance')
        c_grp={  x:S[i] for i,x in enumerate(S_go) }
        if out is None:
            return c_grp
        out['GROUP_ID']=out.GO.apply(lambda x: c_grp[x])
        self.similarity=similarity
        if l_go_selective:
            out.sort_values(['GROUP_ID','GiniIndex','LogP','Enrichment'], ascending=[True,False,True,False], inplace=True)
            out['FirstInGroupByGini']=0
        else:
            out.sort_values(['GROUP_ID','LogP','Enrichment'], ascending=[True,True,False], inplace=True)
        out['FirstInGroupByEnrichment']=0
        out['FirstInGroupByLogP']=0
        iB=iE=0
        n=len(out)
        out.index=list(range(n))
        for i in range(1,n+1):
            if i>=n or out.loc[i,'GROUP_ID']!=out.loc[i-1,'GROUP_ID']:
                iE=i-1
                out.loc[iB:iE,'BestLogPInGroup']=out.loc[iB:iE,'LogP'].min()
                out.loc[iB:iE,'BestEnrichmentInGroup']=out.loc[iB:iE,'Enrichment'].max()
                idx=out.loc[iB:iE,'LogP'].values.argmin()+iB
                out.loc[idx, 'FirstInGroupByLogP']=1
                out.loc[iB, 'FirstInGroupByEnrichment']=1
                if l_go_selective:
                    out.loc[iB:iE,'BestGiniInGroup']=out.loc[iB:iE,'GiniIndex'].max()
                    idx=out.loc[iB:iE,'GiniIndex'].values.argmax()+iB
                    out.loc[idx, 'FirstInGroupByGini']=1
                iB=i
        if l_go_selective:
            out.sort_values(['BestGiniInGroup','BestLogPInGroup','GROUP_ID','FirstInGroupByGini','GiniIndex','LogP','Enrichment'], ascending=[False,True,True,False,False,True,False], inplace=True)
            out.index=list(range(n))
#            out.to_csv('t0.csv', index=False)
            # iteratively pick unique patterns
            S_pattern=util.unique2(out._PATTERN_) # unique but preserve order
            n_pattern=len(S_pattern)
            iB=iE=0
            i_pattern={k:(i+1) for i,k in enumerate(S_pattern)}
            c_pattern={k:0 for k in S_pattern}
            out['NEW_GROUP_ID']=0
            for i in range(1,n+1):
                if i>=n or out.loc[i,'GROUP_ID']!=out.loc[i-1,'GROUP_ID']:
                    iE=i-1
                    s_pat=out.loc[iB, '_PATTERN_']
                    out.loc[iB:iE, 'NEW_GROUP_ID']=c_pattern[s_pat]*n_pattern+i_pattern[s_pat]
                    c_pattern[s_pat]+=1
                    iB=i
            out.sort_values(['NEW_GROUP_ID'], inplace=True)
            out.drop('NEW_GROUP_ID', axis=1, inplace=True)
        else:
            out.sort_values(['BestLogPInGroup','GROUP_ID','FirstInGroupByLogP','LogP','Enrichment'], ascending=[True,True,False,True,False], inplace=True)

        # relabel group id, so that group id are in order of statistical significance
        c_order={}
        cnt=1
        for grp in out.GROUP_ID:
            if grp not in c_order:
                c_order[grp]=cnt
                cnt+=1
        out['GROUP_ID']=out['GROUP_ID'].apply(lambda x: c_order[x])

        out['URL']=''
        out.index=list(range(len(out)))
        S_URL=out.URL.tolist()
        for i in out.index:
            if out.loc[i,'GO'].startswith('M'): #MsigDB
                if re.search(r'\s\(.+\)$', out.loc[i,'Description']):
                    # Notice: description may contain ")" "GSE8515_IL1_VS_IL6_4H_STIM_)MAC_DN"
                    s_key=re.search(r'\s\(.+\)$', out.loc[i,'Description']).group()[2:-1]
                    S_URL[i]='http://http://www.broadinstitute.org/gsea/msigdb/geneset_page.jsp?geneSetName='+s_key
        out['URL']=S_URL
        return out

    @staticmethod
    def representative_rows(t, max_clusters=0):
        t=t[t.FirstInGroupByLogP>0]
        if max_clusters>0:
            t=t[t.GROUP_ID<=max_clusters]
        return set(t.GO)

    @staticmethod
    def sample_rows(t, max_clusters=20, max_members=10, max_nodes=300, l_go_selective=False):
        if max_clusters>0 and max_clusters<t.GROUP_ID.max():
            t=t[t.GROUP_ID<=max_clusters].copy()
        else:
            t=t.copy()

        if (max_nodes>0 or max_members>0):
            S_grp_by_size=[]  # element 5 will contain the 5th GO term from each cluster
            iB=iE=0
            n=len(t)
            #t=t.sort('GROUP_ID')
            if l_go_selective:
                t=t.sort_values(['GROUP_ID','FirstInGroupByGini','GiniIndex','LogP','Enrichment'], ascending=[True,False,False,True,False])
            else:
                t=t.sort_values(['GROUP_ID','FirstInGroupByLogP','LogP','Enrichment'], ascending=[True,False,True,False])
            t.index=list(range(len(t)))
            for i in range(1,n+1):
                if i>=n or t.loc[i,'GROUP_ID']!=t.loc[i-1,'GROUP_ID']:
                    iE=i-1
                    if max_members>0 and iE>iB+max_members-1:
                        iE=iB+max_members-1
                    for j in range(iE-iB+1):
                        if j>=len(S_grp_by_size):
                            S_grp_by_size.append([])
                        S_grp_by_size[j].append(t.loc[iB+j, 'GO'])
                    iB=i
            # need to select a subset of terms for network, b/c a network of too many nodes are hard to visualize
            i_top=i_sz=0
            for i in range(len(S_grp_by_size)):
                i_sz+=len(S_grp_by_size[i])
                if i_sz>max_nodes: break
                i_top+=1
            if i_top==0: i_top=1
            S_node=set()
            for i in range(i_top):
                S_node|=set(S_grp_by_size[i])
            return S_node
        return set(t.GO)

    def network(self, max_clusters=20, max_members=10, max_nodes=300, l_go_selective=False):
        """Construct a GO network, b/c too many nodes will lead to useless networks, we can impose an upper bound
        on the max_nodes size.  We try to take 1 node from each cluster, then 2, then 3, until max_nodes is reached."""
        if len(self.data)==0:
            return None
        if self.DM is None:
            util.error_msg('Please run cluster first!')
        S_node=GO_Cluster.sample_rows(self.t_go, max_clusters=max_clusters, max_members=max_members, max_nodes=max_nodes, l_go_selective=l_go_selective)
        T_node=self.t_go[self.t_go.GO.apply(lambda x: x in S_node)].copy()
        S_go=self.data.header()
        M=self.data.values
        n,m=M.shape
        S_node=set(T_node.GO)
        S_idx=[i for i,x in enumerate(S_go) if x in S_node ]
        S_name=[ S_go[i] for i in S_idx]
        T_node.rename2({'GO':'Gene'})
        s_name='GOCluster'
        if 'Name' in T_node.header():
            s_name=list(T_node.Name)[0]
            T_node.drop('Name', axis=1, inplace=True)
        if 'URL' in T_node.header():
            T_node.drop('URL', axis=1, inplace=True)

        c_has_neighbor={}
        data=[]
        c_cluster={ T_node.loc[i,'Gene']:T_node.loc[i,'GROUP_ID'] for i in T_node.index}
        n2=len(S_idx)
        for _i in range(n2):
            i=S_idx[_i]
            for _j in range(_i+1, n2):
                j=S_idx[_j]
                idx=i*(2*m-i-1)//2+(j-i)-1
                #print (_i, _j, n2, m, i, j, idx, S_name[_i], c_cluster[S_name[_i]], S_name[_j], c_cluster[S_name[_j]], K[idx])
                if self.DM[idx]>=self.similarity:
                    data.append({'Gene_A':S_go[i], 'Gene_B':S_go[j], 'TYPE':'Direct', 'SCORE':self.DM[idx]})
                    c_has_neighbor[S_go[i]]=True
                    c_has_neighbor[S_go[j]]=True
        # keep singletons
        for i in S_idx:
            if S_go[i] not in c_has_neighbor:
                data.append({'Gene_A':S_go[i], 'Gene_B':S_go[i], 'TYPE':'Direct', 'SCORE':1.0})
        if len(data):
            T_edge=pd.DataFrame(data)
        T_node.index=list(range(len(T_node)))
        net=xgmml.Network(T_edge, T_node=T_node, name=s_name)
        return net

if __name__=="__main__":
    Cache.loadL1k()
    exit()
    Cache.load(tax_id=10090, l_use_GPDB=True)
    print(list(Cache.GO_GENE.keys()))
    mydb=db.DB('METASCAPE')
    t=mydb.from_sql('select distinct gid from gid2terms where tax_id=10090')
    S=set(util.iarray2sarray(t.gid.tolist()))
    S2=list(Cache.GENE_GO['GPDB'][10090].keys())
    print(len(S2))
    for x in S2:
        if x not in S:
            print(x)
            exit()

    #print Cache.GO_GENE['LOCAL'][9606].keys()
    Cache.info()
    exit()
    go=GO(tax_id=9606, l_use_GPDB=False)
    print(len(go.GO_GENE))
    exit()
    t=pd.read_csv('/home/RM_Hits.txt')
    S_hit=util.sarray2sarray(t.Gene)
    t_go=go.analysis(S_hit)
    if t_go is not None:
        util.df2sdf(t_go).to_csv('go.csv', index=False)
