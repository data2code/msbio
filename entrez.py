#!/usr/bin/env python
from __future__ import absolute_import
from __future__ import print_function
import sys
import numpy as np
import pandas as pd
import re
import util
import os
import db
from six.moves import range
import setting

class Cache(object):
    DATA_DIR=setting.entrez['DATA_DIR']

    C_TAX_ID={'human':9606, 'mouse':10090, 'rat':10116, 'yeast':4932, 'malaria':5833, 'c. elegans':6239, 'fly':7227, 'zebrafish':7955, 'arabidopsis':3702, 's. pombe':4896}
    C_TAX_NAME={9606:'human', 10090:'mouse', 10116:'rat', 4932:'yeast', 5833:'malaria', 6239:'c. elegans', 7227:'fly', 7955:'zebrafish', 3702:'arabidopsis', 4896:'s. pombe'}
    # in Ensembl, the tax id that has homology data can be diff
    # P. falciparum 3D7: 36329
    # S. cerevisiae: 559292 (baker's yeast)
    # S. pombe: 284812 (fission yeast, strain 972)
    # see http://ensemblgenomes.org/info/genomes

    # share by all tax_id
    C_PATHWAY=None
    C_GENENEW={'LOCAL':{}, 'GPDB':{}}
    C_GENEDESCRIPTION={'LOCAL':{}, 'GPDB':{}}
    C_SPECIES={'LOCAL':{}, 'GPDB':{}}
    C_GENEMAP={'LOCAL':{}, 'GPDB':{}}
    C_CATEGORY_ID={"BP":"19", "MF":"21", "CC":"20"}
    C_CATEGORY={"19":"BP", "20":"CC", "21":"MF"}

    @staticmethod
    def dump():
        def swap(c):
            if 'LOCAL' in c:
                c['LOCAL'], c['GPDB']=c['GPDB'], c['LOCAL']
            return c

        util.dump_object([swap(Cache.C_PATHWAY), swap(Cache.C_GENENEW), swap(Cache.C_GENEDESCRIPTION), \
            swap(Cache.C_SPECIES), swap(Cache.C_GENEMAP), swap(Cache.C_CATEGORY_ID), swap(Cache.C_CATEGORY)],
            os.path.join(Cache.DATA_DIR, "entrez.pickle.gz"))

    @staticmethod
    def get(l_use_GPDB=True, tax_id=9606, user_db=None):
        tax_id=abs(tax_id)
        if tax_id==0: util.error_msg('Tax ID cannot be 0!')
        s_key=Cache.key(l_use_GPDB)
        if not l_use_GPDB:
            X=util.load_object(os.path.join(Cache.DATA_DIR, "entrez.pickle.gz"))
            (Cache.C_PATHWAY, Cache.C_GENENEW, Cache.C_GENEDESCRIPTION, Cache.C_SPECIES, \
            Cache.C_GENEMAP, Cache.C_CATEGORY_ID, Cache.C_CATEGORY) = X
        if tax_id not in Cache.C_GENENEW[s_key]:
            Cache.load(tax_id=tax_id, l_use_GPDB=l_use_GPDB, user_db=user_db)
        return (Cache.C_GENENEW[s_key][tax_id],
            Cache.C_GENEDESCRIPTION[s_key][tax_id],
            Cache.C_SPECIES[s_key],
            Cache.C_GENEMAP[s_key][tax_id],
            Cache.C_PATHWAY
            )

    @staticmethod
    def unload(tax_id, l_use_GPDB):
        tax_id=abs(tax_id)
        s_key=Cache.key(l_use_GPDB)
        if tax_id in Cache.C_GENENEW[s_key]:
            del Cache.C_GENENEW[s_key][tax_id]
            del Cache.C_GENEDESCRIPTION[s_key][tax_id]
            del Cache.C_GENEMAP[s_key][tax_id]
            Cache.C_SPECIES[s_key]={k:v for k,v in Cache.C_SPECIES[s_key].items() if v!=tax_id}

    @staticmethod
    def info():
        if Cache.C_PATHWAY is not None:
            print("C_PATHWAY: %d" % len(Cache.C_PATHWAY))
        for s_key in ('LOCAL','GPDB'):
            print(">Databases: %s" % s_key)
            print("C_SPECIES=%d" % len(Cache.C_SPECIES[s_key]))
            for tax_id in Cache.C_GENENEW[s_key].keys():
                print("TAX_ID=%d (%s)" % (tax_id, Cache.C_TAX_NAME.get(tax_id,"UNKNOWN")))
                print("C_GENENEW: %s" % len(Cache.C_GENENEW[s_key][tax_id]))
                print("C_GENEDESCRIPTION=%d" % len(Cache.C_GENEDESCRIPTION[s_key][tax_id]))
                print("C_GENEMAP=%d" % len(Cache.C_GENEMAP[s_key][tax_id]))
            print("")

    @staticmethod
    def key(l_use_GPDB):
        return 'GPDB' if l_use_GPDB else 'LOCAL'

    @staticmethod
    def load(tax_id=9606, l_use_GPDB=True, user_db=None):
        """tax_id is None, defaults to 9606, if 0, means load all supported species,
        entrez_gene is only used in local mode to accelerate Symbol retrieval"""
        tax_id==abs(tax_id)

        s_key=Cache.key(l_use_GPDB)
        mydb=Cache.get_dbcon(l_use_GPDB, user_db)

        if Cache.C_PATHWAY is None:
            Cache.C_PATHWAY={}
            s_panther=setting.entrez['PANTHER']
            if l_use_GPDB:
                t=mydb.from_sql("SELECT term_id FROM term2term tt WHERE tt.parent_term_id='GO:0008150' and distance=1")
                S_go=["19_"+x for x in t.term_id]
                t=mydb.sql_in("SELECT term_id, term_name FROM term WHERE term_id in (", ")", S_go)
                Cache.C_PATHWAY={ k[3:]:v for k,v in zip(t.term_id.tolist(), t.term_name.tolist()) }
            elif os.path.isfile(s_panther):
                t_panther=pd.read_csv(s_panther, names=['A','B','C','D','E'])
                t_panther.sort_values('C', inplace=True)
                t_panther["RANK"]=range(len(t_panther))
                t_panther["GO"]=""
                for i in t_panther.index:
                    s=t_panther.loc[i,'B']
                    s_name, s_go=re.search(r'^(.+)\s+\((.+)\)$', s).groups()
                    Cache.C_PATHWAY[s_go]=s_name

        if not l_use_GPDB:
            if tax_id not in (0,9606):
                util.error_msg('Local database only supports human!')
            tax_id=9606
            if tax_id in Cache.C_GENENEW[s_key]: return
            S_tax_id=[9606]
        else:
            if tax_id>0 and tax_id in Cache.C_GENENEW[s_key]: return
            if tax_id==0:
                t=mydb.from_sql('SELECT DISTINCT tax_id FROM gid2source_id')
                S_tax_id=[x for x in t.tax_id.astype(int).tolist() if x not in Cache.C_GENENEW[s_key]]
            else:
                S_tax_id=[tax_id]
        if len(S_tax_id)==0: return
        s_tax_id=",".join(util.iarray2sarray(S_tax_id))
        print("Loading Entrez gene table for %s..." % s_tax_id)
        if l_use_GPDB:
            #t_ann=mydb.sql_in("select gs.gid GENE_ID, source_id SYMBOL,a.content DESCRIPTION,gs.tax_id tax_id from gid2source_id gs left join annotation a on gs.gid=a.gid where gs.id_type_id=1 and a.annotation_type_id=4 and gs.tax_id in (", ")", S_tax_id)
            if tax_id==0:
                t1=mydb.from_sql("select gs.gid GENE_ID, source_id SYMBOL,tax_id from gid2source_id gs where gs.id_type_id=1")
                t2=mydb.from_sql("select a.gid GENE_ID, a.content DESCRIPTION, a.tax_id tax_id from annotation a where a.annotation_type_id=4")
            else:
                t1=mydb.sql_in("select gs.gid GENE_ID, source_id SYMBOL,tax_id from gid2source_id gs where gs.id_type_id=1 and gs.tax_id in (", ")", S_tax_id)
                t2=mydb.sql_in("select a.gid GENE_ID, a.content DESCRIPTION, a.tax_id tax_id from annotation a where a.annotation_type_id=4 and a.tax_id in (", ")", S_tax_id)
            t_ann=t1.merge(t2, left_on=['GENE_ID','tax_id'], right_on=['GENE_ID','tax_id'], how='left')
        elif os.path.exists(Cache.DATA_DIR+"N_GENE.csv"):
            t_ann=pd.read_csv(Cache.DATA_DIR+"N_GENE.csv")
            S_name=t_ann.header()
            if "GeneID" in S_name:
                t_ann.rename(columns={'GeneID':'GENE_ID'}, inplace=True)
            if "Symbol" in S_name:
                t_ann.rename(columns={'Symbol':'SYMBOL'}, inplace=True)
            if "Description" in S_name:
                t_ann.rename(columns={'Description':'DESCRIPTION'}, inplace=True)
            if "TAX_ID" in S_name:
                t_ann.rename(columns={'TAX_ID':'tax_id'}, inplace=True)
                t_ann=t_ann[t_ann.tax_id==tax_id].copy()
            else:
                t_ann['tax_id']=tax_id
        t_ann['GENE_ID']=t_ann['GENE_ID'].astype(str)
        X=t_ann.DESCRIPTION.isnull()
        if sum(X):
            t_ann.loc[X, 'DESCRIPTION']=''

        #sw=util.StopWatch('EntrezSW')
        for k in S_tax_id:
            Cache.C_GENENEW[s_key][k]={}
            Cache.C_GENEDESCRIPTION[s_key][k]={}
            Cache.C_GENEMAP[s_key][k]={}

        for k,t_v in t_ann.groupby('tax_id'):
            #t_v=t_v.copy()
            #sw.check('t_ann, '+str(k))
            #C_GENENEW={}
            #C_GENEDESCRIPTION={}
            #C_SPECIES={}
            #for i in t_v.index:
            #    C_GENENEW[t_v.ix[i,'GENE_ID']]=t_v.ix[i,'SYMBOL']
            C_GENENEW=dict(zip(t_v.GENE_ID, t_v.SYMBOL))
            #    C_GENEDESCRIPTION[t_v.ix[i,'GENE_ID']]=t_v.ix[i,'DESCRIPTION']
            C_GENEDESCRIPTION=dict(zip(t_v.GENE_ID, t_v.DESCRIPTION))
            #    C_SPECIES[t_ann.ix[i,'GENE_ID']]=t_ann.ix[i,'tax_id']
            C_SPECIES=dict(zip(t_v.GENE_ID, t_v.tax_id))
            Cache.C_GENENEW[s_key][k]=C_GENENEW
            Cache.C_GENEDESCRIPTION[s_key][k]=C_GENEDESCRIPTION
            Cache.C_SPECIES[s_key].update(C_SPECIES)

        print("Loading mapping table for %s..." % s_tax_id)
        if l_use_GPDB:
            if tax_id==0:
                t_ann=mydb.from_sql("select gid GENE_ID,source_id OLD_GENE_ID,tax_id from gid2source_id where id_type_id=5")
            else:
                t_ann=mydb.sql_in("select gid GENE_ID,source_id OLD_GENE_ID,tax_id from gid2source_id where id_type_id=5 and tax_id in (", ")", S_tax_id)
        elif os.path.exists(Cache.DATA_DIR+"N_GENE_HISTORY.csv"):
            t_ann=pd.read_csv(Cache.DATA_DIR+"N_GENE_HISTORY.csv")
            if "TAX_ID" in S_name:
                t_ann.rename(columns={'TAX_ID':'tax_id'}, inplace=True)
                t_ann=t_ann[t_ann.tax_id==tax_id].copy()
            else:
                t_ann['tax_id']=tax_id
            S_name=t_ann.header()
            if "GeneID" in S_name:
                t_ann.rename(columns={'GeneID':'GENE_ID'}, inplace=True)
            if "OldGeneID" in S_name:
                t_ann.rename(columns={'OldGeneID':'OLD_GENE_ID'}, inplace=True)
        #sw.check('T_ann')
        t_ann['GENE_ID']=t_ann['GENE_ID'].astype(str)
        t_ann['OLD_GENE_ID']=t_ann['OLD_GENE_ID'].astype(str)
        for k,t_v in t_ann.groupby('tax_id'):
            #C_GENEMAP={}
            #for i in t_v.index:
            #    if t_v.ix[i, 'GENE_ID'] not in Cache.C_GENENEW[s_key][k]: continue
            #    C_GENEMAP[t_v.ix[i,'OLD_GENE_ID']]=t_v.ix[i, 'GENE_ID']
            t_v=t_v[t_v.GENE_ID.isin(Cache.C_GENENEW[s_key][k])]
            C_GENEMAP=dict(zip(t_v.OLD_GENE_ID, t_v.GENE_ID))
            Cache.C_GENEMAP[s_key][k]=C_GENEMAP
            #sw.check(str(k)+" ...")

    @staticmethod
    def get_dbcon(l_use_GPDB, user_db=None):
        if l_use_GPDB:
            x=db.DB('METASCAPE')
            return x
        elif user_db is not None:
            return db.open_sqlite(self.user_db)
        elif os.path.exists(Cache.DATA_DIR):
            return None
            #return db.DB('CHEMLIMS')
        util.error_msg('Cannot find database connection!')

class EntrezGene(object):

    VERSION=setting.entrez.get('VERSION',2) # EggNog

    def __init__(self, user_db = None, l_use_GPDB=True, tax_id=None):
        """l_use_GPDB: use Gene Prioritization database"""
        ## for GPDB, if tax_id is provided, only data for that tax_id is pre-loaded and method only works for that tax_id
        ## if GPDB, tax_id will be set to 9606 within connect_db()
        self.GPDB=l_use_GPDB # human only
        self.tax_id=9606 if tax_id is None else tax_id
        self.user_db = user_db
        self.db=Cache.get_dbcon(l_use_GPDB=self.GPDB, user_db=self.user_db)
        (self.C_GENENEW, self.C_GENEDESCRIPTION, self.C_SPECIES, self.C_GENEMAP, self.C_PATHWAY)=Cache.get(l_use_GPDB=self.GPDB, tax_id=self.tax_id, user_db=self.user_db)

    def gene_id_by_id(self, gene_id):
        gene_id=str(gene_id)
        return self.C_GENEMAP.get(gene_id, gene_id)

    def symbol(self, gene_id):
        gene_id=self.gene_id_by_id(gene_id)
        return self.C_GENENEW.get(gene_id, gene_id)

    def species(self, gene_id):
        gene_id=self.gene_id_by_id(gene_id)
        return self.C_SPECIES.get(gene_id, gene_id)

    def description(self, gene_id):
        gene_id=self.gene_id_by_id(gene_id)
        return self.C_GENEDESCRIPTION.get(gene_id, gene_id)

    def GP_annotation_type(self):
        t_map=self.db.from_sql("select t.annotation_type_id ID,t.annotation_type_name NAME,display_name DISPLAY from annotation_type t")
        return t_map

    def GP_annotation(self, S_gene_id, S_ann=None):
        #30: Symbol; 1: Synonyms; 4:full_name; 5,6,7:BP,CC,MF; 8:JAX; 9:OMIM; 10:PubMed; 11:Summary; 17:Disease;
        #12: Tissue; 13:drug target DrugBank; 14:LoF; 18:Drug Traget GeneGo; 36:GOProcess;
        #20: GeneGo Functional Class; 21:Brief Class; 22 Location;
        DEFAULT_ANN=[30,36,20,21]
        S_ann=DEFAULT_ANN if S_ann is None else S_ann
        if not self.GPDB:
            util.error_msg('Only implemented for GPDB')
        if len(S_ann)==0: return
        S_gene_id=util.unique([str(x) for x in S_gene_id])
        t=pd.DataFrame(data={'Gene':S_gene_id})
        t_map=self.db.from_sql("select t.annotation_type_id ID,t.annotation_type_name NAME,display_name DISPLAY from annotation_type t")
        c_map={ k:v for k,v in zip(t_map.ID, t_map.NAME) } #t_map.loc[i,'ID']:t_map.loc[i,'NAME'] for i in t_map.index}
        c_name={ k:v for k,v in zip(t_map.ID, t_map.DISPLAY) } # t_map.loc[i,'ID']:t_map.loc[i,'DISPLAY'] for i in t_map.index}
        for i in S_ann:
            if i not in c_map: continue
            if not c_map[i].startswith('Ontology_'):
                tmp=self.db.sql_in("SELECT gid Gene,content Annotation from annotation t where t.annotation_type_id=? and t.tax_id=? and gid in (", ")", S_gene_id, num=1000, params_before=[i, self.tax_id])
                tmp.rename2({'Annotation': c_name.get(i)})
            else:
                cat_id=c_map[i].replace('Ontology_', '')
                s_like=cat_id+"_%"
                tmp=self.db.sql_in("SELECT gt.term_ids, gt.gid Gene FROM gid2terms gt WHERE gt.term_ids LIKE ? AND tax_id=? and gid in (", ")", S_gene_id, params_before=[s_like, self.tax_id])
                for k in tmp.index:
                    s_gene=tmp.loc[k, 'Gene']
                    S_terms=tmp.loc[k, 'term_ids'].split(",")
                    tmp2=self.db.sql_in("SELECT t.term_id,t.term_name,tg.id_count FROM term t,term2gids tg WHERE t.term_id=tg.term_id AND tg.tax_id=? and tg.term_id IN (", ") ORDER BY tg.id_count ASC LIMIT 0,3", S_terms, params_before=[self.tax_id])
                    tmp.loc[k, 'term_ids']="; ".join(list(tmp2.term_name))
                tmp.rename2({'term_ids': c_name.get(i)})
            tmp['Gene']=tmp['Gene'].astype(str)
            t=t.merge(tmp, left_on=['Gene'], right_on=['Gene'], how='left')
        t.fillna(value='', inplace=True)
        return t

    def GP_gene_sarray_by_affyid_sarray(self, S_affy):
        if not self.GPDB:
            util.error_msg('Only implemented for GPDB')
        S_id=[x for x in S_affy if not pd.isnull(x)]
        t=self.db.sql_in("SELECT gid,source_id from gid2source_id where source_id in (",") and id_type_id in (13,14)", S_id)
        c={ k:str(v) for k,v in zip(t.source_id, t.gid) } #t.ix[i,'source_id']:str(t.ix[i,'gid']) for i in t.index}
        S_gene=[ '' if pd.isnull(x) else c.get(x, '') for x in S_affy ]
        return S_gene

    def _ortholog_version1(self, gene_id, toTaxID):
        """Using homologene"""
        if self.GPDB:
            t=self.db.from_sql("SELECT DISTINCT A.gid GENE_ID FROM homologene A,homologene B WHERE A.homologene_id=B.homologene_id AND B.gid=? AND A.tax_id=? ORDER BY A.gid", [gene_id, toTaxID])
        else:
            t=self.db.from_sql("SELECT DISTINCT A.GENE_ID FROM N_HOMOLOGENE A,N_HOMOLOGENE B WHERE A.HOMOLOGENE_ID=B.HOMOLOGENE_ID AND B.GENE_ID=? AND A.TAX_ID=? ORDER BY A.GENE_ID", [gene_id, toTaxID])
        S=t['GENE_ID'].astype(str).tolist()
        if self.tax_id==toTaxID and str(gene_id) not in S:
            S.append(str(gene_id))
        return list(set(S))

    def ortholog(self, gene_id, toTaxID):
        """Using EggNOG"""
        if not self.GPDB: util.error_msg("Ortholog without GP database Not supported!")
        if EntrezGene.VERSION==1:
            return self._ortholog_version1(gene_id, toTAXID)
        t=self.db.from_sql("SELECT O.gid_B GENE_ID FROM ortholog O WHERE O.gid_A=? AND O.tax_id_B=? ORDER BY O.gid_B", [gene_id, toTaxID])
        S=t['GENE_ID'].astype(str).tolist()
        if self.tax_id==toTaxID and str(gene_id) not in S:
            S.append(str(gene_id))
        return list(set(S))

    def gene_by_RefSeq(self, refseq, toTaxID=None):
        if pd.isnull(refseq) or refseq=='': return None
        toTaxID = toTaxID or self.tax_id
        refseq=re.sub(r'\.\d+', '', refseq)
        refseq=refseq.upper()
        if len(refseq)<=3: return None
        prefix=refseq[:3]
        if prefix in set(["NM_","XM_","NR_","XR_"]):
            if self.GPDB:
                t=self.db.from_sql("SELECT gid GENE_ID,tax_id TAX_ID FROM gid2source_id gs WHERE source_id=? and id_type_id=3", [refseq])
            else:
                t=self.db.from_sql("SELECT GENE_ID,TAX_ID FROM N_GENE_TRANSCRIPT GT WHERE TRANSCRIPT_ID=?", [refseq])
        elif prefix in set(["NP_","XP_","YP_"]):
            if self.GPDB:
                t=self.db.from_sql("SELECT gid GENE_ID,tax_id TAX_ID FROM gid2source_id gs WHERE source_id=? and id_type_id=2", [refseq])
            else:
                t=self.db.from_sql("SELECT GENE_ID,TAX_ID FROM N_GENE_PROTEIN GT WHERE PROTEIN_ID=?", [refseq])
        else:
            return None
        if len(t):
            tmp=t[t.TAX_ID==toTaxID]
            #if toTaxID in t.TAX_ID: always False, as t.TAX_ID is int64 type, does not match toTaxID
            if len(tmp)>0:
                return str(tmp.GENE_ID.iloc[0])
            else:
                S=self.ortholog(t.GENE_ID.iloc[0], toTaxID)
                if len(S):
                    return S[0]
        return None

    def GO_description(self, go_id):
        if self.GPDB:
            t=self.db.from_sql("select t.term_source_id NAME, IF(tc.term_category_id=19, 'BP', IF(tc.term_category_id=20, 'CC', 'MF')) TYPE, t.term_name DESCRIPTION from term_category tc, term t where tc.term_category_id=t.term_category_id and term_source_id=? and tc.ds='GO'", [go_id])
        else:
            t=self.db.from_sql("SELECT NAME,TYPE,DESCRIPTION FROM E_GO WHERE GO_ID=?", [go_id])
        if not len(t): return None
        return {'type':t.loc[0,'TYPE'], 'name':t.loc[0,'NAME'], 'description':t.loc[0,'DESCRIPTION']}

    def all_GO_gene(self, toTaxID=None):
        toTaxID = toTaxID or self.tax_id
        if self.GPDB:
            t=self.db.from_sql("SELECT term_id GO_ID,gids,tax_id from term2gids t where ds='go' and tax_id=?", params=[toTaxID])
            t['GO_ID']=t['GO_ID'].apply(lambda x: re.sub(r'^\d+_', '', x))
            data=[]
            for i in t.index:
                S_gene=t.loc[i,'gids'].split(',')
                data.append(pd.DataFrame({'GO_ID': [t.loc[i,'GO_ID']]*len(S_gene), 'GENE_ID':S_gene, 'TAX_ID':t.loc[i,'tax_id']}))
            t=pd.concat(data, ignore_index=True)
        else:
            t=self.db.from_sql("SELECT DISTINCT GG.PARENT_GO_ID GO_ID,N.GENE_ID FROM E_GO2GO GG,N_GENE_GO N WHERE GG.CHILD_GO_ID=N.GO_ID AND N.TAX_ID=? ORDER BY GO_ID,N.GENE_ID", [toTaxID])
        return t

    def gene_by_GO(self, go_id, toTaxID=None):
        toTaxID = toTaxID or self.tax_id
        if self.GPDB:
            go_id='%'+go_id
            t=self.db.from_sql("SELECT gids,tax_id from term2gids t where t.tax_id=? and ds='go' and term_id like ?", [toTaxID, go_id])
            if len(t):
                return t.loc[0, 'gids'].split(',')
            else:
                return []
        else:
            t=self.db.from_sql("SELECT DISTINCT N.GENE_ID FROM E_GO2GO GG,N_GENE_GO N WHERE GG.CHILD_GO_ID=N.GO_ID AND GG.PARENT_GO_ID=? and N.TAX_ID=? ORDER BY GO_ID,N.GENE_ID", [go_id, toTaxID])
            return list(t.GENE_ID)

    def filter_genes_by_GO(self, go_id, S_gene):
        if self.GPDB:
            S=set(self.gene_by_GO(go_id))
            return [x for x in S_gene if x in S]
        else:
            t=self.db.sql_in("SELECT DISTINCT N.GENE_ID FROM E_GO2GO GG,N_GENE_GO N WHERE GG.CHILD_GO_ID=N.GO_ID AND GG.PARENT_GO_ID=? and N.GENE_ID in (", ")", S_gene, params_before=[go_id])
            return list(t.GENE_ID) if len(t) else []

    def IPRDescription(self, ipr_id):
        if self.GPDB:
            util.error_msg('This method is not implemented for GPDB')
        t=self.db.from_sql("SELECT DESCRIPTION FROM E_IPR WHERE IPR_ID=?", [ipr_id])
        if not len(t): return ""
        return t.loc[0,'DESCRIPTION']

    def all_IPR_gene(self, toTaxID=None):
        if self.GPDB:
            util.error_msg('This method is not implemented for GPDB')
        toTaxID=toTaxID or self.tax_id
        t=self.db.from_sql("SELECT DISTINCT I.PARENT_IPR IPR,N.N_GENE_ID GENE_ID FROM E_IPR2IPR I,E_PROTEIN_IPR P,E_IDMAP M,ENSEMBL_NCBI N WHERE I.CHILD_IPR=P.IPR_ID AND P.TAX_ID=? AND P.PROTEIN_ID=M.PROTEIN_ID AND M.GENE_ID=N.E_GENE_ID ORDER BY IPR,GENE_ID", params=[toTaxID])
        return t

    def GO_descript_by_gene(self, gene_id, s_type="BP", nTop=3):
        GO_ROOT={"BP":"GO:0008150","MF":"GO:0003674","CC":"GO:0005575"}
        s_root=GO_ROOT[s_type]
        if self.GPDB:
            C_ANN_ID={"BP":5, "CC":6, "MF":7}
            s_cat=C_ANN_ID[s_type]
            t=self.db.from_sql("Select content from annotation where gid=? and annotation_type_id=?", [gene_id, s_cat])
        else:
            t=self.db.from_sql("SELECT GO.GO_ID,GO.NAME,GO.TYPE,GG.DISTANCE FROM N_GENE_GO G,E_GO2GO GG,E_GO GO WHERE G.GO_ID=GG.CHILD_GO_ID AND GG.CHILD_GO_ID=GO.GO_ID AND G.GENE_ID=? AND GG.PARENT_GO_ID=? GROUP BY GO.GO_ID ORDER BY DISTANCE DESC,GO.GO_ID LIMIT ?", [gene_id, s_root, nTop])
        S=[]
        for i in range(len(t)):
          S.append(t.loc[i,'GO_ID']+"("+t.loc[i,'TYPE']+"):"+t.loc[i,'NAME'])
        s="; ".join(S)
        if len(s)>200: s=s[:200]
        return s

    def child_GO(self, go_id, l_recursive=False):
        if self.GPDB:
            s_dist='distance>0' if l_recursive else 'distance=1'
            t=self.db.from_sql("SELECT tt.term_id GO_ID,t.term_name NAME,MIN(tt.distance) DISTANCE FROM term2term tt,term t WHERE tt.term_id=t.term_source_id AND tt.parent_term_id=? AND "+s_dist+"  group by tt.term_id,t.term_name order by distance,term_id", [go_id])
        else:
            s_dist='DISTANCE>0' if l_recursive else 'DISTANCE=1'
            t=self.db.from_sql("SELECT GO.GO_ID,GO.NAME,MIN(GG.DISTANCE) DISTANCE FROM E_GO GO,E_GO2GO GG WHERE GG.CHILD_GO_ID=GO.GO_ID and GG.PARENT_GO_ID=? and "+s_dist+" GROUP BY GO.GO_ID,GO.NAME ORDER BY DISTANCE,GO_ID", [go_id])
        t.rename(columns={'GO_ID':'GO'}, inplace=True)
        return t

    def parent_GO(self, go_id, l_recursive=False):
        if self.GPDB:
            s_dist='distance>0' if l_recursive else 'distance=1'
            t=self.db.from_sql("SELECT tt.parent_term_id GO_ID,t.term_name NAME,MIN(tt.distance) DISTANCE FROM term2term tt,term t WHERE tt.parent_term_id=t.term_source_id AND tt.term_id=? AND "+s_dist+" group by tt.parent_term_id,t.term_name order by distance,parent_term_id", [go_id])
        else:
            s_dist='DISTANCE>0' if l_recursive else 'DISTANCE=1'
            t=self.db.from_sql("SELECT GO.GO_ID,GO.NAME,MIN(GG.DISTANCE) DISTANCE FROM E_GO GO,E_GO2GO GG WHERE GG.PARENT_GO_ID=GO.GO_ID and GG.CHILD_GO_ID=? and "+s_dist+" GROUP BY GO.GO_ID,GO.NAME ORDER BY DISTANCE,GO_ID", [go_id])
        t.rename(columns={'GO_ID':'GO'}, inplace=True)
        return t

    def id_conversion(self, S_source_id, target_tax_id=None, Source_tax_id=None, s_source_type=None):
        if not self.GPDB:
            util.error_msg('Only works for GPDB')
        target_tax_id=target_tax_id or self.tax_id
        Source_tax_id=Source_tax_id or [self.tax_id]
        if EntrezGene.VERSION==1:
            t=EntrezGene._id_conversion_version1(S_source_id, target_tax_id=target_tax_id, Source_tax_id=Source_tax_id, s_source_type=s_source_type, con=self.db)
        else:
            t=EntrezGene._id_conversion(S_source_id, target_tax_id=target_tax_id, Source_tax_id=Source_tax_id, s_source_type=s_source_type, con=self.db)
        t=t[['RANK','InputID','gene_id']].copy()
        return t

    @staticmethod
    def _id_conversion_version1(S_source_id, target_tax_id=None, Source_tax_id=None, s_source_type=None, con=None):
        """Use homologene table"""
        # s_source_type has to be one of: Entrez, Symbol, RefSeq
        # symbol, RefSeq_Proteins, RefSeq_RNAs, gene_synonym, Gene_History, uniprot, ensembl, ensembl_protein, ensembl_transcript, ucsc
        def source_id_strip(s):
            p=re.compile('\W+')
            S=[]
            s=s.upper()
            for m in p.finditer(s):
                S.append(s[:m.start()])
            S.append(s)
            return S

        data=[]
        for x in util.unique2(S_source_id):
            for y in source_id_strip(str(x)):
                data.append([x, y])
        # If input is S000001.1:pep-3, it will be cast into ["S000001.1:pep-3", "S000001.1:pep", "S000001.1", "S000001"] for search purpose
        t_id=pd.DataFrame(data, columns=['InputID','source_id'])
        t_id['RANK']=range(len(t_id))

        S_source_type=[]
        if s_source_type is None:
            pass
        elif s_source_type=='Entrez':
            S_source_type=['Gene_History']
        elif s_source_type=='RefSeq':
            S_source_type=['RefSeq_RNAs','RefSeq_Proteins']
        elif s_source_type=='Symbol':
            S_source_type=['Symbol','gene_synonym']
        elif s_source_type=='dbxref':
            S_source_type=['locus_tag','dbxref']
        elif s_source_type=='Ensembl':
            S_source_type=['ensembl_gene_id', 'ensembl_peptide_id', 'ensembl_transcript_id']
        elif type(s_source_type) is not list:
            S_source_type=[s_source_type]
        if Source_tax_id==0: Source_tax_id=None
        if Source_tax_id is not None and type(Source_tax_id) is not list:
            Source_tax_id=[Source_tax_id]
        if S_source_type is not None and len(S_source_type)>0:
            s_id_type=" AND t.id_type_name IN ({0})".format("'"+"','".join(S_source_type)+"'")
        else:
            s_id_type=""

        if con is None:
            con=Cache.get_dbcon(l_use_GPDB=True, user_db=None)

        s_source_tax_id=""
        if Source_tax_id is not None:
            s_source_tax_id="AND h1.tax_id in ({0})".format(",".join([str(x) for x in Source_tax_id]))

        s_target_tax_id=""
        if target_tax_id is not None:
            s_target_tax_id="AND h2.tax_id IN ({0})".format(target_tax_id)

        t=con.sql_in("""
        SELECT ucase(gs.source_id) source_id,
            CONVERT(gs.gid, CHAR) as gid,
            t.id_type_name,
            convert(gs.tax_id, char) as tax_id,
            CAST(h2.gid AS UNSIGNED) AS homologene_gid,
            convert(h2.tax_id,char) AS homologene_tax_id,
            a.content AS priority2,
            t.keep_first_order as priority1,
            b.content as priority3,
            CASE
                WHEN gs.gid = h2.gid  THEN  1
                ELSE 0
            END as priority4
        FROM id_type t, gid2source_id gs
        LEFT JOIN homologene h1 ON h1.gid = gs.gid {0}
        LEFT JOIN homologene h2 ON h1.homologene_id = h2.homologene_id {1}
        LEFT JOIN annotation a ON h2.gid  = a.gid AND a.annotation_type_id = 10
        LEFT JOIN annotation b ON h2.gid  = b.gid AND b.annotation_type_id = 69
        WHERE gs.source_id IN (""".format(s_source_tax_id, s_target_tax_id),
        """)
        {0}
        AND gs.id_type_id = t.id_type_id""".format(s_id_type), t_id.source_id.tolist())

        for i in t.index:
            if pd.isnull(t.loc[i, 'homologene_gid']) and (int(t.loc[i, 'tax_id'])==target_tax_id):
                t.loc[i, 'homologene_gid']=int(t.loc[i, 'gid'])
                t.loc[i, 'homologene_tax_id']=str(target_tax_id)
        t2=t[t['homologene_gid'].apply(lambda x: not pd.isnull(x))].copy()
        if len(t2)==0:
            t2=t[:0].copy()
        #t3=pd.DataFrame({'InputID':S_source_id, 'source_id':[x.upper() for x in S_source_id]})
        t2['priority2']=t2['priority2'].apply(lambda x: 0 if x is None else x)
        t2['priority3']=t2['priority3'].apply(lambda x: 0 if x is None else x)
        t2['homologene_gid']=t2['homologene_gid'].apply(lambda x: util.r2i2s(x))
        t=t_id.merge(t2, left_on='source_id', right_on='source_id')
        # longest source_id first
        t.sort_values(['homologene_tax_id','InputID','source_id', 'priority4','priority1', 'priority3','priority2','homologene_gid'], ascending=[True, True, False, False, True, False, False, True], inplace=True)
        t.drop_duplicates(['InputID','homologene_tax_id'], inplace=True)
        t=t[['RANK','InputID','homologene_gid','homologene_tax_id','tax_id']]
        t.sort_values('RANK', inplace=True)
        #t['homologene_gid']=t['homologene_gid'].apply(lambda x: '' if pd.isnull('homologene_id') else str(x))
        t.rename2({'homologene_gid':'gene_id'})
        return t

    @staticmethod
    def _id_conversion(S_source_id, target_tax_id=None, Source_tax_id=None, s_source_type=None, con=None):
        """Use EggNOG ortholog table"""
        # s_source_type has to be one of: Entrez, Symbol, RefSeq
        # symbol, RefSeq_Proteins, RefSeq_RNAs, gene_synonym, Gene_History, uniprot, ensembl, ensembl_protein, ensembl_transcript, ucsc
        def source_id_strip(s):
            p=re.compile('\W+')
            S=[]
            s=s.upper()
            for m in p.finditer(s):
                S.append(s[:m.start()])
            S.append(s)
            return S

        data=[]
        for x in util.unique2(S_source_id):
            for y in source_id_strip(str(x)):
                data.append([x, y])
        # If input is S000001.1:pep-3, it will be cast into ["S000001.1:pep-3", "S000001.1:pep", "S000001.1", "S000001"] for search purpose
        t_id=pd.DataFrame(data, columns=['InputID','source_id'])
        t_id['RANK']=range(len(t_id))

        S_source_type=[]
        if s_source_type is None:
            pass
        elif s_source_type=='Entrez':
            S_source_type=['Gene_History']
        elif s_source_type=='RefSeq':
            S_source_type=['RefSeq_RNAs','RefSeq_Proteins']
        elif s_source_type=='Symbol':
            S_source_type=['Symbol','gene_synonym']
        elif s_source_type=='dbxref':
            S_source_type=['locus_tag','dbxref']
        elif s_source_type=='Ensembl':
            S_source_type=['ensembl_gene_id', 'ensembl_peptide_id', 'ensembl_transcript_id']
        elif type(s_source_type) is not list:
            S_source_type=[s_source_type]
        if Source_tax_id==0: Source_tax_id=None
        if Source_tax_id is not None and type(Source_tax_id) is not list:
            Source_tax_id=[Source_tax_id]
        if S_source_type is not None and len(S_source_type)>0:
            s_id_type=" AND t.id_type_name IN ({0})".format("'"+"','".join(S_source_type)+"'")
        else:
            s_id_type=""

        if con is None:
            con=Cache.get_dbcon(l_use_GPDB=True, user_db=None)

        s_source_tax_id=""
        if Source_tax_id is not None:
            s_source_tax_id="AND h.tax_id_A in ({0})".format(",".join([str(x) for x in Source_tax_id]))

        s_target_tax_id=""
        if target_tax_id is not None:
            s_target_tax_id="AND h.tax_id_B IN ({0})".format(target_tax_id)

        t=con.sql_in(f"""
        SELECT ucase(gs.source_id) source_id,
            CONVERT(gs.gid, CHAR) as gid,
            t.id_type_name,
            convert(gs.tax_id, char) as tax_id,
            CAST(h.gid_B AS UNSIGNED) AS homologene_gid,
            convert(h.tax_id_B,char) AS homologene_tax_id,
            h.pubmed AS priority2,
            t.keep_first_order as priority1,
            h.rif as priority3,
            0 as priority4
        FROM id_type t, gid2source_id gs
        LEFT JOIN ortholog h ON h.gid_A = gs.gid {s_source_tax_id} {s_target_tax_id}
        WHERE gs.source_id IN (""",
        f""") {s_id_type} AND gs.id_type_id = t.id_type_id""", t_id.source_id.tolist())

        mask=t.homologene_gid.isnull()
        t.loc[mask, 'homologene_gid']=t.loc[mask, 'gid']
        t.loc[mask, 'homologene_tax_id']=t.loc[mask, 'tax_id']
        t.loc[mask, 'priority4']=1
        t['priority1'].fillna(0, inplace=True)
        t['priority2'].fillna(0, inplace=True)
        t['priority3'].fillna(0, inplace=True)
        t['homologene_gid']=t['homologene_gid'].astype(str)
        if target_tax_id is not None:
            t=t[t.homologene_tax_id==str(target_tax_id)].copy()
        t=t_id.merge(t, left_on='source_id', right_on='source_id')
        # longest source_id first
        t.sort_values(['homologene_tax_id','InputID','source_id', 'priority4','priority1', 'priority3','priority2','homologene_gid'], ascending=[True, True, False, False, True, False, False, True], inplace=True)
        t.drop_duplicates(['InputID','homologene_tax_id'], inplace=True)
        t=t[['RANK','InputID','homologene_gid','homologene_tax_id','tax_id']]
        t.sort_values('RANK', inplace=True)
        #t['homologene_gid']=t['homologene_gid'].apply(lambda x: '' if pd.isnull('homologene_id') else str(x))
        t.rename2({'homologene_gid':'gene_id'})
        return t

    def gene_by_symbol(self, s_symbol, toTaxID=None, l_keep_all=False):
        toTaxID=toTaxID or self.tax_id
        #print ">>> ", s_symbol
        if pd.isnull(s_symbol) or s_symbol=='':
            return None
        if type(s_symbol)!=str:
            s_symbol=str(s_symbol)
        s_symbol=s_symbol.upper()
        if re.search(r'^(ENSG|ENSMUSG|ENSRNOG|LRG_)\d+$', s_symbol):
            # ENSEMBL GENE
            s_symbol=re.sub(r'\.\d+$', '', s_symbol)
            if self.GPDB:
                t=self.db.from_sql("SELECT gid GENE_ID,tax_id FROM gid2source_id WHERE source_id=? and tax_id=? and id_type_id in (8,10,11)", [s_symbol, toTaxID])
            else:
                t=self.db.from_sql("SELECT N_GENE_ID GENE_ID,TAX_ID FROM ENSEMBL_NCBI where E_GENE_ID=?", [s_symbol])
                t.rename2({'TAX_ID':'tax_id'})
            if len(t):
                tmp=t[t.tax_id==toTaxID]
                if len(tmp)>0:
                    return str(tmp.GENE_ID.iloc[0]) if not l_keep_all else util.iarray2sarray(tmp.GENE_ID['GENE_ID'])
                else:
                    S=self.ortholog(t.GENE_ID.iloc[0], toTaxID)
                    if len(S):
                        return S[0] if not l_keep_all else S
        elif re.search(r'^(ENST|ENSMUST|ENSRNOT|LRG_)\d+(_t1)?$', s_symbol):
            util.warn_msg('We currenting missing table E_GENE_TRANSCRIPT!!!')
        else:
            if self.GPDB:
                t=self.db.from_sql("SELECT gid GENE_ID FROM gid2source_id WHERE source_id=? and tax_id=? and id_type_id in (1,6,15,16)", [s_symbol, toTaxID])
            else:
                t=self.db.from_sql("SELECT GENE_ID FROM N_GENE_SYNONYM WHERE SYMBOL=? and TAX_ID=?", [s_symbol, toTaxID])
        if len(t)>0: return str(t.loc[0,'GENE_ID']) if not l_keep_all else util.iarray2sarray(t['GENE_ID'])
        return None

    def fix_gene_id(self, gene_id=None, s_refseq=None, s_symbol=None, toTaxID=None, l_ortholog=True):
        toTaxID=toTaxID or self.tax_id
        s_gene_id=None
        gene_id=str(gene_id or 0)
        if re.search(r'^\d+$', gene_id):
            if gene_id in self.C_GENEMAP: gene_id=self.C_GENEMAP[gene_id]
            if gene_id in self.C_GENENEW and self.C_SPECIES[gene_id]==toTaxID: return gene_id
            # valid gene_id
            # make sure gene_id is one of human/mouse gene, so we don't
            # convert genes from other species to human
            if gene_id in self.C_GENENEW and self.C_SPECIES[gene_id]!=toTaxID and l_ortholog:
                S=self.ortholog(gene_id, toTaxID)
                if len(S): return S[0]
        # try to use refseq
        if s_refseq:
            s_gene_id=self.gene_by_RefSeq(s_refseq, toTaxID)
            if s_gene_id: return s_gene_id
        # try to use symbol
        if s_symbol:
            s_gene_id=self.gene_by_symbol(s_symbol, toTaxID)
            if s_gene_id: return s_gene_id
        return None

    def gene_sarray_to_gene_sarray(self, S_gene_id, toTaxID=None, l_ortholog=True):
        toTaxID=toTaxID or self.tax_id
        S_gene_id=util.sarray2sarray(S_gene_id)
        n=len(S_gene_id)
        for i in range(n):
            S_gene_id[i]=self.fix_gene_id(gene_id=S_gene_id[i], s_refseq=None, s_symbol=None, toTaxID=toTaxID, l_ortholog=l_ortholog)
        return S_gene_id

    def gene_sarray_to_ortholog(self, S_gene_id, toTaxID=None):
        toTaxID=toTaxID or self.tax_id
        S_gene_id=util.sarray2sarray(S_gene_id)
        S_target_id=['']*len(S_gene_id)
        for i in range(len(S_gene_id)):
            gene_id=S_gene_id[i]
            if not gene_id: continue
            if gene_id in self.C_GENEMAP: gene_id=self.C_GENEMAP[gene_id]
            if gene_id in self.C_GENENEW and self.C_SPECIES[gene_id]==toTaxID:
                S_target_id[i]=gene_id
                continue
            S=self.ortholog(gene_id, toTaxID)
            if len(S):
                S_target_id[i]=S[0]
        return S_target_id

    def gene_table_to_gene_sarray(self, T, C_columns=None, toTaxID=None, l_ortholog=True):
        C_columns=C_columns or {'GENE_ID':"Gene",'REFSEQ':"RefSeq",'SYMBOL':"Symbol"}
        toTaxID=toTaxID or self.tax_id
        S_gene_id=['']*len(T)
        S_name=T.header()
        gene_idx=refseq_idx=symbol_idx=-1
        if 'GENE_ID' in C_columns: gene_idx=util.index(C_columns['GENE_ID'], S_name)
        if 'REFSEQ' in C_columns: refseq_idx=util.index(C_columns['REFSEQ'], S_name)
        if 'SYMBOL' in C_columns: symbol_idx=util.index(C_columns['SYMBOL'], S_name)
        n=len(T)
        for i in range(n):
            S_gene_id[i]=self.fix_gene_id(
                gene_id=T.iloc[i,gene_idx] if gene_idx>=0 else None,
                s_refseq=T.iloc[i,refseq_idx] if refseq_idx>=0 else None,
                s_symbol=T.iloc[i,symbol_idx] if symbol_idx>=0 else None,
                toTaxID=toTaxID, l_ortholog=l_ortholog)
            #print T.Symbol_.iloc[i], T.RefSeq_.iloc[i], util.info_msg(S_gene_id[i])
        S_gene_id=[ '' if pd.isnull(x) else str(x) for x in S_gene_id]
        return S_gene_id

    def gene_sarray_to_table(self, S_gene, l_description=True):
        n=len(S_gene)
        S_gene=util.sarray2sarray(S_gene)
        T=pd.DataFrame({'Gene':S_gene, 'Symbol':['']*n})
        if l_description: T['Description']=['']*n
        for i in range(len(S_gene)):
            if S_gene[i] and S_gene[i] in self.C_GENENEW:
                T.loc[i,'Symbol']=self.C_GENENEW[S_gene[i]] or ''
                if l_description: T.loc[i,'Description']=self.C_GENEDESCRIPTION[S_gene[i]] or ''
        return T

    def annotate_table_by_gene_id(self, T, s_col="Gene"):
        idx=util.index(s_col, T.header())
        if idx<0: util.error_msg("column "+s_col+" is not found in the table!")
        S=util.sarray2sarray(T.iloc[:,idx])
        T.iloc[:,idx]=S
        t=self.gene_sarray_to_table(S)
        S_name=T.header()
        #if "Symbol" not in S_name: T['Symbol']=util.sarray2sarray(t['Symbol'])
        #if "Description" not in S_name: T['Description']=util.sarray2sarray(t['Description'])
        T['Symbol']=util.sarray2sarray(t['Symbol'])
        T['Description']=util.sarray2sarray(t['Description'])
        return T

    def go_by_gene(self, gene_id, s_keyword=""):
        t=None
        if self.GPDB:
            s_sql_keyword="AND term_name like '%"+s_keyword+"%'" if s_keyword else ""
            t=self.db.from_sql("SELECT term_ids from gid2terms where gid=? and term_category_id in (19,20,21)", params=[gene_id])
            s=",".join(t.term_ids.tolist())
            S_go=util.unique(s.split(','))
            t=self.db.sql_in("SELECT term_id GO_ID, term_name NAME from term where term_id in (", ") "+s_sql_keyword, S_go)
            t['GO_ID']=t['GO_ID'].apply(lambda x: re.sub(r'^\d+_', '', x))
        else:
            s_sql_keyword="AND GO.NAME like '%"+s_keyword+"%'" if s_keyword else ""
            t=self.db.from_sql("SELECT GO.GO_ID,GO.NAME,GO.TYPE,MIN(GG.DISTANCE) DISTANCE FROM N_GENE_GO G,E_GO2GO GG,E_GO GO WHERE G.GO_ID=GG.CHILD_GO_ID AND GG.PARENT_GO_ID=GO.GO_ID AND G.GENE_ID=? {} GROUP BY TYPE,GO.GO_ID ORDER BY TYPE,DISTANCE,GO.GO_ID".format(s_sql_keyword), params=[gene_id])
        return t

    def protein_accession_to_gene_id(self, s_IPI, toTaxID=None, l_ortholog=True):
        if self.GPDB:
            util.error_msg('Method not implemented for GPDB')
        toTaxID=toTaxID or self.tax_id
        t_ann=self.db.from_sql("SELECT DISTINCT EN.N_GENE_ID,EN.TAX_ID FROM ENSEMBL_NCBI EN,E_IDMAP M,E_PROTEIN_ACCESSION P WHERE EN.E_GENE_ID=M.GENE_ID AND M.PROTEIN_ID=P.PROTEIN_ID AND P.ACCESSION=? ORDER BY N_GENE_ID", [s_IPI])
        if len(t_ann):
            t=t_ann[t_ann['TAX_ID']==toTaxID]
            if len(t):
                return str(t['N_GENE_ID'][0])
            else:
                S=self.ortholog(t_ann.loc[0,'N_GENE_ID'], toTaxID)
                return S or None
        return None

    def all_genes(self, taxID=None):
        taxID=taxID or Cache.C_TAX_ID['human']
        if self.GPDB:
            #t_ann=db.from_sql(self.dbcon, "SELECT distinct gid GENE_ID FROM gid2source_id where id_type_id=1 and tax_id=?", [taxID])
            t_ann=self.db.from_sql("SELECT gid GENE_ID FROM annotation WHERE annotation_type_id=3 AND content='protein-coding' and tax_id=?", [taxID])
        else:
            t_ann=self.db.from_sql("SELECT GENE_ID FROM N_GENE WHERE TAX_ID=?", [taxID])
        return util.sarray2sarray(t_ann['GENE_ID'])

    def all_gene_names(self, taxID=None):
        taxID=taxID or self.tax_id
        if self.GPDB:
            t_ann=self.db.from_sql("SELECT gid GENE_ID,source_id SYMBOL FROM gid2source_id WHERE tax_id=? and id_type_id=1", [taxID])
        else:
            t_ann=self.db.from_sql("SELECT GENE_ID,SYMBOL FROM N_GENE WHERE TAX_ID=?", [taxID])
        c_ann={}
        for i in range(len(t_ann)):
            c_ann[str(t_ann.loc[0,'GENE_ID'])]=t_ann.loc[i,'SYMBOL']
        return c_ann

    def pathway_GO(self, s_go):
        tmp=self.parent_GO(s_go, l_recursive=True)
        S_GO=util.unique([self.C_PATHWAY[x] for x in tmp.GO if x in self.C_PATHWAY])
        S_GO.sort()
        return "; ".join(S_GO)

    def pathway_label(self, s_gene):
        s_inGO="'"+"','".join(list(self.C_PATHWAY.keys()))+"'"
        if not s_gene: return ""
        if self.GPDB:
            tmp=self.go_by_gene(s_gene)
        else:
            s_sql="SELECT DISTINCT gg.PARENT_GO_ID GO_ID FROM E_GO2GO gg,N_GENE_GO g WHERE gg.CHILD_GO_ID=g.GO_ID AND gg.PARENT_GO_ID in ("+s_inGO+") AND g.GENE_ID="
            tmp=self.db.from_sql(s_sql+s_gene)
        S_GO=util.unique([self.C_PATHWAY[x] for x in tmp.GO_ID if x in self.C_PATHWAY])
        S_GO.sort()
        return "; ".join(S_GO)

if __name__=='__main__':
    #Cache.load(tax_id=0, l_use_GPDB=True)
    ##Cache.load(tax_id=0, l_use_GPDB=False)
    #Cache.info()
    #exit()
    ez=EntrezGene(tax_id=9606, l_use_GPDB=True)
    ez.id_conversion(["Tlr7","tlr9","KRAS"], target_tax_id=10090, Source_tax_id=None, s_source_type=None).display()
    ez.id_conversion(["Tlr7","tlr9","KRAS"], target_tax_id=9606, Source_tax_id=None, s_source_type=None).display()
    ez.id_conversion([ 51284,54106,3845 ], target_tax_id=10090, Source_tax_id=None, s_source_type=None).display()
    print("OLD VERSION 1")
    EntrezGene.VERSION=1
    ez.id_conversion(["Tlr7","tlr9","KRAS"], target_tax_id=10090, Source_tax_id=None, s_source_type=None).display()
    ez.id_conversion(["Tlr7","tlr9","KRAS"], target_tax_id=9606, Source_tax_id=None, s_source_type=None).display()
    ez.id_conversion([ 51284,54106,3845 ], target_tax_id=10090, Source_tax_id=None, s_source_type=None).display()

    ez=EntrezGene(tax_id=10090, l_use_GPDB=True)
    ez.id_conversion([ 170743,16653,81897 ], target_tax_id=9606, Source_tax_id=None, s_source_type=None).display()
    print("OLD VERSION 1")
    ez.id_conversion([ 170743,16653,81897 ], target_tax_id=9606, Source_tax_id=None, s_source_type=None).display()
    exit()
    ez=EntrezGene(tax_id=9606, l_use_GPDB=True)
    print(ez.gene_by_symbol('Prkci', Cache.C_TAX_ID['rat'], True))
    exit()
    print(ez.child_GO('GO:0007420', False)[:20])
    print(ez.parent_GO('GO:0007420', True)[:20])
    exit()
    print(ez.fix_gene_id(None, 'NM_013867'))
    print(ez.fix_gene_id(None, 'NM_013867', toTaxID=Cache.C_TAX_ID['mouse']))
    exit()
    print(ez.fix_gene_id(None))
    print(ez.fix_gene_id(22995))
    print(ez.fix_gene_id(99100))
    print(ez.fix_gene_id(331241334,"NM_006327"))
    print(ez.fix_gene_id(331241334,"NM_","TLR7"))
