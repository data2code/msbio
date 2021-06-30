#!/usr/bin/env python
from __future__ import absolute_import
from __future__ import print_function
from __future__ import division
import util
import re
import glob
import os
import math
import pandas as pd
import numpy as np
import time
import scipy.stats as ss
import matplotlib
from six.moves import range
matplotlib.use('Agg')
import sklearn.metrics.pairwise as skmp

class DM(object):

    def __init__(self, s_file='', R_DM=None):
        self.file_name = s_file
        self.DM=None
        self.n=0
        self.N=0
        self.dmax=0.0
        if s_file!='':
            with open(s_file) as f:
                self.DM=f.readlines()
            self.n=int(self.DM.pop(0))
            self.DM=util.sarray2rarray(self.DM)
            self.N=len(self.DM)
        elif R_DM is not None:
            self.DM = R_DM
            self.N=len(self.DM)
            self.n=DM.calc_n(self.N)
        else:
            error_msg('DM.__init__: either file name or DM array must be given!')
        self.dmax=np.nanmax(self.DM)
        #self.DM[np.isnan(self.DM)] = self.dmax if fillna<0 else 1
        #if fillna>0: self.dmax=np.max([self.dmax, fillna])

    @staticmethod
    def calc_n(N):
        n=int((1+math.sqrt(1+8*N))//2)
        if DM.calc_N(n) != N: util.error_msg('DM.calc_n: DM size ('+str(N)+') not correct!')
        return n

    @staticmethod
    def calc_N(n):
        return n*(n-1)//2

    def normalize(self, byMax=False):
        if self.dmax!=1:
            self.dmax=self.DM.max()
            if byMax or self.dmax>1:
                self.DM/=self.dmax
                self.dmax=1.0

    def save(self, s_file='', s_format='%.2f'):
        s = s_file if s_file!='' else self.file_name
        # reimplement util.save_list here, as memory may not be large enough to hold the list, when we do self.n+self.DM
        with open(s, 'w') as f:
            f.write(str(self.n)+'\n')
            for r in self.DM:
                f.write('\n' if np.isnan(r) else (s_format % r)+'\n')
        #util.save_list(s, [str(self.n)]+util.rarray2sarray(self.DM, s_format), s_end='\n')

    def get(self, i, j):
        if j<i: (i, j) = (j, i)
        return self.DM[i*(2*n-i-1)//2+(j-i)-1]

    def pair(self, x):
        'Given index in DM, return i, j'
        i=int(math.floor(2*self.n-1-math.sqrt((2*self.n-1)**2-8*x)//2))
        j=int(x-i*(2*self.n-i-1)//2+i+1)
        return (i, j)

import setting
class Clustering(object):
    "Wrapping for running various clustering binary program for hierachical clustering"
    BIN_CWC=setting.cluster['BIN_CWC']
    BIN_HYB=setting.cluster['BIN_HYB']
    BIN_OPT=setting.cluster['BIN_OPT']
    DEFAULT_INPUT_OPT={'ID':'Gene', 'DESCRIPTION':'Description', 'WEIGHT_COL':'WEIGHT', 'GENE_WEIGHT':[], 'DATA_COLS':[], 'GENE_NORMALIZE':False, 'NORMALIZE_METHOD':'Z', 'EXP_WEIGHT':[]}
    # BIN can be HYBRID or CWC
    # set default to CWC, HYBRID has bug
    DEFAULT_CLUSTER_OPT={'BIN':'CWC', 'GENE':True, 'EXP':False, 'DMG':'', 'DME':'', 'GENE_METRICS':'BUILD_IN', 'EXP_METRICS':'BUILD_IN', 'OPTIMIZE':True, 'CLEANUP':True, 'HAS_NULL':True, 'FINGERPRINT':False, 'SKIP_DM':False }
    COUNT_INTERVAL=10000

    def __init__(self, input='', table=None, input_options=None, cluster_options=None, user_hybrid = None):
        if user_hybrid is not None:
            user_hybrid_path = user_hybrid.split()[-1]
            if os.path.exists(user_hybrid_path):
                Clustering.BIN_HYB=user_hybrid
            else:
                util.error_msg("Clustering tool " + user_hybrid_path + " does not exist!")
        self.input_opt={}
        self.cluster_opt={}
        self.input=''
        self.table=None
        if input!='':
            if re.search(r'\.input$', input) is not None:
                self.input=re.sub(r'\.input$', '', s_input)
            else:
                self.input=input
        input_options = input_options or {}
        cluster_options = cluster_options or {}
        self.input_opt=Clustering.DEFAULT_INPUT_OPT.copy()
        self.input_opt.update(input_options)
        if 'EXP_WEIGHT' in cluster_options:
            util.error_msg('Clustering.__init__: EXP_WEIGHT has been moved from cluster_options into input_options!')
        if 'DATA_COLS' in cluster_options:
            util.error_msg('Clustering.__init__: DATA_COLS should be into input_options, not cluster_options!')
        self.cluster_opt=Clustering.DEFAULT_CLUSTER_OPT.copy()
        self.cluster_opt.update(cluster_options)
        if table is not None:
            self.table=table
            if self.table.col_type(self.input_opt['ID'])!='s':
                self.table[self.input_opt['ID']]=self.table[self.input_opt['ID']].astype(str)
        if input=='' and table is None:
            util.error_msg('Clustering.__init__: Missing both input and table!')
        if type(self.input_opt['DESCRIPTION']) is str: self.input_opt['DESCRIPTION']=[self.input_opt['DESCRIPTION']]
        R_w=self.input_opt['EXP_WEIGHT']
        #R_w=input_options['EXP_WEIGHT']
        if R_w is not None and len(R_w)>0:
            R_w=util.sarray2rarray(R_w) # cast to np array
            if np.allclose(R_w, 1.0, atol=1e-5):
                input_options['EXP_WEIGHT']=None
            else:
                input_options['EXP_WEIGHT']=R_w
        R_w=self.input_opt['GENE_WEIGHT']
        #R_w=input_options['GENE_WEIGHT']
        if R_w is not None and len(R_w)>0:
            R_w=util.sarray2rarray(R_w) # cast to np array
            if np.allclose(R_w, 1.0, atol=1e-5):
                input_options['GENE_WEIGHT']=None
            else:
                input_options['GENE_WEIGHT']=R_w

    def set_input_options(self, options):
        self.input_opt.update(options)

    def set_cluster_options(self, options):
        self.cluster_opt.update(options)

    @staticmethod
    def _fix_missing(s_file):
        'This is no longer needed, fixed cwc on May 2, 2013'
        lines = []
        with open(s_file) as f:
            for line in f:
                lines.append(re.sub(r'(?<=\s)340282346638528859811704183484516925440.000', '', line))
        util.save_list(s_file, lines)

    @staticmethod
    def _strip_array_line(s_file):
        with open(s_file) as f:
            lines = f.readlines()
        s=lines[1]
        if s.startswith('AID'):
            del lines[1:2]
            util.save_list(s_file, lines)
            return s
        return ''

    @staticmethod
    def _insert_array_line(s_file, s):
        with open(s_file) as f:
            lines = f.readlines()
        lines[1:0]=[s]
        util.save_list(s_file, lines)

    @staticmethod
    def _remove_temp_files(s_file):
        S_suffix=[".log",".cdt",".gtr",".atr",".jtv",".kdt",".hexp2",".hgene",".dmg",".dme"]
        S_file = glob.glob(s_file+".*")
        for s_file in S_file:
            s_name, s_ext = os.path.splitext(s_file)
            if s_ext in S_suffix:
                os.remove(s_file)

    @staticmethod
    def _remove_extra_files(s_file):
        S_suffix=[".cdt",".gtr",".atr",".kdt",".jtv",".input"]
        S_file = glob.glob(s_file+".*")
        for s_file in S_file:
            s_name, s_ext = os.path.splitext(s_file)
            if s_ext not in S_suffix:
                os.remove(s_file)

    @staticmethod
    def restore_distance(s_file, max_dist=1.0):
        'No longer needed, as we use cwc.new'
        df=pd.read_table(s_file, header=None)
        s_col=util.header(df)[-1]
        R=(-df[s_col]+1.0).abs()*max_dist
        if R.max()<=1.0:
            df[s_col]=util.rarray2sarray((-R+1.0).abs(), s_format='%.3f') # convert to similarity score and output
            df.to_csv(s_file, sep='\t', index=False, header=False)
        else:
            util.warn_msg('Cannot restore distance to similarity, as the max distance '+str(R.max())+' > 1.0, restore skipped!')

    @staticmethod
    def make_JTV(s_file, r_max=2.0):
        s='''<DocumentConfig>
    <UrlExtractor/>
    <ArrayUrlExtractor/>
    <Views>
        <View type="Dendrogram" dock="1">
            <ColorExtractor contrast="%s">
                <ColorSet up="#D8181C" zero="#D8D8D8" down="#3A6C9A" missing="#D8D8D8"/>
            </ColorExtractor>
            <ArrayDrawer/>
            <GlobalXMap current="Fill">
                <FixedMap type="Fixed" scale="7.0"/>
                <FillMap type="Fill"/><NullMap type="Null"/>
            </GlobalXMap>
            <GlobalYMap current="Fill">
                <FixedMap type="Fixed" scale="11.0"/>
                <FillMap type="Fill"/>
                <NullMap type="Null"/>
            </GlobalYMap>
            <ZoomXMap current="Fill">
                <FixedMap type="Fixed"/>
                <FillMap type="Fill"/>
                <NullMap type="Null"/>
            </ZoomXMap>
            <ZoomYMap current="Fill">
                <FixedMap type="Fixed"/>
                <FillMap type="Fill"/>
                <NullMap type="Null"/>
            </ZoomYMap>
            <TextView>
                <TextView>
                    <GeneSummary/>
                </TextView>
                <TextView>
                    <GeneSummary/>
                </TextView>
                <TextView>
                    <GeneSummary/>
                </TextView>
                <TextView>
                    <GeneSummary/>
                </TextView>
            </TextView>
            <ArrayNameView>
                <ArraySummary included="0"/>
            </ArrayNameView>
            <AtrSummary/>
            <GtrSummary/>
        </View>
    </Views></DocumentConfig>''' % (str(r_max))
        util.save_list(s_file+".jtv", s)

    @staticmethod
    def optimize(s_file):
        if re.search(r'\.cdt$', s_file) is not None:
            s_file=re.sub(r'\.cdt$', '', s_file)
        S_cmd=[Clustering.BIN_OPT, "-m O -d P "+s_file, s_file+"Opt"]
        print(util.unix(" ".join(S_cmd)))
        #print S_cmd
        #Clustering._fix_missing(s_file+"Opt.cdt")
        Clustering._remove_temp_files(s_file)

    def kmeans(self, k, s_dm="", i_iteration=1, l_cleanup=False):
        if self.input=='': util.error_msg('Clustering.kmeans: Input file has not been prepared; use make_input() first!')
        S_cmd=[Clustering.BIN_CWC, "-k -E", "-c "+str(k), "-n "+str(i_iteration), "-i "+self.input+".input", "-o "+self.input]
        S_cmd.append("-p" if s_dm=="" else "-d -dm "+s_dm)
        util.unix(" ".join(S_cmd))
        Clustering._strip_array_line(self.input+".kdt")
        #Clustering._fix_missing(self.input+".kdt")

    def get_default_exp_cols(self, S_col=None):
        'Filter S_col for numerical columns that are not in key columns (ID, Description)'
        S_col = S_col or []
        S_key=set([self.input_opt['ID']]+self.input_opt['DESCRIPTION']+[self.input_opt['WEIGHT_COL']])
        if len(S_col)==0: S_col=util.header(self.table)
        # filter out key columns
        S_col=[s for s in S_col if s not in S_key]
        # filter out non-numeric columns
        return [s for s in S_col if (self.table[s].dtype is not np.dtype(object))]

    def make_DM(self, S_col=None, R_weight=None, metrics='PEARSON', by='GENE', l_normalize=True):
        S_col = S_col or []
        if self.table is None: util.error_msg('Clustering.make_input: missing Clustering.table!')
        S_col=self.get_default_exp_cols(S_col)
        n=len(self.table) if by=='GENE' else len(S_col)
        if n==0: util.error_msg('Clustering.make_input: no data record to cluster!')
        if len(S_col)==0:
            util.error_msg('Clustering.make_input: no data column to cluster!')
        if len(S_col)<2 and metrics=='PEARSON':
            util.error_msg('Clustering.make_input: not enough data column for Pearson!')
        N=DM.calc_N(n)
        R_DM=np.empty(N)
        l_gene= by=='GENE'
        M=self.table.reindex(columns=S_col).astype(float).values;
        if not l_gene:
            M=M.T
        cnt=0
        next_cnt=cnt+Clustering.COUNT_INTERVAL
        d_start=time.time()
        if not self.cluster_opt['HAS_NULL']: # and metrics=='PEARSON':
            if metrics=="PEARSON" or metrics=="BUILD_IN":
                if R_weight is None:
                    M_dist=np.maximum((1-np.corrcoef(M, rowvar=1))/2, 0)
                else:
                    R_weight/=np.sum(R_weight)
                    M-=(M*R_weight).sum(axis=1).reshape(n,1)
                    c=np.dot(M*R_weight, M.T)
                    d=np.diag(c)
                    M_dist=np.maximum((1-c/np.sqrt(np.multiply.outer(d,d)))/2, 0)
                R_DM=M_dist[np.triu_indices(n, k=1)]
                #for i in range(n):
                #    R_DM[cnt:(cnt+(n-i-1))]=M_dist[i,i+1:]
                #    cnt+=n-i-1

            elif metrics=="MANHATTAN":
                if R_weight is not None:
                    R_weight/=np.sum(R_weight)
                    M*=R_weight
                R_DM=skmp.manhattan_distances(M)[np.triu_indices(n, k=1)]
                #for i in range(n):
                #    R_DM[cnt:(cnt+(n-i-1))]=np.sum(np.abs(M[i+1:,:]-M[i]), axis=1)
                #    cnt+=n-i-1

            elif metrics=="EUCLIDEAN":
                cnt=0
                if R_weight is not None:
                    R_weight/=np.sum(R_weight)
                    M*=np.sqrt(R_weight)
                R_DM=skmp.euclidean_distances(M)[np.triu_indices(n, k=1)]
                #for i in range(n):
                #    R_DM[cnt:(cnt+(n-i-1))]=np.sqrt(np.sum((M[i+1:,:]-M[i])**2, axis=1))
                #    cnt+=n-i-1
            else:
                util.error_msg('Unsupported Metrics: %s!' % metrics)
        else:
            for i in range(n):
                R1=M[i]
                #print "i=%d" % i
                for j in range(i+1, n):
                    R2=M[j]
                    R_DM[cnt]=util.distance(R1, R2, R_weight=R_weight, metrics=metrics, has_null=self.cluster_opt['HAS_NULL'])
                    cnt+=1
                    if cnt>next_cnt:
                        i_pass=(time.time()-d_start)/60
                        print(" Distance Matrix: "+('%4.2f' % (cnt*100.0/N))+"% at row:"+str(i)+" "+('%.1f' % i_pass)+"min(s) passed, estimate total:"+str('%.1f' % (i_pass*N/cnt))+"min(s)\r");
                        next_cnt+=Clustering.COUNT_INTERVAL

        #print ">>>>>>>>>", time.time()-d_start
        dm=DM(R_DM=R_DM)
        if l_normalize: dm.normalize()
        return dm

    def make_table(self):
        opt=self.input_opt
        self.table=pd.read_table(self.input+".input")
        self.input_opt['EXP_WEIGHT']=self.table.iloc[0][3:].astype(float).values
        self.table=self.table.drop([0], axis=0)
        opt['DATA_COLS']=util.header(self.table)[3:]

    def auto_center(self, S_col=None):
        S_col = S_col or []
        if self.table is None: self.make_table()
        if len(S_col)==0:
            S_col=self.get_default_exp_cols(opt['DATA_COLS']) if len(self.input_opt['DATA_COLS'])==0 else self.input_opt['DATA_COLS']
        t=self.table[S_col].astype(float)
        Rm=t.mean(axis=1)
        Rs=t.std(axis=1)
        self.table[S_col]=t.sub(Rm, axis=0).div(Rs, axis=0)

    def make_input(self, s_file='untitled', options=None):
        if self.table is None: util.error_msg('Clustering.make_input: missing Clustering.table!')
        S=self.table.header()
        S_up=[ s.upper() for s in S]
        opt=self.input_opt
        opt.update(options or {})
        self.input_opt=opt
        S_miss=[s for s in opt['DATA_COLS'] if S.index(s)<0]
        if len(S_miss)>0: util.error_msg('Clustering.make_input: missing data column: '+", ".join(S_miss))
        i_id=util.index(opt['ID'], S)
        if (i_id<0):
            i_id=S_up.index('GENE')
            if i_id<0: util.error_msg('Clustering.make_input: no column is specified as the ID!')
            opt['ID']=S[i_id]
        if type(opt['DESCRIPTION']) is str: opt['DESCRIPTION']=[opt['DESCRIPTION']]
        I_des=[util.index(s, S) for s in opt['DESCRIPTION'] if util.index(s, S)>=0]

        if (len(I_des)==0):
            I_des=[i_id]
            opt['DESCRIPTION']=[opt['ID']]
        else:
            for i in I_des:
                self.table.iloc[:, i]=util.sarray2sarray(self.table.iloc[:,i])
        i_w=util.index(opt['WEIGHT_COL'], S)
        opt['DATA_COLS']=self.get_default_exp_cols(opt['DATA_COLS'])
        n_exp=len(opt['DATA_COLS'])
        if n_exp==0: util.error_msg('Clustering.make_input: no data column is specified!')

        S_out=[]
        S_out.append('Gene\tDescription\tWeight\t'+'\t'.join(opt['DATA_COLS']))
        if opt['EXP_WEIGHT'] is None or len(opt['EXP_WEIGHT'])!=n_exp:
            S_out.append('Exp\t\t'+'\t1'*n_exp)
        else:
            S_out.append('Exp\t\t\t'+'\t'.join(util.rarray2sarray(opt['EXP_WEIGHT'], s_format='%g', s_null=1.0)))
        #df.fillna('', inplace=True)
        i_cols=[S.index(s) for s in opt['DATA_COLS']]
        if opt['GENE_WEIGHT'] is not None and len(opt['GENE_WEIGHT'])==len(self.table):
            if opt['WEIGHT_COL']=='':
                opt['WEIGHT_COL']='WEIGHT'
            self.table[opt['WEIGHT_COL']]=opt['GENE_WEIGHT']
        for i in range(len(self.table)):
            s=str(self.table.iloc[i, i_id])+'\t'+":".join(self.table.iloc[i, I_des])+'\t'+str(self.table.iloc[i, i_w] if i_w>=0 else 1)
            R=np.array([x for x in self.table.iloc[i,i_cols]])
            if opt['GENE_NORMALIZE'] and opt['NORMALIZE_METHOD']=='Z':
                valid=util.no_nan(R)
                if len(valid)>1:
                    R=(R-np.mean(valid))/np.std(R, ddof=1)
            s+='\t'+'\t'.join(['' if pd.isnull(x) else str(x) for x in R])
            S_out.append(s)
        if re.search(r'\.input$', s_file) is not None:
            s_file=re.sub(r'\.input$', '', s_file)
        util.save_list(s_file+".input", S_out, s_end='\n')
        self.input=s_file

    def hint(self):
        if self.input!='':
            print("Input file: "+self.input+".input")
        else:
            print("Please make input file before continue, use make_input() first!")
            return
        print("Cluster Genes? "+('Y' if self.cluster_opt['GENE'] else 'N'))
        print('Distance: '+(self.cluster_opt['DMG'] if self.cluster_opt['DMG']!='' else self.cluster_opt['GENE_METRICS']))
        print('Cluster Experiments? '+('Y' if self.cluster_opt['EXP'] else 'N'))
        print('Distance: '+(self.cluster_opt['DME'] if self.cluster_opt['DME']!='' else self.cluster_opt['EXP_METRICS']))
        print('Optimization? '+('Y' if self.cluster_opt['OPTIMIZE'] else 'N'))
        print('Clean up? '+('Y' if self.cluster_opt['CLEANUP'] else 'N'))
        print('EXP_WEIGHT? '+('Y' if self.input_opt['EXP_WEIGHT'] else 'N'))
        print('GENE_WEIGHT? '+('Y' if self.input_opt['GENE_WEIGHT'] else 'N'))

    def hierarchical(self, options=None):
        if self.input=='': util.error_msg('Clustering.hierachical: Input file has not been prepared; use make_input() first!')
        if self.table is None: self.make_table()
        opt=self.cluster_opt
        opt.update(options or {})
        self.cluster_opt=opt
        l_CWC=self.cluster_opt['BIN'] == 'CWC'
        if self.cluster_opt['FINGERPRINT'] and l_CWC:
            util.error_msg('Clustering.hierachical: fingerprint mode has to be used with hybrid binary, not CWC!')
        #l_CWC=False
        #XXXXXXXXXXXXXXX
        if l_CWC:
            S_cmd=[Clustering.BIN_CWC, "-h -a -E -P", "-i "+self.input+".input", "-o "+self.input]
        else:
            S_cmd=[Clustering.BIN_HYB, "-eis", "-i "+self.input+".input", "-o "+self.input]
            if self.cluster_opt['SKIP_DM']:
                S_cmd.append('-ctr')
        s_dme=opt['DME']
        s_dmg=opt['DMG']
        r_maxe=1
        r_maxg=1
        iopt=self.input_opt
        d_start=time.time()
        if opt['GENE']:
            if opt['GENE_METRICS']=='BUILD_IN' and opt['DMG']=='' and not opt['HAS_NULL']:
                S_cmd.append("-p")
            else:
                if opt['GENE_METRICS']=='BUILD_IN':
                    opt['GENE_METRICS']='PEARSON'
                if opt['DMG']=='':
                    R_w=self.input_opt['EXP_WEIGHT']
                    #R_w=R_w+np.random.randn(len(R_w))*0.001
                    if R_w is not None and np.allclose(R_w, 1.0, atol=1e-5): R_w=None
                    #print R_w
                    dmg=self.make_DM(S_col=iopt['DATA_COLS'], metrics=opt['GENE_METRICS'], R_weight=R_w, by='GENE')
                    dmg.save(s_file=self.input+'.dmg', s_format='%.2f')
                    opt['DMG']=self.input+'.dmg'
                else:
                    dmg=DM(s_file=opt['DMG'])
                r_maxg=dmg.dmax
                del dmg
                if l_CWC:
                    S_cmd.append("-dmg "+opt['DMG'])
                else:
                    if self.cluster_opt['FINGERPRINT']:
                        S_cmd.append("-f "+opt['DMG'])
                    else:
                        S_cmd.append("-d "+opt['DMG'])
        if opt['EXP']:
            if not l_CWC:
                #util.warn_msg('Clustering.hierachical: experiment clustering currently is only supported by CWC!')
                if opt['EXP_METRICS']=='BUILD_IN':
                    opt['EXP_METRICS']='PEARSON'
                if opt['DME']=='':
                    R_w=None
                    dme=self.make_DM(S_col=iopt['DATA_COLS'], metrics=opt['EXP_METRICS'], R_weight=R_w, by='EXP')
                    dme.save(s_file=self.input+'.dme', s_format='%.2f')
                    opt['DME']=self.input+'.dme'
                else:
                    dme=DM(s_file=opt['DME'])
                r_maxe=dme.dmax
                del dme
                S_cmd.append("-de "+opt['DME'])
            else:
                S_cmd.append("-eg" if opt['GENE'] else '-e')
                if opt['EXP_METRICS']=='BUILD_IN' and opt['DME']=='' and not opt['HAS_NULL']:
                    if "-p" not in S_cmd: S_cmd.append("-p")
                else:
                    if opt['EXP_METRICS']=='BUILD_IN':
                        opt['EXP_METRICS']='PEARSON'
                    if opt['DME']=='':
                        R_w=None
                        if (iopt['WEIGHT_COL']!='' and util.index(iopt['WEIGHT_COL'], self.table.header())>=0):
                            R_w=self.table[iopt['WEIGHT_COL']].values
                        if R_w is not None and np.allclose(R_w, 1, atol=1e-5): R_w=None
                        dme=self.make_DM(S_col=iopt['DATA_COLS'], metrics=opt['EXP_METRICS'], R_weight=R_w, by='EXP')
                        dme.save(s_file=self.input+'.dme', s_format='%.2f')
                        opt['DME']=self.input+'.dme'
                    else:
                        dme=DM(s_file=opt['DME'])
                    r_maxe=dme.dmax
                    del dme
                    S_cmd.append("-dme "+opt['DME'])
        # cwc sends standard message to error channel
        util.unix(" ".join(S_cmd), l_error=False, l_print=False)
        #### ZZZ
        print(" ".join(S_cmd))
        #Clustering._fix_missing(self.input+".cdt")
        #if opt['RESTORE_DISTANCE']:
        #    if opt['GENE'] and opt['GENE_METRICS']!='BUILD_IN': Clustering.restore_distance(self.input+".gtr", max_dist=r_maxg)
        #    if opt['EXP'] and opt['EXP_METRICS']!='BUILD_IN': Clustering.restore_distance(self.input+".atr", max_dist=r_maxe)

        if not opt['EXP']:
            # old CWC version will generate an AID row
            s_array=Clustering._strip_array_line(self.input+".cdt")
        if (opt['OPTIMIZE'] and opt['GENE']):
            Clustering.optimize(self.input)
            # optimization can handle Array line
            #if opt['EXP']: Clustering._insert_array_line(self.input+"Opt.cdt", s_array)
            Clustering.make_JTV(self.input+"Opt")
        else:
            Clustering.make_JTV(self.input)
        if opt['CLEANUP']:
            if opt['OPTIMIZE']:
                Clustering._remove_extra_files(self.input+"Opt")
            else:
                Clustering._remove_extra_files(self.input)

    def make_guide_DM(self, S_data_cols, S_guide_cols, data_weight=1.0, guide_weight=1.0, data_metrics='PEARSON', guide_metrics='PEARSON'):
        dm_data=self.make_DM(S_data_cols, metrics=data_metrics)
        dm_guide=self.make_DM(S_guide_cols, metrics=guide_metrics)
        return DM(R_DM=(dm_data.DM*data_weight+dm_guide.DM*guide_weight)/(data_weight+guide_weight))

    def hierarchical_guide(self, S_data_cols, S_guide_cols, data_weight=1.0, guide_weight=1.0, data_metrics='PEARSON', guide_metrics='PEARSON'):
        dm=self.make_guide_DM(S_data_cols, S_guide_cols, data_weight, guide_weight, data_metrics, guide_metrics)
        dm.save(self.input+".dmg")
        self.hierachical({'GENE':True, 'DMG':(self.input+".dmg")})

import scipy.cluster.hierarchy as clst
import fastcluster
class FastCluster:

    def __init__(self, data, S_col=None, S_row=None, S_description=None, Zr=None, Zc=None):
        """data is DataFrame or numpy 2d-array"""
        if type(data) is pd.DataFrame:
            if S_col is not None:
                data=data.loc[:, S_col].values
            else:
                S_col=data.header()
                data=data.values
            if S_row is None:
                S_row=[str(x) for x in range(data.shape[0])]
        self.data=data
        self.S_row=S_row
        self.S_description=self.S_row if S_description is None else S_description
        self.S_col=S_col
        self.Zr=Zr
        self.Zc=Zc

    def cluster(self, method='average', metric='euclidean', l_row=True, l_col=True):
        """
        https://github.com/cran/fastcluster/blob/master/src/python/fastcluster.py
        metric: euclidean, minkowski, cityblock, seuclidean, sqecuclidean
            cosine, hamming, jaccard, chebychev, canberra, braycurtis,
            mahalanobis, yule, matching, sokalmichener, dice, rogerstanimoto
            russelrao, sokasneath, kulsinski, USER

            correlation
        """
        if l_row:
            self.Zr=fastcluster.linkage(self.data, method=method, metric=metric, preserve_input=True)
            #left_dendrogram=clst.dendrogram(Zr, orientation='left')
        if l_col:
            self.Zc=fastcluster.linkage(self.data.T, method=method, metric=metric, preserve_input=True)
            #top_dendrogram=clst.dendrogram(Zc, orientation='top')

    def cluster_rc(self, method_r='average', metric_r='euclidean', method_c='average', metric_c='euclidean', l_row=True, l_col=True):
        if l_row:
            self.Zr=fastcluster.linkage(self.data, method=method_r, metric=metric_r, preserve_input=True)
            #left_dendrogram=clst.dendrogram(Zr, orientation='left')
        if l_col:
            self.Zc=fastcluster.linkage(self.data.T, method=method_c, metric=metric_c, preserve_input=True)


    @staticmethod
    def quick_plot(data, s_out, S_row, S_col, method='average', metric='euclidean', l_row=True, l_col=True, l_norm_row=False, l_pdf=False):
        fc=FastCluster(data, S_col=S_col, S_row=S_row, S_description=S_row)
        fc.cluster(method=method, metric=metric, l_row=True, l_col=True)
        s_out, s_ext=os.path.splitext(s_out)
        fc.plot(s_out+".png", colormap=None, l_pdf=l_pdf)
        fc.save(s_out, l_norm_row=l_norm_row)

    def plot(self, s_imgfile, colormap=None, row_labels_size=0, col_labels_size=14, l_pdf=False, l_normalize_for_color=True, l_legend_pvalue=False):
        import pydendroheatmap as pdh
        heatmap=pdh.DendroHeatMap(heat_map_data=self.data, row_labels=self.S_description, col_labels=self.S_col, left_dendrogram=self.Zr, top_dendrogram=self.Zc, row_labels_size=row_labels_size, col_labels_size=col_labels_size, l_normalize_for_color=l_normalize_for_color, l_legend_pvalue=l_legend_pvalue)
        if colormap is None:
            heatmap.colormap=heatmap.color_brewer(brewer_name='Oranges', map_type='sequential', number=3, reverse=False)
        else:
            heatmap.colormap=colormap
        heatmap.export(s_imgfile, l_pdf)

    def save(self, s_cdtfile, l_norm_row=False, r_max=2.0):
        """r_max controls the max matrix value to be color saturated"""
        X=self.data
        s_cdtfile, s_ext=os.path.splitext(s_cdtfile)
        S_row=self.S_row
        if self.Zr is not None:
            #den_r=clst.dendrogram(self.Zr)
            den_r=FastCluster.linkage2order(self.Zr)
            #X=X[den_r['leaves'], :]
            X=X[den_r, :]
            S=[] #"NODEID\tLEFT\RIGHT\tCORRELATION"]
            r_dist=max(self.Zr[:, 2].max(), 1.0) if len(self.Zr)>0 else 0.0
            n=X.shape[0]
            node_cnt=0
            S_gene=["GENE%dX" % (x+1) for x in den_r] #['leaves']]
            S_row=[self.S_row[i] for i in den_r] #['leaves']]
            S_description=[self.S_description[i] for i in den_r]#['leaves']]
            for i,R in enumerate(self.Zr):
                node_cnt+=1
                s_left="GENE%dX" % int(R[0]+1) if int(R[0])<n else "NODE%dX" % (int(R[0]-n+1))
                s_right="GENE%dX" % int(R[1]+1) if int(R[1])<n else "NODE%dX" % (int(R[1]-n+1))
                S.append("NODE%dX\t%s\t%s\t%.4f" % (node_cnt, s_left, s_right, max(1.0-R[2]/r_dist, 0.0)))
                util.save_list(s_cdtfile+'.gtr', S, s_end="\n")
        S_col=self.S_col
        if self.Zc is not None:
            #den_c=clst.dendrogram(self.Zc)
            den_c=FastCluster.linkage2order(self.Zc)
            X=X[:, den_c]#['leaves']]
            S=[]
            r_dist=max(self.Zc[:, 2].max(), 1.0) if len(self.Zc)>0 else 0.0
            n=X.shape[1]
            node_cnt=0
            S_array=["ARRY%dX" % (x+1) for x in den_c]#['leaves']]
            S_col=[self.S_col[i] for i in den_c]#['leaves']]
            for i,R in enumerate(self.Zc):
                node_cnt+=1
                s_left="ARRY%dX" % int(R[0]+1) if int(R[0])<n else "NODE%dX" % (int(R[0]-n+1))
                s_right="ARRY%dX" % int(R[1]+1) if int(R[1])<n else "NODE%dX" % (int(R[1]-n+1))
                S.append("NODE%dX\t%s\t%s\t%.4f" % (node_cnt, s_left, s_right, max(1.0-R[2]/r_dist, 0.0)))
                util.save_list(s_cdtfile+'.atr', S, s_end="\n")
        n_exp=len(S_col)
        S=["GID\tGENE\tNAME\tGWEIGHT\t"+"\t".join(S_col)]
        if self.Zc is not None and len(self.Zc):
            S.append("AID\t\t\t\t"+"\t".join(S_array))
        S.append("EWEIGHT\t\t\t"+"\t1"*n_exp)
        for i,R in enumerate(X):
            if l_norm_row:
                R=(R-R.mean())/R.std()
            S.append(S_gene[i]+"\t"+S_row[i]+"\t"+S_description[i]+"\t1\t"+"\t".join(util.rarray2sarray(R, s_format="%.3f")))
        util.save_list(s_cdtfile+".cdt", S, s_end="\n")
        import cluster
        Clustering.make_JTV(s_cdtfile, r_max=r_max)

    @staticmethod
    def linkage2order(Z, M=None):
        """Convert n-1 X 4 linkage matrix to row order, if M is provided, rows with larger sum are place first"""
        r,c=Z.shape
        if r==0: return [0] #only one row
        n=r+1
        X={}
        if M is None:
            for i in range(r):
                left=int(Z[i,0])
                right=int(Z[i,1])
                if left>=n:
                    left=X.pop(left)
                else:
                    left=[left]
                if right>=n:
                    right=X.pop(right)
                else:
                    right=[right]
                left.extend(right)
                X[(n+i)]=left
            return X[n+r-1]
        else:
            # sort so that large-sum rows on left
            for i in range(r):
                left=int(Z[i,0])
                right=int(Z[i,1])
                if left>=n:
                    left=X.pop(left)
                else:
                    left=[1, M[left,:].sum(), [left]]
                if right>=n:
                    right=X.pop(right)
                else:
                    right=[1, M[right,:].sum(), [right]]
                if (left[1]*1.0/left[0])>(right[1]*1.0/right[0]):
                    left[2].extend(right[2])
                    left[0]+=right[0]
                    left[1]+=right[1]
                    X[(n+i)]=left
                else:
                    right[2].extend(left[2])
                    right[0]+=left[0]
                    right[1]+=left[1]
                    X[(n+i)]=right
            return X[n+r-1][2]

def cluster_array_to_k_groups(R, k):
    Z=fastcluster.linkage(R, method='average', metric='euclidean', preserve_input=True)
    import tree
    tr=tree.Tree(Z=Z)
    X=tr.representatives(n_picks=k, l_keep_members=True)
    return X

def ostu(R):
    # https://github.com/scikit-image/scikit-image/blob/v0.13.1/skimage/filters/thresholding.py#L231
    # class probabilities for all possible thresholds
    IDX=np.argsort(R)
    R=R[IDX]
    n=len(R)
    weight1 = np.arange(1., n)
    weight2 = np.arange(n-1, 0, -1)
    tot=np.sum(R)
    # class means for all possible thresholds
    mean1 = np.cumsum(R)[:-1]
    mean2 = tot-mean1
    mean1/=weight1
    mean2/=weight2
    # Clip ends to align class 1 and class 2 variables:
    # The last value of `weight1`/`mean1` should pair with zero values in
    # `weight2`/`mean2`, which do not exist.
    variance12 = weight1*weight2 * (mean1-mean2) ** 2
    idx = np.argmax(variance12)
    threshold = (R[idx]+R[idx+1])/2
    # a bool mask, True foreground, False background
    mask=np.array(IDX>idx)
    return (mask,threshold)

if __name__ == '__main__':
    R=np.random.rand(10,1)
    print(R)
    cluster_array_to_k_groups(R, 2)
    exit()
    s_input="test/clustering.input"
    c=Clustering(input=s_input)
    c.hierachical(options={'EXP':True, 'GENE_METRICS':'PEARSON', 'OPTIMIZE':True, 'CLEANUP':False})

    df=pd.read_csv('test/heatmap.csv')
    c=Clustering(table=df)
    c.make_input(s_file='test/heatmap')
    c.hint()
    c.hierarchical(options={'EXP':True, 'GENE_METRICS':'PEARSON', 'EXP_METRICS':'PEARSON', 'OPTIMIZE':False, 'CLEANUP':False})
    #Clustering.hierachical_table(df, s_file='test/fromT')
