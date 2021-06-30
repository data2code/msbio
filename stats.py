#!/usr/bin/env python
from __future__ import absolute_import
from __future__ import print_function
import scipy.stats as ss
import numpy as np
import pandas as pd
import util
import scipy.linalg as la
import math
import parallel
from six.moves import range

def z2p(z, two_tail=True):
    x=2.0 if two_tail else 1.0
    return ss.norm.sf(z)*x

def p2z(x, two_tail=True):
    if two_tail: x/=2
    return ss.norm.isf(x)

def lnfactorial(n):
    return math.lgamma(n+1)

def lnbinomial(n, k):
    return math.lgamma(n+1)-math.lgamma(n-k+1)-math.lgamma(k+1)

def hyper(n, N, n1, n2, tolerance=1e-300, n_chunk=1000):
    '''N M_total: total number of objects in bin
    n1 n_white: total number of white objects in bin
    n2 N_pick: number of draws without replacement
    n x_white: x out of N_pick are white'''
    n_chunk=1000 if n_chunk<=0 else n_chunk
    min_idx=min(n1,n2)
    l_left=n*1.0/n1 < n2*1.0/N and n<min_idx-n+1
    term=1.0
    P=0.0 if l_left else 1.0 #when l_left, do not include pvalue2(N,n1,n2,n) itself
    if l_left:
        for x in range(n-1,-1,-n_chunk):
            ## vectorize in chunks of 1000
            ## in case N is huge, we stop when the remaining area is too small
            if term*(x+1)<tolerance: break # no need to run, too small already
            X=np.arange(x, max(-1, x-n_chunk), -1.0)
            X=(X+1)*(N-n1-n2+1+X)/(n1-X)/(n2-X)
            X=X.cumprod()*term
            term=X[-1]
            P+=X.sum()
    else:
        for x in range(n+1, min_idx+1, n_chunk):
            if term*(min_idx-x+1)<tolerance: break
            X=np.arange(x, min(min_idx+1, x+n_chunk), 1.0)
            X=(1+n1-X)*(1+n2-X)/X/(X+N-n1-n2)
            X=X.cumprod()*term
            term=X[-1]
            P+=X.sum()
    P*=math.exp(lnbinomial(n2,n)+lnbinomial(N-n2,n1-n)-lnbinomial(N,n1))
    return 1.0-P if l_left else P

def hyper_previous(n, N, n1, n2):
    '''N M_total: total number of objects in bin
    n1 n_white: total number of white objects in bin
    n2 N_pick: number of draws without replacement
    n x_white: x out of N_pick are white'''
    min_idx=min(n1,n2)
    l_left=n*1.0/n1 < n2*1.0/N and n<min_idx-n+1
    term=1.0
    P=0.0 if l_left else 1.0 #when l_left, do not include pvalue2(N,n1,n2,n) itself
    if l_left:
        for x in range(n-1,-1,-1):
            term*=(x+1.0)*(N-n2-n1+x+1.0)/(n1-x)/(n2-x)
            P+=term
    else:
        for x in range(n+1, min_idx+1):
            term*=(n1-x+1.0)*(n2-x+1.0)/x/(N-n2-n1+x)
            P+=term
    P*=math.exp(lnbinomial(n2,n)+lnbinomial(N-n2,n1-n)-lnbinomial(N,n1))
    return 1.0-P if l_left else P

# too slow, replace by my own implementation above
def hyper_(x_white, M_total, n_white, N_pick):
    '''M_total: total number of objects in bin
    n_white: total number of white objects in bin
    N_pick: number of draws without replacement
    x_white: x out of N_pick are white'''
    return ss.hypergeom.sf(x_white-1, M_total, n_white, N_pick)

def ZScore_GeneGo(n, N, n1, n2):
    """Each subnetwork is associated with a Z-score which ranks the subnetworks according to saturation with the objects from the initial gene list. The Z-score ranks the subnetworks of the analyze network algorithm with regard to their saturation with genes from the experiment. A high Z-score means the network is highly saturated with genes from the experiment. The formula for the Z-score is:
    n (white in picked), N (total), n1 (white), n2 (picked, my gene list)
    The standard deviation of hypergeometric distribution comes from
    https://en.wikipedia.org/wiki/Hypergeometric_distribution
    where Ki->n1, N->N, Xi=n, n->n2, var=n1/N*(1-n1/N)*n2*(N-n2)/(N-1), Sigma=Sqrt(var)
    Z = (n-n1*n2/N)/Sqrt(n1*n2/(1-n2/N)(1-n1/N)/(N-1)))
    Z = (n*N-n1*n2)/Sqrt(n1*n2*(N-n1)*(N-n2)/(N-1))
    notice this formula is symmetrical for n1 and n2"""
    r=math.sqrt(n1*n2*(N-n1)*(N-n2)/(N-1))
    if r<1e-100: return 0.0
    return (n*N*1.0-n1*n2)/r

def adjust_p(R_p, N=None, method="BH"):
    """Calculate FDR for multiple test. N is the total # of tests run, if not given, set to len(R_p)
    R_p: an array of p-values
    N: int, total number of tests run
    method: currently fixed to Benjamini and Hochberg method.
    Output has been validated with adjust.p in R"""
    l_old=False # old implementation, slower, keep in case there is bug, new code has been tested
    N=len(R_p) if N is None else N
    if method.upper()=="BONFERRONI":
        return np.clip(np.array(R_p)*N, 0.0, 1.0)
    elif method.upper()=="HOLM":
        n=len(R_p)
        t=pd.DataFrame({'p':R_p, 'q':R_p, 'I':list(range(len(R_p)))})
        t.sort_values('p', ascending=True, inplace=True)
        t.index=range(n)
        if l_old:
            q=0.0
            for i in range(n):
                q=t.loc[i, 'q']=min(max(q, t.loc[i, 'p']*(N-i)),1)
        else:
            q=np.clip(t.p.values*(N-np.arange(n)), 0.0, 1.0)
            q=np.maximum.accumulate(q)
            t['q']=q
        t.sort_values('I', inplace=True)
        return t.q.values
    elif method.upper() in ("BH","FDR"):
        n=len(R_p)
        t=pd.DataFrame({'p':R_p, 'q':R_p, 'I':list(range(len(R_p)))})
        t.sort_values('p', ascending=False, inplace=True)
        t.index=range(n)
        if l_old:
            q=1.0
            for i in range(n):
                q=t.loc[i, 'q']=min(q, t.loc[i, 'p']*N*1.0/(len(t)-i))
        else:
            q=np.clip(t.p.values*N*1.0/(n-np.arange(n)), 0.0, 1.0)
            q=np.minimum.accumulate(q)
            t['q']=q
        t.sort_values('I', inplace=True)
        return t.q.values
    else:
        util.error_msg('Unsupported method: %s' % method)

def t_test_mean(R, m=0, two_tail=True):
    dict={'t':0, 'p':1}
    try:
        if not two_tail:
            if np.mean(R)>m:
                out=ss.ttest_1samp(R, m)
                dict={'t':out[0], 'p':out[1]/2}
        else:
            out=ss.ttest_1samp(R, m)
            dict={'t':out[0], 'p':out[1]}
    except Exception as e:
        print(e)
    return dict['p']

def empty_df(m, n):
    """ create an empty two-way df
        elements should be assigned np array later
    """
    return pd.DataFrame(np.zeros([m,n], dtype=object))

def rm_nan_df(t):
    """Remove nan from two-way dataframe"""
    m,n=t.shape
    for i in range(m):
        for j in range(n):
            R=t.iloc[i, j]
            t.iloc[i,j]=R[~np.isnan(R)]
    return t

def rm_nan_pair(R1, R2):
    """Only keep common non-nan elements, for paired test"""
    I=~np.isnan(R1) & ~np.isnan(R2)
    return (R1[I], R2[I])

def t_test(R1, R2, two_tail=True, paired=False):
    if paired:
        R=ss.ttest_rel(R1, R2)
    else:
        R=ss.ttest_ind(R1, R2)
    dict={'t':R[0], 'p':R[1]}
    if not two_tail: dict['p']/=2
    return dict['p']

def log10_P(p, MIN=1e-100):
    return math.log(np.maximum(p, MIN), 10)

def pearson(R1, R2):
    return np.corrcoef(R1, R2)[0,1]

def wilcoxon(R1, R2, two_tail=True): # paired non-parametric
    dict={'T':0, 'p':1}
    try:
        R=ss.wilcoxon(R1, R2)
        dict={'T':R[0], 'p':R[1]}
        if not two_tail: dict['p']/=2
    except Exception as e:
        print(e)
    return dict['p']

def anova(*args):
    dict={'F':0, 'p':1}
    try:
        R=ss.f_oneway(*args)
        dict={'F':R[0], 'p':R[1]}
    except Exception as e:
        print(e)
    return dict['p']

def kruskal(*args):
    dict={'H':0, 'p':1}
    try:
        R=ss.kruskal(*args)
        dict={'H':R[0], 'p':R[1]}
    except Exception as e:
        print(e)
    return dict['p']

def mannwhitney(R1, R2, two_tail=True):
    dict={'u':0, 'p':1}
    try:
        R=ss.mannwhitneyu(R1,R2)
        dict={'u':R[0], 'p':R[1]}
        if two_tail: dict['p']*=2
    except Exception as e:
        print(e)
    return dict['p']

def chi2(M):
    '''M is a numpy 2D array'''
    R=ss.chi2_contingency(M)
    return R[1]/2.0

def fisher_exact(M, alternative='greater'):
    """First row, foreground data (Y,N), second row, background data.
    greater mean foreground %Y > background %Y
    E.g., when alternative is 'greater'
    fisher_exact([[900, 1000], [500, 1000]]) returns 7.6486698436085204e-17
    fisher_exact([[10, 1000], [500, 1000]]) returns 0.99999999999819367
    stats.fisher_exact([[10, 1000], [500, 1000]], alternative='two-sided')
    """

    R=ss.fisher_exact(M, alternative=alternative)
    return R[1]

def make_contingency(S_row, S_col, S_all):
    '''     Col=Y  Col=N
    Row=Y
    Row=N
    '''
    M=np.zeros([2,2], dtype=int)
    s1=set(S_row)
    s2=set(S_col)
    s3=set(S_all)
    M[0,0]=len(s1.intersection(s2))
    M[0,1]=len(s1.difference(s2))
    M[1,0]=len(s2.difference(s1))
    M[1,1]=len(s3.difference(s1.union(s2)))
    return M

def chi2_by_matrix(M):
    s=M.sum()
    #print M[0,0], M[0,1], M[1,0], M[1,1]
    f1=(M[0,0]+M[0,1])*1.0/s
    f2=(M[0,0]+M[1,0])*1.0/s
    if f1<=f2 or M[0,0]==0:
        p_fisher=p_chi2=1.0
    else:
        p_fisher=fisher_exact(M)
        p_chi2=chi2(M)
    return {'p':np.min([p_fisher, p_chi2]), 'p_fisher':p_fisher, 'p_chi2':p_chi2, 'f1':f1, 'f2':f2, 'enrichment':f1/max(f2,0.5), 'M':M}

def chi2_by_lists(S_row, S_col, S_all):
    '''method is fisher or chi2'''
    M=make_contingency(S_row, S_col, S_all)
    return chi2_by_matrix(M)

def rank(R):
    return pd.core.algorithms.rank(R)

def Z_norm(R):
    R2=(rank(R)-0.5)/float(len(R))
    return ss.norm.ppf(R2)

def Z_norm2(X, window):
    """Parametic Z-norm, using a sliding window, X-mean/stdv
    For non-parametric, use Z_norm"""
    n=len(X)
    delta=window/2
    X2=X*X
    sumX=np.zeros(n)
    sumX2=np.zeros(n)
    N=np.zeros(n)
    iB0=iE0=0
    sumX[0]=x=X[0]
    sumX2[0]=x2=X[0]*X[0]
    N[0]=1
    for i in range(n):
        #print ">>", i
        iB=max(i-delta,0)
        iE=min(i+delta, n-1, iB+delta*2)
        x-=np.sum(X[iB0:iB])
        x+=np.sum(X[iE0+1:iE+1])
        x2-=np.sum(X2[iB0:iB])
        x2+=np.sum(X2[iE0+1:iE+1])
        #print ">>><<", iE-iB+1, "==", iE0-iB0+1-(iB-iB0)+(iE-iE0), x==np.sum(X[iB:iE+1])
        N[i]=iE-iB+1
        iB0=iB
        iE0=iE
        sumX[i]=x
        sumX2[i]=x2
    m=sumX/N
    std=np.sqrt(sumX2/N-m*m)
    return (X-m)/std

def quantile_norm(M):
    """Input a list of arrays of the same size"""
    from scipy.stats import rankdata
    m, n = np.size(M,0), np.size(M,1)
    M_sorted=np.zeros_like(M)
    M_rank=np.zeros_like(M,dtype='int')
    for i in xrange(n):
        M_sorted[:,i]=sorted(M[:,i])
        M_rank[:,i]=rankdata(M[:,i],method='min') # minimum rank for ties
    Quantiles=M_sorted.mean(axis=1)
    return Quantiles[M_rank-1]

def unique(R):
    return pd.core.algorithms.unique(R)

def rmsd(R):
    return np.std(R, ddof=1)

def anova2(T, l_interaction=False):
    '''Sokal 95 Chapter 11. Two-way analysis of variance'''
    b=len(T)
    a=len(util.header(T))
    n=len(T.values[0,0])
    if n==1 and l_interaction: return None;
    N=a*b*n;
    Ya=np.zeros(a);
    Yb=np.zeros(b);
    for i in range(a):
        for j in range(b):
            Ya[i]+=np.sum(T.values[j,i])
        Ya[i]/=n*b;

    for i in range(b):
        for j in range(a):
            Yb[i]+=np.sum(T.values[i,j])
        Yb[i]/=n*a

    Yab=np.mean(Ya)
    SSa=n*b*np.var(Ya)*a
    SSb=n*a*np.var(Yb)*b
    SSab=0.0
    for j in range(a):
        for i in range(b):
            SSab+=(np.mean(T.values[i,j])-Ya[j]-Yb[i]+Yab)**2
    SSab*=n;
    SSwn=0.0;
    if n>1:
        for j in range(a):
            for i in range(b):
                SSwn+=np.var(T.values[i,j])*n

    dfa=a-1.0;
    MSa=SSa/dfa;
    if l_interaction:
        df=a*b*(n-1.0)
        F=max(MSa/(SSwn/df+np.finfo(float).eps),0.0)
    else:
        df=(a-1.0)*(b-1.0)+a*b*(n-1.0)
        F=max(MSa/((SSwn+SSab)/df+np.finfo(float).eps),0.0)

    return ss.beta.cdf(df/(df+dfa*F),df*0.5,dfa*0.5)

def anova2_unbalanced(T):
    '''Unbalanced two-way ANOVA
    #First two way model
    # Y = B0+a1*A1+a2*A2+...+aa-1*Aa-1+b1*B1+b2*B2+...+bb-1*Bb-1
    # J.D. Jobson 91, page 458, Applied multivariate data analysis'''
    a=len(util.header(T))
    b=len(T);
    Y=list()
    for i in range(b):
        for j in range(a):
            Y.append(T.values[i,j])
    Y=np.hstack(Y)
    N=len(Y)
    Y=Y.reshape(N,1)
    X=np.zeros([N,a+b-1])
    cnt=0;
    for i in range(b):
        for j in range(a):
            for k in range(len(T.values[i,j])):
                X[cnt,0]=1 #B0
                if j<a-1: X[cnt, j+1]=1 #Aj
                if i<b-1: X[cnt, a+i]=1 #Bi
                cnt+=1
    # Jobson page 226
    C2=la.inv(np.dot(X.T, X))
    Xb=np.dot(X,np.dot(C2,np.dot(X.T,Y)))
    SSE1=np.dot((Y-Xb).T,(Y-Xb)).sum()
    dfE1=N-(a+b-1)
    # Now try drop the A factor from the model
    # Y=B0+b1*B1+b2*B2+...+bb-1*Bb-1;
    X=np.zeros([N,b])
    cnt=0
    for i in range(b):
        for j in range(a):
            for k in range(len(T.values[i,j])):
                X[cnt,0]=1 #B0
                if i<b-1: X[cnt,i+1]=1 #Bi
                cnt+=1
    C2=la.inv(np.dot(X.T,X))
    Xb=np.dot(X,np.dot(C2,(np.dot(X.T,Y))))
    SSE2=np.dot((Y-Xb).T,(Y-Xb)).sum()
    dfE2=N-b
    df=dfE2-dfE1
    F=max(SSE2-SSE1,0.0)/df/((SSE1+np.finfo(float).eps)/dfE1)
    return ss.beta.cdf(dfE1/(dfE1+df*F),dfE1*0.5,df*0.5)

def mackskillings(T):
    '''Non-parametric multi-group test with multiple measurements
    # T is a table, each row is a block (a probe)
    # each column is a treatment
    # each cell is a Rarray()
    # See Hollander & Wolfe, Nonparametric statistical methods'''
    n=len(T);
    c=len(T.values[0,0]); # nof measurements, assume all share the same # of measurements
    k=len(util.header(T))
    S=np.zeros(k)
    N=n*c*k;
    T2=T.copy();

    Rc=None
    for i in range(n):
        Rr=rank(np.hstack(T.values[i,:]))
        for j in range(k):
            T2.values[i,j]=Rr[j*c : (j+1)*c]
    for j in range(k):
        for i in range(n):
            S[j]+=np.sum(T2.values[i,j])/c
    MS=max(np.sum((S-(N+n)/2.0)**2)*12.0/k/(N+n), 0.0)
    return ss.gamma.sf(MS/2, (k-1.0)/2.0)

def mackskillings2(T):
    '''Two-way Mack-Skillings
    # T is a table, each row is a block (a probe)
    # each column is a treatment
    # each cell is a Rarray()
    # See Hollander & Wolfe, Nonparametric statistical methods'''
    n=len(T)
    k=len(util.header(T))
    T2=T.copy()
    for i in range(n):
        Rr=rank(np.hstack(T.values[i,:]))
        cnt=0
        for j in range(k):
            T2.values[i,j]=Rr[cnt:cnt+len(T.values[i,j])];
            cnt+=len(T.values[i,j]);
    q=np.zeros(n)
    for i in range(n):
        for j in range(k):
            q[i]+=len(T2.values[i,j])

    V0=np.zeros([k-1,k-1])
    for s in range(k-1):
        for t in range(s, k-1):
            if s==t:
                for i in range(n):
                    cis=len(T2.values[i,s]);
                    V0[s,t]+=cis*(q[i]-cis)*(q[i]+1)/12.0/q[i]/q[i];
            else:
                for i in range(n):
                    cis=len(T2.values[i,s])
                    cit=len(T2.values[i,t])
                    V0[s,t]-=cis*cit*(q[i]+1)/12.0/q[i]/q[i];
                V0[t,s]=V0[s,t];

    V=np.zeros([1, k-1])
    for j in range(k-1):
        V[0,j]=0.0;
        for i in range(n):
            V[0,j]+=np.sum(T2.values[i,j])/q[i];
            V[0,j]-=len(T2.values[i,j])*(q[i]+1)/2.0/q[i];
    MS=np.dot(np.dot(V, la.inv(V0)), V.T).sum();
    return ss.gamma.sf(MS/2, (k-1.0)/2.0)

def RSA_rank(R_score, I_index):
    # I_index contains index for postives, when R_score contains ties, I_index is not the same as I_rank
    # assume R_score has been sorted, so that I_rank will be set to the last score, if ties are found
    # I_index is zero-based, not starting from one
    I_rank=I_index[:]
    I_rank.sort()
    n=len(I_rank)
    N=len(R_score)
    for i in range(n-1,-1,-1):
        idx0=I_rank[i];
        idx=idx0+1
        while (idx<N and R_score[idx]==R_score[idx0] and (i==n-1 or idx<I_rank[i+1])):
            # make sure I_Rank[i[<I_Rank[i+1]
            I_rank[i]=idx
            idx+=1
    return I_rank

def RSA_score(I_rank,N,i_min=None,i_max=None,l_BonferroniCorrection=False,tolerance=1e-100, l_keep_most=False, p_cutoff=0.01):
    """l_keep_most, as long as corrected p-value is < p_cutoff, we aim to keep the longest hit list
    p_cutoff, only used when l_keep_most is True"""
    #I_rank, zero-based
    # i_min, i_max also zero-based
    cutoff=0
    logP_min=1.0
    n=len(I_rank);
    if i_max is None: i_max=n-1
    if i_min is None: i_min=0
    i_scale=(i_max-i_min+1) if l_BonferroniCorrection else 1
    for i in range(i_max, i_min-1, -1):
    #for i in xrange(i_min,i_max+1):
        #print i, N, n, I_rank[i], i_scale
        #print "%g" % (hyper(i+1, N, n, I_rank[i]+1)*i_scale)
        #print ">>>>>>>", hyper(i+1, N, n, I_rank[i])*i_scale
        # P can be > 1 after BonferroniCorrection
        if (i<i_max) and (I_rank[i]==I_rank[i+1]): continue
        # a rare bug, where I_rank={0,23,23}, N=24, i_min=1, i_max=3, i=2
        # we will call hyper(3, 24, 2, 24)
        #print ">>>", i, i+1, N, n, I_rank[i]+1
        # logP=min(math.log(max(hyper(i+1, N, n, I_rank[i]+1)*i_scale, tolerance), 10), 0) YaZ 11132019: clip after changing cutoff to fix wrong #hitWell
        logP=math.log(max(hyper(i+1, N, n, I_rank[i]+1)*i_scale, tolerance), 10)
        if (logP < logP_min):
            logP_min=logP
            cutoff=i
        logP=min(logP, 0)
        if l_keep_most and logP_min<=np.log10(p_cutoff):
            break # good enough
    return {'logP':logP_min, 'cutoff':cutoff}

def RSA(T, s_gene="GeneID", s_score="Score", l_reverse=False, LB=0.2, UB=0.8, l_randomize=False, l_BonferroniCorrection=False, N_total=None):
    """N_total: total number of genes in the screen, defaults to len(T)
        This is introduced, beause in HTS screen, we may have 1 million compounds, N=1000000, but we only need say 10000 active compounds
        and those inactive compounds (say 20000) that are analogs to actives in our T, so len(T) is 30000.
        By providing N_total=1000000, we do not have to provide all the entries in T"""
    t=T.copy()
    S=util.header(t)
    if t[s_gene].dtype is not np.dtype(object):
        t[s_gene]=t[s_gene].astype(str)
    t=t[ (pd.notnull(t[s_gene])) & (t[s_gene]!="") & (pd.notnull(t[s_score])) ]
    N=len(t)
    if N_total is None: N_total=N
    R_logP=np.zeros(N)
    R_bestActivity=np.zeros(N)
    I_hit=np.zeros(N).astype(int)
    I_totWell=np.zeros(N).astype(int)
    I_hitWell=np.zeros(N).astype(int)
    if l_randomize:
        R=t[s_score].values
        R=R[np.random.permutation(len(R))]
        t[s_score]=R
    t.sort_values(s_score, ascending=(not l_reverse), inplace=True)
    c_gene=dict()
    c_rank=dict()
    # we need to hash the max rank of a given score.
    # if t is a membership matrix, there are lots of ties, obtaining
    # c_rank can be the bottleneck
    c_score=dict()
    R_score=t[s_score].values
    for i in range(N):
        c_score[R_score[i]]=i

    for i in range(len(t)):
        s=t.iloc[i][s_gene]
        if s not in c_gene:
            c_gene[s]=[] # store the exact index for this gene
            c_rank[s]=[] # modify the rank, in case there are ties
        # updated on 10/19/2012
        c_gene[s].append(i)
        # the following can be the slowest part, if there are lots of ties
        #for j in xrange(i+1, N):
        #    if t[s_score].iloc[j]!=t[s_score].iloc[i]: break
        #c_rank[s].append(j-1)
        c_rank[s].append(c_score[R_score[i]])
    for s in c_gene:
        #if s!='1221200': continue
        I_rank=c_rank[s]
        i_max=None
        i_min=None
        for k in range(len(I_rank)):
            if l_reverse:
                if R_score[I_rank[k]]>=LB: i_max=k
                if R_score[I_rank[k]]>=UB: i_min=k
                if (R_score[I_rank[k]]<LB and i_max is None): i_max=k-1
            else:
                if R_score[I_rank[k]]<=UB: i_max=k
                if R_score[I_rank[k]]<=LB: i_min=k
                if (R_score[I_rank[k]]>UB and i_max is None): i_max=k-1
        #print I_rank, N, i_min, i_max, l_BonferroniCorrection
        rslt=RSA_score(I_rank,N_total,i_min,i_max,l_BonferroniCorrection=l_BonferroniCorrection)
        #print rslt
        logP=rslt['logP']
        cutoff=rslt['cutoff']
        I_idx=c_gene[s]
        for k in range(len(I_idx)):
            R_logP[I_idx[k]]=logP
            R_bestActivity[I_idx[k]]=R_score[I_idx[0]]
            I_hitWell[I_idx[k]]=cutoff+1
            I_totWell[I_idx[k]]=len(I_idx)
            if (k<=cutoff): I_hit[I_idx[k]]=1

    t["LogP"]=R_logP
    t["BestActivity"]=R_bestActivity
    t["RSA_Hit"]=I_hit
    t["#hitWell"]=I_hitWell
    t["#totalWell"]=I_totWell
    # q-value
    t_p=t.drop_duplicates([s_gene,'LogP'])
    Rq=np.log10(adjust_p(np.power(10, t_p.LogP.clip(-np.inf, 0))))
    t_q=pd.DataFrame({s_gene: t_p[s_gene], 'Log_q':Rq})
    t=t.merge(t_q, left_on=[s_gene], right_on=[s_gene])
    ###
    t.sort_values(['LogP',s_gene,s_score], ascending=[True, True, not l_reverse], inplace=True)
    #t["LogP"]=util.rarray2sarray(t["LogP"], s_format='%.4f')
    t["RSA_Rank"]=np.zeros(N).astype(int)+999999
    cnt=0
    for k in range(N):
        if t["RSA_Hit"].values[k]>0:
            cnt+=1
            t["RSA_Rank"].values[k]=cnt
    return t

def kappa_stat(M, n_CPU=0):
    """M: 2d-np array, n x m, see http://david.abcc.ncifcrf.gov/helps/linear_search.html
    rows: genes (categories), cols: GOs (judges)
    M: membership int array of 1/0
    see Jacob Cohen, A coefficient of agreement for nominal scales, Educational & Psychological Measurement, Vol. 20 No. 1, 1960:37-46
    Note: chi2 tests the null hypothesis with regard to association, not agreement.
    When chi2 is significant, two judges are not associated, but it does not say they tend to agree or disagree.
        Oab: # cases they agree
        Aab: # cases they agree by chance
        K=(Oab-Aab)/(N-Aab) # portion of missing disagreement cases
    return D, similarity Kappa of m*(m-1)/2"""

    n,m=M.shape
    N=M.sum(axis=0)
    if m<=1: return []

    def K(ab):
        # ab: tuple of range lb, ub index
        out=[]
        for i in range(ab[0], ab[1]):
            Oab=(M[:, list(range(i+1,m))]-M[:,i:i+1]==0).sum(axis=0)
            Aab=(N[i]*N[i+1:]+(n-N[i])*(n-N[i+1:]))*1.0/n
            Kab=(Oab-Aab+1e-10)/(n-Aab+1e-10)
            out.append(Kab)
        return np.concatenate(out)

    if n_CPU<=1:
        return K((0, m-1))
    else: # split rows into approximately n_CPU chunks
        n_tot=m*(m-1)/2
        n_size=n_tot*1.0/n_CPU
        L=[]
        iB=cnt=0
        for i in range(m-1):
            cnt+=m-1-i
            if cnt>=n_size or i==m-2:
                L.append((iB, i+1))
                iB=i+1
                cnt=0
        mp=parallel.MP()
        mp.start(f=K, n_CPU=n_CPU)
        out=mp.map(L)
        return np.concatenate(out)
    return K

def calc_maxpct(s_ab,moa2cpd,i):
    Kab=[]
    for ma,mb in s_ab:
        a=moa2cpd[ma]
        b=moa2cpd[mb]
        Oab=len(set(a).intersection(set(b)))
        if Oab == 0:
            Kab.append(0)
            continue
        na=len(set(a))
        nb=len(set(b))
        Kab.append(max(Oab/na,Oab/nb))
    return (Kab,i)

def DM_maxpct(moa2cpd, n_cpus=1,mp=None,l_quit=True):
    """Compute maximum pairwise overlap ratio.
       moa2cpd: dict of list name to values
    """
    import itertools
    pairs=list(itertools.combinations(moa2cpd.index.tolist(),2))
    if n_cpus==1:
        return np.array(calc_maxpct(pairs,moa2cpd,0)[0])
    if mp is None:
        mp=parallel.parprep(n_CPU=n_cpus)
    L=[(calc_maxpct,x,moa2cpd,i) for i,x in enumerate(util.split(pairs,n_cpus))]
    out=mp.map(L,l_quit=l_quit)
    print(len(out),out[0][1])
    out={i:v for v,i in out}
    DM=np.concatenate([out[i] for i in range(n_cpus) if i in out])
    return DM

def find_linear_comb(M, I=None, I_DEL=None, c_comb=None, n_sample=0):
    """http://www.inside-r.org/packages/cran/caret/docs/findLinearCombo
    M: input 2d array
    I: labels for columns, defaults to X0, X1, ...
    I_DEL: column labels to be removed (previously found to be deleted columns, should be None if you call this from outside.
    c_comb: previously found linear dependencies, should be None, if you call from outside
    n_sample: sample n rows, if zero, keep all rows. This is to speed things up with a subsample
    return (lst, c_comb): list of column labels to be removed, dict of linear dependencies"""

    m,n=M.shape
    if n_sample>0 and n_sample<m:
        M=M[np.random.permutation(m)[:n_sample], :]
    I=["X%d" % i for i in np.arange(n)] if I is None else I
    I_DEL=[] if I_DEL is None else I_DEL
    c_comb={} if c_comb is None else c_comb
    U,s,V=np.linalg.svd(M)
    s/=np.absolute(s).max()
    #print U.shape, V.shape, M.shape, s.shape
    I_found=np.arange(n)[np.absolute(s)<=1e-6]
    if len(s)<len(V[:,0]): # extra NULL space
        I_found=np.concatenate([I_found, list(range(len(s), len(V[:,0])))])
    V=V.T
    del_cols=set()
    for i in I_found:
        j=np.argmax(np.absolute(V[:,i]))
        V[:,i]/=V[j,i]
        del_cols.add(j)
        c_comb[I[j]]="".join("".join(["%s%.3g*%s" % ("+" if V[x,i]>0 else "-", abs(V[x,i]), I[x]) for x in range(n) if abs(V[x,i])>1e-4]))+"=0"
    if len(del_cols):
        del_cols=np.sort(-np.array(list(del_cols)))
        for i in del_cols:
            I_DEL.append(I[-i])
            I=np.delete(I, -i)
            M=np.delete(M, -i, 1)
        I_DEL, c_comb = find_linear_comb(M, I, I_DEL, c_comb)
    return (list(np.sort(I_DEL)), c_comb)

def Z_factor(R_p, R_n, l_robust=True):
    """Calculate Z factor for given Positive and Negative numpy arrays"""
    # QC warning
    if l_robust:
        m_p=np.median(R_p)
        m_n=np.median(R_n)
    else:
        m_p=np.mean(R_p)
        m_n=np.mean(R_n)
    std_p=np.std(R_p)
    std_n=np.std(R_n)
    Z=1.0-3.0*(std_p+std_n)/max(abs(m_p-m_n), 1e-10)
    return Z

def Otsu_threshold(R):
    """Split R into two groups using Otsu's algorithm
    Own implementation that deal with nparray, rather than image array
    adopted from http://en.m.wikipedia.org/wiki/Otsu%27s_method, JS implementation"""
    R=np.copy(R)
    R.sort()
    N=len(R)
    sumA=t1=t2=0.0
    sumB=R.sum()
    cutoff1=cutoff2=None
    max_btw=-np.inf
    for i,r in enumerate(R):
        sumA+=r
        sumB-=r
        if i<len(R)-1 and R[i]==R[i+1]: continue
        mA=sumA/(i+1.0)
        mB=sumB/(N-i-1.0)
        btw=(sumA-sumB)*(sumA-sumB)*i*(N-i-1)/N/N
        if (btw>=max_btw):
            r=(R[i]+R[i+1])/2 if i<len(R)-1 else R[i]+1e-100
            if (btw>max_btw):
                cutoff1=cutoff2=r
                max_btw=btw
            else: # equal
                cutoff2=r
    return (cutoff1+cutoff2)/2

def gini(R):
    """Finally use http://www3.nccu.edu.tw/~jthuang/Gini.pdf"""
    R2 = np.sort(R) #sorted(list_of_values)
    n=len(R2)
    return (n+1-np.sum(np.linspace(n,1,n)*R2)*2.0/np.sum(R2))/n

if __name__ == '__main__':
    def df2str(df, s_format='%.4g', s_null=''):
        df2=df.copy()
        for s in util.header(df2):
            df2[s]=[",".join([s_format % x for x in X]) for X in df[s].values]
        return df2

    def str2df(df):
        df2=df.copy()
        for s in util.header(df2):
            df2[s]=[np.array(util.sarray2rarray(x.split(","))) for x in df[s].values]
        return df2

    S1=["a","b","c","d","f","g"]
    S2=["a","c","e","f","g"]
    S_all=["a","b","c","d","e","f","g","h","i","j","k","p","q","r","s","t"]
    #print(chi2_by_lists(S1, S2, S_all))
    #exit()
    import os.path
    if not os.path.isfile("test/test.csv"):
        dict={'c'+str(i):pd.Series(np.array([np.random.randn(3*10)+i/2.0]).reshape(10,3).tolist()) for i in range(3)}
        df=pd.DataFrame(dict)
        df2=df2str(df)
        print(df2)
        df2.to_csv('test/test.csv', index=False)
    df=str2df(pd.read_csv('test/test2.csv'))
    for i in range(2):
        df.values[i,1]+=7.0
    #print mackskillings(df)
    #print mackskillings2(df)
    #print anova2(df)
    #print anova2_unbalanced(df)
    exit()
    #R1=np.array([6.4,6.8,7.2,8.3,8.4,9.1,9.4,9.7])
    #R2=np.array([2.5,3.7,4.9,5.4,5.9,8.1,8.2])
    #R3=np.array([1.3,4.1,4.9,5.2,5.5,8.2])
    #print kruskal(R1, R2, R3)

    n=10
    M=np.zeros([20,n])
    M[:,0]=1
    M[:,1]=np.random.randn(20)
    M[:,2]=np.random.randn(20)
    M[:,3]=np.random.randn(20)
    M[:,4]=0.5*M[:,1]-0.25*M[:,2]-0.25*M[:,3]
    M[:5,5]=1
    M[5:10,6]=1
    M[10:20,7]=1
    I_del, c = find_linear_comb(M, I=["A","B","C","D","E","F","G","H","I","J"])
    print(I_del)
    print(c)

def mad(data, axis=0):
    """median absolute deviation, robust estimation of standard deviation
        http://www.programcreek.com/python/example/10046/numpy.median
        http://en.wikipedia.org/wiki/Median_absolute_deviation
    """
    #return np.median(np.abs(R-np.median(R)))/0.67449
    if axis == 0:
        demeaned = data - np.median(data, axis=0)
        return np.median(np.abs(demeaned), axis=0)*1.4826
    else:
        demeaned = data-np.median(data, axis=1).reshape(-1,1)
        return np.median(np.abs(demeaned), axis=1)*1.4826


