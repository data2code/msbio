#!/usr/bin/env python3

from __future__ import absolute_import
from __future__ import print_function
import numpy as np
import util
import pandas as pd
import sklearn as skl
from sklearn.model_selection import StratifiedKFold
import scipy.optimize as optz
from sklearn.metrics import auc, roc_curve
#from sklearn.metrics import average_precision_score
from sklearn.metrics import precision_recall_curve
import sys
import os
from six.moves import range
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

#sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)

class FindWeight():

    def __init__(self):
        self.W=None

    @staticmethod
    def norm_weight(Weights=None):
        if Weights is None:
            return (0.5, 0.5)
        else:
            tot=sum(Weights)
            w1=Weights[0]*1.0/tot
            w2=Weights[1]*1.0/tot
            return (w1, w2)

    @staticmethod
    def auto_weight(y):
        w1=sum(y)
        w0=len(y)-w1
        return FindWeight.norm_weight([1.0/w1, 1.0/w0])

    @staticmethod
    def f(W, X, y, penalty, Weights=None):
        """If Weights [Weight(y==1), Weight(y==1)] array is None, samples are equally weight,
        otherwise, the LR prob is no longer the real probability of the sample, needs to be corrected
        by passing the same Weights array to predict()"""
        #l_weighted=False # if we weight, then the logistic formula no longer predicts the probability
        W2=np.sum(W[1:]*W[1:])
        W=np.reshape(W,(-1,1))
        q=np.exp(np.clip(np.dot(X, W), -100, 100))
        q=np.clip(q/(1.0+q), 1e-15, 1-1e-15)
        # cross-entropy
        w1, w2 = FindWeight.norm_weight(Weights)
        return w1*np.sum(-np.log(q[y==1]))+w2*np.sum(-np.log(1-q[y==0]))+penalty*W2

    @staticmethod
    def accuracy(W, X, y):
        q=np.exp(np.dot(X, W))
        q=q/(1.0+q)
        y_pred=np.array(q>=0.5, dtype=int)
        n=len(y)
        pos=sum(y)
        neg=n-pos
        R_w=y*0.5/pos+(1-y)*0.5/neg
        #print pos, neg, sum(y_pred == y)*1.0/n, R_w
        r=1-np.sum((y_pred != y)*R_w)
        return r

    @staticmethod
    def F1(W, X, y):
        q=np.exp(np.dot(X, W))
        q=q/(1.0+q)
        y_pred=np.array(q>=0.5, dtype=int)
        n=len(y)
        pos=sum(y)
        neg=n-pos
        tp=sum(y[y_pred>0.5])
        precision=tp*1.0/(sum(y_pred)+1e-5)
        recall=tp*1.0/(pos+1e-5)
        f1=2*precision*recall/(precision+recall+1e-5)
        return f1

    @staticmethod
    def metrics(W, X, y):
        q=np.exp(np.dot(X, W))
        q=q/(1.0+q)
        y_pred=np.array(q>=0.5, dtype=int)
        n=len(y)
        P=sum(y)
        N=n-P
        TP=sum(y[y_pred>=0.5])
        FP=sum(y_pred>=0.5)-TP
        TN=sum(1-y[y_pred<0.5])
        FN=sum(y_pred<0.5)-TN
        precision=TP*1.0/(TP+FP+1e-5)
        recall=TP*1.0/(P+1e-5)
        accuracy=(TP+TN)/(n+1e-5)
        avg_accuracy=(TP/(P+1e-5)+TN/(N+1e-5))/2
        F1=2*precision*recall/(precision+recall+1e-5)
        MCC=(TP*TN-FP*FN)/np.sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)+1e-5)
        return {'accuracy':accuracy, 'avg_accuracy':avg_accuracy, 'F1':F1, 'MCC':MCC}

    def fit(self, X, y, penalty=0.0, signs=1, Weights=None):
        n,m=X.shape
        kwargs={'maxiter':1000, 'ftol':1e-6}
        bounds=[(None,None)]
        if type(signs) is int:
            signs=[signs]*(m-1)
        for x in signs:
            if x>0.01:
                bounds.append((0, None)) # must be >=0
            elif x<-0.01:
                bounds.append((None, 0)) # must be <=0
            else:
                bounds.append((None, None)) # no constrain
        #W0=np.zeros(m)
        # set initial vector to point from the center of 0 to 1
        m0=X[y==0].mean(axis=0)
        m1=X[y==1].mean(axis=0)
        W0=m1-m0
        W0+=(np.random.rand(* W0.shape)-0.5)*1e-5 # avoid zero W0
        W0/=np.sqrt(sum(W0*W0)+1e-5)
        for j,x in enumerate(signs):
            if x>0:
                W0[j+1]=abs(W0[j+1])
            else:
                W0[j+1]=-abs(W0[j+1])
        res=optz.minimize(FindWeight.f, W0, (X, y, penalty, Weights), method='L-BFGS-B', bounds=bounds, options=kwargs)
        best_score=res.fun
        self.W=res.x

    def predict(self, X):
        q=np.exp(np.dot(X, self.W))
        q=q/(1.0+q)
        return q

    def auc(self, y, y_pred):
        fpr, tpr, thresholds=roc_curve(y, y_pred, pos_label=1)
        return auc(fpr, tpr)

    @staticmethod
    def evidence_weight(X, y, signs=1, folds=5, lb_penalty=-5, ub_penalty=5, num_penalty=11, Weights=None):
        """Give X (n*m) and y (n*1), returns weights (n+1) elements, and estimated weighted accuracy
        Weights is for sample weight, for unbalanced case"""
        fw=FindWeight()
        R_penalty=np.logspace(lb_penalty, ub_penalty, num=num_penalty)
        out=[]
        n,m=X.shape
        #print(n,m)
        X=np.hstack([np.ones([n,1]), X]) # add 1 as a dummie column for bias
        kf=StratifiedKFold(folds, shuffle=True)
        for k_train, k_test in kf.split(X, y):
            X_train=X[k_train]
            y_train=y[k_train]
            X_test=X[k_test]
            y_test=y[k_test]
            for penalty in R_penalty:
                fw.fit(X_train, y_train, penalty=penalty, signs=signs, Weights=Weights)
                #y_pred=fw.predict(X_test)
                #score=fw.auc(y[k_test], y_pred)
                score_cross_entropy=FindWeight.f(fw.W, X_test, y_test, 0, Weights)
                #score_error_rate=1-FindWeight.accuracy(fw.W, X_test, y_test)
                #f1=FindWeight.F1(fw.W, X_test, y_test)
                c=FindWeight.metrics(fw.W, X_test, y_test)
                #print(c)
                out.append({'Penalty':penalty, 'ErrorRate':(1-c['avg_accuracy']), 'CrossEntropy':score_cross_entropy, 'F1':c['F1'], 'MCC':c['MCC']})
        t_score=pd.DataFrame(out)
        out=[]
        scores=t_score.groupby('Penalty').mean().reset_index()
        scores.sort_values(['MCC','F1','ErrorRate','CrossEntropy'], ascending=[False, False, True, True], inplace=True)
        scores.index=list(range(len(scores)))
        print(scores[:])
        idx=0
        #print scores, idx
        best_penalty=scores.Penalty[idx]
        best_score=scores.ErrorRate[idx]
        fw.fit(X, y, penalty=best_penalty, signs=signs, Weights=Weights)
        #print myregress.W
        return (fw.W, 1.0-best_score)

    @staticmethod
    def score(X, R_w, Weights=None):
        y_pred=(X*R_w[1:]).sum(axis=1)+R_w[0]
        weights=FindWeight.norm_weight(Weights)
        if weights[0]/weights[1]==1:
            R_prob=1.0/(1.0+np.exp(-y_pred))
        else:
            R_prob=1.0/(1.0+np.exp(-y_pred)*weights[0]/(weights[1]+1e-10))
        return (y_pred, R_prob)

    @staticmethod
    def evidence_weight_wrapper(X, y, signs=1, folds=5, lb_penalty=-5, ub_penalty=5, num_penalty=11, l_auto_weight=True, X_test=None):
        """We automatically balance the sample classes, run prediction if X_test is not None"""
        Weights=None
        if l_auto_weight:
            Weights=FindWeight.auto_weight(y)
        R_w, accuracy = FindWeight.evidence_weight(X, y, signs=signs, folds=folds, lb_penalty=lb_penalty, ub_penalty=ub_penalty, num_penalty=num_penalty, Weights=Weights)
        if X_test is None:
            y_pred=R_prob=None
        else:
            y_pred, R_prob = FindWeight.score(X_test, R_w, Weights=Weights)
        return (R_w, accuracy, y_pred, R_prob)

    @staticmethod
    def nested_loop(X, y, signs=1, folds=5, lb_penalty=-5, ub_penalty=5, num_penalty=11, test_folds=5, Weights=None):
        """First reserve 1/test_folds portion of data as test set, use the remaining set for training (including validation), use the test set to evaluation performance.
        The model produce by evidence_weight() is the optimal model, however, the performance there is an overestimation, as
        the performance was tuned with all training data.  The better performance measure is to reserve some test data for
        evaluation purpose.  We get better performance accuracy by extra 5-fold computing time."""
        sw=util.StopWatch("FindWeight::nested_loop")
        kf=StratifiedKFold(test_folds, shuffle=True)
        R_accuracy=[]
        for k_train_validate, k_test in kf.split(X, y):
            X_train_validate=X[k_train_validate]
            y_train_validate=y[k_train_validate]
            X_test=X[k_test]
            n,m=X_test.shape
            X_test=np.hstack([np.ones([n,1]), X_test]) # add 1 as a dummie column for bias
            y_test=y[k_test]
            #print(len(y_train_validate), len(y_test))
            (R_w, score)=FindWeight.evidence_weight(X_train_validate, y_train_validate, signs=signs, folds=folds, lb_penalty=lb_penalty, ub_penalty=ub_penalty, num_penalty=num_penalty, Weights=Weights)
            R_accuracy.append(FindWeight.accuracy(R_w, X_test, y_test))
        R=np.array(R_accuracy)
        sw.check('Nested Loop')
        return (R.mean(), R.std())

    @staticmethod
    def plot(S_evidence, R_weight, y, R_score, R_prob, s_out, Weights=None):
        """S_evidence and R_weight are a list of evidence and their weights
        R_score, R_prob are evidence scores and their predicted probabilities
        y is the corresponding truth 1/0"""

        plt.figure(figsize=(4,4))
        plt.subplots(ncols=2, nrows=2)

        precision, recall, _ = precision_recall_curve(y, R_prob)
        # buggy over-estimate, as it use (1.0, 0.0) point incorrectly
        #average_precision = average_precision_score(y, R_prob)
        average_precision=np.sum(-np.diff(recall)*precision[:-1])

        ax=plt.subplot(221)
        plt.step(recall, precision, color='#3182bd', where='post')
        plt.fill_between(recall, precision, step='post', color='#9ecae1')
        plt.xlabel('Recall')
        plt.ylabel('Precision')
        plt.ylim([0.0, 1.0])
        plt.xlim([0.0, 1.0])
        ax.text(0.7, 0.9, 'AUC={0:0.2f}'.format(average_precision), ha='left', va='center', transform=ax.transAxes, color='#2c7bb6', fontsize=8)
        #plt.title('Average Precision={0:0.2f}'.format(average_precision))

        ax=plt.subplot(222)
        tpr=np.cumsum(y[np.argsort(-R_score)])
        #print(tpr[:50])
        tpr=tpr*1.0/np.sum(y)
        n=len(y)
        plt.plot([1,n+1], [0, 1], '-', c='#9ecae1')
        plt.step(np.arange(1, n+1), tpr, where='post', color='#3182bd')
        plt.xlabel('Rank')
        plt.ylabel('TPR')
        plt.ylim([0.0, 1.0])

        ax=plt.subplot(223)
        t2=pd.DataFrame(data={'Evidence':S_evidence, 'Weight':R_weight[1:]})
        t2['Sign']=np.sign(t2.Weight.values)
        t2['Sign']=t2.Sign.astype(int)
        t2['Weight']=np.abs(t2.Weight.values)
        t2['Weight']/=t2['Weight'].max()
        t2.sort_values(['Weight','Sign'], ascending=[False, False], inplace=True)
        t2=t2[t2.Weight >0].copy()
        max_evidence=10
        t2=t2[:max_evidence]
        t2.index=range(len(t2))
        n=len(t2)
        rotation=45 if len(t2) <= 4 else 90
        if n>0:
            X=np.arange(n)
            S_clr=['#9ecae1' if v > 0 else '#fdae6b' for v in t2.Sign ]
            plt.bar(X, t2.Weight.values, color=S_clr)
            plt.ylabel('Weight')
            ax.set_xticks(range(n))
            ax.set_xticklabels(t2.Evidence.tolist(), rotation=rotation)

        ax=plt.subplot(224)
        xmin=min(-3, R_score.min())
        xmax=max(3, R_score.max())
        X=np.linspace(xmin, xmax, num=50)
        weights=FindWeight.norm_weight(Weights)
        Y=1.0/(1.0+np.exp(-X)*weights[0]/weights[1])
        plt.plot(X, Y, '-', color='#9ecae1')
        X=R_score[ y>0 ]
        Y=R_prob[ y>0 ]
        n=len(X)
        plt.plot(X, Y+np.random.rand(n)*0.1-0.05, '+', color='#e6550d', alpha=0.8)
        plt.xlabel('Evidence Score')
        plt.ylabel('Probability')
        plt.ylim([-0.05, 1.05])

        plt.tight_layout()
        if type(s_out) is not list:
            s_out=[s_out]
        for s_file in s_out:
            plt.savefig(s_file)


if __name__ == "__main__":
    X=np.random.rand(100,4)
    y=np.sum(X*np.array([1,2,3,4]), axis=1)+np.random.rand(100)*0.3
    y=np.array(y> np.median(y), dtype=int)
    #X2=np.ones([100,5])
    #X2[:,1:]=X
    fw=FindWeight()
    #R_w, acc=fw.evidence_weight(X, y, signs=1, folds=5, lb_penalty=-5, ub_penalty=5, num_penalty=11)
    R_w, acc, y_pred, prob=fw.evidence_weight_wrapper(X, y, signs=1, folds=5, lb_penalty=-5, ub_penalty=5, num_penalty=11, l_auto_weight=True, X_test=None)
    print(R_w/np.max(R_w)*4, acc)

    S_evidence=['A','B','C','D']
    R_score, R_prob=fw.score(X, R_w)
    R_prob=1.0/(1+np.exp(-R_score))
    fw.plot(S_evidence, R_w, y, R_score, R_prob, ['example.png','example.pdf'])

