#!/usr/bin/env python
from __future__ import absolute_import
import re
import pandas as pd
import util
import os
import xgmml
from six.moves import range
from six.moves import zip

class Node(object):

    def __init__(self, id, label='', left=None, right=None, similarity=1.0):
        self.id, self.label, self.left, self.right, self.similarity = id, label, left, right, similarity
        if self.label=='': self.label=self.id

    def has_child(self):
        return not(self.left is None and self.right is None)

    def all_children(self, l_keep_node=False, l_include_self=True):
        # need to rewrite to not use recursion, otherwise, won't work for deep trees
        children=[]
        queue=[self]
        base=self.id
        while len(queue)>0:
            x=queue.pop(0)
            if x.has_child():
                if l_keep_node:
                    children.append(x.id)
                if x.left is not None:
                    queue.append(x.left)
                if x.right is not None:
                    queue.append(x.right)
            elif not l_keep_node:
                children.append(x.id)
        if not l_include_self:
            children=[x for x in children if x!=base]
        return children

    #def all_children(self, l_keep_node=False, l_include_self=True):
    #    # need to rewrite to not use recursion, otherwise, won't work for deep trees
    #    if l_keep_node:
    #        return self.all_children_nodes(l_include_self=l_include_self)
    #    children=[]
    #    if self.left is not None:
    #        if self.left.has_child():
    #            children.extend(self.left.all_children())
    #        else:
    #            ch ildren.append(self.left.id)
    #    if self.right is not None:
    #        if self.right.has_child():
    #            children.extend(self.right.all_children())
    #        else:
    #            children.append(self.right.id)
    #    return children

    def all_children_nodes(self, l_include_self=True):
        """Only keep children nodes, not leaves."""
        return self.all_children(l_keep_node=True, l_include_self=l_include_self)
        #children=[]
        #if l_include_self and self.has_child(): children.append(self.id)
        #if self.left is not None and self.left.has_child():
        #    children.extend(self.left.all_children_nodes(l_include_self=True))
        #if self.right is not None and self.right.has_child():
        #    children.extend(self.right.all_children_nodes(l_include_self=True))
        #return children

    # return the ID of the most representative node, and the # of nodes it represents
    #def representative(self):
    #    if not self.has_child():
    #        return (self.id, 1)
    #    else:
    #        (ln, i_l)=self.left.representative()
    #        (rn, i_r)=self.right.representative()
    #        if (i_l<i_r):
    #            return (rn, i_l+i_r)
    #        else:
    #            return (ln, i_l+i_r)

    def all_nof_leaves(self):
        """compute the total number of leaves under each node"""
        c_cnt={}
        q=[self]
        seen=[]
        while len(q):
            k=q.pop(0)
            seen.append(k)
            if k.left is not None:
                q.append(k.left)
            if k.right is not None:
                q.append(k.right)
        seen.reverse()
        for x in seen:
            if not x.has_child():
                c_cnt[x.id]=1
            else:
                c_cnt[x.id]=0
                if x.left is not None:
                    c_cnt[x.id]+=c_cnt[x.left.id]
                if x.right is not None:
                    c_cnt[x.id]+=c_cnt[x.right.id]
        return c_cnt

    def representative(self):
        """return the ID of the most representative node, and the # of nodes it represents
        Rewrite to avoid recursion"""
        # first compute the total number of leaves under each node
        c_cnt=self.all_nof_leaves()
        # now pick the representative gene from the larger tree branch
        k=self
        while True:
            if not k.has_child():
                return (k.id, c_cnt[self.id])
            else:
                if c_cnt[k.left.id]<c_cnt[k.right.id]:
                    k=k.right
                else:
                    k=k.left

    #def node_similarities(self):
    #    """Return list of (node, similarity) and sort them from small to large"""
    #    if self.has_child():
    #        out=self.left.node_similarities()+self.right.node_similarities()
    #        out.append((self, self.similarity))
    #        return sorted(out, key=lambda(x): x[1])
    #    else:
    #        return []

    def node_similarities(self):
        """Return list of (node, similarity) and sort them from small to large
        Rewrite to avoid recursion."""
        q=[self]
        out=[]
        while len(q):
            k=q.pop(0)
            if k.left is not None:
                q.append(k.left)
            if k.right is not None:
                q.append(k.right)
            if k.has_child():
                out.append((k, k.similarity))
        return sorted(out, key=lambda x: x[1])

    ## n_picks most representative nodes
    #def representatives(self, n_picks=1):
    #    if n_picks==1 or not self.has_child():
    #        return [self.representative()]
    #    else:
    #        out=[]
    #        # take the n_picks most representative subtrees
    #        L_nodes=self.node_similarities()[:n_picks-1]
    #        c_nodes={ n.id:True for n,s in L_nodes }
    #        for node,s in L_nodes:
    #            #print node.id, "<<<", node.left.id, ">>>", node.right.id
    #            if node.left.id not in c_nodes:
    #                out.append(node.left.representative())
    #            if node.right.id not in c_nodes:
    #                out.append(node.right.representative())
    #        return sorted(out, key=lambda(x): -x[1])

    # n_picks most representative nodes, each node must represent at least min_size nodes
    def representatives(self, n_picks=1, min_size=1, l_keep_members=False):

        def cut_grps(L_nodes, i_cut):
            L=L_nodes[:i_cut]
            c_nodes={ n.id:True for n,s in L }
            out=[]
            for node,s in L:
                #print node.id, "<<<", node.left.id, ">>>", node.right.id
                if node.left.id not in c_nodes and c_cnt[node.left.id]>=min_size:
                    X=node.left.representative()
                    if l_keep_members: X=(X[0], X[1], node.left.all_children())
                    out.append(X)
                if node.right.id not in c_nodes and c_cnt[node.right.id]>=min_size:
                    X=node.right.representative()
                    if l_keep_members: X=(X[0], X[1], node.right.all_children())
                    out.append(X)
            return (len(out), out, )

        c_cnt=self.all_nof_leaves()
        if n_picks==1 or not self.has_child():
            if c_cnt[self.id]>=min_size:
                X=self.representative()
                if l_keep_members: X=(X[0], X[1], self.all_children())
                return [X]
            else:
                return []
        else:
            # take the n_picks most representative subtrees
            L_nodes=self.node_similarities()
            if min_size>1:
                L_nodes=[(n,s) for n,s in L_nodes if c_cnt[n.id]>=min_size]
            out=[]
            # have not found a clever way, so just try different cutoffs, until we get n_picks groups
            # if min_size==1, it should get it right the first time
            (best_i, best_n, best_out)=(1, 1, [])
            for i in range(n_picks-1, len(L_nodes)+1):
                n, out=cut_grps(L_nodes, i)
                #print ">>", i, n, n_picks
                if n>=n_picks:
                    return sorted(out, key=lambda x: -x[1])
                if abs(n-n_picks)<abs(best_n-n_picks) or i==n_picks-1:
                    (best_i, best_n, best_out)=(i, n, out)
            #print ">>>", best_i, best_n
            return sorted(best_out, key=lambda x: -x[1])

    #def cut_(self, min_similarity=0.8, l_keep_node=False):
    #    """Replaced by cut() below. If l_keep_node, only output the nodes, instead of genes."""
    #    if self.similarity>=min_similarity:
    #        if self.has_child():
    #            return [self.all_children(l_keep_node=l_keep_node)]
    #        else:
    #            return [] if l_keep_node else [[self.id]]
    #    else:
    #        out=[]
    #        if self.left is not None:
    #            out_left=self.left.cut(min_similarity=min_similarity, l_keep_node=l_keep_node)
    #            if out_left: out.extend(out_left)
    #        if self.right is not None:
    #            out_right=self.right.cut(min_similarity=min_similarity, l_keep_node=l_keep_node)
    #            if out_right: out.extend(out_right)
    #        return out

    def cut(self, min_similarity=0.8, l_keep_node=False):
        """If l_keep_node, only output the nodes, instead of genes. Rewrite it to avoid recursion, so it works for flat trees."""
        q=[self]
        out=[]
        while len(q):
            k=q.pop(0)
            if k.similarity>=min_similarity:
                if k.has_child():
                    out.append(k.all_children(l_keep_node=l_keep_node))
                elif not l_keep_node:
                    out.append([k.id])
            else:
                if k.left is not None:
                    q.append(k.left)
                if k.right is not None:
                    q.append(k.right)
        return out

    #def bicut(self, high_similarity=0.8, low_similarity=0.6):
    #    if self.similarity<low_similarity:
    #        out=[]
    #        #print self.id, "too low"
    #        if self.left is not None:
    #            #if self.left.has_child():
    #                out_left=self.left.bicut(high_similarity=high_similarity, low_similarity=low_similarity)
    #                if out_left: out.extend(out_left)
    #        if self.right is not None:
    #            #if self.right.has_child():
    #                out_right=self.right.bicut(high_similarity=high_similarity, low_similarity=low_similarity)
    #                if out_right: out.extend(out_right)
    #    elif self.similarity<high_similarity:
    #        #print self.id, "middle"
    #        out=[[]]
    #        if self.left is not None:
    #            #if self.left.has_child():
    #                out_left=self.left.bicut(high_similarity=high_similarity, low_similarity=low_similarity)
    #                out[0].extend(out_left[0])
    #            #else:
    #            #    out[0].extend([self.left.id])
    #        if self.right is not None:
    #            #if self.right.has_child():
    #                out_right=self.right.bicut(high_similarity=high_similarity, low_similarity=low_similarity)
    #                out[0].extend(out_right[0])
    #            #else:
    #            #    out[0].extend([self.right.id])
    #        if not self.has_child(): out.append([self.id])
    #    else: # >=high_similarity
    #        #print self.id, "high"
    #        if not self.has_child():
    #            out=[[[self.id]]]
    #        else:
    #            out=[[self.all_children()]]
    #    #print "return >>>", self.id, out
    #    return out

    def bicut(self, high_similarity=0.8, low_similarity=0.6):
        """Rewrite into non-recursive version"""
        q=[self]
        out=[]
        while(len(q)):
            k=q.pop(0)
            if k.similarity<low_similarity:
                if k.left is not None:
                    q.append(k.left)
                if k.right is not None:
                    q.append(k.right)
            elif k.similarity<high_similarity:
                out2=[]
                if k.left is not None:
                    out2.extend(k.left.cut(min_similarity=high_similarity, l_keep_node=False))
                if k.right is not None:
                    out2.extend(k.right.cut(min_similarity=high_similarity, l_keep_node=False))
                if not k.has_child():
                    out2.extend([k.id])
                out.append(out2)
            else: # >=high_similarity
                if not k.has_child():
                    out.append([[k.id]])
                else:
                    #print k.all_children()
                    out.append([k.all_children()])
        return out

    def __str__(self, level=0):
        s=''
        if self.has_child():
            if self.left is not None: s+=self.left.__str__(level+1)
            if self.right is not None: s+=self.right.__str__(level+1)
        else:
            s='  '*level+self.id+':'+self.label+'\n'
        return s

class Tree(object):

    def __init__(self, s_file='', Z=None, l_gene_tree=True):
        """Z: linkage matrix, if None, assume s_file is not empty"""
        self.l_gene_tree=l_gene_tree
        self.root=Node('ROOT')
        self.l_gene_tree=l_gene_tree    # gene tree or array tree
        self.c_name={}
        self.c_node={}
        self.size=0
        self.parent={} # track the parent node for each node
        self.tree_file=None

        if Z is not None:
            self.l_gene_tree=True
            r,c=Z.shape
            n=r+1
            r_dist=max(Z[:, 2].max(), 1.0)
            for i in range(r):
                id_l=str(int(Z[i, 0]))
                id_r=str(int(Z[i, 1]))
                id_n=str(n+i)
                r=max(1.0-Z[i, 2]/r_dist, 0.0)
                self.new_node(id_n, label=self.c_name.get(id_n, ''), left=self.new_node(id_l), right=self.new_node(id_r), similarity=r)
                self.parent[id_l]=id_n
                self.parent[id_r]=id_n
            self.root=self.get_node(id_n)
            self.size=n-1
        else:
            self.l_gene_tree=l_gene_tree
            if re.search(r'\.[ag]tr$', s_file):
                if re.search(r'\.atr$', s_file):
                    l_gene_tree=False
                s_file=re.sub(r'\.[ag]tr$', '', s_file)
            self.root=Node('ROOT')
            self.l_gene_tree=l_gene_tree    # gene tree or array tree
            self.c_name={}
            self.c_node={}
            self.size=0
            self.parent={} # track the parent node for each node
            if not os.path.exists(s_file+".cdt"):
                util.error_msg("File not exist: "+s_file+".cdt!")
            f=open(s_file+'.cdt')
            S_header=f.readline().strip().split("\t")
            if not l_gene_tree:
                while True:
                    line=f.readline()
                    if not line: break
                    if line.startswith("AID\t"):
                        S_AID=line.strip().split("\t")
                        self.c_name={s:x for s,x in zip(S_AID, S_header) if str(s).startswith('ARRY')}
                        break
            else:
                s_col='GENE'
                if s_col not in S_header and 'NAME' in S_header:
                    s_col='NAME'
                i_GID=util.index('GID', S_header)
                i_NAME=util.index(s_col, S_header)
                while True:
                    line=f.readline()
                    if not line: break
                    if line.startswith('AID') or line.startswith('EWEIGHT'):
                        continue
                    S=line.strip().split("\t")
                    self.c_name[S[i_GID]]=S[i_NAME]
            f.close()
            self.size=len(self.c_name)
            if self.size==0: error_msg("Tree:__init_: No node is found to build the tree!")
            s_filename=s_file+('.gtr' if l_gene_tree else '.atr')
            # check if file has column header
            self.tree_file=s_filename
            df=Tree.read_tree_file(s_filename)
            self.parse(df)

    def get_node(self, id):
        return self.c_node.get(id, None)

    #def nof_nodes(self):
    #    # this includes both NODE* and GENE*
    #    return len(self.c_node)

    def nof_leaves(self):
        # c_node contains n leaves and n-1 nodes
        return (len(self.c_node)+1)/2

    def get_node_by_name(self, s_name):
        for k,v in self.c_name.items():
            if v==s_name:
                return self.get_node(k)
        return None

    def get_parents(self, id):
        parents=[]
        while True:
            p=self.parent.get(id, None)
            if p is None: break
            parents.append(p)
            id=p
        return parents

    # find the nearest node that contains all ids in id_list as child
    # this is useful to reverse find the subtree using subtree members
    # e.g., find the tree contains a group made by cut
    # warning, if id_list has only one gene, it will return itself
    def get_umbrella_node(self, id_list):
        if len(id_list)==0:
            return None
        if len(id_list)==1:
            return self.get_node(id_list[0])
        reference_gene=None # the node that descends fastest in tree
        reference_depth=0
        paths={}
        for id in id_list:
            p=self.get_parents(id)
            paths[id]=p
            if reference_gene is None or len(p)<reference_depth:
                reference_gene=id
                reference_depth=len(p)
        #print paths
        #print reference_gene, reference_depth, paths[reference_gene]
        for p in paths[reference_gene]:
            #print ">>>>>>>> "+p
            l_all=True
            for k,v in paths.items():
                if p not in v:
                    l_all=False
                    break
            if l_all: return self.get_node(p)
        return None

    def new_node(self, id, label='', left=None, right=None, similarity=1.0):
        if self.get_node(id) is None:
            self.c_node[id]=Node(id, label, left=left, right=right, similarity=similarity)
        return self.get_node(id)

    def parse(self, df):
        n=len(df)
        for i in range(n):
            id_l=df.loc[i, 'left']
            id_r=df.loc[i, 'right']
            id_n=df.loc[i, 'node']
            r=float(df.loc[i, 'similarity'])
            #print id_n, id_l, id_r, r
            self.new_node(id_n, label=self.c_name.get(id_n, ''), left=self.new_node(id_l), right=self.new_node(id_r), similarity=r)
            self.parent[id_l]=id_n
            self.parent[id_r]=id_n
        self.root=self.get_node(id_n)
        self.size=n

    def representatives(self, n_picks=1, min_size=1, l_keep_members=False):
        return self.root.representatives(n_picks=n_picks, min_size=min_size, l_keep_members=l_keep_members)

    def cut(self, min_similarity=0.8, l_keep_node=False):
        return self.root.cut(min_similarity=min_similarity, l_keep_node=l_keep_node)

    def bicut(self, high_similarity=0.8, low_similarity=0.6):
        return self.root.bicut(high_similarity=high_similarity, low_similarity=low_similarity)

    @staticmethod
    def read_tree_file(s_filename):
        f=open(s_filename, "r")
        s=f.readline()
        l_has_header = s.startswith('NODEID')
        f.close()
        if l_has_header:
            df=pd.read_table(s_filename)
        else:
            df=pd.read_table(s_filename, header=None)
        S=['node','left','right','similarity']
        if len(df.header())==5:
            S.append('color')
        df.columns=S
        return df

    @staticmethod
    def color_map(nodes, cm=None):
        """create a dictionary of node-hex color mapping
        nodes:
            (1) dict of {node_id: color}, nodes not found in the dict will be colored as black
                value color:
                (a) matplotlib.colors, tuple of floats, or hex string #FF0000
                (b) int: index into cm (a list of colors)
                (c) float or int, when cm is not a list, values will be normalized into
                    [0,1] and colormap cm is used to translate the value into a color
            (2) list of list, [['node1','node5'],['node2']], each sublist is assigned one color
                In this case cm should be a list of colors of the same length ['red','blue']
        cm: matplotlib.mcolors.LinearSegmentedColormap or a list of colors
            if None, we use rainbow colormap
        Examples:
            Tree.color_map({'node1':'#ff0000', 'node3':'#0000ff'})
            Tree.color_map({'node1':0, 'node3':1}, ['#ff0000','#0000ff'])
            Tree.color_map({'node1':0, 'node3',1}) # cm is set to matplotlib.cm.gist_rainbow
            Tree.color_map([['node1'], ['node2','node3']], ['#ff0000','#0000ff'])
        return dict of {node: hex_color}
        """
        import matplotlib
        import matplotlib.cm
        import matplotlib.colors as mcolors
        if cm is None:
            cm=matplotlib.cm.gist_rainbow
        if type(nodes) is dict:
            R=[]
            for k,v in nodes.items():
                if type(v) in (int, float):
                    R.append(v)
            if len(R):
                r_min,r_max=min(R), max(R)
            for k,v in nodes.items():
                if type(v) is int:
                    if type(cm) is list:
                        nodes[k]=list[v % len(cm)]
                    else: # cm must be a colormap
                        nodes[k]=cm((v-r_min)/(r_max-r_min))
                elif type(v) is float:
                    nodes[k]=cm((v-r_min)/(r_max-r_min))
        else: # nodes must be a list
            c={}
            n=len(nodes)
            for i,X in enumerate(nodes):
                if type(cm) is list:
                    clr=cm[i%n]
                else:
                    clr=cm(i*1.0/(n-1)) if n>1 else cm(1.0)
                for x in X:
                    c[x]=clr
            nodes=c
        for k,v in nodes.items():
            if type(v) is tuple:
                v=[ min(max(int(x*255),0), 255) for x in v]
                v='#%02x%02x%02x' % tuple(v[:3])
                nodes[k]=v
        return nodes

    def color(self, nodes, l_name_to_id=True, cm=None):
        """Color tree nodes
        nodes: dict or list of nodes lists
        nodes and colormap cm are combined to passed to Tree.color_map (see document)
        l_name_to_id: bool, default True, use leave name or node_id
            Warning: this method color nodes, so if leave name is provided and two leaves
            under the same node has different colors, only one color is used.
        """
        df=Tree.read_tree_file(self.tree_file)
        # ['node','left','right','similarity', 'color']
        df['color']='#000000'
        c_id2node={}
        c_name2id={}
        if l_name_to_id:
            for k,v in self.c_name.items():
                c_name2id[v]=k
            for i in df.index:
                if not df.loc[i, 'left'].startswith('NODE'):
                    c_id2node[df.loc[i, 'left']]=df.loc[i, 'node']
                if not df.loc[i, 'right'].startswith('NODE'):
                    c_id2node[df.loc[i, 'right']]=df.loc[i, 'node']
        t=df.loc[:, ['node','color']]
        t.set_index('node', inplace=True)
        c=Tree.color_map(nodes, cm)
        for k,v in c.items():
            if type(v) is tuple:
                v=[ min(max(int(x*255),0), 255) for x in v]
                v='#%02x%02x%02x' % tuple(v[:3])
            if l_name_to_id:
                k=c_name2id.get(k,k)
                k=c_id2node.get(k,k)
            if k in t.index:
                t.loc[k, 'color']=v
            else:
                util.warn_msg('Node not in tree: '+k)
        df['color']=list(t.color)
        df.columns=['NODEID','LEFT','RIGHT','CORRELATION','NODECOLOR']
        util.df2sdf(df, s_format="%.4f").to_csv(self.tree_file, index=False, sep="\t")

    def color_cut(self, min_similarity=0.8, S_COLOR=None):
        #colorbrewer2.org, qualitative, 5-class set
        #S_COLOR=['#E41A1C', '#277EB8', '#4DAF4A', '#984EA3', '#FF7F00']
        c_cut=self.cut(min_similarity=min_similarity, l_keep_node=True)
        self.color(c_cut, False, S_COLOR)

    # tries to pick up to n_picks most representative genes within each group
    def cut2table(self, c_cut, n_picks=1):
        rows=[]
        for i,grp in enumerate(c_cut):
            tr=self.get_umbrella_node(grp)
            L_rep=tr.representatives(n_picks)
            c_rep={ n:c for n,c in L_rep }
            sz=len(grp)
            for g in grp:
                rows.append({'GroupID':i+1, 'GroupSize':sz, 'Entry':self.c_name.get(g, g), 'RepresentCounts':c_rep.get(g, 0)})
        df=pd.DataFrame(rows)
        df=df.sort_values(['GroupSize', 'GroupID'], ascending=[False, True])
        return df

    def bicut2table(self, c_bicut, n_opt_picks=1, n_ok_picks=0):
        rows=[]
        for i,ok_grp in enumerate(c_bicut):
            ok_g=[]
            for opt_g in ok_grp:
                ok_g.extend(opt_g)
            ok_sz=len(ok_g)
            c_ok_rep={}
            if n_ok_picks:
                tr=self.get_umbrella_node(ok_g)
                L_rep=tr.representatives(n_ok_picks)
                c_ok_rep={ n:c for n,c in L_rep }
            for j,opt_g in enumerate(ok_grp):
                opt_sz=len(opt_g)
                c_opt_rep={}
                if n_opt_picks:
                    tr=self.get_umbrella_node(opt_g)
                    L_rep=tr.representatives(n_opt_picks)
                    c_opt_rep={ n:c for n,c in L_rep }
                for g in opt_g:
                    one={'OkayGroupID':i+1, 'OkayGroupSize':ok_sz, 'OptimalGroupID':j+1, 'OptimalGroupSize':opt_sz, 'Entry':self.c_name.get(g, g)}
                    if n_ok_picks:
                        one['OkayRepresentCounts']=c_ok_rep.get(g,0)
                    if n_opt_picks:
                        one['OptimalRepresentCounts']=c_opt_rep.get(g,0)
                    rows.append(one)
        df=pd.DataFrame(rows)
        df=df.sort_values(['OkayGroupSize', 'OptimalGroupSize', 'OkayGroupID', 'OptimalGroupID'], ascending=[False, False, True, True])
        return df

    def __str__(self):
        return str(self.root)

    def to_network_(self, l_digraph=False):
        out=[]
        for x in self.c_node.values():
            is_node=1 if x.has_child() else 0
            if x==self.root:
                is_node=2
            out.append({'Gene':x.id, 'Symbol': self.c_name.get(x.id, x.id), 'IsNode': is_node})
        t_node=pd.DataFrame(out)
        queue=[self.root]
        out=[]
        while len(queue)>0:
            x=queue.pop(0)
            #print ">>>", x.id, x.left.id, x.right.id
            if x.left is not None:
                #print x.left, x.left.similarity
                dist=abs(x.left.similarity - x.similarity)
                out.append({'Gene_A': x.id, 'Gene_B': x.left.id, 'Length': dist, 'Similarity':x.similarity, 'Length_inv':abs(1.0-dist), 'Similarity_inv':(1-x.similarity)})
                if x.left.has_child():
                    queue.append(x.left)
            if x.right is not None:
                #print x.right, x.right.similarity
                dist=abs(x.right.similarity - x.similarity)
                out.append({'Gene_A': x.id, 'Gene_B': x.right.id, 'Length': dist, 'Similarity':x.similarity, 'Length_inv':abs(1.0-dist), 'Similarity_inv':(1-x.similarity)})
                if x.right.has_child():
                    queue.append(x.right)
        t_edge=pd.DataFrame(out)
        #print t_edge[:]
        if l_digraph:
            import digraph
            net=digraph.Digraph(t_edge, name='Untitled', T_node=t_node, s_noa=None)
        else:
            net=xgmml.Network(t_edge, allow_indirect=False, name='Untitled', T_node=t_node, s_noa=None)
        return net

    def to_text(self):
        # tuple is (ident_prefix, is_left, tree)
        queue=[("*-", False, self.root)]
        out=[]
        while len(queue)>0:
            s_ident, is_left, x=queue.pop(0)
            out.append(s_ident+x.id+":"+self.c_name.get(x.id,'')+(" (%.4f)" % x.similarity))
            s_ident=s_ident[:-2]+"| +-" if is_left else s_ident[:-2]+"  +-"
            if x.right is not None:
                queue.insert(0, (s_ident, False, x.right))
            if x.left is not None:
                queue.insert(0, (s_ident, True, x.left))
        return "\n".join(out)

    def to_network(self):
        return self.to_network_(False)

    def to_digraph(self):
        return self.to_network_(True)

def read_cdt(s_file):
    if not s_file.endswith('.cdt'):
        s_file+='.cdt'
    if not os.path.exists(s_file):
        util.error_msg("File not exist: "+s_file+"!")
    f=open(s_file)
    S_header=f.readline().strip().split("\t")
    i_w=util.index("GWEIGHT", S_header)
    i_gene=util.index('GENE', S_header)
    i_name=util.index('NAME', S_header)
    l_start=False
    R_exp=[]
    R_gene=[]
    data=[]
    offset=0
    while True:
        line=f.readline()
        if not line: break
        S=line.strip().split("\t")
        if S[0]=='EWEIGHT':
            for i in range(1,len(S)):
                if S[i]!="":
                    offset=i
                    break
            tmp=[]
            if i_gene>=0: tmp.append(S_header[i_gene])
            if i_name>=0: tmp.append(S_header[i_name])
            S_header=tmp+S_header[offset:]
            R_exp=util.sarray2rarray(S[offset:])
            if i_w<0: i_w=offset-1
            l_start=True
        elif l_start:
            one=[]
            if i_gene>=0: one.append(S[i_gene])
            if i_name>=0: one.append(S[i_name])
            one.extend(util.sarray2rarray(S[offset:]))
            data.append(one)
            R_gene.append(float(S[i_w]))
    f.close()
    t=pd.DataFrame(data, columns=S_header)
    return (t, R_exp, R_gene)

def add_column(s_file, R, s_name, l_separator=True):
    """Add an extra column using value R array to an existing heat map.
    s_file: str, file name without extension, it will modify .cdt and .atr
    R: array(int/float), values to add
    s_name: str, column name
    l_separator: bool, default True. If True, add a column of blank value to separate the new column from existing ones."""
    if re.search('\.\w{3}$', s_file):
        s_file=s_file[:-4]
    if not os.path.exists(s_file+'.cdt'):
        util.error_msg("File not exist: "+s_file+".cdt!")
    f=open(s_file+'.cdt')
    S=[]
    cnt=0

    while True:
        line=f.readline()
        if not line: break
        SS=line.strip().split("\t")
        if SS[0].startswith('GENE'):
            if l_separator:
                SS.append('')
            SS.append('%.2f' % R[cnt])
            cnt+=1
        elif SS[0]=='GID':
            if l_separator:
                SS.append('separator')
            SS.append(s_name)
        elif SS[0]=='AID':
            X=[ int(re.sub(r'\D', '', x)) for x in SS if x.startswith('ARRY') ]
            n_array=max(X)+1
            SS.append('ARRY%dX' % n_array)
            if l_separator:
                SS.append('ARRY%dX' % (n_array+1))
        elif SS[0]=='EWEIGHT':
            if l_separator:
                SS.append('0')
            SS.append('0')
        S.append(SS)
    f.close()
    S=["\t".join(X) for X in S]
    util.save_list(s_file+'.cdt', S, s_end="\n")

    if os.path.exists(s_file+'.atr'):
        S=util.read_list(s_file+'.atr')
        SS=S[-1].split("\t")
        n_node=int(re.sub(r'\D', '', SS[0]))+1
        S.append('NODE%dX\tNODE%dX\tARRY%dX\t0' % (n_node, n_node-1, n_array))
        if l_separator:
            S.append('NODE%dX\tNODE%dX\tARRY%dX\t0' % (n_node+1, n_node, n_array+1))
        util.save_list(s_file+'.atr', S, s_end="\n")

def color_cdt(s_file, exps=None, exp_bgcolor=None, genes=None, gene_bgcolor=None):
    if not s_file.endswith('.cdt'):
        s_file+='.cdt'
    if not os.path.exists(s_file):
        util.error_msg("File not exist: "+s_file+"!")
    BG='#ffffff'
    f=open(s_file)
    S=[]
    c_first={}
    i=0
    while True:
        line=f.readline()
        if not line: break
        SS=line.strip().split("\t")
        c_first[SS[0]]=i
        i+=1
        S.append(SS)
    f.close()
    S_header=S[0]
    i_gene=util.index('GENE', S_header)
    i_name=util.index('NAME', S_header)
    i_gid=util.index('GID', S_header)
    i_w=util.index("GWEIGHT", S_header)
    offset=max([i_gene, i_name, i_gid, i_w])+1
    n_exp=len(S_header)-offset
    if 'EWEIGHT' not in c_first:
        # add EWEIGHT ROW
        i_w=max([c_first.get('GID',-1), c_first.get('AID',-1)])+1
        S.insert(i_w, ['EWEIGHT']+['']*(offset-1)+['1.000']*n_exp)
        c_first['EWEIGHT']=i_w
    i_w=util.index("GWEIGHT", S_header)
    if i_w<0: # add GWEIGHT column
        i_w=offset
        S_header.insert(i_w,'GWEIGHT')
        for i in range(1,len(S)):
            if i<=c_first['EWEIGHT']:
                S[i].insert(i_w,'')
            else:
                S[i].insert(i_w,'1.000')
        offset+=1
    i_gene_color=util.index('BGCOLOR', S_header)
    if i_gene_color<0 and genes is not None:
        i_gene_color=offset-1
        S_header.insert(i_gene_color, 'BGCOLOR')
        offset+=1
        for i in range(1, len(S)):
            if i<=c_first['EWEIGHT']:
                S[i].insert(i_gene_color,'')
            else:
                S[i].insert(i_gene_color, BG)
    i_exp_color=c_first.get('BGCOLOR', -1)
    if i_exp_color<0 and exps is not None:
        i_exp_color=c_first['EWEIGHT']
        S.insert(i_exp_color, ['BGCOLOR']+['']*(offset-1)+[BG]*n_exp)
        c_first['EWEIGHT']+=1
    if genes is not None:
        c_m=Tree.color_map(genes, gene_bgcolor)
        idx=i_gene if i_gene>=0 else i_name
        for i in range(c_first['EWEIGHT']+1, len(S)):
            S[i][i_gene_color]=c_m.get(S[i][idx], BG)
    if exps is not None:
        c_m=Tree.color_map(exps, exp_bgcolor)
        SS=S[c_first['EWEIGHT']-1]
        for i in range(offset, len(SS)):
            SS[i]=c_m.get(S_header[i], BG)
    S=["\t".join(X) for X in S]
    util.save_list(s_file, S, s_end="\n")

if __name__ == '__main__':
    #tr=Tree(s_file="test/clusteringOpt", l_gene_tree=False)
    #print tr.root.all_children()
    tr=Tree(s_file="test/clusteringOpt", l_gene_tree=True)
    #print(tr.to_text())
    grp=tr.cut(min_similarity=0.8)
    df=tr.cut2table(grp)
    df.to_csv('test/grp.csv', index=False)
    #tr=Tree(s_file="test/clusteringOpt", l_gene_tree=False)
    #print str(tr)
    grp=tr.bicut(high_similarity=0.8, low_similarity=0.4)
    df=tr.bicut2table(grp)
    df.to_csv('test/grp2.csv', index=False)
    #print grp

#  function ToCytoscape(s_node, S_noa=Sarray(), t_sif=Table())
#    local node;
#    node=data[s_node]
#    if node==null then
#      #S_noa//=s_node;
#    else
#      S_noa//=node.name;
#      if Type(node.left)=="string" then
#        #S_noa//=node.left;
#        t_sif//={Gene_A=node.name, Gene_B=node.left};
#      else
#        ToCytoscape(node.left.name, S_noa, t_sif);
#        t_sif//={Gene_A=node.name, Gene_B=node.left.name};
#      endif
#      if Type(node.right)=="string" then
#        #S_noa//=node.right;
#        t_sif//={Gene_A=node.name, Gene_B=node.right};
#      else
#        ToCytoscape(node.right.name, S_noa, t_sif);
#        t_sif//={Gene_A=node.name, Gene_B=node.right.name};
#      endif
#    endif
#  endfunction
#endclass
