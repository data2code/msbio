#!/usr/bin/env python
from __future__ import absolute_import
from __future__ import print_function
import pandas as pd
import re
import util
import os
import glob
import tempfile
import brewer2mpl
import numpy as np
from six.moves import range
import ms.mssetting
#from gpsetting import Setting

# Install circos-0.67-7
# modify etc/housekeeping.conf, so that it supports chromosome name with space
# Use TAB as delimiters in our conf files
# file_delim = \t
# put metascape.ideogram.conf and metascape.ticks.conf into etc subfolder

class Circos(object):

    def __init__(self, BIN=None, TEMPLATE=None):
        self.BIN=BIN if BIN is not None else ms.mssetting.circos['CIRCOS_BIN']
        DEFAULT_TEMPLATE=os.path.join(os.path.dirname(__file__), "circos", "circos.conf.template")
        self.TEMPLATE=TEMPLATE if TEMPLATE is not None else DEFAULT_TEMPLATE
        #if not os.path.exists(self.BIN):
        #    util.error_msg("Circos tool " + self.BIN + " does not exist!")
        if not os.path.exists(self.TEMPLATE):
            util.error_msg("Circos template " + self.TEMPLATE + " does not exist!")

    @staticmethod
    def hex_to_rgb(value):
        value = value.lstrip('#')
        lv = len(value)
        return ",".join([str(int(value[i:i + lv // 3], 16)) for i in range(0, lv, lv // 3)])

    @staticmethod
    def get_qualitative_colors(n):
        if n>21: n=21
        if n<=9:
            C=brewer2mpl.get_map("Set1", 'qualitative', max(n,3)).colors[:]
        elif n<=12:
            C=brewer2mpl.get_map("Paired", 'qualitative', max(n,3)).colors[:]
        else:
            C=brewer2mpl.get_map("Set1", 'qualitative', 9).colors[:]
            C+=brewer2mpl.get_map("Set3", 'qualitative', n-9).colors[:]
        C=C[:n]
        if len(C)<n:
            C=[C[i%len(C)] for i in range(n)]
        C=[",".join([str(x) for x in X]) for X in C]
        return C

    def plot_membership(self, t_membership, t_go=None, outputdir=None, outputfile="CircosPlot", s_symbol="Symbol", S_label=None, S_color=None):
        """t_mem is a table, where index is Gene ID, contains an optional column for Symbol,
        a list of columns starts with _MEMBER_, which is 1/0 indicating membership in certain hit list.  This method is for membership table generated within biolist.py.
        rename is a dict that maps existing list name (without _MEMBER_ prefix) into new name. For Circos, should not contain space.
        SIDE EFFECT: if t_go is specified, t_membership will be modified, 0.5 added to where GO-indirect membership is found!!!
        """
        #sw=util.StopWatch()
        MAX_NOF_LISTS=100
        D=ms.mssetting.circos['DELIMITER']
        t_mem=t_membership.copy()
        S_col=[x for x in t_mem.header() if x.startswith('_MEMBER_')]
        if len(S_col)>MAX_NOF_LISTS:
            S_col=S_col[:MAX_NOF_LISTS]
            util.warn_msg('Too many lists, only show the first %d in Circos plot!' % MAX_NOF_LISTS)
        if len(S_col)<=1:
            util.warn_msg('Need at least two lists to plot Circos plot, only %d found!' % len(S_col))
            return self.plot(karyotype="", symbol="", links=None, hits="", outputdir=outputdir, outputfile=outputfile)
        if S_label is None:
            S_label=[re.sub(r'[.;]', ' ', x[8:]) for x in S_col]
        S_chr=[str(i+1) for i in range(len(S_col))]
        c_rename={ x: S_chr[i] for i,x in enumerate(S_col) }
        t_mem.rename2(c_rename)
        n=len(S_chr)
        l_symbol=s_symbol is not None and (s_symbol in t_mem.header())
        if not l_symbol:
            t_mem['Symbol']=t_mem.index.astype(str)
        t_mem['_SUM_']=t_mem.loc[:, S_chr].sum(axis=1)
        if S_color is None:
            S_color=Circos.get_qualitative_colors(n)
        else:
            # convert color from hex to RGB
            S_color=[Circos.hex_to_rgb(x) for x in S_color][:n]
        c_coord={}
        s_karyotype=""
        s_symbol="" if l_symbol else None
        s_hits=""
        #sw.check('Prepare Membership')
        for i,x in enumerate(S_chr):
            s_karyotype+="chr"+D+"-"+D+S_chr[i]+D+S_label[i]+D+"0"+D+str(int(t_mem.loc[:,x].sum()))+D+S_color[i%len(S_color)]+"\n"
            tmp=t_mem[t_mem.loc[:, x]>0].copy()
            #tmp["_PATTERN_"]=tmp["_PATTERN_"].apply(lambda x: x[i:]+x[:i])
            tmp.sort_values(['_SUM_', '_PATTERN_', 'Symbol'], ascending=[False, False, True], inplace=True)
            j=0
            for idx in tmp.index:
                j+=1
                if l_symbol:
                    s_symbol+=S_chr[i]+D+str(j-1)+D+str(j)+D+t_mem.loc[idx, 'Symbol']+"\n"
                if idx not in c_coord:
                    c_coord[idx]={}
                c_coord[idx][S_chr[i]]=j
                s_hits+=S_chr[i]+D+str(j-1)+D+str(j)+D+str(int(t_mem.loc[idx, '_SUM_']))+"\n"
        #sw.check('prepare Hit')
        tmp=t_mem[t_mem._SUM_ > 1].copy()
        tmp.sort_values('_SUM_', inplace=True)
        link_cnt=0
        out=[]
        for s_gene,v in c_coord.items():
            S_on_chr=list(v.keys())
            for i in range(len(S_on_chr)):
                chr1=S_on_chr[i]
                i_from=c_coord[s_gene][chr1]
                for j in range(i+1, len(S_on_chr)):
                    chr2=S_on_chr[j]
                    i_to=c_coord[s_gene][chr2]
                    out.append("e"+str(link_cnt)+D+chr1+D+str(i_from-1)+D+str(i_from))
                    out.append("e"+str(link_cnt)+D+chr2+D+str(i_to-1)+D+str(i_to))
                    #s_link+=str(link_cnt)+D+chr1+D+str(i_from-1)+D+str(i_from)+D+str(link_cnt)+D+chr2+D+str(i_to-1)+D+str(i_to)+"\n"
                    link_cnt+=1
        s_link="\n".join(out)
        S_link=[s_link]
        #sw.check('prepare link')

        def go_links(t_go, l_sample=False):
            link_cnt=0
            c_plot={}
            out=[]
            MAX_EDGES=10000
            MAX_GENES=20
            s_link=None
            for _,r in t_go.iterrows():
                s=r['GeneID']
                #print '>>>>', r['GO']
                #if r['GO']!='GO:0007156': continue
                S=s.split("|")
                if link_cnt>MAX_EDGES: # we start to get aggressive
                    if not l_sample: return None # the caller should call us again with l_sample=True
                    S=np.random.choice(S, MAX_GENES)
                for i in range(len(S)):
                    g1=S[i]
                    for j in range(i+1, len(S)):
                        g2=S[j]
                        #if r['GO']=='GO:0007156':
                        #    print g1+":"+g2, g1+":"+g2 in c_plot, g2+":"+g1 in c_plot
                        if g1+":"+g2 in c_plot or g2+":"+g1 in c_plot: continue
                        c_plot[g1+":"+g2]=True
                        if g1 not in c_coord or g2 not in c_coord: continue # if nof lists exceed MAX_OF_LISTS, c_coord may not contain all hits
                        for chr1,i_from in c_coord[g1].items():
                            #if r['GO']=='GO:0007156': print ">>> "+str(chr1)+" "+str(i_from)
                            for chr2,i_to in c_coord[g2].items():
                                #if r['GO']=='GO:0007156': print "   <<< "+str(chr2)+" "+str(i_to)
                                if chr1==chr2: continue # same chromosome
                                out.append("e"+str(link_cnt)+D+chr1+D+str(i_from-1)+D+str(i_from))
                                out.append("e"+str(link_cnt)+D+chr2+D+str(i_to-1)+D+str(i_to))
                                #str(link_cnt)+D+chr1+D+str(i_from-1)+D+str(i_from)+D+str(link_cnt)+D+chr2+D+str(i_to-1)+D+str(i_to))
                                link_cnt+=1
                                # use 0.4 instead 0.5 to avoid degeneracy in Euclidean distance calculation
                                if t_membership.loc[g1, S_col[int(chr2)-1]]==0: t_membership.loc[g1, S_col[int(chr2)-1]]=0.333
                                if t_membership.loc[g2, S_col[int(chr1)-1]]==0: t_membership.loc[g2, S_col[int(chr1)-1]]=0.333
            s_link="\n".join(out)
            return s_link

        if t_go is not None:
            for x in S_col:
                t_membership[x]=t_membership[x].astype(float)
            s_link=go_links(t_go, False)
            if s_link is None:
                s_link=go_links(t_go, True)
            S_link.append(s_link)
            #sw.check('GO links')
        return self.plot(karyotype=s_karyotype, symbol=s_symbol, links=S_link, hits=s_hits, outputdir=outputdir, outputfile=outputfile)

    def plot(self, karyotype="", symbol=None, links=None, hits=None, outputdir=None, outputfile="CircosPlot"):
        #sw=util.StopWatch()
        outputdir=outputdir if outputdir is not None else '/tmp'
        for ext in [".png" , ".svg"]:
            s_file=os.path.join(outputdir, outputfile+ext)
            if os.path.exists(s_file):
                os.remove(s_file)
        if links is None:
            util.warn_msg('No link to plot, simply ignore')
            #return
        tmp=tempfile.NamedTemporaryFile(dir=outputdir, delete=False, prefix="CIRCOS_", suffix=".txt")
        conf_file=tmp.name
        S_tmp_file=[conf_file]

        s_conf=util.read_string(self.TEMPLATE)
        kary_file=re.sub('CIRCOS_', 'CIRCOS_KARYOTYPE_', conf_file)
        util.save_string(kary_file, karyotype)
        S_tmp_file.append(kary_file)
        s_conf=re.sub(r'@KARYOTYPE@', kary_file, s_conf)

        r0=0.90
        s_plot=""
        if hits is None:
            hits=[]
        elif type(hits) is not list:
            hits=[hits]
        for i,s_hit in enumerate(hits):
            hit_file=re.sub('CIRCOS_', 'CIRCOS_HIT_%d_' % i, conf_file)
            util.save_list(hit_file, s_hit)
            S_tmp_file.append(hit_file)
            s_plot+="<plot>\n"
            s_plot+="file = "+hit_file+"\n"
            s_plot+="r0 = "+('%.3f' % r0)+"r\n"
            s_plot+="r1 = "+('%.3f' % r0)+"r+70p\n"
            s_plot+="stroke_thickness = 0\n"
            s_plot+="min = 0\n"
            s_plot+="max = 2\n"
            s_plot+="color = oranges-3-seq\n"
            s_plot+="</plot>\n\n"
            r0-=0.05;

        s_conf=re.sub(r'@PLOTS@', s_plot, s_conf)
        #t_chr=pd.read_csv(os.path.join(Circos.HOME, "karyotype_"+pid+".tmp"), sep=r'\s+', header=None)
        #s_conf=re.sub(r'@CHROMOSOMES@', ";".join(t_chr[2]), s_conf)
        #avoid using Pandas, so that this script can be used in CGI on ldweb server, where numpy is not installed correctly
        S=karyotype.split("\n")
        S_chr=[]
        for s in S:
            if s.strip()=='': break
            S_chr.append(re.split(ms.mssetting.circos['DELIMITER'], s)[2])
        s_conf=re.sub(r'@CHROMOSOMES@', ";".join(S_chr), s_conf)
        s_symbol=""
        if symbol is not None:
            symbol_file=re.sub('CIRCOS_', 'CIRCOS_SYMBOL_', conf_file)
            util.save_string(symbol_file, symbol)
            S_tmp_file.append(symbol_file)
            s_symbol+="<plot>\n"
            s_symbol+="type = text\n"
            s_symbol+="color = black\n"
            s_symbol+="file = "+symbol_file+"\n"
            s_symbol+="r0=1.02r\n"
            s_symbol+="r1=1.2r\n"
            s_symbol+="label_size = 12p\n"
            s_symbol+="label_font = condensed\n"
            s_symbol+="padding = 0p\n"
            s_symbol+="rpadding = 0p\n"
            s_symbol+="</plot>\n"
        s_conf=re.sub(r'@SYMBOL@', s_symbol, s_conf)

        #S_color=['107,174,214,0.85', '116,196,118,0.85', '106,81,163,0.85']
        S_color=['107,174,214', '116,196,118', '106,81,163']
        s_link=""
        MAX_EDGES=10000 # Circos does not seem to work well with too many edges, it will not draw edges after maybe 20000-ish
        if links is not None:
            if type(links) is str: links=[links]
            for i in range(len(links)-1, -1, -1):
                link_file=re.sub('CIRCOS_', 'CIRCOS_LINK%02d_' % (i+1), conf_file)
                S_tmp_file.append(link_file)
                S_edges=links[i].strip().split("\n")
                n_edge=len(S_edges)/2
                if n_edge>MAX_EDGES:
                #if True:
                    # we need to make sure circos.conf.template has record_limit set to > MAX_EDGES*2
                    # otherwise, we should always shuffle links to avoid the dropping of edges appear later in the file
                    # randomly sample a subset
                    IDX=np.repeat(np.random.permutation(list(range(0,len(S_edges),2)))[:MAX_EDGES], 2)
                    IDX[list(range(1,len(IDX),2))]+=1
                    #S_in=[x for x in IDX if  x >=len(S_edges) ]
                    S_edges=pd.Series(S_edges)[IDX].astype(str)
                    links[i]="\n".join(S_edges)
                util.save_string(link_file, links[i])
                s_link+="<link link"+str(i+1)+">\n"
                s_link+="show = yes\n"
                s_link+="color = "+S_color[(i+len(S_color)-1)%len(S_color)]+"\n"
                s_link+="file = "+link_file+"\n"
                s_link+="</link>\n\n"
        s_conf=re.sub(r'@LINKS@', s_link, s_conf)
        #print s_conf
        util.save_string(conf_file, s_conf)
        #print s_conf
        ## run Circos
        #s_cmd = "cd "+os.path.join(os.path.dirname(__file__), "circos")+"; "
        s_cmd=self.BIN+" -conf "+conf_file
        s_cmd+=' -outputdir '+outputdir
        s_cmd+=' -outputfile '+outputfile
        #sw.check('prepare conf file')
        print(s_cmd)
        util.unix(s_cmd, l_print=False, l_error=True)
        l_remove_temp=True
        if l_remove_temp:
            for f in S_tmp_file:
                os.remove(f)
        s_file=os.path.join(outputdir, outputfile+".png")
        #sw.check('make circos image')
        if os.path.exists(s_file):
            return s_file
        else:
            return None

    @staticmethod
    def clean_tmp(outputdir="/tmp", l_remove_plot=False):
        for f in glob.glob(outputdir+"/CIRCOS_*.txt"):
            os.remove(f)
        if l_remove_plot:
            for f in glob.glob(outputdir+"/CIRCOS_*.png"):
                os.remove(f)
            for f in glob.glob(outputdir+"/CIRCOS_*.svg"):
                os.remove(f)

if __name__=='__main__':
    c=Circos(BIN="~/circos-0.67-7/bin/circos", TEMPLATE="~/circos/circos.conf.template")
    s_dir="~/python/lib/circos/example"
    s_karyotype=util.read_string(s_dir+"/karyotype.txt")
    s_symbol=util.read_string(s_dir+"/symbol.txt")
    s_links=util.read_string(s_dir+"/link.txt")
    c.plot(karyotype=s_karyotype, symbol=s_symbol, links=s_links, outputdir="/tmp", outputfile="CircosPlot")
   # Circos.clean_tmp()

