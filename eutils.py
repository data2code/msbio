#!/usr/bin/env python
from __future__ import absolute_import
from __future__ import print_function
import os
import requests
import pytz
import datetime
import time
import util
import xml.etree.ElementTree as ET
import re
from six.moves import range
import setting
#import progress_bar as pb

class EUtils(object):
    url_base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils'

    ### code from here to _fetch method are taken and modified based on eutils-0.0.9 package at
    # https://pypi.python.org/pypi/eutils
    # The package has good implementation of throttling, Bio.Entrez has no throttling at all
    # However, the package misses batch processing feature and do not make use of history
    # Therefore, we need to write a new package on our own

    def __init__(self, def_args=None, debug=False):
        # we set sehistory=yes by default
        self.def_args={'retmode': 'xml', 'usehistory': 'y', 'retmax': 250, 'email': setting.eutils['EMAIL'], 'tool':'Metascape eUtils Client'}
        if def_args is not None:
            self.def_args=def_args
        self.eastern_tz = pytz.timezone('US/Eastern')
        self._last_request_clock = 0
        self._request_count = 0
        self.debug=debug
        self._giveup=3

    def request_interval(self, utc_dt=None):
        # From https://www.ncbi.nlm.nih.gov/books/NBK25497/:
        # "In order not to overload the E-utility servers, NCBI recommends that
        # users post no more than three URL requests per second and limit
        # large jobs to either weekends or between 9:00 PM and 5:00 AM Eastern
        # time during weekdays."
        # Translation: Weekdays 0500-2100 => 0.333s between requests; no throttle otherwise
        default_request_interval = 0.333
        if utc_dt is None:
            utc_dt = datetime.datetime.utcnow()
        eastern_dt = self.eastern_tz.fromutc( utc_dt )
        return default_request_interval if (0<=eastern_dt.weekday()<=4 and 5<=eastern_dt.hour<21) else 0

    # Please check https://www.ncbi.nlm.nih.gov/books/NBK25499/
    # for syntax of args dictionary
    #def efetch(self,args={}):
    #    return self._fetch('/efetch.fcgi',args)

    def egquery(self,args={}):
        return self._fetch('/egquery.fcgi',args)

    def einfo(self,args={}):
        return self._fetch('/einfo.fcgi',args)

    # replace elink with batch mode
    #def elink(self,args={}):
    #    return self._fetch('/elink.fcgi',args)

    def epost(self,args={}):
        return self._fetch('/epost.fcgi',args)

    # repleace esearch with batch mode
    #def esearch(self,args={}):
    #    return self._fetch('/esearch.fcgi',args)

    # replace esummary with batch mode
    #def esummary(self,args={}):
    #    return self._fetch('/esummary.fcgi',args)

    def espell(self,args={}):
        return self._fetch('/espell.fcgi',args)

    def ecitmatch(self,args={}):
        return self._fetch('/ecitmatch.cgi',args)

    ############################################################################
    ## Internals
    def _fetch(self,path,args={}):
        """return results for a NCBI query, possibly from the cache

        :param: path: relative query path (e.g., 'einfo.fcgi')
        :param: args: dictionary of query args
        :rtype: xml string

        The args are joined with args required by NCBI (tool and email
        address) and with the default args declared when instantiating
        the client.
        """

        url = EUtils.url_base + path
        if type(args) is dict:
            args = dict( list(self.def_args.items()) + list(args.items()) )
        trial=1
        while True:
            if self.debug:
                print("trial %d" % trial)
            # else args is str, pass as it is
            req_int = self.request_interval()
            sleep_time = req_int - (time.time()-self._last_request_clock)
            #print "Sleep: ", sleep_time
            if sleep_time > 0:
                if self.debug:
                    print("sleep_time %d" % sleep_time)
                time.sleep(sleep_time)
            r = requests.post(url, args)
            self._last_request_clock = time.time()
            self._request_count += 1
            if self.debug:
                print(r.text)

            if not r.ok:
                if trial==self._giveup:
                    if any(bad_word in r.text for bad_word in ['<error>','<ERROR>']):
                        xml = ET.fromstring(r.text.encode('utf-8'))
                        util.error_msg('{r.reason} ({r.status_code}): {error}'.format(
                            r=r, error=xml.find('ERROR').text))
                    else:
                        util.error_msg('{r.reason} ({r.status_code}): {r.error}'.format(
                            r=r, error=r.text))
                else:
                    time.sleep(10)
            else:
                return r.content#.encode('utf-8')
            trial+=1

    def _check_interval(self):
        """Delay a request as needed, if it's too frequent"""
        req_int = self.request_interval()
        sleep_time = req_int - (time.time()-self._last_request_clock)
        if sleep_time > 0:
            time.sleep(sleep_time)

    def _flag_check(self):
        """Flag upon finish of a request"""
        self._last_request_clock = time.time()
        self._request_count += 1

    def _fetch_pcassay(self, args):
        """PubChem BioAssay cannot be obtained via efetch, use web URL to get one at a time instead"""
        if 'id' in args:
            s_id=args['id']
        else:
            # retrieve id list from WebEnv+query_key
            s_id=self.esearch(args)
        if type(s_id) is str:
            s_id=s_id.split(",")
        out=[]
        for id in s_id:
            url="https://pubchem.ncbi.nlm.nih.gov/assay/assay.cgi?aid=%s&version=1.2&q=expdesc_xmldisplay" % str(id)
            trial=1
            while True:
                self._check_interval()
                r = requests.get(url)
                self._flag_check()
                if not r.ok:
                    if trail==self._giveup:
                        util.error_msg('{r.reason} ({r.status_code}): {r.error}'.format(
                            r=r, error=r.text))
                    else:
                        time.sleep(10)
                else:
                    x=re.sub(r'<PC-AssayDescription\s.*?>', '<dummy_tag><PC-AssayDescription>', re.sub(r'\n', '', r.content))
                    out.append(x+'</dummy_tag>')
                    break
                trial+=1
        return out

    def _fetch_pccompound(self, args):
        """PubChem Compound cannot be obtained via efetch, use web URL to get one at a time instead"""
        if 'id' in args:
            s_id=args['id']
        else:
            # retrieve id list from WebEnv+query_key
            s_id=self.esearch(args)
        if type(s_id) is str:
            s_id=s_id.split(",")
        out=[]
        for id in s_id:
            url="https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/%s/XML/?response_type=display" % str(id)
            trial=1
            while True:
                self._check_interval()
                r = requests.get(url)
                self._flag_check()
                if not r.ok:
                    if trail==self._giveup:
                        util.error_msg('{r.reason} ({r.status_code}): {r.error}'.format(
                            r=r, error=r.text))
                    else:
                        time.sleep(10)
                else:
                    x=re.sub(r'<Record\s.*?>', '<dummy_tag><Record>', re.sub(r'\n', '', r.content))
                    out.append(x+'</dummy_tag>')
                    break
                trial+=1
        return out

    def esearch(self, args, max_hits=0):
        """Aims to handle arbitrary return size
        Return:
        list(str) for id list
        dict: {WebEnv, query_key} for later queries"""
        #https://www.ncbi.nlm.nih.gov/books/NBK25499/
        # esearch can also be used to retrieve all IDs for a given WebEnv+query_key
        retmax=100000 # maximum allowed
        retstart=0
        # should not escape, as we are using post
        #if 'term' in args:
        #    args['term']=args['term'].replace(' ', '+')
        out=self._fetch('/esearch.fcgi',
            dict(list(args.items())+[('retmax',retmax),('retstart',retstart)]))
        xml=ET.fromstring(out)
        webenv=xml.find('WebEnv').text
        query_key=xml.find('QueryKey').text
        count=int(xml.find('Count').text)
        max_hits=min(count, max_hits) if max_hits>0 else count
        retmax=int(xml.find('RetMax').text)
        retstart=int(xml.find('RetStart').text)
        ids=[x.text for x in xml.findall('./IdList/Id')]
        print('Total>', count, webenv, query_key)
        while (retstart+retmax < max_hits):
            retstart+=retmax
            out=self._fetch('/esearch.fcgi',
                dict([('db',args['db']),('WebEnv',webenv),('query_key',query_key),('retmax',retmax),('retstart',retstart)]))
            xml=ET.fromstring(out)
            retmax=int(xml.find('RetMax').text)
            retstart=int(xml.find('RetStart').text)
            ids.extend([x.text for x in xml.findall('./IdList/Id')])
        if len(ids)>max_hits:
            ids=ids[:max_hits]
        return ids, {'WebEnv':webenv, 'query_key':query_key, 'count':count}

    def _format_ids(self, s_id):
        if type(s_id) is not str:
            count=len(s_id)
            s_id=",".join(util.iarray2sarray(s_id))
        else:
            s_id=re.sub(r',+', ',', re.sub(r'\D', ',', s_id))
            count=len(s_id.split(","))
        return count, s_id

    def _batch_retrieve(self, action, args, count=0, func=None, retmax=10000):
        """Return dict {WebEnv, query_key}"""
        # action can be efetch or esummary
        # according to https://www.ncbi.nlm.nih.gov/books/NBK25499/
        # maximumly allowed retmax is 10000 for efetch and esummary
        retstart=0
        if 'id' in args:
            count, s_id=self._format_ids(args['id'])
            args['id']=s_id
            args=dict(list(self.def_args.items())+list(args.items()))
            out=self.epost(args)
            del args['id']
            xml=ET.fromstring(out)
            args['WebEnv']=xml.find('WebEnv').text
            args['query_key']=xml.find('QueryKey').text
        else:
            if 'WebEnv' not in args and 'query_key' not in args:
                util.error_msg('Missing id, WebEnv, query_key!')
        #if action=='elink' and 'cmd' in args and args['cmd']=='neighbor_history':
        #    # not very meaningful for our use, as we separate id by id=&id=...,
        #    # it will return one query_key per input id
        #    del args['cmd']
        out=self._fetch("/"+action+'.fcgi', dict(list(args.items())+[('retmax',retmax),('retstart',0)]))
        #if 'cmd' in args and args['cmd']=='neighbor_history':
        #    xml=ET.fromstring(out)
        #    webenv=xml.find('./LinkSet/WebEnv').text
        #    query_key=xml.find('./LinkSet/LinkSetDbHistory/QueryKey').text
        #    return {'WebEnv': webenv, 'query_key': query_key}

        if func is not None and callable(func):
            out=func(out)
        S_xml=[out]
        # if there are more entries than retmax, we need to make additional trip
        while (count>0 and retstart+retmax < count):
            retstart+=retmax
            print("Fetching batch: %d ..." % (retstart))
            out=self._fetch("/"+action+'.fcgi', dict(list(args.items())+[('retmax',retmax),('retstart',retstart)]))
            if func is not None and callable(func):
                out=func(out)
            S_xml.append(out)
        return S_xml

    def efetch(self, args, count=0, func=None):
        """Return a list of XML strings, which can be parsed by Entity.FetchList()"""
        if 'db' in args and args['db']=='pcassay':
            return self._fetch_pcassay(args)
        elif 'db' in args and args['db']=='pccompound':
            return self._fetch_pccompound(args)
        elif 'db' in args and args['db'] in ('mesh','books','structure','pcsubstance'):
            #  these database do not return meaninful XML content, use esummary instead of efetch
            return self.esummary(args, count=count, func=func)
        return self._batch_retrieve(action='efetch', args=args, count=count, func=func, retmax=10000)

    def esummary(self, args, count=0, func=None):
        """Return a list of XML strings, which can be parsed by Entity.SummaryList()"""
        if args['db'] in ('mesh','pubmed'):
            args['version']=2.0
        return self._batch_retrieve(action='esummary', args=args, count=count, func=func, retmax=10000)

    def elink(self, args, count=0):
        """
        All link names are in https://eutils.ncbi.nlm.nih.gov/entrez/query/static/entrezlinks.html
        if cmd=neighbor_history, current code will fail, as WebEnv and query_key are not in the standard location.
        However, the mapping between source and destination IDs are lost, it basically only returns a list of destinate IDs
        Return
        We should not use cmd=neighbor_history
        """
        def func(out):
            try:
                xml=ET.fromstring(out)
            except:
                print(out)
                xml=ET.fromstring(out)

            c_map={}
            for x in xml.findall('./LinkSet'):
                id=x.find('./IdList/Id').text
                Y=x.findall('./LinkSetDb')
                if Y is None: continue # no mapping
                for y in Y:
                    linkname=y.find('./LinkName').text
                    if linkname not in c_map: c_map[linkname]={}
                    c_map[linkname][id]=[ z.text for z in y.findall('./Link/Id')]
            return c_map

        if 'cmd' in args and args['cmd']=='neighbor_history':
            util.error_msg('Should not use cmd=neighbor_history for elink, as it returns one query_key per input id (id=&id=...')

        if 'id' not in args: # let's first get id from WebEnv
            S_id, tmp=self.esearch({'WebEnv':args['WebEnv'], 'query_key':args['query_key'], 'db':args['dbfrom']})
        else:
            S_id=args['id']
        if type(S_id) is str:
            S_id=S_id.split(',')
        count=len(S_id)

        retmax=5000
        retmax=50
        data={}
        args=dict(list(self.def_args.items())+list(args.items()))
        if 'id' in args: del args['id']
        if 'WebEnv' in args: del args['WebEnv']
        if 'query_key' in args: del args['query_key']
        if 'count' in args: del args['count']
        if 'retstart' in args: del args['retstart']
        if 'retmax' in args: del args['retmax']

        args="&".join([ k+"="+str(v) for k,v in args.items() ])

        for retstart in range(0, count, retmax):
            print(">>>>", retstart)
            my_args=args+'&retstart='+str(retstart)+'&'
            my_args+="&".join([ "id="+str(x) for x in S_id[retstart:(retstart+retmax)]])
            #print my_args
            out=self._fetch('/elink.fcgi', my_args)
            out=func(out)
            for k,v in out.items():
                #print ">>>", len(v)
                if k not in data:
                    data[k]=v
                else:
                    data[k].update(v)
        return data

if __name__=='__main__':

    # try some examples in https://www.ncbi.nlm.nih.gov/books/NBK25497/
    # remember it is much slower to run efetch than esummyar, as efetch contains tons of extra data
    # so if your attributes are in esummary, use that instead
    # If you want to fetch the full list of ids from a histry entry, use esearch+webenv+query_key
    # Below are some examples showing the usage
    eu=EUtils()

    if False:
        # https://www.ncbi.nlm.nih.gov/books/NBK25498/#chapter3.ESearch__ESummaryEFetch
        pubmed_ids, args=eu.esearch({'db':'pubmed', 'term':'asthma[mesh] AND leukotrienes[mesh] AND 2009[pdat]'})
        import Entity.PubMed as pm
        args['db']='pubmed'
        out=eu.esummary(args)
        X=pm.SummaryList(out)
        print(X.to_list(['pubmed_id','title']))
        out=eu.efetch(args)
        X=pm.FetchList(out)
        print(X.to_list(['pubmed_id','abstract']))

    if False:
        # https://www.ncbi.nlm.nih.gov/books/NBK25498/#chapter3.EPost__ESummaryEFetch
        args={'db':'protein', 'id':'194680922,50978626,28558982,9507199,6678417'}
        import Entity.GI as gi
        out=eu.esummary(args)
        X=gi.SummaryList(out)
        print(X.to_list(['gi','title']))
        out=eu.efetch(args)
        X=gi.FetchList(out)
        print(X.to_list(['gi','sequence']))

    if False:
        # https://www.ncbi.nlm.nih.gov/books/NBK25498/#chapter3.ESearch__ELink__ESummaryEFetch
        pubmed_ids, args=eu.esearch({'db':'pubmed', 'term':'asthma[mesh] AND leukotrienes[mesh] AND 2009[pdat]'})
        import Entity.PubMed as pm
        import Entity.GI as gi
        args['db']='protein'
        args['dbfrom']='pubmed'
        args['linkname']='pubmed_protein'
        #args['cmd']='neighbor_history'
        c_map=eu.elink(args)
        print(c_map)
        args['id']=[]
        for v in c_map['pubmed_protein'].values():
            args['id'].extend(v)
        print(args['id'])
        args['db']='protein'
        out=eu.esummary(args)
        X=gi.SummaryList(out)
        print(X.to_list(['gi','title']))
        out=eu.efetch(args)
        X=gi.FetchList(out)
        print(X.to_list(['gi','description']))

    if False:
        # https://www.ncbi.nlm.nih.gov/books/NBK25498/#chapter3.ELink__ESearch
        args={}
        args['db']='gene'
        args['dbfrom']='protein'
        args['linkname']='protein_gene'
        #args['cmd']='neighbor_history'
        args['id']='148596974,42544182,187937179,4557377,6678417'
        c_map=eu.elink(args)
        print(c_map)
        args={}
        args['id']=[]
        for v in c_map['protein_gene'].values():
            args['id'].extend(v)
        args['db']='gene'
        args['term']='human[orgn] AND x[chr]'
        id, args=eu.esearch(args)
        import Entity.Gene as gene
        out=eu.esummary(args)
        X=gene.SummaryList(out)
        print(X.to_list(['gene_id','description','tax_id']))

    if True:
        # Example of fetching Scripps publications in 2014
        id,args=eu.esearch({'db':'pubmed', 'term':'Scripps[Affiliation] AND Research Institute[Affiliation] AND 2014[pdat]'})
        print(args)

        import Entity.PubMed as pm
        args['db']='pubmed'
        out=eu.esummary(args)
        x=pm.SummaryList(out)
        S_id=[y['pubmed_id'] for y in x.to_list(['pubmed_id'])]
        print(S_id)
        # find out what gene targets are studied at Scripps
        c_map=eu.elink({'db':'gene','dbfrom':'pubmed','id':S_id, 'linkname':'pubmed_gene'})
        print(c_map)

        import Entity.Gene as g
        for k,v in c_map['pubmed_gene'].items():
            print(pm.SummaryList(eu.esummary({'db':'pubmed','id':k})).to_list(['title']))
            print(g.SummaryList(eu.esummary({'db':'gene', 'id':v})).to_list(['gene_id','symbol', 'description']))
            print("\n")
        exit()
        x.data=x.data[:5] # use 5 publications as an example
        print(x.to_list(['pubmed_id', 'title', 'journal', 'author']))

        # Example of fetching gene summary text
        id='6774,7089,7124'
        out=eu.esummary({'db':'gene', 'id':id})
        x=g.SummaryList(out)
        print(x.to_list(['gene_id', 'tax_id', 'symbol', 'description', 'summary']))


