from __future__ import absolute_import
import re
import socket

import numpy as np
import pandas as pd
import util
import ms.biolists as bl
import ms.meta_engine as me
import db
import os
from pprint import pprint

def _report(myopt, dbcon=None):
    #mylib
    DEFAULT_OPT={
        #'URL_BASE':'',
        'FOLDER_BASE': '.',
        'LISTS':1,
        'ANNOTATION':[],
        'MEMBERSHIP': None,
        'MEMBER_SOURCE':[],
        'IMG_MEMBERSHIP':None,
        'ENRICHMENT':None,
        'DATAFRAME_LISTS':None,
        'DATABASE_UPDATE':'1/1/2016',
        'IMG_CIRCO_BY_GENE':None,
        'IMG_CIRCO_BY_GO':None,
        'GO_SOURCE':[],
        'MIN_COUNT':3,
        'P_CUTOFF':0.01,
        'ENRICH_CUTOFF':3,
        'IMG_GO_HEATMAP':None,
        'IMG_GO_HEATMAP_TOP100':None,
        'IMG_GO_HEATMAP_PARENT':None,
        'IMG_GO_NETWORK_CLUSTER':None,
        'IMG_GO_NETWORK_PVALUE':None,
        'IMG_GO_NETWORK_COUNT':None,
        'IMG_QC_HEATMAP':[],
        'TARGET_SPECIES':None,
        'GPEC':False,
        'GPEC_BENCHMARK':0,
        'BACKGROUND_GENE_COUNT':0,
        'engine':None
        #'l_WEB_MODE':True
    }

    REF_METASCAPE='Zhou et al., Metascape provides a biologist-oriented resource for the analysis of systems-level datasets. Nature Communications (2019) 10(1):1523.'
    URL_METASCAPE='http://metascape.org'
    REF_CIRCOS='Krzywinski M. et al. Circos: an Information Aesthetic for Comparative Genomics. Genome Res (2009) 19:1639-1645'
    URL_CIRCOS='http://circos.ca'
    REF_JTREEVIEW='Saldanha AJ. Java Treeview - extensible visualization of microarray data. Bioinformatics (2004) 20:3246-3248'
    URL_JTREEVIEW='http://jtreeview.sourceforge.net'
    REF_CYTOSCAPE='Shannon P. et al., Cytoscape: a software environment for integrated models of biomolecular interaction networks. Genome Res (2003) 11:2498-2504.'
    URL_CYTOSCAPE='http://www.cytoscape.org'
    REF_DAVID='Huang DW, et al. Systematic and integrative analysis of large gene lists using DAVID Bioinformatics Resources. Nature Protoc. (2009) 4:44-57.'
    REF_PVALUE='Zar, J.H. Biostatistical Analysis 1999 4th edn., NJ Prentice Hall, pp. 523'
    REF_QVALUE='Hochberg Y., Benjamini Y. More powerful procedures for multiple significance testing. Statistics in Medicine (1990) 9:811-818.'
    REF_KAPPA='Cohen, J. A coefficient of agreement for nominal scales. Educ. Psychol. Meas. (1960) 20:27-46.'
    URL_DAVID='https://david.ncifcrf.gov'
    URL_CYTOSCAPE_JS_CLIENT = ''
    REF_MCODE='Bader, G.D. et al. An automated method for finding molecular complexes in large protein interaction networks. BMC bioinformatics (2003) 4:2.'
    REF_BIOGRID='Stark C. et al. BioGRID: a general repository for interaction datasets. Nucleic Acids Res. (2006) 34:D535-539.'
    REF_INWEB_IM='Li T. et al. A scored human protein-protein interaction network to catalyze genomic interpretation. Nat. Methods. (2017) 14:61-64.'
    REF_OMNIPATH='Turei D. et al. A scored human protein-protein interaction network to catalyze genomic interpretation. Nat. Methods. (2016) 13:966-967.'
    REF_ROC='https://en.wikipedia.org/wiki/Sensitivity_and_specificity'
    REF_QC={
        'TTRUST':'Han H, et al. TRRUST v2: an expanded reference database of human and mouse transcriptional regulatory interactions. Nucleic acids research 46, D380-D386 (2018).',
        'DisGeNET':'Pinero J, et al. DisGeNET: a comprehensive platform integrating information on human disease-associated genes and variants. Nucleic acids research 45, D833-D839 (2017).',
        'PaGenBase':'Pan JB, et al. PaGenBase: a pattern gene database for the global and dynamic understanding of gene function. PLoS One 8, e80747 (2013).',
        'L1000':'Subramanian A, et al. A Next Generation Connectivity Map: L1000 Platform and the First 1,000,000 Profiles. Cell 171, 1437-1452 e1417 (2017).',
        'MSigDB':'Subramanian A, et al. Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles. Proc Natl Acad Sci U S A 102, 15545-15550 (2005).',  
        'COVID':'https://metascape.org/COVID.',
    }
    REF_STRING='Szklarczyk D. et al. STRING v11: protein-protein association networks with increased coverage, supporting functional discovery in genome-wide experimental datasets. Nucleic Acids Res. (2019) 47:D607-613.'

    opt=DEFAULT_OPT.copy()
    opt.update(myopt)
    i_tbl=0
    i_fig=0
    i_ref=0
    S_ref=[]

    engine=opt['engine'] #None

    # list of annotation columns
    S_ann=opt['ANNOTATION']
    s_mem=opt['MEMBERSHIP']

    s_title="""<div class="container">
        <h2>Metascape Gene List Analysis Report</h2><p/>"""
    i_ref+=1
    s_title+='  <a href="'+URL_METASCAPE+'">metascape.org</a><sup>%d</sup><p/>\n' % (i_ref)
    S_ref.append(REF_METASCAPE)

    s_body=s_sum=s_ovlp=s_ann=s_mem=s_enrh=s_cytoscape_js_client=""

    t_lists=opt['DATAFRAME_LISTS']
    if t_lists is None or len(t_lists)==0:
        return ""
    S_name=t_lists['Name'].tolist()
    tmp=t_lists[ t_lists.Unique > 3000 ].copy()

    def skip_lists(tmp):
        if len(tmp):
            S_name=list(tmp['Name'])
            if len(S_name)>=2:
                S_name[-2]=S_name[-2]+" and "+S_name[-1]
                S_name=S_name[:-1]
            return "Note: gene lists exceeding the size limit of 3000 are ignored in pathway and process enrichment analysis; this applies to the following list(s): %s.<p/>" % (", ".join(S_name))
        else:
            return ""

    if opt['IMG_GO_HEATMAP']:
        go_heatmap_link='''
                <a href='{0}' title='download PDF file'>
                <img class='link' src='{1}' >
                </a>
                '''.format("@URL_BASE@/"+opt['IMG_GO_HEATMAP'].replace('.png','.pdf'), '@ICON_URL@/PDF48.png')
        i_fig+=1
        graph_type="Heatmap" if opt['LISTS']>1 else "Bar Graph"

        s_extra=""
        if opt['IMG_GO_HEATMAP_TOP100']:
            s_name=opt['IMG_GO_HEATMAP_TOP100'].replace('.png', '')
            s_extra='''
                <tr><td>Metascape only visualizes the top 20 clusters. Up to 100 enriched clusters can be viewed here.<br>
                <a href="{0}"><img src="{0}" style="height:auto; width:auto; max-height:100px;"></a>
                <a href='{1}' title='download PDF file'>
                <img class='link' src='{2}' >
                </a></td></tr>
                '''.format("@URL_BASE@/"+s_name+'.png', "@URL_BASE@/"+s_name+'.pdf', '@ICON_URL@/PDF48.png')
        if opt['IMG_GO_HEATMAP_PARENT']:
            s_name=opt['IMG_GO_HEATMAP_PARENT'].replace('.png', '')
            s_extra+='''
                <tr><td>The top-level Gene Ontology biological processes can be viewed here.<br>
                <a href="{0}"><img src="{0}" style="height:auto; width:auto; max-height:100px;"></a>
                <a href='{1}' title='download PDF file'>
                <img class='link' src='{2}' >
                </a></td></tr>
                '''.format("@URL_BASE@/"+s_name+'.png', "@URL_BASE@/"+s_name+'.pdf', '@ICON_URL@/PDF48.png')

        s_sum+="""
              <h3>%s Summary</h3>
              <div class="panel panel-info">
                <div class="panel-heading">Figure %d. %s of enriched terms across input gene lists, colored by p-values.</div>
                <div class="panel-body">

                <table>
                    <tr>
                        <td>
                            <img src="%s" style="width:1000px;">
                        </td>
                    </tr>
                    <tr>
                        <td align='center'>
                            %s
                        </td>
                    </tr>
                    %s
                </table>
                </div>
              </div><p/>""" % (graph_type, i_fig, graph_type.capitalize(),
                               "@URL_BASE@/"+opt['IMG_GO_HEATMAP'],
                               go_heatmap_link,
                               s_extra
                               )
        if opt['LISTS']>1:
            i_ref+=1
            S_ref.append(REF_JTREEVIEW)
            s_sum+="""The heatmap can be interactively viewed using <a href="%s">JTreeView</a><sup>%d</sup> (.cdt, .gtr and .atr files can be found in the Zip package).<p/>""" % (URL_JTREEVIEW, i_ref)

        #print skip_lists(tmp)
        s_sum+=skip_lists(tmp)
    else:
        if len(tmp):
            s_sum+="<h3>Enrichment Summary</h3>"+skip_lists(tmp)

    i_tbl+=1
    #ZHOUZHOU
    s_sum+="""
      <h3>Gene Lists</h3>
      User-provided gene identifiers are first converted into their corresponding %s Entrez gene IDs using the latest version of the database (last updated on %s).  If multiple identifiers correspond to the same Entrez gene ID, they will be considered as a single Entrez gene ID in downstream analyses.   """ \
           % (opt['TARGET_SPECIES'],opt['DATABASE_UPDATE'])
    if opt['LISTS']>1: # multiple list
        s_sum+="""Each gene list is assigned a unique color, which is used throughout the analysis.   """
    s_sum+="""The gene lists are summarized in Table %d.<p/>"""  % (i_tbl)
    S_color = engine.get('S_color')
    if opt['LISTS']>1: # generate color coding
        S=['<div style="background-color:%s;width:80px;height:20px;"></div>' % x for x in S_color]
        opt['DATAFRAME_LISTS']['Color Code']=S
        if '_Color_' in opt['DATAFRAME_LISTS']:
            opt['DATAFRAME_LISTS'].drop('_Color_', axis=1, inplace=True)

    s_sum+="""
          <div class="panel panel-info">
            <div class="panel-heading">Table %d. Statistics of input gene lists.</div>
            <div class="panel-body">%s</div>
          </div>""" % ( i_tbl, util.html(opt['DATAFRAME_LISTS'], tag_table={'class':'table'}, colors=None, tag_th={'class':'info'}))
    if opt['LISTS']>1: # multiple list
        if opt['IMG_CIRCO_BY_GENE']:
            i_fig+=1
            i_ref+=1
            S_ref.append(REF_CIRCOS)
            l_more=opt['IMG_CIRCO_BY_GO'] is not None
            s_ovlp=""" The overlaps between these lists are shown in a <a href="%s">Circos</a><sup>%d</sup> plot (Figure %d%s).""" % (URL_CIRCOS, i_ref, i_fig, ".a" if l_more else "")


            circo_by_gene_links = """
                                <a href='{0}' title='download SVG file'>
                                    <img class='link' src='{1}' height="48" width="48" >
                                </a>""".format("@URL_BASE@/"+opt['IMG_CIRCO_BY_GENE'].replace('.png','.svg'), '@ICON_URL@/SVG48.png')

            if l_more:
                circo_by_go_links = """
                                    <a href='{0}' title='download SVG file'>
                                        <img class='link' src='{1}' height="48" width="48" >
                                    </a>""".format("@URL_BASE@/"+opt['IMG_CIRCO_BY_GO'].replace('.png','.svg'), '@ICON_URL@/SVG48.png')
                s_ovlp+="  Another useful representation is to overlap genes based on their functions or shared pathways.  The overlaps between gene lists can be significantly improved by considering overlaps between genes sharing the same enriched ontology term(s)(Figure %d.b).  Only ontology terms that contain less than 100 genes are used to calculate functional overlaps to avoid linking genes using very general annotation. (We do not want to link all genes, only genes that belong to specific biological processes.)"
                s_ovlp+="""
                      <div class="panel panel-info">
                        <div class="panel-heading">Figure %d. Overlap between gene lists: (a) only at the gene level, where purple curves link identical genes; (b) including the shared term level, where blue curves link genes that belong to the same enriched ontology term.  The inner circle represents gene lists, where hits are arranged along the arc.  Genes that hit multiple lists are colored in dark orange, and genes unique to a list are shown in light orange.  The publication-quality version of the figures is included in the Zip package as a .svg file under the Overlap_circos folder (readable by popular web browsers and Adobe Illustrator).</div>
                        <div class="panel-body">

                        <table>
                            <tr>
                                <td>
                                    <img src="%s" style="width:500px;">
                                </td>
                                <td>
                                    <img src="%s" style="width:500px;">
                                </td>
                            </tr>
                            <tr>
                                <td align='center'>
                                    %s
                                </td>
                                <td align='center'>
                                    %s
                                </td>
                            </tr>
                        </table>
                        </div>
                      </div><p/>""" % (i_fig,
                                       "@URL_BASE@/"+opt['IMG_CIRCO_BY_GENE'],
                                       "@URL_BASE@/"+opt['IMG_CIRCO_BY_GO'],
                                       circo_by_gene_links,
                                       circo_by_go_links
                                       )
            else:
                s_ovlp+="""<p/>
                          <div class="panel panel-info">
                            <div class="panel-heading">Figure %d. Overlap between gene lists.</div>
                            <div class="panel-body">
                            <table>
                                <tr>
                                    <td>
                                        <img src="%s" style="width:500px;">
                                    </td>
                                </tr>
                                <tr>
                                    <td align='center'>
                                        %s
                                    </td>
                                </tr>
                            </table>

                                </div>
                          </div><p/>""" % (i_fig,
                                           "@URL_BASE@/"+opt['IMG_CIRCO_BY_GENE'],
                                           circo_by_gene_links
                                           )

    if S_ann is not None and len(S_ann)>0:
        i_tbl+=1
        s_ann="""
              <h3>Gene Annotation</h3>
              The following are the list of annotations retrieved from the latest version of the database (last updated on %s) (Table %d).<p/>""" % (opt['DATABASE_UPDATE'], i_tbl)
        t_ann=db.sql_in("select display_name Name, group_name Type, display_description Description from annotation_type where annotation_type_name in (", ") order by display_rank", S_ann, con=dbcon)
        s_ann+="""
              <div class="panel panel-info">
                <div class="panel-heading">Table %d. Gene annotations extracted</div>
                <div class="panel-body">%s</div>
              </div><p/>""" % (i_tbl, util.html(t_ann, tag_table={'class':'table'}, colors=None, tag_th={'class':'info'}))

    if opt['MEMBERSHIP']:
        t=db.sql_in("select category_name Name from term_category where term_category_id in (", ") order by display_order", opt['MEMBER_SOURCE'], con=dbcon)
        s_mem="""
              <h3>Membership Search</h3>
              The search term "%s" has been applied to the following list of ontologies: %s, and the binary 1/0 results have been added to the summary spreadsheet as the following columns: %s.  The Excel file metascape_result.xlsx is included in the Zip package.<p/>""" % (opt['MEMBERSHIP'], ", ".join(list(t['Name'])), '"Membership: '+opt['MEMBERSHIP']+'", "Membership Terms: '+opt['MEMBERSHIP']+'"')
        if opt['IMG_MEMBERSHIP']:
            membership_links='''
                <a href='{0}' title='download PDF file'>
                <img class='link' src='{1}' >
                </a>
                '''.format("@URL_BASE@/"+opt['IMG_MEMBERSHIP'].replace('.png','.pdf'), '@ICON_URL@/PDF48.png')
            i_fig+=1
            s_mem+="""
                      <div class="panel panel-info">
                        <div class="panel-heading">Figure %d. Enrichment of genes matching membership term: %s.  The outer pie shows the number and the percentage of genes in the background that are associated with the membership (in black); the inner pie shows the number and the percentage of genes in the individual input gene list that are associated with the membership.  The p-value indicates whether the membership is statistically significantly enriched in the list.</div>
                        <div class="panel-body">

                        <table>
                            <tr>
                                <td>
                                    <img src="%s" style="width:500px;">
                                </td>
                            </tr>
                            <tr>
                                <td align='center'>
                                    %s
                                </td>
                            </tr>
                        </table>
                        </div>
                      </div><p/>""" % (i_fig,
                                       opt['MEMBERSHIP'],
                                       "@URL_BASE@/"+opt['IMG_MEMBERSHIP'],
                                       membership_links
                                       )


    def get_links(png_url, gene_list_name_ppi=None):
        isPPI = 'Enrichment_PPI/' in png_url;
        linksHtml = """<a href='{0}' title='download PDF file'>
                        <img class='link' src='{1}' >
                        </a>
                     """.format(png_url.replace('.png', '.pdf'), '@ICON_URL@/PDF48.png')
        import re
        #file_url = re.sub('[^=]*$', '', png_url)

        linksHtml += """
                        <a href='{0}' title='download CYS file'>
                        <img class='link' src='{1}' >
                        </a>
                    """.format("@URL_BASE@/"+ ('Enrichment_PPI/MCODE_PPI.cys' if isPPI else 'Enrichment_GO/GONetwork.cys'), '@ICON_URL@/CYS48.png')
        #m=re.search('session_id=([^&]*)&', png_url)
        #if m is not None:
        #    session_id = m.group(1)
        ##elif opt['l_WEB_MODE']:
        ##    print("Missing session_id in URL!")
        ##    return ""
        #else:
        #    session_id=""
        gene_list_name_in_html = ''
        if isPPI:
            isMCODE = '_MCODE_ALL' in png_url
            gene_list_name_in_html = '&isPPI=' + gene_list_name + ',' + ('_MCODE_ALL' if isMCODE else '');
        IS_METASCAPE_DOT_ORG = "metascape.org" in str(socket.getfqdn()).lower()

        cytoscape_js_url=""
        #cytoscape_js_url = (
        url_param=""
        s_network=os.path.splitext(os.path.basename(png_url))[0]
        s_style=re.sub(r'.*_', '', s_network)
        if isPPI:
            if '_MCODE_ALL_' not in s_network:
                s_style+="NoLabel"
            s_network=re.sub(r'ColorByCounts$', 'ColorByCluster', s_network)
            url_param="Network={}&Style={}&isPPI=True".format(s_network, s_style)
        else:
            # style has the same name as the network
            url_param="Network=GONetwork&Style={}".format(s_style)

        cytoscape_js_url = ('/gp/' if IS_METASCAPE_DOT_ORG  else '') + 'Content/CyJS/index.html?session_id=@SESSION_ID@'
        linksHtml += """
                        @ONLINE_BEGIN@
                        <a href='{0}' title='interactive cytoscape' target='_blank'>
                        <img class='link' src='{1}' >
                        </a>
                        @ONLINE_END@
                    """.format(cytoscape_js_url+"&"+url_param,
                               #gene_list_name_in_html,
                               '@ICON_URL@/WEB_CYS48.png'
                               )
        # for offline viewing
        #s_network=re.sub(r'ColorByCounts$', 'ColorByCluster', os.path.splitext(os.path.basename(png_url))[0])
        cytoscape_js_url = "Enrichment_PPI/PPINetwork.html" if isPPI else "Enrichment_GO/GONetwork.html"
        linksHtml += """
                        @OFFLINE_BEGIN@
                        <a href='{0}' title='interactive cytoscape' target='_blank'>
                        <img class='link' src='{1}' >
                        </a>
                        @OFFLINE_END@
                    """.format(cytoscape_js_url+"?"+url_param,
                               '@ICON_URL@/WEB_CYS48.png'
                               )

        return linksHtml

    def get_clustergram(cluster_id):
        #cluster_js_url=('/gp/' if IS_METASCAPE_DOT_ORG  else '') + '#/Content/heatmaponecluster/tkrn1_m2m/'
        cluster_js_url=('/gp/#/heatmaponecluster/@SESSION_ID@/')
        linksHtml = """<a href='{0}{1}' title='view clustergram'>
                        <img class='link' src='{2}' >
                       </a>
                    """.format(cluster_js_url, cluster_id, '@ICON_URL@/WEB_HM48.png')
        return linksHtml

    t_go=opt['ENRICHMENT']

    if t_go is not None and len(t_go)>0:
        t=db.sql_in("select category_name Name from term_category where term_category_id in (", ") order by display_order", opt['GO_SOURCE'], con=dbcon)
        S_cat=list(t['Name'])
        if len(S_cat)>=2:
            S_cat[-2]=S_cat[-2]+" and "+S_cat[-1]
            S_cat=S_cat[:-1]
        i_ref+=3
        S_ref.append(REF_PVALUE)
        S_ref.append(REF_QVALUE)
        S_ref.append(REF_KAPPA)
        s_enrh="""
              <h3>Pathway and Process Enrichment Analysis</h3>
              For each given gene list, pathway and process enrichment analysis has been carried out with the following ontology sources: %s.  %s have been used as the enrichment background. Terms with a p-value &lt; %.2f, a minimum count of %d, and an enrichment factor &gt; %.1f (the enrichment factor is the ratio between the observed counts and the counts expected by chance) are collected and grouped into clusters based on their membership similarities.  More specifically, p-values are calculated based on the accumulative hypergeometric distribution<sup>%d</sup>, and q-values are calculated using the Banjamini-Hochberg procedure to account for multiple testings<sup>%d</sup>.  Kappa scores<sup>%d</sup> are used as the similarity metric when performing hierachical clustering on the enriched terms, and sub-trees with a similarity of > 0.3 are considered a cluster.  The most statistically significant term within a cluster is chosen to represent the cluster.<p/>""" \
                           % (", ".join(S_cat),
                              'A user-supplied list of {0} genes'.format(opt['BACKGROUND_GENE_COUNT'])
                                  if opt['BACKGROUND_GENE_COUNT']
                                  else  'All genes in the genome',
                              opt['P_CUTOFF'],
                              opt['MIN_COUNT'],
                              opt['ENRICH_CUTOFF'],
                              i_ref-2,
                              i_ref-1,
                              i_ref)
        if opt['LISTS']>1:
            s_enrh+="""
              When multiple gene lists are provided, all lists are merged into one list called "_FINAL".  A term may be found enriched in several individual gene lists and/or in the _FINAL gene list, and the best p-value among them is chosen as the final p-value.  The pathway/process clusters that are found to be of interest (either shared or unique based on specific list enrichment) are used to prioritize the genes that fall into those clusters (membership is presented as 1/0 binary columns in the Excel spreadsheet).  Note that individual gene lists containing more than 3000 genes are ignored during the enrichment analysis to avoid superficial terms; this is because long gene lists are often not random and generally trigger too many terms that are not of direct relevance to the biology under study.<p/>"""
        t2=t_go[(t_go['FirstInGroupByLogP']>0) & (t_go['GROUP_ID']<=20)].copy()
        t2['Clustergram']=t2.GROUP_ID.apply(lambda x: get_clustergram(x))
        # Gene Lists in _PATTERN_ are alphabetically sorted, which may not match the order in S_name
        S_order=[re.sub(r'^_MEMBER_', '', x) for x in t2.header() if x.startswith('_MEMBER_')]
        c_order={ x: i for i,x in enumerate(S_order) }
        if opt['LISTS']>1:
            t2=t2[['_PATTERN_','GO','Category','Description','#GeneInGOAndHitList','%InGO','LogP','Log(q-value)','Clustergram']].copy()
        else:
            t2=t2[['GO','Category','Description','#GeneInGOAndHitList','%InGO','LogP','Log(q-value)','Clustergram']].copy()
        # disable Clustergram for now
        t2.drop('Clustergram', axis=1, inplace=True)
        t2.rename2({'#GeneInGOAndHitList':'Count', '%InGO':'%', 'LogP':'Log10(P)','Log(q-value)':'Log10(q)'})
        if opt['LISTS']>1:
            for idx in t2.index:
                s_pat=t2.loc[idx, '_PATTERN_'][1:]
                t2.loc[idx, '_PATTERN_']=" ".join(['<div style="background-color:%s;width:20px;height:20px;float:left;"></div>' % (S_color[i] if ((x in c_order) and (s_pat[c_order[x]]=="1")) else "#FFFFFF") for i,x in enumerate(S_name)])
        i_tbl+=1
        s="""
          <div class="panel panel-info">
          <div class="panel-heading">Table %d. Top %d clusters with their representative enriched terms (one per cluster).  "Count" is the number of genes in the user-provided lists with membership in the given ontology term. "%%" is the percentage of all of the user-provided genes that are found in the given ontology term (only input genes with at least one ontology term annotation are included in the calculation).  "Log10(P)" is the p-value in log base 10.  "Log10(q)" is the multi-test adjusted p-value in log base 10."""
        if opt['LISTS']>1:
            s+="  __PATTERN__ shows the color code used for the gene lists where the term is found statistically significant, i.e., multiple colors indicate a pathway/process that is shared across multiple lists."
        s+="""</div>
                <div class="panel-body">%s</div>
              </div><p/>"""
        s_enrh+= s % (i_tbl, len(t2), util.html(util.df2sdf(t2, s_format='%.2f'), tag_table={'class':'table'}, colors=None, tag_th={'class':'info'}))


        if opt['IMG_GO_NETWORK_CLUSTER']:
            go_network_cluster_links = get_links("@URL_BASE@/"+opt['IMG_GO_NETWORK_CLUSTER'])
            go_network_pvalue_links =get_links("@URL_BASE@/"+opt['IMG_GO_NETWORK_PVALUE'])
            i_ref+=1
            S_ref.append(REF_CYTOSCAPE)
            i_fig+=1
            s_enrh+="""
                  To further capture the relationships between the terms, a subset of enriched terms have been selected and rendered as a network plot, where terms with a similarity &gt; 0.3 are connected by edges.  We select the terms with the best p-values from each of the 20 clusters, with the constraint that there are no more than 15 terms per cluster and no more than 250 terms in total.  The network is visualized using <a href="%s">Cytoscape</a><sup>%d</sup>, where each node represents an enriched term and is colored first by its cluster ID (Figure %d.a) and then by its p-value (Figure %d.b).  These networks can be interactively viewed in Cytoscape through the .cys files (contained in the Zip package, which also contains a publication-quality version as a PDF) or within a browser by clicking on the web icon.  For clarity, term labels are only shown for one term per cluster, so it is recommended to use Cytoscape or a browser to visualize the network in order to inspect all node labels.  We can also export the network into a PDF file within Cytoscape, and then edit the labels using Adobe Illustrator for publication purposes.  To switch off all labels, delete the "Label" mapping under the "Style" tab within Cytoscape, and then export the network view.<p/>""" % (URL_CYTOSCAPE, i_ref, i_fig, i_fig)
            s_enrh+="""
                  <div class="panel panel-info">
                    <div class="panel-heading">Figure %d. Network of enriched terms: (a) colored by cluster ID, where nodes that share the same cluster ID are typically close to each other; (b) colored by p-value, where terms containing more genes tend to have a more significant p-value.</div>
                    <div class="panel-body">
                    <table>
                        <tr>
                            <td>
                                <img src="%s" style="width:500px;">
                            </td>
                            <td>
                                <img src="%s" style="width:500px;">
                            </td>
                        </tr>
                        <tr>
                            <td align='center'>
                                %s
                            </td>
                            <td align='center'>
                                %s
                            </td>
                        </tr>
                    </table>
                    </div>

                  </div><p/>""" % (i_fig,
                                   "@URL_BASE@/"+opt['IMG_GO_NETWORK_CLUSTER'],
                                   "@URL_BASE@/"+opt['IMG_GO_NETWORK_PVALUE'],
                                   go_network_cluster_links,
                                   go_network_pvalue_links
                                   )

        if opt['LISTS']>1 and opt['IMG_GO_NETWORK_COUNT']:
            go_network_count_links = get_links("@URL_BASE@/"+opt['IMG_GO_NETWORK_COUNT'])
            i_fig+=1
            s_enrh+="""
               In the case of when multiple gene lists are provided, the nodes are represented as pie charts, where the size of a pie is proportional to the total number of hits that fall into that specific term.  The pie charts are color-coded based on the gene list identities, where the size of a slice represents the percentage of genes under the term that originated from the corresponding gene list.  This plot is particularly useful for visualizing whether the terms are shared by multiple lists or unique to a specific list, as well as for understanding how these terms associate with each other within the biological context of the meta study (Figure %d).<p/>""" % (i_fig)
            s_enrh+="""
                      <div class="panel panel-info">
                        <div class="panel-heading">Figure %d. Network of enriched terms represented as pie charts, where pies are color-coded based on the identities of the gene lists.</div>
                        <div class="panel-body">
                        <table>
                            <tr>
                                <td>
                                    <img src="%s" style="width:1000px;">
                                </td>
                            </tr>
                            <tr>
                                <td align='center'>
                                    %s
                                </td>
                            </tr>
                        </table>
                        </div>
                      </div><p/>""" % (i_fig,
                                       "@URL_BASE@/"+opt['IMG_GO_NETWORK_COUNT'],
                                       go_network_count_links
                                       )

    #{ZHOUZHOU}
    s_ppi = ''
    if ('PPI' in opt) and (not opt['PPI']['disablePPI']): #user ppi checked:
        s_cat=opt['PPI']['DATASOURCE']
        #if type(S_cat) is not list: S_cat=[S_cat]
        S_cat=["STRING", "BioGrid", "OmniPath", "InWeb_IM"]
        S_ds=[]
        for x in S_cat:
            i_ref+=1
            S_ds.append("%s<sup>%d</sup>" % (x, i_ref))
            if x == 'BioGrid': S_ref.append(REF_BIOGRID)
            if x == 'InWeb_IM': S_ref.append(REF_INWEB_IM)
            if x == 'OmniPath': S_ref.append(REF_OMNIPATH)
            if x == 'STRING': S_ref.append(REF_STRING)
        i_ref+=1
        S_ref.append(REF_MCODE)
        s_db=""
        s_ppi_blog_url="http://metascape.org/blog/?p=219"
        if s_cat=='PHYSICAL_CORE':
            s_db=f'Only physical interactions in STRING (physical score &gt 0.132) and BioGrid are used (<a href="{s_ppi_blog_url}">details</a>).'
        elif s_cat=='PHYSICAL_ALL':
            s_db=f'Only physical interactions in STRING and BioGrid are used (<a href="{s_ppi_blog_url}">details</a>).'
        elif s_cat=='COMBINED_CORE':
            s_db=f'All interactions in STRING (combined score &gt 0.187) in STRING are used (<a href="{s_ppi_blog_url}">details</a>).'
        elif s_cat=='COMBINED_ALL':
            s_db=f'All interactions in STRING are used (<a href="{s_ppi_blog_url}">details</a>).'


        s_ppi+="<h3>Protein-protein Interaction Enrichment Analysis</h3>\n"
        s_ppi+=("For each given gene list, protein-protein interaction enrichment analysis has been carried out with the following databases: %s.%s  " +\
                "The resultant network contains the subset of proteins that form physical interactions with at least one other member in the list.  "+\
                "If the network contains between %s and %s proteins, the Molecular Complex Detection (MCODE) algorithm<sup>%d</sup> has been applied to identify densely connected network components.  ")\
            % (", ".join(S_ds), s_db, opt['PPI']['minSize'], opt['PPI']['maxSize'], i_ref)

        i_nof_ppi=len(opt['PPI']['img_fn'])
        if i_nof_ppi>0:
            i_fig+=1
            s_ppi +="""The MCODE networks identified for individual gene lists have been gathered and are shown in Figure %d.<p/>\n""" \
                    %(i_fig)
#            if i_nof_ppi>1:
#                s_ppi+="""
#  In the case of multiple user-provided gene lists, all gene lists were also merged and then analyzed as a combined gene list.  MCODE algorithm was applied to identify densely connected network components within this combined network as well, the results is shown in Figure %d.<p/>\n""" \
#                    % (i_fig)
            s_ppi+="""
                    Pathway and process enrichment analysis has been applied to each MCODE component independently, and the three best-scoring terms by p-value have been retained as the functional description of the corresponding components, """+\
                    """shown in the tables underneath corresponding network plots within Figure %d.<p/>\n""" \
                    % (i_fig)
            df_top3 = None
            if os.path.exists(opt['PPI']['GO_MCODE_Top3_fn']):
                df_top3 = pd.read_csv(opt['PPI']['GO_MCODE_Top3_fn'])
            df_mcode = None
            if os.path.exists(opt['PPI']['MCODE_CSV_fn']):
                df_mcode = pd.read_csv(opt['PPI']['MCODE_CSV_fn'])
            df_final_mcode = None
            if os.path.exists(opt['PPI']['_FINAL_MCODE_CSV_fn']):
                df_final_mcode = pd.read_csv(opt['PPI']['_FINAL_MCODE_CSV_fn'])

            def get_top3(list_name,isMCODE):
                df_mcode_which = df_final_mcode if '_FINAL' in list_name else df_mcode
                if '_FINAL' in list_name:
                    list_name = '_FINAL'
                if df_top3 is None:
                    return ''
                rows=[]
                if isMCODE:
                    df1 = df_top3[df_top3['Network'].apply(lambda x: x.startswith(list_name) and re.search(r'_SUB\d+_MCODE', x) is not None)]
                else:
                    df1 = df_top3[df_top3['Network'] == list_name]


                S_color = ["#BCBDDC"]
                if df_mcode_which is not None:
                    n_mcode = df_mcode_which[df_mcode_which['Network'].apply(lambda x: x.startswith(list_name))]['Cluster'].max()
                    if n_mcode is not np.NAN:
                        S_color=bl.CytoPlot.get_qualitative_colors(n_mcode)

                for i, r in df1.iterrows():
                    if re.search(r'_SUB\d+_MCODE_', r['Network']):
                        mcode_index = re.search('_(\d+)$', r['Network']).group(1)
                        color_bar = S_color[int(mcode_index)-1]
                    else:
                        color_bar = S_color[0]
                    for ann in r['Annotation'].split(';'):
                        t = ann.split('|')
                        color = '<div style="background-color:{0};width:40px;height:20px;"></div>'.format(color_bar)
                        rows.append({
                            'Color': color,
                            'MCODE': re.sub(r'.+_SUB\d+_', '', r['Network']),
                            'GO':t[0],
                            'Description':t[1],
                            'Log10(P)':t[2],
                        })
                if len(rows)>0:
                    if isMCODE:
                        df2 = pd.DataFrame(rows,columns=['Color','MCODE','GO','Description','Log10(P)'])
                    else:
                        df2 = pd.DataFrame(rows,columns=['GO','Description','Log10(P)'])
                    return util.html(df2, tag_table={'class':'table'}, colors=None, tag_th={'class':'info'})
                return  ''
            all_network = ''
            #split all network to individual and merged
            all_gene_list_1 = [ x.replace('_MCODE_ALL','').replace('_PPIColorByCluster.png','').replace('_PPIColorByCounts.png','') for x in opt['PPI']['img_fn']]
            all_gene_list = util.unique(all_gene_list_1)
            all_gene_list.sort()
            all_gene_list = util.minus(all_gene_list, ['_FINAL'])
            if '_FINAL_PPIColorByCounts.png' in opt['PPI']['img_fn']:
                all_gene_list.append('_FINAL_Count')
            if '_FINAL_PPIColorByCluster.png' in opt['PPI']['img_fn']:
                all_gene_list.append('_FINAL_Cluster')
            if opt['LISTS']<=1 and ('_FINAL' in all_gene_list_1):
                all_gene_list = ['_FINAL']

            for gene_list_name in all_gene_list:
                first_line = ''
                second_line = ''
                third_line = ''
                nws = [ x for x in opt['PPI']['img_fn'] if x.startswith(gene_list_name+'_')]
                if gene_list_name == '_FINAL_Count':
                    nws = []
                    if '_FINAL_PPIColorByCounts.png' in opt['PPI']['img_fn']:
                        nws.append('_FINAL_PPIColorByCounts.png')
                    if '_FINAL_MCODE_ALL_PPIColorByCounts.png' in opt['PPI']['img_fn']:
                        nws.append('_FINAL_MCODE_ALL_PPIColorByCounts.png')
                if gene_list_name == '_FINAL_Cluster':
                    nws = []
                    if '_FINAL_PPIColorByCluster.png' in opt['PPI']['img_fn']:
                        nws.append('_FINAL_PPIColorByCluster.png')
                    if '_FINAL_MCODE_ALL_PPIColorByCluster.png' in opt['PPI']['img_fn']:
                        nws.append('_FINAL_MCODE_ALL_PPIColorByCluster.png')


                nws.sort(reverse=True)
                for i, nw in enumerate(nws):
                    isMCODE = '_MCODE_ALL' in nw
                    if i>0:
                        first_line += '<td width="10px"></td>'
                        second_line +='<td width="10px"></td>'
                        third_line += '<td width="10px"></td>'
                    figure_name = ''
                    if opt['LISTS']>1:
                        gene_list_name_html = gene_list_name
                        colored_by_html = ''
                        if '_FINAL' in gene_list_name:
                            gene_list_name_html = 'All lists merged '
                            if 'ColorByCounts' in nw:
                                colored_by_html = 'Colored by Counts'
                            else:
                                colored_by_html = 'Colored by Cluster'
                        colored_by_html += '(Keep MCODE Nodes Only)' if isMCODE else '(Full Connection)'
                        figure_name =  '<h4>{0} {1} </h4>'.format(gene_list_name_html,colored_by_html)


                    first_line += '''<td><img src="{0}" style="width:500px;"></td>'''.format("@URL_BASE@/Enrichment_PPI/"+nw)
                    second_line += '''<td align="center">{0}</td>'''.format(get_links("@URL_BASE@/Enrichment_PPI/"+nw),gene_list_name)
                    third_line += '''<td align="center" valign="top">{0}{1}</td>'''.format(
                                    figure_name,
                                    get_top3(gene_list_name,isMCODE))

                all_network +='''
                <table border-collapse='collapse'>
                <tr>{0}</tr>
                <tr>{1}</tr>
                <tr>{2}</tr>
                </table>'''.format(first_line, second_line, third_line)
            s_ppi+="""
      <div class="panel panel-info">
        <div class="panel-heading">Figure %d. Protein-protein interaction network and MCODE components identified in the gene lists.</div>
        <div class="panel-body">
        %s
        </div>
      </div><p/>""" % (i_fig, all_network)

    s_gpec=''
    S_gpec=[]
    if opt['GPEC']:
        S_gpec.append("<h3>Gene Prioritization by Evidence Counting (GPEC) (beta feature)</h3>\n")
        S_gpec.append('<p>GPEC is an effective way to identify a subset of genes that are more likely to be of higher quality hits.  As we are still working on a GPEC publication, please be advised of the risk when using GPEC analysis results.  A gene receives an evidence token whenever it is identified as a hit within an input gene list, falls into at least one enriched pathway (pathway size is no more than 100) derived from an input gene list, or is part of the protein-protein interaction network formed by a given input gene list.  The evidence count of a gene is the total number of evidence token it receives.  Given <i>n</i> input gene lists, the maximum possible evidence count is <i>3n</i>.  Our research suggests that the likelihood of a gene being a true biological hit increases as its total evidence count increases.  Therefore, all input genes can be ranked in descending order by their evidence counts, and then various evidence count cutoffs can be applied to generate new gene lists of varying quality.  Shorter gene lists containing the most of the top-ranked genes are of higher quality (i.e., higher precision and lower false discovery rate), and longer gene lists obtained under relaxed cutoffs are more comphrehesive (i.e., higher recall and sensitivity).  Term enrichment analysis is carried out for different cutoffs, and the best p-value per term is retained.  For protein network analysis, an evidence cutoff that yields a _FINAL network of approximately 250 protein nodes is selected, so that the network remains rich yet visually interpretable.  Conceptually, the processes and network components identified by the GPEC algorithm is more robust, as it is based on the subset of genes that have higher evidence counts, compared to relying on the list in which all gene lists are merged into one.  Additional information regarding GPEC can be found on the <a ui-sref="menu.gpec">menu page</a>.<p/>')
    if opt['GPEC'] and opt['GPEC_BENCHMARK'] is not None and opt['GPEC_BENCHMARK']>0:
        n_pos=engine.get('benchmark_pos')
        n_neg=engine.get('benchmark_neg')
        acc=engine.get('ML_accuracy')
        if acc is None:
            S_gpec.append('<p>A _BENCHMARK gene list consisting of {0} genes is provided as "True" hits, based on which {1} positives and {2} negatives are assigned.  Machine learning process has not been carried out in this case.</p>'.format(opt['GPEC_BENCHMARK'], str(n_pos), str(n_neg)))
        else:
            S_gpec.append('<p>A _BENCHMARK gene list consisting of {0} genes is provided as "True" hits, based on which {1} positives and {2} negatives are assigned.  A machine learning algorithm has been employed to identify the optimal weightings of each evidence line to maximize the prediction accuracy ({3}) based on cross validation.  The probability of each gene being a "True" hit is then predicted based on a logistic regression model.  The probability is then used to replace the integer evidence counts, and the same GPEC prioritization algorithm is iteratively applied on the various cutoffs to yield refined term and protein network enrichment analysis results.<p/>'.format(opt['GPEC_BENCHMARK'], str(n_pos), str(n_neg), "%.1f%%" % (acc*100)))

            R_GPEC_Weight=engine.get('R_GPEC_Weight')[1:]
            R_abs=np.abs(R_GPEC_Weight)
            R_GPEC_Weight/=np.max(R_abs)
            S_evi=engine.get('S_evi')
            top3=R_abs.argsort()[-10:][::-1]
            S_top3=[ '%s (%.3f)' % (S_evi[i], R_GPEC_Weight[i]) for i in top3 ]
            s_gpec+='The top %d most weighted evidence lines are: %s.</p>' % (len(S_top3), ", ".join(S_top3))
            i_fig+=1
            i_ref+=1
            S_ref.append(REF_ROC)
            S_gpec.append("""<p>The machine learning results are summarized in Figure %d.  In the precision-recall plot, the precision is defined as TP/(TP+FP) <sup>%d</sup> and recall, a.k.a. TPR, is defined as TP/(TP+FN) <sup>%d</sup>.  If the evidence matrix is effective, we expect that precision decreases as recall increases, which implies that the hits ranked at the top are more likely to be true positives (higher precision).  The AUC (area under curve) in the plot represents the average precision value.  The precision value expected for a random hit prioritization algorithm can be read from the last data point of the step curve, where all genes were classified as hits (recall = 1).  The second plot shows how the number of true positives increases as we include more genes.  A good result would show the TPR curve to be significantly above the diagonal line, where the diagonal line represents the result achieved by random hit prioritization.  The third plot is a bar graph showing the relative weights assigned to each evidence line (up to the top ten lines).  If an evidence line contributes negatively, it will be colored in orange; otherwise, it will be blue.  The last plot shows the predicted probability of benchmark genes.  Since true positives are rare in reality, a mere probability of 20%% is already a fairly impressive number (i.e., one out of five hits is expected to be validated).  A good result would place more benchmark genes (shown as orange "+") towards the right portion of the curve (with probability above, say, 20%%), but the exact cutoff is determined by the validation logistic capacity.  If the majoirty of benchmark genes are below 10%%, it generally implies that the true hit status cannot be effectively predicted from the underlying evidence lines.  If these plots strongly support the effectiveness of evidence lines in retrospectively recapturing benchmark genes, the predicted evidence score and probability can then be used to prioritize gene candidates.</p>""" % (i_fig, i_ref, i_ref))

            evidence_link='''
                <a href='{0}' title='download PDF file'>
                <img class='link' src='{1}' >
                </a>
                '''.format("@URL_BASE@/Evidence/EvidenceWeight.pdf", '@ICON_URL@/PDF48.png')
            S_gpec.append("""
              <div class="panel panel-info">
                <div class="panel-heading">Figure %d. Summary of machine learning results.</div>
                <div class="panel-body">

                <table>
                    <tr>
                        <td>
                            <img src="%s">
                        </td>
                    </tr>
                    <tr>
                        <td align='center'>
                            %s
                        </td>
                    </tr>
                </table>

                </div>
              </div>""" % (i_fig,
                               "@URL_BASE@/Evidence/EvidenceWeight.png",
                               evidence_link
                               ))

        s_gpec="\n".join(S_gpec)

    s_qc=''
    s_covid=os.path.join(opt['FOLDER_BASE'], 'COVID_map.xlsx')
    l_covid=os.path.exists(s_covid)

    if len(opt['IMG_QC_HEATMAP']) or l_covid:
        S_go_category=[re.sub(r'(^.*HeatmapSelectedGO_|\..*$)', '', x) for x in opt['IMG_QC_HEATMAP']]

        S_qc=["<h3>Quality Control and Association Analysis</h3>\n"]
        if l_covid:
            S_qc.append("""Input gene list(s) have been cross compared to all COVID reference gene lists. The gene-level overlap map can be <a href="@URL_BASE@/COVID_map.xlsx"><b style="color:#e6550dd">downloaded here</b></a>.""")
            if "COVID" in S_go_category:
                S_qc.append("Enrichment analysis results are presented in the COVID section below.")
            S_qc.append("<p/>")

        S_fig=[]
        i_fig_start=i_fig+1
        i_L1000_ref=0
        for i,s_cat in enumerate(S_go_category):
            i_fig+=1

            go_heatmap_link='''
                <a href='{0}' title='download PDF file'>
                <img class='link' src='{1}' >
                </a>
                '''.format("@URL_BASE@/"+opt['IMG_QC_HEATMAP'][i].replace('.png','.pdf'), '@ICON_URL@/PDF48.png')

            if s_cat in REF_QC:
                i_ref+=1
                s_ref="<sup>{}</sup>".format(i_ref)
                S_ref.append(REF_QC[s_cat])
            elif s_cat.startswith('L1000'):
                if i_L1000_ref==0:
                    i_ref+=1
                    i_L1000_ref=i_ref
                    S_ref.append(REF_QC['L1000'])
                s_ref="<sup>{}</sup>".format(i_L1000_ref)
            elif s_cat=='Transcription_Factor_Targets':
                i_ref+=1
                s_ref="<sup>{}</sup>".format(i_ref)
                S_ref.append(REF_QC['MSigDB'])
            else:
                s_ref=""

            s_cat2=s_cat.replace(' ', '_')
            s_qc_file=os.path.join(opt['FOLDER_BASE'], 'Enrichment_QC/GO_{}.csv'.format(s_cat2))
            if not os.path.exists(s_qc_file):
                util.warn_msg('Missing QC file: '+s_qc_file)
                continue
            t_go=pd.read_csv(os.path.join(opt['FOLDER_BASE'], 'Enrichment_QC/GO_{}.csv'.format(s_cat2)))
            if 'FirstInGroupByLogP' in t_go.header():
                t2=t_go[(t_go['FirstInGroupByLogP']>0) & (t_go['GROUP_ID']<=20)].copy()
            else:
                t2=t_go[:20].copy()
            # Gene Lists in _PATTERN_ are alphabetically sorted, which may not match the order in S_name
            S_order=[re.sub(r'^_MEMBER_', '', x) for x in t2.header() if x.startswith('_MEMBER_')]
            c_order={ x: i for i,x in enumerate(S_order) }
            if opt['LISTS']>1:
                t2=t2[['_PATTERN_','GO','Description','#GeneInGOAndHitList','%InGO','LogP','Log(q-value)']].copy()
            else:
                t2=t2[['GO','Description','#GeneInGOAndHitList','%InGO','LogP','Log(q-value)']].copy()
            t2.rename2({'#GeneInGOAndHitList':'Count', '%InGO':'%', 'LogP':'Log10(P)','Log(q-value)':'Log10(q)'})
            s_note=''
            if s_cat=='TTRUST':
                S=t2['GO']
                t2['GO']=[ '<a href="https://www.grnpedia.org/trrust/result.php?gene={}&species={}" target="_TERM">{}</a>'.format(des.replace('Regulated by: ',''), opt['TARGET_SPECIES'].lower(), S[i]) for i,des in enumerate(t2.Description.tolist()) ]
            elif s_cat=='DisGeNET':
                t2['GO']=t2['GO'].apply(lambda x: '<a href="http://www.disgenet.org/browser/0/0/3/0/diseaseid__{0}-source__ALL/_b./" target="_TERM">{0}</a>'.format(x))
            elif s_cat=='COVID':
                t2['GO']=t2['GO'].apply(lambda x: '<a href="https://metascape.org/COVID#200_{0}" target="_TERM">{0}</a>'.format(x))
            elif s_cat.startswith('L1000'):
                if s_cat=='L1000_Cpd':
                    s_note="Note: some perturbagens (e.g., many Broad compounds BRD-K* are not accessible via ConnectivityMap web site, although data were released into GEO database)."
                t2['GO']=t2['Description'].apply(lambda x: '<a href="https://clue.io/command?q={0}" target="_TERM">{0}</a>'.format(re.sub(r'(^\w+:\s*|,.*)', '', x)))
            elif s_cat=='Transcription_Factor_Targets':
                pass
            if opt['LISTS']>1:
                for idx in t2.index:
                    s_pat=t2.loc[idx, '_PATTERN_'][1:]
                    t2.loc[idx, '_PATTERN_']=" ".join(['<div style="background-color:%s;width:20px;height:20px;float:left;"></div>' % (S_color[i] if ((x in c_order) and (s_pat[c_order[x]]=="1")) else "#FFFFFF") for i,x in enumerate(S_name)])

            s_table="<br>"+s_note+"<br>"+util.html(util.df2sdf(t2, s_format='%.2f'), tag_table={'class':'table'}, colors=None, tag_th={'class':'info'})

            S_fig.append("""
              <div class="panel panel-info">
                <div class="panel-heading">Figure %d. Summary of enrichment analysis in %s%s.</div>
                <div class="panel-body">

                <table>
                    <tr>
                        <td>
                            <img src="%s" style="width:1000px;">
                        </td>
                    </tr>
                    <tr>
                        <td align='center'>
                            %s
                        </td>
                    </tr>
                    <tr>
                        <td>
                            %s
                        </td>
                    </tr>
                </table>

                </div>
              </div>""" % (i_fig,
                            s_cat,
                            s_ref,
                            "@URL_BASE@/"+opt['IMG_QC_HEATMAP'][i],
                            go_heatmap_link,
                            s_table
                          ))

        S_qc.append("""<p>Gene list enrichments are identified in the following ontology categories: %s. %s have been used as the enrichment background. Terms with a p-value &lt; %.2f, a minimum count of %d, and an enrichment factor &gt; %.1f (the enrichment factor is the ratio between the observed counts and the counts expected by chance) are collected and grouped into clusters based on their membership similarities.  The top few enriched clusters (one term per cluster) are shown in the Figure %s.  The algorithm used here is the same as that is used for pathway and process enrichment analysis.<p/>"""
                           % (", ".join(S_go_category),
                              'A user-supplied list of {0} genes'.format(opt['BACKGROUND_GENE_COUNT'])
                                  if opt['BACKGROUND_GENE_COUNT']
                                  else  'All genes in the genome',
                              opt['P_CUTOFF'],
                              opt['MIN_COUNT'],
                              opt['ENRICH_CUTOFF'],
                                str(i_fig_start) if len(S_fig)==1 else "{}-{}".format(i_fig_start, i_fig)))
        S_qc.append("\n".join(S_fig))
        s_qc="\n".join(S_qc)

    s_body="\n".join([s_sum, s_ovlp, s_ann, s_mem, s_enrh, s_ppi, s_gpec, s_qc, s_cytoscape_js_client])

    s_ref=""
    if len(S_ref):
        s_ref="""    <h3>Reference</h3></p>
    <ol style="list-style: decimal inside;">\n"""
        for i,x in enumerate(S_ref):
            s_ref+="      <li>"+x+"</li>\n"
        s_ref+="</ol>\n"
    s_end="</div>\n"
    return "".join([s_title, s_body, s_ref, s_end])

def report(myopt, dbcon=None):
    s=_report(myopt, dbcon)
    URL_BASE=myopt['URL_BASE']
    m=re.search('session_id=([^&]*)&', URL_BASE)
    session_id=""
    if m is not None:
        session_id = m.group(1)
    x=re.compile(r'@OFFLINE_BEGIN@.*?@OFFLINE_END@', re.MULTILINE|re.DOTALL)
    s_online=re.sub('@SESSION_ID@', session_id, re.sub('@ICON_URL@', 'Content/Images', re.sub(r'@ONLINE_(BEGIN|END)@', '', re.sub('@URL_BASE@', URL_BASE, s))))
    s_online=x.sub('', s_online)

    x=re.compile(r'@ONLINE_BEGIN@.*?@ONLINE_END@', re.MULTILINE|re.DOTALL) # DOTALL is import for . to match \n
    s_offline=re.sub(r'@OFFLINE_(BEGIN|END)@', '', re.sub('@ICON_URL@', 'icon', x.sub('', re.sub('@URL_BASE@', '.', s))))
    s_offline="""
<!DOCTYPE html>
<html lang="en">
  <head>
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css">
    <script src="https://ajax.googleapis.com/ajax/libs/jquery/3.2.1/jquery.min.js"></script>
    <script src="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/js/bootstrap.min.js"></script>
  </head>
<body ng-app="app">
"""+s_offline+"\n</body></html>"
    return (s_online, s_offline)

if __name__=="__main__":
    import ms.msobject
    s_dir='/projects/LDWeb/public_html/tmp/output'
    t_lists=ms.msobject.MSObject.load(os.path.join(s_dir, 't_lists.pickle'))
    #pd.DataFrame(data={'Name':['Brass','Karlas','Konig'], 'Total':[231,128,297], 'Unique':[231,128,297]})
    t_go=util.read_csv(os.path.join(s_dir, 't_go.csv'))
    # pickle file seems broken, only has GeneID column?
    #ms.msobject.MSObject.load(os.path.join(s_dir, 't_go.pickle'))

    S_ann=[30, 11, 5, 8]
    opt={
            'URL_BASE':'http://localhost/tmp/output',
            'LISTS':3,
            'ANNOTATION':S_ann,
            'MEMBERSHIP':'invasion',
            'MEMBER_SOURCE':[11,15],
            'ENRICHMENT':t_go,
            'DATAFRAME_LISTS':t_lists,
            'DATABASE_UPDATE':'1/20/2016',
            'IMG_CIRCO_BY_GENE':'Overlap_circos/CircosOverlapByGene.png',
            'IMG_CIRCO_BY_GO':'Overlap_circos/CircosOverlapByGO.png',
            'GO_SOURCE':[11,15,19],
            'P_CUTOFF':0.01,
            'ENRICH_CUTOFF':3,
            'IMG_GO_HEATMAP':'Enrichment_heatmap/HeatmapSelectedGO.png',
            'IMG_GO_NETWORK_CLUSTER':'Enrichment_GO/ColorByCluster.png',
            'IMG_GO_NETWORK_PVALUE':'Enrichment_GO/ColorByPValue.png',
            'IMG_GO_NETWORK_COUNT':'Enrichment_GO/ColorByCounts.png'

        }
    x=db.DB('LDDB')
    s_online, s_offline=report(opt, x.con())
    s0="""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="utf-8" />
    <meta http-equiv="X-UA-Compatible" content="IE=EDGE" />
    <meta name="viewport" content="initial-scale=1, maximum-scale=1, user-scalable=no, width=device-width">

    <meta name="keywords" content="Metascape, free bioinformatics resource, gene annotation, target discovery, gene prioritization, gene enrichment, DAVID replacement, gene identifier conversion">
    <meta http-equiv="Pragma" content="no-cache">
    <meta http-equiv="CACHE-CONTROL" content="NO-CACHE">
    <meta http-equiv="EXPIRES" content="Mon, 22 Jul 2002 11:12:01 GMT">

    <title>Metascape</title>
    <link rel="shortcut icon" href="Content/Layout/MetascapeLogo_24.png">
    <link href="Content/Vender/bootstrap/bootstrap.min.css" rel="stylesheet" />
    <link href="Content/Vender/kendo/kendo.common.min.css" rel="stylesheet" />
    <link href="Content/Vender/kendo/kendo.bootstrap.min.css" rel="stylesheet" />
    <link href="Content/Vender/jquery.qtip.css" rel="stylesheet" />
    <link href="Content/Layout/site.css" rel="stylesheet" />


</head>
<body ng-app="app">
    <div ui-view="top"></div>
    <div ui-view="body">@REPLACE@</div>


    <script src="Content/Vender/jquery/jquery-1.10.2.min.js"></script>
    <script src="http://maxcdn.bootstrapcdn.com/bootstrap/3.3.5/js/bootstrap.min.js"></script>

    <script src="Content/Vender/jquery.qtip.js"></script>
    <script src="Content/Vender/angular.js"></script>
    <script src="Content/Vender/angular-ui-router/angular-ui-router.js"></script>
    <script src="Content/Vender/angular-sanitize.js"></script>
    <script src="Content/Vender/bootstrap/ui-bootstrap-tpls-0.14.2.min.js"></script>
    <script src="Content/Vender/kendo/kendo.all.min.js"></script>
    <script src="Content/Vender/lodash.js"></script>
    <script src="Content/Vender/linq.js"></script>
    <script src="Content/MenuPages/gp_description_stats.js?v=version1027"></script>
    <script src="Content/MenuPages/gp_stats.js?v=version1027"></script>
    <script src="Content/Main/setting.js?v=version1027"></script>

    <script src="Content/Main/app.js?v=version1027"></script>
    <script src="Content/Main/util.js?v=version1027"></script>
    <script src="Content/Main/addinLibrary.js?v=version1027"></script>
    <script src="Content/Layout/header.js?v=version1027"></script>
    <script src="Content/Layout/side_bar.js?v=version1027"></script>
    <script src="Content/Main/main_step1.js?v=version1027"></script>
    <script src="Content/Main/main_step3_controller.js?v=version1027"></script>
    <script src="Content/Main/main_step3_directive.js?v=version1027"></script>
    <script src="Content/Main/main_step3_service.js?v=version1027"></script>
    <script src="Content/IdMapping/idmapping.js?v=version1027"></script>
    <script src="Content/Annotation/annotation.js?v=version1027"></script>
    <script src="Content/Membership/membership.js?v=version1027"></script>
    <script src="Content/Enrichment/enrichment.js?v=version1027"></script>
    <script src="Content/Report/report.js?v=version1027"></script>
    <script src="Content/Report/report_final.js?v=version1027"></script>
    <script src="Content/MenuPages/manual.js?v=version1027"></script>
    <script src="Content/MenuPages/data_source_logging.js?v=version1027"></script>
    <script src="Content/Tools/table_join.js?v=version1027"></script>
    <script src="Content/Layout/news.js?v=version1027"></script>
    <script src="Content/Enrichment/background_dialog.js?v=version1027"></script>
</body>
</html>"""
    util.save_string('/projects/LDWeb/public_html/tmp/report.html', s0.replace('@REPLACE@', s_online).replace('"Content/','"http://localhost/gp/Content/'))
