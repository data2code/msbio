#
# The MIT License (MIT)
#
# Copyright (c) 2015 Matthew Antalek Jr
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

### ZHOU
from __future__ import absolute_import
from __future__ import print_function
from __future__ import division
import matplotlib as mpl
#see https://stackoverflow.com/questions/5956182/cannot-edit-text-in-chart-exported-by-matplotlib-and-opened-in-illustrator
mpl.rcParams['pdf.fonttype'] = 42
from six.moves import range
import util
mpl.use('Agg')
import matplotlib.pyplot as pylab
import scipy.cluster.hierarchy as sch
import numpy as np
import math

class DendroHeatMap(object):
    """
    Class for quickly and easily plotting heatmaps with dendrograms on the side, as seen in
    http://code.activestate.com/recipes/578175-hierarchical-clustering-heatmap-python/
    """

    def __init__(self, heat_map_data=None, left_dendrogram=None,top_dendrogram=None,
                 window_height=None, window_width = None, color_bar_width = 0.015,
                 left_dendro_x=0.05,left_dendro_y=0.22,left_dendro_width=0.2,left_dendro_height=0.6, left_dendro_x_distance_to_row_cb=0.004, left_dendro_y_distance_to_col_cb=0.004,
                 top_dendro_x=0.273, top_dendro_y=0.843, top_dendro_width=0.5, top_dendro_height=0.117,
                 row_cb_x=0.254,row_cb_y=0.22,row_cb_width=0.015,row_cb_height=0.6,row_cb_on=True,
                 col_cb_x = 0.273, col_cb_y=0.824, col_cb_width=0.5, col_cb_height=0.015, col_cb_on=True,
                 heat_x=0.273, heat_y=0.22,heat_width=0.5,heat_height=0.6,
                 #color_legend_x=0.07,color_legend_y=0.88, color_legend_width=0.2,color_legend_height=0.05, color_legend_ticks=7,
                 color_legend_x=0.80,color_legend_y=0.88, color_legend_width=0.191,color_legend_height=0.05, color_legend_ticks=7,
                 row_labels=None, max_row_labels=100, row_labels_size=8,
                 col_labels=None, max_col_labels=100, col_labels_size=8,
                 l_normalize_for_color=True, ### ZHOU
                 l_legend_pvalue=False, ### ZHOU
                 verbose=False):

        import warnings
        warnings.simplefilter("ignore")

        if left_dendrogram is not None and len(left_dendrogram)==0: left_dendrogram=None
        if top_dendrogram is not None and len(top_dendrogram)==0: top_dendrogram=None

        self.figure = None
        self.verbose= verbose

        ## ZHOU
        self.l_normalize_for_color=l_normalize_for_color
        self.l_legend_pvalue=l_legend_pvalue

        #set the default behaviors
        if window_height is None or window_width is None:
            n,m=heat_map_data.shape
            tile_sz_x=max(0.30, 4/m)
            tile_sz_y=0.30
            heat_height=(tile_sz_y*n) # add extra for axis labels
            heat_width=(tile_sz_x*m) # add extra for axis labels
            # for 20x4 heatmap_data, the heatmap size is 6x4
            # dendrogram should have size of 2, but we need the relative size
            window_width=heat_width+4 # total size: heatmap + dendrogram + label
            window_height=heat_height+4
            unit_x=1.0/window_width
            unit_y=1.0/window_height
            left_dendro_x=0.05*unit_x
            left_dendro_y=2*unit_y
            left_dendro_width=1.9*unit_x
            left_dendro_height=heat_height*unit_y
            left_dendro_x_distance_to_row_cb=0.05*unit_x
            left_dendro_y_distance_to_col_cb=0.05*unit_y
            top_dendro_x=2*unit_x
            top_dendro_y=(2+heat_height+0.05)*unit_y
            top_dendro_width=heat_width*unit_x
            top_dendro_height=1.9*unit_y
            row_cb_x=2.54*unit_x
            row_cb_y=2.2*unit_y
            row_cb_width=0.15*unit_x
            row_cb_height=heat_height*unit_y
            row_cb_on=True,
            col_cb_x = 2.73*unit_x
            col_cb_y=(2+heat_height)*unit_y
            col_cb_width=0.5
            col_cb_height=0.15*unit_y
            col_cb_on=True
            color_legend_x=(2+heat_width+0.05)*unit_x
            color_legend_y=(2+heat_height+1.0)*unit_y
            color_legend_width=1.9*unit_x
            color_legend_height=0.5*unit_y
            color_legend_ticks=7
            heat_x=2*unit_x
            heat_y=2*unit_y
            heat_width*=unit_x
            heat_height*=unit_y
            #print(left_dendro_x,left_dendro_y,left_dendro_width,left_dendro_height)
            #print(top_dendro_x,top_dendro_y,top_dendro_width,top_dendro_height)
            #print(heat_x,heat_y,heat_width,heat_height)

        # print 'should be moving into setter land....'
        self.heat_map_data = heat_map_data
        self.top_dendrogram = top_dendrogram
        self.left_dendrogram = left_dendrogram

        self.window_height=window_height
        self.window_width=window_width
        self.color_bar_width=color_bar_width

        self.left_dendro_x=left_dendro_x
        self.left_dendro_y=left_dendro_y
        self.left_dendro_width=left_dendro_width
        self.left_dendro_height=left_dendro_height
        self.left_dendro_x_distance_to_row_cb=left_dendro_x_distance_to_row_cb
        self.left_dendro_y_distance_to_col_cb=left_dendro_y_distance_to_col_cb


        self.top_dendro_x=top_dendro_x
        self.top_dendro_y=top_dendro_y
        self.top_dendro_width = top_dendro_width
        self.top_dendro_height=top_dendro_height

        self.cluster_cb_colors = mpl.colors.ListedColormap(['r', 'g', 'b', 'y', 'w', 'k', 'm'])

        self.row_cb_x=row_cb_x
        self.row_cb_y = row_cb_y
        self.row_cb_width=row_cb_width
        self.row_cb_height=row_cb_height
        self.row_cb_on=row_cb_on

        self.col_cb_x=col_cb_x
        self.col_cb_y=col_cb_y
        self.col_cb_width=col_cb_width
        self.col_cb_height=col_cb_height
        self.col_cb_on=col_cb_on

        self.heat_x=heat_x
        self.heat_y=heat_y
        self.heat_width=heat_width
        self.heat_height=heat_height

        self.color_legend_x=color_legend_x
        self.color_legend_y=color_legend_y
        self.color_legend_width=color_legend_width
        self.color_legend_height=color_legend_height
        self.color_legend_ticks = color_legend_ticks

        self.row_labels=row_labels
        self.row_labels_size=row_labels_size
        self.max_row_labels=max_row_labels

        self.col_labels=col_labels
        self.col_labels_size=col_labels_size
        self.max_col_labels=max_col_labels

        self.redBlackBlue=self.__RedBlackBlue()
        self.redBlackSkyBlue=self.__RedBlackSkyBlue()
        self.redBlackGreen=self.__RedBlackGreen()
        self.yellowBlackBlue=self.__YellowBlackBlue()
        self.colormap=self.redBlackGreen

        self.left_dendro_title = ''
        self.top_dendro_title = ''
        self.title = ''
        self.color_legend_title = ''
        self.plotRendered = False
        self.exportDPI = 100

    def render_plot(self,showFrames=False):
        self.resetPlot()

        if(self.verbose):
            print('Rendering plot...')

        self.figure = pylab.figure(figsize=[self.window_width, self.window_height])

        #plot the top dendrogram
        if(not self.top_dendrogram is None):
            self.top_dendro_axes = self.figure.add_axes([self.top_dendro_x, self.top_dendro_y, self.top_dendro_width, self.top_dendro_height], frame_on=showFrames)
            ### ZHOU set link color to black instead of default blue
            self.top_dendro_plot = sch.dendrogram(self.top_dendrogram, link_color_func=lambda k: 'black')
            self.top_dendro_axes.set_xticks([])
            self.top_dendro_axes.set_yticks([])
            self.top_dendro_axes.set_title(self.top_dendro_title)
            ### ZHOU
            self.__heat_map_data=self.__heat_map_data[:, self.top_dendro_plot['leaves']]
            self.col_labels=[ self.col_labels[x] for x in self.top_dendro_plot['leaves']]

        #plot the left dendrogram
        if(not self.left_dendrogram is None):
            self.left_dendro_axes = self.figure.add_axes([self.left_dendro_x, self.left_dendro_y, self.left_dendro_width, self.left_dendro_height], frame_on=showFrames)
            ### ZHOU
            self.left_dendro_plot = sch.dendrogram(self.left_dendrogram,orientation='left', link_color_func=lambda k: 'black')
            self.left_dendro_axes.set_xticks([])
            self.left_dendro_axes.set_yticks([])
            self.left_dendro_axes.set_title(self.left_dendro_title,rotation='vertical')
            ### ZHOU
            self.__heat_map_data=self.__heat_map_data[self.left_dendro_plot['leaves'], :]
            self.row_labels=[ self.row_labels[x] for x in self.left_dendro_plot['leaves']]

        #plot the heat map
        if(not self.heat_map_data is None):
            self.heat_map_axes = self.figure.add_axes([self.heat_x, self.heat_y, self.heat_width, self.heat_height], frame_on=showFrames)
            if self.cmap_norm is None:
                self.heat_map_plot = self.heat_map_axes.matshow(self.heat_map_data, aspect='auto', origin='lower', cmap=self.colormap, vmin=0, vmax=1)
            else:
                self.heat_map_plot = self.heat_map_axes.matshow(self.heat_map_data, aspect='auto', origin='lower', cmap=self.colormap, norm=self.cmap_norm)
            self.heat_map_axes.set_xticks([])
            self.heat_map_axes.set_yticks([])
            self.heat_map_rows = self.heat_map_data.shape[0]
            self.heat_map_cols = self.heat_map_data.shape[1]

            w,h= self.get_ax_size(self.heat_map_axes)
            #add the from the labels to the figure
            # print len(self.row_labels)
            if self.row_labels_size==0:
                # with smallest font size 8.0, it can fit at most h/8.1 rows of labels
                i_step=max(int(math.ceil(self.heat_map_rows/(h/8.1))), 1)
                scale=0.7 #if util.is_python3() else 1.0  #maybe it has something to do with the installation
                self.row_labels_size=int(math.floor(max(8, h/(self.heat_map_rows//i_step)*scale)))
            else:
                i_step=1
            #print(h, i_step, self.heat_map_rows, self.row_labels_size)
            #print type(self.row_labels), self.row_labels[:5], len(self.row_labels)
            from matplotlib.font_manager import FontProperties
            import os
            s_file=os.path.join(os.path.dirname(__file__), "ms", "report", "arial.ttf")
            prop=FontProperties(fname=s_file)

            #if(self.row_labels is not None and (len(self.row_labels)*1.0/i_step) < self.max_row_labels):
            for i in range(0, self.heat_map_rows, i_step):
                    #if(self.row_labels):
                    #if(len(self.row_labels) < self.max_row_labels):
                    #print(i-0.25, self.row_labels_size)
                self.heat_map_axes.text(self.heat_map_cols-0.5, i, ' '+self.row_labels[i], fontproperties=prop, size=self.row_labels_size, verticalalignment='center', horizontalalignment='left')
            if self.col_labels_size==0:
                # with smallest font size 8.0, it can fit at most h/8.0 rows of labels
                i_step=int(math.ceil(self.heat_map_rows/(w/8.0)))
                self.col_labels_size=int(math.floor(max(8, w/(self.heat_map_cols//i_step))))
            else:
                i_step=1

            #if (self.col_labels is not None and (len(self.col_labels)*1.0//i_step) < self.max_col_labels):
            for i in range(0, self.heat_map_cols, i_step):
                    #if(self.col_labels):
                    #if(len(self.col_labels) < self.max_col_labels):
                self.heat_map_axes.text(i, -0.5, ' '+self.col_labels[i], size=self.col_labels_size, rotation=270,verticalalignment='top', horizontalalignment='center')
                        #self.heat_map_axes.text(i+0.05, self.heat_map_rows-self.heat_map_rows-0.5, ' '+self.col_labels[i], size=self.col_labels_size, rotation=270,verticalalignment='top')

       # #plot the column colorbar
       # if(not self.top_dendrogram is None):
       #     self.col_cb_axes = self.figure.add_axes([self.col_cb_x, self.col_cb_y, self.col_cb_width, self.col_cb_height], frame_on=True)
       #     # print self.top_colorbar_labels.shape
       #     # print 'Col cb'
       #     # print [self.col_cb_x, self.col_cb_y, self.col_cb_width, self.col_cb_height]
       #     self.col_cb_plot = self.col_cb_axes.matshow(self.top_colorbar_labels,aspect='auto',origin='lower',cmap=self.cluster_cb_colors)
       #     self.col_cb_axes.set_xticks([])
       #     self.col_cb_axes.set_yticks([])

       # #plot the row colorbar
       # if(not self.left_dendrogram is None):
       #     self.row_cb_axes = self.figure.add_axes([self.row_cb_x, self.row_cb_y, self.row_cb_width, self.row_cb_height], frame_on=True)
       #     # print self.left_colorbar_labels.shape
       #     # print 'Row cb'
       #     # print [self.row_cb_x, self.row_cb_y, self.row_cb_width, self.row_cb_height]
       #     self.row_cb_plot = self.row_cb_axes.matshow(self.left_colorbar_labels, aspect='auto',origin='lower',cmap=self.cluster_cb_colors)
       #     self.row_cb_axes.set_xticks([])
       #     self.row_cb_axes.set_yticks([])

        #plot the color legend
        if(not self.heat_map_data is None):
            self.color_legend_axes = self.figure.add_axes([self.color_legend_x, self.color_legend_y, self.color_legend_width, self.color_legend_height], frame_on=showFrames)
            self.color_legend_plot = mpl.colorbar.ColorbarBase(self.color_legend_axes, cmap=self.colormap, norm=self.cmap_norm,orientation='horizontal')
            ### ZHOU
            if self.l_legend_pvalue:
                #def pval(x, pos):
                #    print ">>>>>>>>>>>", x, pos
                #    return "%d" % x*20
                self.color_legend_plot.set_ticks([0,0.1,0.15,0.2,0.3,0.5,1.0])
                self.color_legend_plot.set_ticklabels(['0', '2', '3', '4', '6', '10', '20'])
                #set_major_formatter(mpl.ticker.FuncFormatter(pval))
                self.color_legend_title='-log10(P)' # TeX is not installed, r'$-log_{10}P$'
            else:
                tl=mpl.ticker.MaxNLocator(nbins=self.color_legend_ticks)
                self.color_legend_plot.locator = tl
            self.color_legend_plot.update_ticks()
            self.color_legend_axes.set_title(self.color_legend_title)
            self.heat_map_axes.format_coord = self.__formatCoords

        self.figure.suptitle(self.title)

        self.plotRendered = True

        if(self.verbose):
            print('Plot rendered...')


    def show(self):
        self.resetPlot()
        self.render_plot()
        pylab.show()


    def export(self,filename, l_pdf=False):
        self.resetPlot()
        if('.' not in filename):
            filename += '.png'
        else:
            filename = filename[:-4] + '.png'

        if(self.verbose):
            print('Saving plot to: ', filename)
        self.render_plot()
        #if filename.lower().endswith('.pdf'):
        #    #does not work, not knowing why
        #    pylab.savefig(filename)
        #else:
        pylab.savefig(filename, bbox_inches='tight')
        if l_pdf:
            pylab.savefig(filename.replace('.png', '.pdf'), bbox_inches='tight')
        pylab.close()

    @property
    def heat_map_data(self):
        return self.__heat_map_data

    @heat_map_data.setter
    def heat_map_data(self, heat_map_data):
        # print 'In the setter...'
        self.__heat_map_data=heat_map_data
        self.resetPlot()
        # print type(heat_map_data)
        if self.l_normalize_for_color:
            if((isinstance(heat_map_data,np.ndarray)) | (isinstance(heat_map_data,np.matrix))):
                hm_min = heat_map_data.min()
                hm_max = heat_map_data.max()
                self.cmap_norm = mpl.colors.Normalize(hm_min,hm_max)
            else:
                raise TypeError('Data for the heatmap must be a numpy.ndarray or numpy.matrix object!')
        else:
            self.cmap_norm=None


    def resetPlot(self):
        self.plotRendered = False
        if(self.figure):
            pylab.close(self.figure)
            self.figure = None
        else:
            self.figure = None

    @property
    def figure(self):
        return self.__figure

    @figure.setter
    def figure(self,figure):
        self.__figure = figure
        if((not isinstance(figure, pylab.Figure)) & (isinstance(figure,object))):
            #this force's the figure to either be "None" type or a pylab.Figure object
            self.__figure = None


    @property
    def row_labels(self):
        return self.__row_labels

    @row_labels.setter
    def row_labels(self, row_labels):
        if(not isinstance(self.heat_map_data,np.ndarray) or not isinstance(self.heat_map_data, np.matrix)):
            if(self.verbose):
                print("""Warning: data for heat map not yet specified, be sure that the number of elements in row_labels
                is equal to the number of rows in heat_map_data.
                """)
            self.__row_labels = row_labels
        else:
            if(len(row_labels) != self.heat_map_data.shape[0]):
                print("""Invalid entry for row_labels. Please be sure that the number of elements in row_labels is equal
                to the number of rows in heat_map_data.""")
                self.__row_labels = None
            else:
                self.__row_labels = row_labels


    @property
    def col_labels(self):
        return self.__col_labels

    @col_labels.setter
    def col_labels(self, col_labels):
        if(not isinstance(self.heat_map_data,np.ndarray) or not isinstance(self.heat_map_data, np.matrix)):
            if(self.verbose):
                print("""Warning: data for heat map not yet specified, be sure that the number of elements in col_labels
                is equal to the number of columns in heat_map_data.
                """)
            self.__col_labels = col_labels
        else:
            if(len(col_labels) != self.heat_map_data.shape[0]):
                print("""Invalid entry for col_labels. Please be sure that the number of elements in col_labels is equal
                to the number of columns in heat_map_data.""")
                self.__col_labels = None
            else:
                self.__col_labels = col_labels


    @property
    def colormap(self):
        return self.__colormap

    @colormap.setter
    def colormap(self, colormap):
        self.__colormap = colormap
        self.resetPlot()


    @property
    def top_dendrogram(self):
        return self.__top_dendrogram

    @top_dendrogram.setter
    def top_dendrogram(self,top_dendrogram):
        if(isinstance(top_dendrogram,np.ndarray)):
            self.__top_dendrogram = top_dendrogram
            self.resetPlot()
            self.top_colorbar_labels = np.array(sch.fcluster(top_dendrogram,0.7*max(top_dendrogram[:,2]),'distance'),dtype=int)
            self.top_colorbar_labels.shape = (1,len(self.top_colorbar_labels))
            temp_dendro = sch.dendrogram(top_dendrogram,no_plot=True)
            self.top_colorbar_labels = self.top_colorbar_labels[:,temp_dendro['leaves']]
        elif top_dendrogram is None:
            self.__top_dendrogram = top_dendrogram
            self.resetPlot()
        else:
            raise TypeError('Dendrograms must be a n-1 x 4 numpy.ndarray as per the scipy.cluster.hierarchy implementation!')

    @property
    def left_dendrogram(self):
        return self.__left_dendrogram

    @left_dendrogram.setter
    def left_dendrogram(self,left_dendrogram):

        if isinstance(left_dendrogram,np.ndarray):
            self.__left_dendrogram = left_dendrogram
            self.resetPlot()
            self.left_colorbar_labels = np.array(sch.fcluster(left_dendrogram,0.7 * max(left_dendrogram[:,2]),'distance'), dtype=int)
            self.left_colorbar_labels.shape = (len(self.left_colorbar_labels),1)
            temp_dendro = sch.dendrogram(left_dendrogram,no_plot=True)
            self.left_colorbar_labels = self.left_colorbar_labels[temp_dendro['leaves'],:]
        elif left_dendrogram is None:
            self.__left_dendrogram = left_dendrogram
            self.resetPlot()

        else:
            raise TypeError('Dendrograms must be a n-1 x 4 numpy.ndarray as per the scipy.cluster.hierarchy implementation!')


    def __RedBlackSkyBlue(self):
        cdict = {'red':   ((0.0, 0.0, 0.0),
                           (0.5, 0.0, 0.1),
                           (1.0, 1.0, 1.0)),

                 'green': ((0.0, 0.0, 0.9),
                           (0.5, 0.1, 0.0),
                           (1.0, 0.0, 0.0)),

                 'blue':  ((0.0, 0.0, 1.0),
                           (0.5, 0.1, 0.0),
                           (1.0, 0.0, 0.0))
                }

        my_cmap = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,256)
        return my_cmap

    def __RedBlackBlue(self):
        cdict = {'red':   ((0.0, 0.0, 0.0),
                           (0.5, 0.0, 0.1),
                           (1.0, 1.0, 1.0)),

                 'green': ((0.0, 0.0, 0.0),
                           (1.0, 0.0, 0.0)),

                 'blue':  ((0.0, 0.0, 1.0),
                           (0.5, 0.1, 0.0),
                           (1.0, 0.0, 0.0))
                }

        my_cmap = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,256)
        return my_cmap


    def __RedBlackGreen(self):
        cdict = {'red':   ((0.0, 0.0, 0.0),
                           (0.5, 0.0, 0.1),
                           (1.0, 1.0, 1.0)),

                 'blue': ((0.0, 0.0, 0.0),
                           (1.0, 0.0, 0.0)),

                 'green':  ((0.0, 0.0, 1.0),
                           (0.5, 0.1, 0.0),
                           (1.0, 0.0, 0.0))
                }

        my_cmap = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,256)
        return my_cmap

    def __YellowBlackBlue(self):
        cdict = {'red':   ((0.0, 0.0, 0.0),
                           (0.5, 0.0, 0.1),
                           (1.0, 1.0, 1.0)),

                 'green': ((0.0, 0.0, 0.8),
                           (0.5, 0.1, 0.0),
                           (1.0, 1.0, 1.0)),

                 'blue':  ((0.0, 0.0, 1.0),
                           (0.5, 0.1, 0.0),
                           (1.0, 0.0, 0.0))
                }
        ### yellow is created by adding y = 1 to RedBlackSkyBlue green last tuple
        ### modulate between blue and cyan using the last y var in the first green tuple
        my_cmap = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,256)
        return my_cmap

    ### ZHOU
    def BlueGrayRed(self):
        cdict = {'red':   ((0.0, 0.0, 0.227451),
                           (0.5, 0.847059, 0.847059),
                           (1.0, 0.847059, 0.847059)),

                 'green': ((0.0, 0.0, 0.423529),
                           (0.5, 0.847059, 0.847059),
                           (1.0, 0.094118, 0.094118)),

                 'blue':  ((0.0, 0.0, 0.603922),
                           (0.5, 0.847059, 0.847059),
                           (1.0, 0.109804, 0.109804))
                }
        ### yellow is created by adding y = 1 to RedBlackSkyBlue green last tuple
        ### modulate between blue and cyan using the last y var in the first green tuple
        my_cmap = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,256)
        return my_cmap

    ### ZHOU
    def color_brewer(self, brewer_name='Blues', map_type='sequential', number=6, reverse=False):
        import  brewer2mpl
        C=brewer2mpl.get_map(brewer_name, map_type, number, reverse=reverse).colors[:]
        def color2array(I_color):
            n=len(I_color)
            data=[]
            for i in range(n):
                X= [min(max(i/(n-1.0), 0.0), 1.0), 0.0, 0.0]
                X[1]=X[2]=I_color[i]/255.
                data.append(X)
            return data

        cdict={
            'red': color2array([x[0] for x in C]),
            'green': color2array([x[1] for x in C]),
            'blue': color2array([x[2] for x in C]),
        }

        ### yellow is created by adding y = 1 to RedBlackSkyBlue green last tuple
        ### modulate between blue and cyan using the last y var in the first green tuple
        my_cmap = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,256)
        return my_cmap

    @staticmethod
    def color_by_pvalue():
        import  brewer2mpl
        P=[0,0.1,0.15,0.2,0.3,0.5,1.0]
        C=brewer2mpl.get_map("YlOrBr", 'sequential', len(P)-1, reverse=False).colors[:]
        #P.insert(0, 0.0)
        C.insert(0, (217,217,217))
        #C.append(C[-1])

        def color2array(Value, I_color):
            n=len(Value)
            data=[]
            for i in range(n):
                X= [Value[i], 0.0, 0.0]
                X[1]=X[2]=I_color[i]/255.
                data.append(X)
            return data

        cdict={
            'red': color2array(P, [x[0] for x in C]),
            'green': color2array(P, [x[1] for x in C]),
            'blue': color2array(P, [x[2] for x in C]),
        }

        ### yellow is created by adding y = 1 to RedBlackSkyBlue green last tuple
        ### modulate between blue and cyan using the last y var in the first green tuple
        my_cmap = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,256)
        return my_cmap

    def get_ax_size(self, ax):
        """http://stackoverflow.com/questions/19306510/determine-matplotlib-axis-size-in-pixels"""
        bbox = ax.get_window_extent().transformed(self.figure.dpi_scale_trans.inverted())
        width, height = bbox.width, bbox.height
        width *= self.figure.dpi
        height *= self.figure.dpi
        return width, height

    def __formatCoords(self, x,y):
        col = int(x+0.5)
        row = int(y+0.5)
        if col>=0 and col<self.heat_map_cols and row>=0 and row<self.heat_map_rows:
            z = self.heat_map_data[row,col]
            return 'x=%1.4f, y=%1.4f, z=%1.4f'%(x, y, z)
        else:
            return 'x=%1.4f, y=%1.4f'%(x, y)

if __name__=="__main__":
    n_rows=20
    n_cols=6
    data_dist=np.random.rand(n_rows, n_cols)
    S_go=["Description text: whatever ... {}".format(i+1) for i in range(n_rows)]
    S_label=["Dataset #{}".format(i+1) for i in range(n_cols)]
    cm=DendroHeatMap.color_by_pvalue()
    import fastcluster
    Zr=fastcluster.linkage(data_dist, method='average', metric='euclidean', preserve_input=True)
    Zc=fastcluster.linkage(data_dist.T, method='average', metric='euclidean', preserve_input=True)

    import cluster
    den_r=cluster.FastCluster.linkage2order(Zr)
    fc=cluster.FastCluster(data_dist, S_col=S_label, S_row=S_go, S_description=None, Zr=Zr, Zc=Zc)
    fc.plot("heatmap.png", colormap=cm, l_normalize_for_color=False, l_legend_pvalue=True, l_pdf=False)

