import plotly.graph_objs as go
from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot

import numpy as np
import pandas as pd
import os


from matplotlib import pyplot as plt

from matplotlib.collections import BrokenBarHCollection

import collections

def recursively_default_dict():
        return collections.defaultdict(recursively_default_dict)


def read_focus(index_file):
    indxs = []
    
    Input = open(index_file,'r')
    for line in Input:
        line = line.split()
        indxs.append(line[0])
    
    Input.close()
    
    return indxs



def display_ideo(select_column,CHR,ID= 'none'):
    house= 'Ideos/'
    
    if select_column== 0:
        select_column= 'all'
    else: 
        select_column= str(select_column - 1)

    filename= house + 'Ideo_id.' + ID + '_label.' + select_column +'_CHR' + str(CHR).zfill(2) + '.png'
    img_width = 1900
    img_height = 1300
    scale_factor = .5

    layout = go.Layout(
        xaxis = go.layout.XAxis(
            visible = False,
            range = [0, img_width*scale_factor]),
        yaxis = go.layout.YAxis(
            visible=False,
            range = [0, img_height*scale_factor]),
            # the scaleanchor attribute ensures that the aspect ratio stays constant
            #scaleanchor = 'x'),
        width = img_width*scale_factor,
        height = img_height*scale_factor,
        margin = {'l': 0, 'r': 0, 't': 0, 'b': 0},
    images = [go.layout.Image(
        x=0,
        sizex=img_width*scale_factor,
        y=img_height*scale_factor,
        sizey=img_height*scale_factor,
        xref="x",
        yref="y",
        opacity=1.0,
        layer="below",
        sizing="stretch",
        source=filename)]
    )
    
    # we add a scatter trace with data points in opposite corners to give the Autoscale feature a reference point
    fig = go.Figure(data=[{
        'x': [0, img_width*scale_factor], 
        'y': [0, img_height*scale_factor], 
        'mode': 'markers',
        'marker': {'opacity': 0}}],layout = layout)
    iplot(fig)


def target_ideogram(gp,Coordinates,Focus,Colors= {}, Chr= 1,background= True,height_chrom= .5,height= 10,width= 5):
    
    Chromo_coord= Coordinates[Coordinates.chrom == Chr]
    
    #Where= Chromo_coord[Chromo_coord.label == gp-1]
    
    ##
    trace0 = go.Scatter(
        x=[0, max(Chromo_coord.end)],
        y=[0, height_chrom * len(Focus)],
        mode='text',
    )

    data = [trace0]
    layout = {
        'xaxis': {
            'showgrid': False
        },
        'yaxis': {
            'showgrid': False
        },
        'shapes': [{
            'type': 'rect',
            'y0': 0,
            'x0': 0,
            'y1': height_chrom * len(Focus) + 1,
            'x1': max(Chromo_coord.end) + 1,
            'line': {
                'color': 'rgba(1, 1, 1, 1)',
                'width': 2,
            }
        }]
    }

    for row in range(Chromo_coord.shape[0]):
            # filled Rectangle
            CHR= Chromo_coord.iloc[row,:].chrom
            start= Chromo_coord.iloc[row,:].start
            end= Chromo_coord.iloc[row,:].end
            trigger= [float(x) for x in Chromo_coord.iloc[row,:].members.split('.')]
            label= Chromo_coord.iloc[row,:].label + 1
            color= Colors[-1]
            
            if label in gp:
                color= Colors[label]
            else:
                if not background:
                    continue
            
            for v in trigger:
                v= len(Focus) - v - 1
                rekt= {
                    'type': 'rect',
                    'y0': v * height_chrom,
                    'x0': start,
                    'y1': (v + 1)*height_chrom,
                    'x1': end,
                    'line': {
                        'color': color,
                        'width': .2,
                    },
                    'fillcolor': color,
                }

                layout['shapes'].append(rekt)



    fig = {
        'data': data,
        'layout': layout,
    }
    
    iplot(fig)



def chromosome_collections(df, y_positions, height,  **kwargs):
    """
    Yields BrokenBarHCollection of features that can be added to an Axes
    object.
    Parameters
    ----------
    df : pandas.DataFrame
        Must at least have columns ['chrom', 'start', 'end', 'color']. If no
        column 'width', it will be calculated from start/end.
    y_positions : dict
        Keys are chromosomes, values are y-value at which to anchor the
        BrokenBarHCollection
    height : float
        Height of each BrokenBarHCollection
    Additional kwargs are passed to BrokenBarHCollection
    """
    del_width = False
    if 'width' not in df.columns:
        del_width = True
        df['width'] = df['end'] - df['start']
    for chrom, group in df.groupby('chrom'):
        
        yrange = (y_positions[chrom], height)
        xranges = group[['start', 'width']].values
        yield BrokenBarHCollection(
            xranges, yrange, facecolors=group['colors'], **kwargs)
    if del_width:
        del df['width']



def mpl_target_ideo(gp,
                    Coordinates,
                    Focus,
                    order= [],
                    background= True,
                    Chr= 1,
                    ideo_height= 1,
                    ideo_spacing= 0,
                    xticks= 1e6,
                    fig_save= True,
                    fig_id= 'ideo_target',
                    Colors= {},
                    height= 5,
                    width= 15,
                    Home= ''):
    
    Chromo_coord= Coordinates[Coordinates.chrom == Chr]
    
    if not order:
        order= list(range(len(Focus)))
    
    Ideo = []
    chromosome_list= ['chr'+str(Chr)+ '_' + Focus[x] for x in order]
    
    for row in range(Chromo_coord.shape[0]):
            # filled Rectangle
            CHR= Chromo_coord.iloc[row,:].chrom
            start= Chromo_coord.iloc[row,:].start
            end= Chromo_coord.iloc[row,:].end
            trigger= [int(x) for x in Chromo_coord.iloc[row,:].members.split('.')]
            
            label= Chromo_coord.iloc[row,:].label + 1
            color= Colors[-1]
            
            if label in gp:
                color= Colors[label]
            else:
                if not background:
                    continue
            
            
            color= [round(y / float(255),1) for y in color]
            for v in trigger:
                
                Subject = Focus[v]
                leine= ['chr'+str(CHR)+ '_' + Subject,start,end,color]
                Ideo.append(leine)


    print(len(Ideo))
    if Ideo:

        # Height of each ideogram
        chrom_height = ideo_height

        # Spacing between consecutive ideograms
        chrom_spacing = ideo_spacing

        # Height of the gene track. Should be smaller than `chrom_spacing` in order to
        # fit correctly
        gene_height = 0.0

        # Padding between the top of a gene track and its corresponding ideogram
        gene_padding = 0.0


        # Keep track of the y positions for ideograms and genes for each chromosome,
        # and the center of each ideogram (which is where we'll put the ytick labels)
        ybase = 0
        chrom_ybase = {}
        gene_ybase = {}
        chrom_centers = {}

        # Iterate in reverse so that items in the beginning of `chromosome_list` will
        # appear at the top of the plot
        for chrom in chromosome_list[::-1]:
            chrom_ybase[chrom] = ybase
            chrom_centers[chrom] = ybase + chrom_height / 2.
            gene_ybase[chrom] = ybase - gene_height - gene_padding
            ybase += chrom_height + chrom_spacing

        # Read in ideogram.txt, downloaded from UCSC Table Browser
        ideo = pd.DataFrame(Ideo,columns = ['chrom', 'start', 'end', 'colors'])

        # Filter out chromosomes not in our list
        ideo = ideo[ideo.chrom.apply(lambda x: x in chromosome_list)]

        # Add a new column for width
        ideo['width'] = ideo.end - ideo.start

        # Width, height (in inches)
        figsize = (width, height)

        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111)

        # Now all we have to do is call our function for the ideogram data...
        print("adding ideograms...")
        for collection in chromosome_collections(ideo, chrom_ybase, chrom_height, edgecolors=None, linewidths= 0):
            ax.add_collection(collection)


        # Axes tweaking

        ax.set_xlim(min(Chromo_coord.start),max(Chromo_coord.end))
        ax.set_xticks(list(range(min(Chromo_coord.start),max(Chromo_coord.end),int(xticks))))
        plt.xticks(fontsize = 10,rotation = 90)
        ax.tick_params(axis = 'x',pad = 10)

        ax.tick_params(axis='y', which='major', pad=30)
        ax.set_yticks([chrom_centers[i] for i in chromosome_list])
        ax.set_yticklabels(chromosome_list, fontsize = 5)
        
        if fig_save:
            if Home:
                Home= Home + '/'
                filename= Home+ 'Ideo_id.' + fig_id + '_label.' + '.'.join([str(x) for x in gp]) +'_CHR' + str(Chr).zfill(2)+'.png'
                os.makedirs(os.path.dirname(filename), exist_ok=True)
            else:
                filename= Home+ 'Ideo_id.' + fig_id + '_label.' + '.'.join([str(x) for x in gp]) +'_CHR' + str(Chr).zfill(2)+'.png'
                plt.savefig(filename,bbox_inches = 'tight')
        
        plt.show()
        plt.clf()
        
        
