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


def paint_bars(selected,cluster_pca= {}):
    clusters= cluster_pca
    if selected == 0:
        color_range = ['rgb(158,202,225)'] * len(list(set(clusters.iloc[:,0])))
    else:
        color_range = ['rgb(158,202,225)'] * len(list(set(clusters.iloc[:,0])))
        color_range[selected - 1] = 'rgb(204,0,0)'
    
    whom= sorted(list(set(clusters.iloc[:,0])))
    nb= [round(len([x for x in clusters.iloc[:,0] if x == y]) / float(len(clusters)),3) for y in whom]
    nc= [str(x + 1) for x in whom]
    trace = [go.Bar(
    x= nc,
    y= nb,
    text= nb,
    marker=dict(
        color= color_range,
        line=dict(
            color='rgb(8,48,107)',
            width=1.5),
    ),
    opacity= .6
    )]
    layout= go.Layout(
    title= 'cluster proportions'
    )
    fig= go.Figure(data=trace,layout=layout)
    iplot(fig)


def plot_clusters(selected,pc1,pc2,cluster_pca= {},pca_var= [],):
    scheme = [x for x in list(set(cluster_pca.label))]
    if selected == 0:
        opac_range = [.9]*len(list(set(cluster_pca.label)))
    else:
        opac_range = [.3]*len(list(set(cluster_pca.label)))
        opac_range[selected-1] = 1
    fig_data= [go.Scatter(
        x = cluster_pca.iloc[[x for x in range(len(cluster_pca)) if cluster_pca.iloc[x,0] == i],pc1],
        y = cluster_pca.iloc[[x for x in range(len(cluster_pca)) if cluster_pca.iloc[x,0] == i],pc2],
        mode= "markers",
        marker= {
            'line': {'width': 0},
            'size': 5,
            'symbol': 'circle',
            'opacity': opac_range[i]
          },
          name = i + 1
        ) for i in scheme]
    
    layout = go.Layout(
        yaxis= dict(
            title= 'PC{}, var: {}'.format(str(pc2),str(round(pca_var[pc2-1],4)))
        ),
        xaxis= dict(
            title= 'PC{}, var: {}'.format(str(pc1),str(round(pca_var[pc1-1],4)))
        ),
        showlegend= True
        )
    
    fig = go.Figure(data=fig_data, layout=layout)
    iplot(fig)


def plot_accessions3D(selected_column,df= {},orderCore= {},vectors= [],allow_sbgp= [],allow_geo= [],color_ref= []):
    
    layout= go.Layout(
        margin=dict(
            l=0,
            r=0,
            b=0,
            t=0
        ),
        xaxis= dict(
            title= 'PC1',
        ),
        showlegend= True
        )
    names_index = [[f for f in orderCore.ID].index(x) for x in [str(y) for y in df.iloc[:,1]]]
    opac= .8
    if selected_column == -2:
        scheme= [orderCore.loc[x,'COUNTRY'] for x in names_index]
        allow= allow_geo
        scheme= [['Other',x][int(x in allow)] for x in scheme]
        coords = {y:[x for x in range(len(scheme)) if scheme[x] == y and x in names_index] for y in list(set(scheme))}
        
        fig= [go.Scatter3d(
        x = df.iloc[coords[i],2],
        y = df.iloc[coords[i],3],
        z = df.iloc[coords[i],4],
        mode= "markers",
        text= orderCore.iloc[[names_index[x] for x in coords[i]],:][["ID","NAME","COUNTRY","Initial_subpop"]].apply(lambda lbgf: (
      "<b>{}</b><br>Name: {}<br>Country: {}<br>{}".format(lbgf[0],lbgf[1],lbgf[2],lbgf[3])),
        axis= 1),
        marker= {
        'line': {'width': 0},
        'size': 4,
        'symbol': 'circle',
      "opacity": opac
      },
      name= i
    ) for i in coords.keys() if coords[i]]
    if selected_column == -1:
        scheme= [orderCore.loc[x,'Initial_subpop'] for x in names_index]
        
        allow= allow_sbgp
        scheme= [['Other',x][int(x in allow)] for x in scheme]
        
        coords = {y:[x for x in range(len(scheme)) if scheme[x] == y and x in names_index] for y in list(set(scheme))}
        pop_refs= ["Indica","cAus","Japonica","GAP","cBasmati","Admix"]
        color_here= color_ref
        
        fig= [go.Scatter3d(
        x = df.iloc[coords[i],2],
        y = df.iloc[coords[i],3],
        z = df.iloc[coords[i],6],
        mode= "markers",
        text= orderCore.iloc[[names_index[x] for x in coords[i]],:][["ID","NAME","COUNTRY","Initial_subpop"]].apply(lambda lbgf: (
      "<b>{}</b><br>Name: {}<br>Country: {}<br>{}".format(lbgf[0],lbgf[1],lbgf[2],lbgf[3])),
        axis= 1),
        marker= {
        'line': {'width': 0},
        'size': 4,
        'symbol': 'circle',
      "opacity": opac
      },
      name= i
    ) for i in coords.keys() if coords[i]]
    if selected_column == 0:
        scheme = [orderCore.loc[x,'sNMF_K3'] for x in names_index]
        coords = {y:[x for x in range(len(scheme)) if scheme[x] == y and x in names_index] for y in list(set(scheme))}
        
        pop_refs= ["Indica","cAus","Japonica","GAP","cBasmati","Admix"]
        color_here= color_ref
        
        fig= [go.Scatter3d(
        x = df.iloc[coords[i],2],
        y = df.iloc[coords[i],3],
        z = df.iloc[coords[i],6],
        mode= "markers",
        text= orderCore.iloc[[names_index[x] for x in coords[i]],:][["ID","NAME","COUNTRY","Initial_subpop"]].apply(lambda lbgf: (
      "<b>{}</b><br>Name: {}<br>Country: {}<br>{}".format(lbgf[0],lbgf[1],lbgf[2],lbgf[3])),
        axis= 1),
        marker= {
#        'color': scheme,
        'color': color_here[i],
        'line': {'width': 0},
        'size': 4,
        'symbol': 'circle',
      "opacity": opac
      },
      name= pop_refs[i]
    ) for i in list(set(scheme)) if coords[i]]
    
    if selected_column > 0:
        fig= [go.Scatter3d(
            x = df.iloc[:,2],
            y = df.iloc[:,3],
            z = df.iloc[:,6],
            mode= "markers",
            text= orderCore.iloc[names_index,:][["ID","NAME","COUNTRY","Initial_subpop"]].apply(lambda lbgf: (
          "<b>{}</b><br>Name: {}<br>Country: {}<br>{}".format(lbgf[0],lbgf[1],lbgf[2],lbgf[3])),
            axis= 1),
        marker= {
            'color': vectors.iloc[:,selected_column - 1],
                'colorbar': go.scatter3d.marker.ColorBar(
                    title= 'ColorBar'
                ),
                'colorscale': 'Viridis',
            'line': {'width': 0},
            'size': 4,
            'symbol': 'circle',
          "opacity": opac
          }
        )]
            
    fig = go.Figure(data=fig,layout= layout)
    iplot(fig)

def plot_accessions(selected_column,pc1,pc2,df= {},orderCore= {},vectors= [],allow_sbgp= [],allow_geo= [],color_ref= [],pca_var= []):
    layout= go.Layout(
        title= 'cluster: {}'.format(selected_column),
        yaxis= dict(
            title= 'PC{}, var: {}'.format(str(pc2),str(round(pca_var[pc2-1],4)))
        ),
        xaxis= dict(
            title= 'PC{}, var: {}'.format(str(pc1),str(round(pca_var[pc1-1],4)))
        ),
        showlegend= True
        )
    names_index = [[f for f in orderCore.ID].index(x) for x in [str(y) for y in df.iloc[:,1]]]
    opac= .8
    soiz= 8
    if selected_column == -2:
        scheme= [orderCore.loc[x,'COUNTRY'] for x in names_index]
        allow= allow_geo
        scheme= [['Other',x][int(x in allow)] for x in scheme]
        coords = {y:[x for x in range(len(scheme)) if scheme[x] == y] for y in list(set(scheme))}
        
        fig= [go.Scatter(
        x = df.iloc[coords[i],pc1+1],
        y = df.iloc[coords[i],pc2+1],
        mode= "markers",
        text= orderCore.iloc[[names_index[x] for x in coords[i]],:][["ID","NAME","COUNTRY","Initial_subpop"]].apply(lambda lbgf: (
      "<b>{}</b><br>Name: {}<br>Country: {}<br>{}".format(lbgf[0],lbgf[1],lbgf[2],lbgf[3])),
        axis= 1),
        marker= {
        'line': {'width': 0},
        'size': soiz,
        'symbol': 'circle',
      "opacity": opac
      },
      name= i
    ) for i in coords.keys() if coords[i]]
    if selected_column == -1:
        scheme= [orderCore.loc[x,'Initial_subpop'] for x in names_index]
        allow= allow_sbgp
        scheme= [['Other',x][int(x in allow)] for x in scheme]
        
        coords = {y:[x for x in range(len(scheme)) if scheme[x] == y] for y in list(set(scheme))}
        pop_refs= ["Indica","cAus","Japonica","GAP","cBasmati","Admix"]
        color_here= color_ref
        
        fig= [go.Scatter(
        x = df.iloc[coords[i],pc1+1],
        y = df.iloc[coords[i],pc2+1],
        mode= "markers",
        text= orderCore.iloc[[names_index[x] for x in coords[i]],:][["ID","NAME","COUNTRY","Initial_subpop"]].apply(lambda lbgf: (
      "<b>{}</b><br>Name: {}<br>Country: {}<br>{}".format(lbgf[0],lbgf[1],lbgf[2],lbgf[3])),
        axis= 1),
        marker= {
        'line': {'width': 0},
        'size': soiz,
        'symbol': 'circle',
      "opacity": opac
      },
      name= i
    ) for i in coords.keys()]
    if selected_column == 0:
        scheme = [int(orderCore.sNMF_K3[x]) for x in names_index]
        coords = {y:[x for x in range(len(scheme)) if scheme[x] == y and x in names_index] for y in list(set(scheme))}
        
        pop_refs= ["Indica","cAus","Japonica","GAP","cBasmati","Admix"]
        color_here= color_ref
        
        fig= [go.Scatter(
        x = df.iloc[coords[i],pc1+1],
        y = df.iloc[coords[i],pc2+1],
        mode= "markers",
        text= orderCore.iloc[[names_index[x] for x in coords[i]],:][["ID","NAME","COUNTRY","Initial_subpop"]].apply(lambda lbgf: (
      "<b>{}</b><br>Name: {}<br>Country: {}<br>{}".format(lbgf[0],lbgf[1],lbgf[2],lbgf[3])),
        axis= 1),
        marker= {
        'color': color_here[i],
        'line': {'width': 0},
        'size': soiz,
        'symbol': 'circle',
      "opacity": opac
      },
      name= pop_refs[i]
    ) for i in list(set(scheme)) if coords[i]]
    
    if selected_column > 0:
        fig= [go.Scatter(
            x = df.iloc[:,pc1+1],
            y = df.iloc[:,pc2+1],
            mode= "markers",
            text= orderCore.iloc[names_index,:][["ID","NAME","COUNTRY","Initial_subpop"]].apply(lambda lbgf: (
          "<b>{}</b><br>Name: {}<br>Country: {}<br>{}".format(lbgf[0],lbgf[1],lbgf[2],lbgf[3])),
            axis= 1),
        marker= {
            'color': vectors.iloc[:,(selected_column - 1)],
            'colorbar': go.scatter.marker.ColorBar(
                title= 'ColorBar'
            ),
            'colorscale': 'Viridis',
            'line': {'width': 0},
            'size': soiz,
            'symbol': 'circle',
          "opacity": opac
          }
        )]
            
    fig = go.Figure(data=fig,layout= layout)
    iplot(fig)


def plot_countries(selected_column,threshold= 0.05,df= {},vectors= {},orderCore= {},allow_geo= [],scope= 'asia',change_name= {}):

    world_codes= pd.read_csv('https://raw.githubusercontent.com/plotly/datasets/master/2014_world_gdp_with_codes.csv')

    names_index = [[f for f in orderCore.ID].index(x) for x in [str(y) for y in df.iloc[:,1]]]
    scheme= [orderCore.loc[x,'COUNTRY'] for x in names_index]
    allow= allow_geo
    scheme= [['Other',x][int(x in allow)] for x in scheme]

    coords = {y:[x for x in range(len(scheme)) if scheme[x] == y and x in names_index] for y in list(set(scheme))}
    
    
    lik_country= {
        z: np.mean([x for x in vectors.iloc[coords[z],(selected_column - 1)] if x >= threshold]) for z in allow_geo
    }
    
    if change_name:
        for ctry in change_name.keys():
            if ctry in lik_country.keys():
                lik_country[change_name[ctry]]= lik_country[ctry]
    

    scl = [[0.0, 'rgb(242,240,247)'],[0.05, 'rgb(218,218,235)'],[0.1, 'rgb(188,189,220)'],\
                [0.15, 'rgb(158,154,200)'],[0.2, 'rgb(117,107,177)'],[0.25, 'rgb(84,39,143)']]

    file_CTR= [str(x) for x in world_codes.COUNTRY]
    countries_geo= [x for x in lik_country.keys() if x in file_CTR]

    countries_code= [file_CTR.index(x) for x in countries_geo] 
    countries_code= [world_codes.CODE[x] for x in countries_code]

    data = [go.Choropleth(
        colorscale = 'Viridis',
        autocolorscale = False,
        locations = countries_code,
        z = [float(lik_country[x]) for x in countries_geo],
        #locationmode = 'USA-states',
        #text = countries_geo,
        marker = go.choropleth.Marker(
            line = go.choropleth.marker.Line(
                color = 'rgb(255,255,255)',
                width = 2
            )),
        colorbar = go.choropleth.ColorBar(
            title = "average likelihood")
    )]

    layout = go.Layout(
        title= 'geo_plot: {}; cl: {}'.format(scope,selected_column),
        geo = go.layout.Geo(
            scope = scope,
            projection = go.layout.geo.Projection(),
            showlakes = True,
            lakecolor = 'rgb(255, 255, 255)'),
    )

    fig = go.Figure(data = data, layout = layout)
    iplot(fig)
