# -*- coding: utf-8 -*-
"""
Created on Sat Nov 18 22:09:46 2017

@author: jgarcia
"""

import dash
from dash.dependencies import Input, Output
import dash_html_components as html
import dash_core_components as dcc
import plotly.graph_objs as go

import plotly.figure_factory as ff

import numpy as np
import pandas as pd

import base64

### read data sets

ID= "cBasmati_Jap"

Where= "12"

ref= "3"

pop_refs= ['unlabelled', '1', '3', '4']

###

Ordercore_file= 'Order_core.txt'

df = pd.read_csv('DIM_private_'+ ref +'_request_CHR'+ Where +'.'+ID+'.txt',sep= '\t',header= None)

cluster_pca = pd.read_csv('DIM_private_'+ref+'_comp_CHR'+Where+'.'+ID+'.txt',sep= '\t')

pca_var= [str(round(float(x),3)) for x in  cluster_pca.columns[1:(cluster_pca.shape[1]-1)]]

cluster_pca.columns = ['label','pc1','pc2','pc3','pc4','pc5','Na']

orderCore= pd.read_csv(Ordercore_file,sep = '\t')

vectors = pd.read_csv("Profile_"+ref+"_CHR"+Where+"."+ID+".txt",sep= '\t')

Coordinates= pd.read_csv('Profile_coordinates_' + ref + '_CHR'+Where+"."+ID+".txt",sep= '\t')

color_ref= "1"


#### code to label:

code_set= orderCore[['code','label']].drop_duplicates()


#### Prepqre Dash apps
app = dash.Dash(__name__)
server = app.server

app.css.append_css({'external_url': 'https://codepen.io/chriddyp/pen/dZVMbK.css'})


app.layout = html.Div([
    
    html.Div([
    html.H4(
    id= "header1",className= "six columns",children= "Feature space"
    ),
    html.H4(
    id= "header2",className= "six columns",children= "selected accessions"
    )    
    ],className= "row"),
    
    
    html.Div([
    dcc.Graph(id='local_pca', className='six columns',animate = True),
    html.Div(id='table-content', className='six columns')
    ],className='row'),
    
    html.Div([
    html.Div(dcc.Markdown(children="""**Fig. 1** relative distances among accessions given cluster profiles selected and analysed. 
    In fact, loadings plot\n of PCA run on the former."""),className='six columns'),
    html.Div(dcc.Markdown(children= """**Table 1** Passport information on accessions shown in Fig. 1. If cluster cloud is selected
    below, then only accessions in red (updated plot) are shown."""),className= 'six columns')
    ],className= "row"),
    
    html.Hr(),
    
    html.Div([
    html.H5(
    id= "header1",className= "six columns",children= "opacity"
    ),
    html.H5(
    id= "header2",className= "six columns",children= "Likelihood threshold"
    )    
    ],className= "row"),
    
    html.Div([

    dcc.Slider(
    updatemode= "drag",
    className='six columns',
    id = 'opacity',
    min = .1,
    max = 1,
    value = .8,
    step = .1,
    marks = {.1*x:str(round(.1*x,1)) for x in range(3,9)}
    ),
    dcc.Slider(
    updatemode= "drag",
    className= "six columns",
    id= "threshold",
    value= .1,
    min= .05,
    max= 1,
    marks= {i:str(round(i,2)) for i in np.arange(.05,1,.05)},
    step= .025
    )
    ],className='row'),

    
    html.Hr(),
    
    html.Div([
    html.H5(children= 'Chose a color')
    ],className= "row"),
    
    html.Div(dcc.Markdown(children= """Clusters profiles in **Fig. 3** were grouped by colour and given a number each (hover over the points to see number).
    Proportion of different cluster types is plotted below to help you chose interesting clusters to analyse.
    Proximity among cluster profiles correlates with individual contribution patterns. When a group is chosen, the density of mean individual likelihoods across
    that group's profiles is plotted below (**Fig. 5**). Accessions with mean likelihoods above the *Lilelihood threshold* selected above will appear in red in **Fig. 1**."""),className= 'row'),
    
    
    html.Div([
    dcc.Dropdown(
    className='six columns',
    id = 'chose_color',
    value = 0,
    options = [{"label":x,"value": x} for x in range(len(cluster_pca.label.unique())+1)]
    )
    ],className= "row"),
    
    html.Hr(),
    
    html.Div([
    dcc.Graph(id = "clusters",className="six columns"),
    dcc.Graph(id= "density_plot",className= "six columns")
    ]),

    html.Div([html.Div(dcc.Markdown(children= """**Fig.3** Relation among cluster profiles selected. In fact the distribution in feature space 
    of these profiles\n following principal component analysis."""),className= 'six columns'),
    html.Div(dcc.Markdown(children= """**Fig. 5** Density plot of average Likelihood by accession across cluster cloud selected."""),className= 'six columns')],
    className= 'row'),

    html.Div([
    dcc.Graph(id= "bars",className= "six columns")
    ]),
    
    html.Div([html.Div(dcc.Markdown(children= """**Fig. 5** proportion of cluster profiles by cloud (read *cluster*) in **clusters - observations**"""),
    className= 'six columns'),
    html.Hr(),
    html.H4(
    id= "header2",className= "row",children= "selected accessions"
    )]),
    
    html.Div([
    dcc.Dropdown(
    id= 'View',
    className= 'three columns',
    value= 1,
    options= [{'label':'ideogram','value':0},
                {'label': 'Image','value':1}]
    ),
    dcc.Dropdown(
    id= 'chromosome',
    className= 'three columns',
    value= int(Where),
    options= [{'label':x,'value':x} for x in Coordinates.chrom.unique()]
    ),
    dcc.Dropdown(
    id= 'profile',
    className= 'three columns',
    value= 0,
    options = [{"label":x,"value": int(x)} for x in range(len(list(set(cluster_pca.label)))+1)]
    )
    ],className= 'row'),
    
    html.Div(
    id= "ideogram"
    ),
    html.Hr(),
    
    html.Div(id = 'Intermediate_nothing',children='high',style= {'display': 'none'})

])




def generate_table(dataframe):
    return html.Table(
        # Header
        [html.Tr([html.Th(col) for col in dataframe.columns])] +

        # Body
        [html.Tr([
            html.Td(dataframe.iloc[i][col]) for col in dataframe.columns
        ]) for i in range(len(dataframe))]
    )



@app.callback(
    Output('ideogram','children'),
    [Input('View','value'),
    Input('chromosome','value'),
    Input('profile','value')]
)
def return_Ideogram(View,chrom,selected_group):    
    if selected_group == 0:
        filter= 'all'
    else:
        filter= str(selected_group - 1)
    image_filename = 'Ideos/Ideo_id.{0}_label.{1}_CHR{2}.png'.format(ID,filter,str(chrom).zfill(2))
    encoded_image = base64.b64encode(open(image_filename, 'rb').read())
    return [html.Img(id= 'spore',src='data:image/png;base64,{}'.format(encoded_image.decode()))]


@app.callback(
   Output('density_plot','figure'),
    [Input('chose_color','value')])
def update_density(selected_group):
    if selected_group != 0:
        dense_plot=  ff.create_distplot([vectors.iloc[:,selected_group-1]], [selected_group-1])
        dense_plot['layout'].update(title='<b>likelihood density</b>')
        return dense_plot


@app.callback(
    Output('table-content','children'),
    [Input('chose_color','value'),
     Input('threshold','value')])
def update_table(selected_group,threshold):
    names_index = [[f for f in orderCore.ID].index(x) for x in [str(y) for y in df.iloc[:,1]]]
    if selected_group== 0:
        show_table = [names_index[x] for x in range(len(vectors))]
    else:
        show_table = [names_index[x] for x in range(len(vectors)) if vectors.iloc[x,selected_group-1] >= threshold]
    return html.Div(
        children= generate_table(orderCore.iloc[show_table,:]),
        style={
            'overflowX': 'scroll',
            'overflowY': 'scroll',
            'height': '450px',
            'display': 'block',
            'paddingLeft': '15px'
        }
)

@app.callback(
    Output('bars','figure'),
    [Input('chose_color','value')])
def cluster_bars(selected):
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
    return fig


@app.callback(
    Output("clusters","figure"),
    [Input('chose_color','value')])
def update_secondFigure(selected):
    scheme = [x for x in list(set(cluster_pca.label))]
    if selected == 0:
        opac_range = [.8]*len(list(set(cluster_pca.label)))
    else:
        opac_range = [.3]*len(list(set(cluster_pca.label)))
        opac_range[selected-1] = 1
    return {'data': [go.Scatter3d(
        x = cluster_pca.iloc[[x for x in range(len(cluster_pca)) if cluster_pca.iloc[x,0] == i],:].pc1,
        y = cluster_pca.iloc[[x for x in range(len(cluster_pca)) if cluster_pca.iloc[x,0] == i],:].pc2,
        z = cluster_pca.iloc[[x for x in range(len(cluster_pca)) if cluster_pca.iloc[x,0] == i],:].pc3,
        type='scatter3d',
        mode= "markers",
        marker= {
            'line': {'width': 0},
            'size': 4,
            'symbol': 'circle',
            'opacity': opac_range[i]
          },
          name = i + 1
        ) for i in scheme],
        'layout': {
      "autosize": True, 
      "hovermode": "closest",
      "legend": {
        "x": 0.873529411765, 
        "y": 0.877829326396, 
        "borderwidth": 1, 
        "font": {"size": 13}
      },
      "scene": {
        "aspectmode": "auto", 
        "aspectratio": {
          "x": 1.02391505715, 
          "y": 0.737436541286, 
          "z": 1.3243763495
        }, 
        "camera": {
          "center": {
            "x": 0, 
            "y": 0, 
            "z": 0
          }, 
          "eye": {
            "x": 1.80578427889, 
            "y": 1.17729688569, 
            "z": 0.201532084509
          }, 
          "up": {
            "x": 0, 
            "y": 0, 
            "z": 1
          }
        }, 
        "xaxis": {
          "title": "PC1", 
          "type": "linear"
        }, 
        "yaxis": {
          "title": "PC2", 
          "type": "linear"
        }, 
        "zaxis": {
          "title": "PC3", 
          "type": "linear"
        }
      }, 
      "showlegend": True, 
      "title": "<b>clusters - observations</b>", 
      "xaxis": {"title": "V3"}, 
      "yaxis": {"title": "V2"}
    }
    }



@app.callback(
    Output(component_id= 'local_pca',component_property = 'figure'),
    [Input(component_id= 'chose_color',component_property = 'value'),
     Input(component_id= 'opacity',component_property= 'value'),
    Input(component_id= "threshold",component_property= "value")])
def update_figure(selected_column,opac,threshold):
    names_index = [[f for f in orderCore.ID].index(x) for x in [str(y) for y in df.iloc[:,1]]]
    fig = {'layout': {
  "autosize": True, 
  "hovermode": "closest",
  "legend": {
  "x": 0.798366013072, 
  "y": 0.786064620514, 
  "borderwidth": 1, 
  "font": {"size": 13}
      },
  "margin": {
    "r": 40, 
    "t": 50, 
    "b": 20, 
    "l": 30, 
    "pad": 0
  }, 
  "scene": {
    "aspectmode": "auto", 
    "aspectratio": {
      "x": 1.02391505715, 
      "y": 0.737436541286, 
      "z": 1.3243763495
    }, 
    "camera": {
      "center": {
        "x": 0, 
        "y": 0, 
        "z": 0
      }, 
      "eye": {
        "x": 1.80578427889, 
        "y": 1.17729688569, 
        "z": 0.201532084509
      }, 
      "up": {
        "x": 0, 
        "y": 0, 
        "z": 1
      }
    }, 
    "xaxis": {
      "title": "PC1: " + pca_var[0], 
      "type": "linear"
    }, 
    "yaxis": {
      "title": "PC2: " + pca_var[1], 
      "type": "linear"
    }, 
    "zaxis": {
      "title": "PC3: " + pca_var[2], 
      "type": "linear"
    }
  }, 
  "showlegend": True, 
  "title": "<b>Accessions</b>", 
  "xaxis": {"title": "V3"}, 
  "yaxis": {"title": "V2"}
    }}
    if selected_column == 0:
        scheme = [orderCore.code[x] for x in names_index]
        
        coords = {y:[x for x in range(len(scheme)) if scheme[x] == y and x in names_index] for y in list(set(scheme))}
        
        color_here= color_ref
        
        fig['data'] = [go.Scatter3d(
        x = df.iloc[coords[i],2],
        y = df.iloc[coords[i],3],
        z = df.iloc[coords[i],4],
        type='scatter3d',
        mode= "markers",
        text= orderCore.iloc[[names_index[x] for x in coords[i]],:][["ID"]].apply(lambda lbgf: (
            "<b>ID: {}</b>".format(lbgf[0])),
        axis= 1),
        marker= {
        'line': {'width': 0},
        'size': 4,
        'symbol': 'circle',
      "opacity": opac
      },
      name= code_set.loc[code_set['code'] == i, 'label'].iloc[0]
    ) for i in list(set(scheme))]
    
    else:
        fig['data']= [go.Scatter3d(
            x = df.iloc[:,2],
            y = df.iloc[:,3],
            z = df.iloc[:,4],
            type='scatter3d',
            mode= "markers",
            text= orderCore.iloc[names_index,:][["ID"]].apply(lambda lbgf: (
          "<b>ID: {}</b>".format(lbgf[0])),
            axis= 1),
        marker= {
            'color': vectors.iloc[:,selected_column - 1],
            'colorbar': go.ColorBar(
                title= 'ColorBar'
            ),
            'colorscale': 'Viridis',
            'line': {'width': 0},
            'size': 4,
            'symbol': 'circle',
          "opacity": opac
          }
        )]
    return fig 


if __name__ == '__main__':
    app.run_server(debug=True)
    