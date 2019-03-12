import plotly.graph_objs as go
from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot

import numpy as np
import pandas as pd
import itertools as it

import os

from sklearn.neighbors import KernelDensity
from sklearn.decomposition import PCA
from sklearn.model_selection import GridSearchCV
from sklearn.cluster import estimate_bandwidth

from matplotlib import pyplot as plt

from matplotlib.collections import BrokenBarHCollection

import collections

def recursively_default_dict():
        return collections.defaultdict(recursively_default_dict)


def plot_clusters(selected,clusters= {},clusters_labels= {},ID= ''):
    scheme = [x for x in list(set(clusters_labels))]
    if selected == 0:
        opac_range = [.8]*len(list(set(clusters_labels)))
    else:
        opac_range = [.3]*len(list(set(clusters_labels)))
        opac_range[selected-1] = 1
    fig_data= [go.Scatter(
        y = clusters.iloc[[x for x in range(len(clusters_labels)) if clusters_labels[x] == i],:].pc1,
        x = clusters.iloc[[x for x in range(len(clusters_labels)) if clusters_labels[x] == i],:].pc2,
        #z = clusters.iloc[[x for x in range(len(clusters_labels)) if clusters_labels[x] == i],:].pc3,
        #type='scatter',
        mode= "markers",
        marker= {
            'line': {'width': 0},
            'size': 4,
            'symbol': 'circle',
            'opacity': opac_range[i]
          },
          name = str(i + 1)
        ) for i in scheme]
    
    layout = go.Layout(
        title= '{} Distance to references PCA'.format(ID),
        xaxis= dict(
            title= 'PC2'
        ),
        yaxis= dict(
            title= 'PC1'
        ),
        showlegend= True
        )
    
    fig = go.Figure(data=fig_data, layout=layout)
    iplot(fig)


def plot_accessions(selected_column, 
                    df= {},
                    vectors= {},
                    orderCore= {},
                    allow_geo= [],
                    allow_sbgp= [],
                    color_ref= [],
                    opac= .8,
                    cols_text= ["ID","NAME","COUNTRY","Initial_subpop"]):
    
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
    names_index = [[f for f in orderCore.ID].index(x) for x in [str(y) for y in df.iloc[:,0]]]
    opac= .8
    if selected_column == -2:
        scheme= [orderCore.loc[x,'COUNTRY'] for x in names_index]
        allow= allow_geo
        scheme= [['Other',x][int(x in allow)] for x in scheme]
        coords = {y:[x for x in range(len(scheme)) if scheme[x] == y and x in names_index] for y in list(set(scheme))}
        
        fig= [go.Scatter3d(
        x = df.iloc[coords[i],1],
        y = df.iloc[coords[i],2],
        z = df.iloc[coords[i],3],
        #type='scatter3d',
        mode= "markers",
        text= orderCore.iloc[[names_index[x] for x in coords[i]],:][cols_text].apply(lambda lbgf: (
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
        
        color_here= color_ref
        
        fig= [go.Scatter3d(
        x = df.iloc[coords[i],1],
        y = df.iloc[coords[i],2],
        z = df.iloc[coords[i],3],
        #type='scatter3d',
        mode= "markers",
        text= orderCore.iloc[[names_index[x] for x in coords[i]],:][cols_text].apply(lambda lbgf: (
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
        
        scheme = [int(orderCore.sNMF_K3[x]) for x in names_index]
        coords = {y:[x for x in range(len(scheme)) if scheme[x] == y and x in names_index] for y in list(set(scheme))}
        
        pop_refs= ["Indica","cAus","Japonica","GAP","cBasmati","Admix"]
        color_here= color_ref
        
        fig= [go.Scatter3d(
        x = df.iloc[coords[i],1],
        y = df.iloc[coords[i],2],
        z = df.iloc[coords[i],3],
        #type='scatter3d',
        mode= "markers",
        text= orderCore.iloc[[names_index[x] for x in coords[i]],:][cols_text].apply(lambda lbgf: (
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
    ) for i in coords.keys() if coords[i]]
    
    if selected_column > 0:
        fig= [go.Scatter3d(
            x = df.iloc[:,1],
            y = df.iloc[:,2],
            z = df.iloc[:,3],
            #type='scatter3d',
            mode= "markers",
            text= orderCore.iloc[names_index,:][["ID","NAME","COUNTRY","Initial_subpop"]].apply(lambda lbgf: (
          "<b>{}</b><br>Name: {}<br>Country: {}<br>{}".format(lbgf[0],lbgf[1],lbgf[2],lbgf[3])),
            axis= 1),
        marker= {
            'color': vectors.iloc[:,selected_column - 1],
                'colorbar': go.ColorBar(
                    title= 'ColorBar'
                ),
                'colorscale':[[0.0, 'rgb(165,0,38)'], [0.1111111111111111, 'rgb(215,48,39)'], [0.2222222222222222, 'rgb(244,109,67)'], [0.3333333333333333, 'rgb(253,174,97)'], [0.4444444444444444, 'rgb(254,224,144)'], [0.5555555555555556, 'rgb(224,243,248)'], [0.6666666666666666, 'rgb(171,217,233)'], [0.7777777777777778, 'rgb(116,173,209)'], [0.8888888888888888, 'rgb(69,117,180)'], [1.0, 'rgb(49,54,149)']],
            'line': {'width': 0},
            'size': 4,
            'symbol': 'circle',
          "opacity": opac
          }
        )]
            
    fig = go.Figure(data=fig,layout= layout)
    iplot(fig)


def return_densities(label,Distances= {},data= {},Centre_distances= {},Dr= ''):
    select_l1= [label]
    selected1= [x for x in range(Distances.shape[0]) if data[Dr]['labels_l1'][x] + 1 in select_l1]    

    x= np.mean(Distances[selected1,:],axis= 0)
    y= np.mean(Centre_distances[selected1,:],axis= 0)

    fig= [go.Scatter(
        x= x,
        y= y,
        mode= 'markers'  
    )]

    colorscale = ['#7A4579', '#D56073', 'rgb(236,158,105)', (1, 1, 0.2), (0.98,0.98,0.98)]

    trace1 = go.Scatter(
        x=x, y=y, mode='markers', name='points',
        marker=dict(color='rgb(102,0,0)', size=5, opacity=0.4)
    )
    trace2 = go.Histogram2dContour(
        x=x, y=y, name='density', ncontours=20,
        colorscale='Hot', reversescale=True, showscale=False
    )
    trace3 = go.Histogram(
        x=x, name='x density',
        marker=dict(color='rgb(102,0,0)'),
        yaxis='y2'
    )
    trace4 = go.Histogram(
        y=y, name='y density', marker=dict(color='rgb(102,0,0)'),
        xaxis='x2'
    )
    fig= {'data':[trace1, trace2, trace3, trace4]}

    fig['layout'] = go.Layout(
        showlegend=False,
        autosize=True,
        #width=1200,
        #height=550,
        xaxis=dict(
            domain=[0, 0.85],
            showgrid=False,
            zeroline=False
        ),
        yaxis=dict(
            domain=[0, 0.85],
            showgrid=False,
            zeroline=False
        ),
        margin=dict(
            t=50
        ),
        hovermode='closest',
        bargap=0,
        xaxis2=dict(
            domain=[0.85, 1],
            showgrid=False,
            zeroline=False
        ),
        yaxis2=dict(
            domain=[0.85, 1],
            showgrid=False,
            zeroline=False
        )
    )
    

    iplot(fig)

def return_cluster_refs(label,
                Z,
                df= {},
                orderCore= {},
               Ref_stats= {},
                Distances= {},
                Centre_distances= {},
                Slice= {},
                data= {},
                coeff_list= {},
                const_list= {},
                color_ref= [],
                Dr= '',
                ID= ''):
    
    select_l1= [label]
    trim= [x for x in range(Ref_stats.shape[0]) if Ref_stats[x,1] > 0.1 and Slice.Nsnps[x] > 50]
    selected1= [x for x in range(Distances.shape[0]) if data[Dr]['labels_l1'][x] + 1 in select_l1 and x in trim]    
    
    names_index = [[f for f in orderCore.ID].index(x) for x in [str(y) for y in df.iloc[:,0]]]

    ###
    scheme = [int(orderCore.sNMF_K3[x]) for x in names_index]
    scheme= [x for x in scheme if x != 5]
    coords = {y:[x for x in range(len(scheme)) if scheme[x] == y and x in names_index] for y in list(set(scheme))}
    
    print(coords.keys())
    pop_refs= ["Indica","cAus","Japonica","GAP","cBasmati","Admix"]
    color_here= color_ref

    ###
    
    ref_names= ['Indica','cAus','Japonica','cBasmati','Admix']

    #scheme= [orderCore.loc[x,'sNMF_K3'] for x in names_index]
    #scheme= {z:[x for x in range(len(scheme)) if scheme[x] == z] for z in list(set(scheme))}
    if Z== 0:
        x= np.mean(Distances[selected1,:],axis= 0)
        print(len(scheme))
        y= np.mean(Centre_distances[selected1,:],axis= 0)
        
        max_dists= [max(Centre_distances[x]) for x in selected1]
        min_centre= [min(Distances[x]) for x in selected1]
        
        Normalized_mean= np.mean(max_dists)
        Normalized_centre= np.mean(min_centre)
        
        range_raw_X= [-4,11]
        range_raw_Y= [-2,8.5]
        x_lab= 'PCA'
        y_lab= 'PCA'
        
    if Z== 1:
        meansVars= Ref_stats[selected1,:]
        trim= [x for x in range(meansVars.shape[0]) if meansVars[x,1] > 0.05]
        
        max_dists= [max(Centre_distances[x]) for x in selected1]
        min_centre= [min(Distances[x]) for x in selected1]
        meansVars= meansVars[trim,:]
        select_trim= [selected1[x] for x in trim]
        
        Normalized_mean = [(max_dists[x] - meansVars[x,0]) / meansVars[x,1] for x in range(meansVars.shape[0])]
        Normalized_mean= np.array(Normalized_mean)
        #Normalized_mean= [x for x in Normalized_mean if np.isnan()]
        Normalized_mean[Normalized_mean > 20] = 10
        Normalized_mean[Normalized_mean < -20] = -10
        Normalized_mean= np.mean(Normalized_mean)
        
        Normalized_centre = [(min_centre[x] - meansVars[x,0]) / meansVars[x,1] for x in range(meansVars.shape[0])]
        Normalized_centre= np.array(Normalized_centre)
        #Normalized_mean= [x for x in Normalized_mean if np.isnan()]
        Normalized_centre[Normalized_centre > 20] = 10
        Normalized_centre[Normalized_centre < -20] = -10
        Normalized_centre= np.mean(Normalized_centre)
        
        sel_d= Distances[select_trim,:]
        sel_d= [(sel_d[x] - meansVars[x,0]) / meansVars[x,1] for x in range(len(select_trim))]
        x= np.mean(sel_d,axis=0)
        
        sel_refd=Centre_distances[select_trim,:]
        sel_refd= [(sel_refd[x] - meansVars[x,0]) / meansVars[x,1] for x in range(len(select_trim))]
        y= np.mean(sel_refd,axis=0)
        
        range_raw_X= [-4,11]
        range_raw_Y= [-2,8.5]
        
        x_lab= 'scaled PCA '
        y_lab= 'scaled PCA'

    if Z== 2:
        coeffs_select= [coeff_list[x] for x in selected1]
        const_select= [const_list[x] for x in selected1]
        
        X= Distances[selected1,:]
        
        Y= Centre_distances[selected1,:]
        
        X= [np.exp(coeffs_select[x]*np.log(X[x]) + const_select[x]) for x in range(len(X))]
        Y= [np.exp(coeffs_select[x]*np.log(Y[x]) + const_select[x]) for x in range(len(Y))]
        
        x= np.mean(X,axis= 0)
        y= np.mean(Y,axis= 0)
        
        max_dists= [max(Centre_distances[x]) for x in selected1]
        min_centre= [min(Distances[x]) for x in selected1]
        
        range_raw_X= [-.02,1]
        range_raw_Y= [-.02,1]
                
        Normalized_mean= np.mean(max_dists)
        Normalized_centre= np.mean(min_centre)
        x_lab= 'Fst'
        y_lab= 'Fst'

    
    trace1 = [go.Scatter(
        x= [x[z] for z in coords[i]],
        y= [y[z] for z in coords[i]],
        mode= 'markers',
        name= ref_names[i],
        text= orderCore.iloc[[names_index[x] for x in coords[i]],:][["ID","NAME","COUNTRY","Initial_subpop"]].apply(lambda lbgf: (
      "<b>{}</b><br>Name: {}<br>Country: {}<br>{}".format(lbgf[0],lbgf[1],lbgf[2],lbgf[3])),
        axis= 1),
        marker=dict(color=color_ref[i], size=7, opacity=0.6)
    ) for i in coords.keys()]# if coords[i] and i != 5]


    trace1.append(go.Scatter(
        x= [Normalized_centre],
        y= [Normalized_mean],
        mode= 'markers',
        marker=dict(color='black', size=10, opacity=0.6)
    ))

    layout = go.Layout(
        title= '{} query, {} distances label {}; Target mean: {}'.format(ID,x_lab,select_l1[0],round(Normalized_mean,2)),
        showlegend=False,
        autosize=True,
        xaxis= dict(
        range= range_raw_X,
        title= '{} euc distances'.format(y_lab)
        ),
        yaxis= dict(
        range= range_raw_Y,
        title= '{} distances to center'.format(y_lab))
    )

    fig= go.Figure(data=trace1,layout= layout)

    iplot(fig)



def return_refs(t1,t2,registered= 'inlier',Z= 2,Wsnps= 50,df= {},Ref_stats= {},Slice= {},orderCore= {},Centre_distances= {},Distances= {},coeff_list= [],const_list= [],color_ref= [],ID= ''):
    
    trim= [x for x in range(Ref_stats.shape[0]) if Ref_stats[x,1] > 0.1 and Slice.Nsnps[x] > 60]   
    
    
    Normalized_mean = [(Ref_stats[x,2] - Ref_stats[x,0]) / Ref_stats[x,1] for x in trim]
    
    Dists_mean= np.mean(Normalized_mean)
    Dists_std= np.std(Normalized_mean)
    Dists_median= np.median(Normalized_mean)

    from scipy.stats import norm

    Dists_endpoints= norm.interval(0.95, loc=Dists_mean, scale=Dists_std)
    
    Normalized_vec= Normalized_mean

    if registered == 'inlier':
        t1, t2= Dists_endpoints
        #Normalized_mean= Normalized_vec
     
    Normalized_mean= np.nan_to_num(Normalized_mean)
    Normalized_mean[Normalized_mean > 20] = 10
    Normalized_mean[Normalized_mean < -20] = -10
    
    coords_threshold= [int(x > t1 and x < t2) for x in Normalized_mean]
    coords_threshold= {z:[x for x in range(len(coords_threshold)) if coords_threshold[x] == z]  for z in list(set(coords_threshold))}
    
    names_index = [[f for f in orderCore.ID].index(x) for x in [str(y) for y in df.iloc[:,0]]]
    
    average_selected= round(np.mean([Normalized_mean[x] for x in coords_threshold[1] if Normalized_vec[x] > - 20 and Normalized_vec[x] < 20]),2)
    
    prop_thresh= round(len(coords_threshold[1]) / float(len(Normalized_mean)) * 100,2)
    
    selected1= coords_threshold[1]
    
    ###
    scheme = [int(orderCore.sNMF_K3[x]) for x in names_index]
    scheme= [x for x in scheme]
    coords = {y:[x for x in range(len(scheme)) if scheme[x] == y and x in names_index] for y in list(set(scheme))}
    
    print(coords.keys())
    pop_refs= ["Indica","cAus","Japonica","GAP","cBasmati","Admix"]
    color_here= color_ref

    ###
    
    ref_names= ['Indica','cAus','Japonica','cBasmati','Admix']

    #scheme= [orderCore.loc[x,'sNMF_K3'] for x in names_index]
    #scheme= {z:[x for x in range(len(scheme)) if scheme[x] == z] for z in list(set(scheme))}
    if Z== 'raw':

        X= np.mean(Distances[selected1,:],axis= 0)
        print(len(scheme))
        Y= np.mean(Centre_distances[selected1,:],axis= 0)
        range_raw_X= [min(X)-1,max(X)+1]
        range_raw_Y= [min(Y)-1,max(Y)+1]
        title_raw= 'PCA distances'
        
        max_dists= [max(Centre_distances[x]) for x in selected1]
        min_centre= [min(Distances[x]) for x in selected1]
        
        Normalized_mean= np.mean(max_dists)
        Normalized_centre= np.mean(min_centre)
        
    if Z== 'scaled':
        meansVars= Ref_stats[selected1,:]
        trim= [x for x in range(meansVars.shape[0]) if meansVars[x,1] > 0.05]
        
        meansVars= meansVars[trim,:]
        select_trim= [selected1[x] for x in trim]
        
        max_dists= [max(Centre_distances[x]) for x in select_trim]
        min_centre= [min(Distances[x]) for x in select_trim]
        
        average_selected= round(np.mean(meansVars[:,2]),2)
        
        sel_d= Distances[select_trim,:]
        sel_d= [(sel_d[x] - meansVars[x,0]) / meansVars[x,1] for x in range(len(select_trim))]
        X= np.mean(sel_d,axis=0)
        
        sel_refd=Centre_distances[select_trim,:]
        sel_refd= [(sel_refd[x] - meansVars[x,0]) / meansVars[x,1] for x in range(len(select_trim))]
        Y= np.mean(sel_refd,axis=0)
        
        range_raw_X= [min(X)-1,max(X)+1]
        range_raw_Y= [min(Y)-1,max(Y)+1]
        title_raw= 'PCA distances'

    if Z== 'Fst':
        selected1= [x for x in selected1 if Slice.Nsnps[x] >= Wsnps]
        coeffs_select= [coeff_list[x] for x in selected1]
        const_select= [const_list[x] for x in selected1]
        
        X= Distances[selected1,:]
        print(len(scheme))
        Y= Centre_distances[selected1,:]
        
        X= [np.exp(coeffs_select[x]*np.log(X[x]) + const_select[x]) for x in range(len(X))]
        Y= [np.exp(coeffs_select[x]*np.log(Y[x]) + const_select[x]) for x in range(len(Y))]
        
        X= np.mean(X,axis= 0)
        Y= np.mean(Y,axis= 0)
        
        max_dists= [max(Centre_distances[x]) for x in selected1]
        min_centre= [min(Distances[x]) for x in selected1]
        
        range_raw_X= [-.02,1]
        range_raw_Y= [-.02,1]
        
        title_raw= 'Fst'
    
    trace1 = [go.Scatter(
        visible = True,
        x= [X[z] for z in coords[i]],
        y= [Y[z] for z in coords[i]],
        mode= 'markers',
        name= ref_names[i],
        text= orderCore.iloc[[names_index[x] for x in coords[i]],:][["ID","NAME","COUNTRY","Initial_subpop"]].apply(lambda lbgf: (
      "<b>{}</b><br>Name: {}<br>Country: {}<br>{}".format(lbgf[0],lbgf[1],lbgf[2],lbgf[3])),
        axis= 1),
        marker=dict(color=color_ref[i], size=7, opacity=0.6)
    ) for i in coords.keys() if coords[i]]
        
    layout = go.Layout(
        title= '{} distances; % {}; {} < Tz. < {} local sd; sel: {}'.format(ID,prop_thresh,round(t1,2),round(t2,2),average_selected),
        showlegend=False,
        autosize=True,
        xaxis= dict(
        range= range_raw_X,
        title= '{}, median= {}'.format(title_raw,round(Dists_median,3))
        ),
        yaxis= dict(
        range= range_raw_Y,
        title= '{} to center'.format(title_raw))
    )

    fig= go.Figure(data=trace1,layout= layout)

    iplot(fig)


def dist_Centre_Fst(Slice= {},Ref_stats= {},Wsnps= 50,coeff_list= {},const_list= {},ID= ''):
    selected1= [x for x in range(Ref_stats.shape[0]) if Ref_stats[x,1] > 0.1 and Slice.Nsnps[x] > Wsnps]   

    ###

    Means= Ref_stats[:,0]
    Fst_means= [coeff_list[x] * np.log(Means[x]) + const_list[x] for x in range(Ref_stats.shape[0])]
    Fst_means= [np.exp(x) for x in Fst_means]
    Fst_means= [Fst_means[x] for x in selected1]
    Varsity= Ref_stats[selected1,1]

    target_mean= Ref_stats[:,2]
    Fst_target= [coeff_list[x] * np.log(target_mean[x]) + const_list[x] for x in range(Ref_stats.shape[0])]
    Fst_target= [np.exp(x) for x in Fst_target]
    Fst_target= [Fst_target[x] for x in selected1]
    target_var= Ref_stats[selected1,3]


    ###
    ###

    X_plot = np.linspace(0, 1, 1000)

    bandwidth = estimate_bandwidth(np.array(Fst_means).reshape(-1,1), quantile=0.05, n_samples=1000)

    kde = KernelDensity(kernel='gaussian', bandwidth=bandwidth).fit(np.array(Fst_means).reshape(-1,1))

    log_dens = kde.score_samples(X_plot.reshape(-1,1))

    fig_roost_dens= [go.Scatter(x=X_plot, y=np.exp(log_dens), 
                                mode='lines', fill='tozeroy', name= 'Ref mean',
                                line=dict(color='blue', width=2))]

    ## Change means and variances to those of selected clusters


    ##

    X_plot = np.linspace(0, 1, 1000)

    bandwidth = estimate_bandwidth(np.array(Fst_target).reshape(-1,1), quantile=0.05, n_samples=5000)


    kde = KernelDensity(kernel='gaussian', bandwidth=bandwidth).fit(np.array(Fst_target).reshape(-1,1))

    log_dens = kde.score_samples(X_plot.reshape(-1,1))

    fig_roost_dens.append(go.Scatter(x=X_plot, y=np.exp(log_dens), 
                                mode='lines', fill='tozeroy', name= 'Target mean',
                                line=dict(color='red', width=2)))


    layout= go.Layout(
        title= '{}, Ztarget and normal Dists.'.format(ID),
        xaxis= dict(
        range= [-.1,1]
        )
    )

    fig = go.Figure(data=fig_roost_dens, layout= layout)
    iplot(fig)



def plot_clust_dist_vectors(select_l1,N_view= 15,N= 400, P= 40,
                           Distances= {}, Centre_distances= {}, 
                            coeff_list= {},
                            const_list= {},
                            data= {},
                            Dr= '',
                            Fst= False):
    #select_l2= [1,7]

    selected1= [x for x in range(Distances.shape[0]) if data[Dr]['labels_l1'][x] + 1 in select_l1]    
    
    Focus_dist= np.array(Distances[selected1])
    Focus_center= np.array(Centre_distances[selected1])
    
    if Fst:
        Focus_dist= [coeff_list[selected1[x]] * np.log(Focus_dist[x]) + const_list[selected1[x]] for x in range(len(Focus_dist))]
        Focus_dist= np.nan_to_num(np.array(Focus_dist))
        Focus_center= [coeff_list[selected1[x]] * np.log(Focus_center[x]) + const_list[selected1[x]] for x in range(len(Focus_dist))]
        Focus_center= np.nan_to_num(np.array(Focus_center))
        
        print(Focus_center.shape)

    
    #Centre_distances= (Centre_distances - np.mean(Focus_center,axis= 0)) / np.std(Focus_center,axis= 0)

    print('{} of clusters selected'.format(len(selected1)))

    Distances_mean= np.mean(Distances[selected1],axis= 0)


    range_distances= [np.percentile(Focus_dist,1),np.percentile(Focus_dist,99),N]
    range_cdist= [np.percentile(Focus_center,1),np.percentile(Focus_center,99),N]

    params = {'bandwidth': np.linspace(.05,.3,10)}
    grid = GridSearchCV(KernelDensity(algorithm = "ball_tree",breadth_first = False), params,verbose=0)
    distances_dens= []
    Cdist_dens= []

    i_coords, j_coords = np.meshgrid(np.linspace(range_distances[0],range_distances[1],P),
                          np.linspace(range_cdist[0],range_cdist[1],P),indexing= 'ij')
    traces= [x for x in it.product(range(P),range(P))]

    params_unique = {'bandwidth': np.linspace(.1, .3,20)}
    grid_unique = GridSearchCV(KernelDensity(algorithm = "ball_tree",breadth_first = False), params_unique,verbose=0)

    params_dens= {'bandwidth': np.linspace(.4, .8,20)}
    grid_dens = GridSearchCV(KernelDensity(algorithm = "ball_tree",breadth_first = False), params_dens,verbose=0)

    traces= [x for x in it.product(range(P),range(P))]

    background= np.array([i_coords, j_coords])

    background= [background[:,c[0],c[1]] for c in traces]
    background=np.array(background)

    for karl in np.random.choice(list(range(len(Focus_center))), N_view, replace= False):

        """
        kde = KernelDensity(kernel='gaussian', bandwidth=0.2).fit(np.array(Distances[karl,:])[:,np.newaxis])
        scores= kde.score_samples(np.linspace(*range_distances)[:,np.newaxis])


        distances_dens.append(scores)

        kde = KernelDensity(kernel='gaussian', bandwidth=0.2).fit(np.array(Centre_distances[karl,:])[:,np.newaxis])
        scores= kde.score_samples(np.linspace(*range_cdist)[:,np.newaxis])

        Cdist_dens.append(scores)


        ### Density measure
        datum = np.array([[Distances[karl,x],Centre_distances[karl,x]] for x in range(Distances.shape[1])])
        grid_dens.fit(datum)
        kde = grid_dens.best_estimator_

        P_dist= kde.score_samples(datum)
        scores= kde.score_samples(background)

        scores= np.exp(scores)
        scores= np.array([x for x in scipy.stats.norm(np.mean(scores),np.std(scores)).cdf(scores)])

        """
        ### haplotypes measure
        dotum= np.array([[Focus_dist[karl,x],Focus_center[karl,x]] for x in range(Focus_dist.shape[1])])

        datum= np.unique(dotum,axis= 0)
        if len(datum) < 3:
            datum= dotum

        grid_unique.fit(datum)
        kde = grid_unique.best_estimator_

        P_dist= kde.score_samples(datum)
        scores_haps= kde.score_samples(background)
        scores_haps= np.exp(scores_haps)

        #scores= scores_haps
        scores= np.hstack((scores_haps))
        #
        Cdist_dens.append(scores)


    #distances_dens= np.array(distances_dens)
    Cdist_dens= np.array(Cdist_dens)
    #coords= {i:[x for x in range(Distances.shape[0]) if data['labels_l1'][x] == i] for i in list(set(data['labels_l1']))}

    ### Cdist_dens must have been calculated.


    i_coords, j_coords = np.meshgrid(np.linspace(range_distances[0],range_distances[1],P),
                          np.linspace(range_cdist[0],range_cdist[1],P),indexing= 'ij')


    traces= [x for x in it.product(range(P),range(P))]

    background= np.array([i_coords, j_coords])

    background= [background[:,c[0],c[1]] for c in traces]
    background=np.array(background)

    Xs= []
    Ys= []
    Zs= []
    scores= []

    for target in range(N_view):
        Xs.extend(background[:,0])
        Ys.extend(background[:,1])
        scores.extend(Cdist_dens[target,:P**2])
        Zs.extend([target] * background.shape[0])


    thresho= .005
    tie= [x for x in range(len((scores))) if scores[x] < thresho]


    win= np.array([Xs,Ys,Zs,scores]).T
    unmasked= win[[x for x in range(win.shape[0]) if x not in tie],:]

    fig= [go.Scatter3d(
        x= unmasked[:,0],
        y= unmasked[:,1],
        z= unmasked[:,2],
        mode= 'markers',
        marker= {
            'color':unmasked[:,3],
            'colorbar': go.ColorBar(
                title= 'ColorBar'
            ),
            'colorscale':'Viridis',
            'line': {'width': 0},
            'size': 4,
            'symbol': 'circle',
          "opacity": .8
          }
    )]


    fig = go.Figure(data=fig)
    iplot(fig)
