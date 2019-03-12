import numpy as np
import collections
import StructE_tools as Ste
from scipy.stats import beta
from sklearn.decomposition import PCA
from sklearn.metrics.pairwise import pairwise_distances
from sklearn.metrics.pairwise import euclidean_distances
from sklearn.preprocessing import scale

from IPython.display import clear_output

def recursively_default_dict():
        return collections.defaultdict(recursively_default_dict)

import plotly.graph_objs as go
from plotly.offline import iplot

### Converting distances to fst's.
def Euc_to_fst(vector_lib,n_comp= 5,pop_max= 8,Iter= 20,bias_range= [20,300],Eigen= False, Scale= False,Centre= True,ploidy= 1):
    ### Select pre and post processing measures. 
        
    length_haps= vector_lib.shape[1]
    
    Iter= 20 # repeats
    
    #### Predict
    predicted= []

    #def controled_fsts(vector_lib,Eigen,length_haps,Scale,Center,N_pops,n_comp,Iter,N_sims,MixL,MixP,Pairs):
    lengths_vector= []

    ### store distances between centroids
    biased_pairwise= []

    ### store PC projection:
    dist_PC_corrected= {x:[] for x in range(n_comp)}

    ### store fsts
    fst_store= []


    ### proceed.

    for rep in range(Iter):
        clear_output()
        
        N_pops= np.random.choice(range(3,pop_max),1,replace= False)[0]
        
        ## Population Sizes and labels
        bias_scheme= np.random.choice(range(bias_range[0],bias_range[1]),N_pops,replace= False)
        
        bias_labels= np.repeat(np.array([x for x in range(N_pops)]),bias_scheme)
        
        ### triangular matrices extract.
        iu1= np.triu_indices(N_pops,1) # for centroid comparison

        iu_bias= np.triu_indices(sum(bias_scheme),1)

        iu_control= np.triu_indices(2,1)

        Pops= np.random.choice(vector_lib.shape[0],N_pops,replace= False)
        print('Iter: {}, vectors selected: {}, hap length: {}'.format(rep,Pops,length_haps))
        ########## FST

        freqs_selected= vector_lib[Pops,:length_haps]
        Pairwise= Ste.return_fsts2(freqs_selected)

        #fsts_compare = scale(Pairwise.fst)
        fsts_compare= Pairwise.fst
        
        fst_store.extend(fsts_compare)

        ## lengths
        lengths_vector.extend([length_haps] * len(fsts_compare))
        
        #### generate data and perform PCA
        data= []

        for k in range(N_pops):

            probs= vector_lib[Pops[k],:]

            m= bias_scheme[k]
            Haps= [[np.random.choice([ploidy,0],p= [1-probs[x],probs[x]]) for x in range(length_haps)] for acc in range(m)]

            data.extend(Haps)

        data2= np.array(data)

        if Scale:
            data2= scale(data2)

        pca = PCA(n_components=n_comp, whiten=False,svd_solver='randomized').fit(data2)

        feat_bias= pca.transform(data2)

        if Eigen:
            feat_bias= feat_bias * pca.explained_variance_ratio_

        #### Centroid distances
        
        bias_centroids= [np.mean(feat_bias[[y for y in range(feat_bias.shape[0]) if bias_labels[y] == z],:],axis= 0) for z in range(N_pops)]
        bias_centroids= np.array(bias_centroids)

        bias_pair_dist= pairwise_distances(bias_centroids,metric= 'euclidean')
        bias_pair_dist= bias_pair_dist[iu1]
        #bias_pair_dist= scale(bias_pair_dist)
        
        biased_pairwise.extend(bias_pair_dist)

    
    Size= length_haps
    fst_lm_range= [0,.3]
    
    Lindexes= [x for x in range(len(lengths_vector)) if lengths_vector[x] == Size and fst_store[x] >= fst_lm_range[0] and fst_store[x] <= fst_lm_range[1]]
    y_true= [np.log(biased_pairwise[x]) for x in Lindexes]
    fst_x= [np.log(fst_store[x]) for x in Lindexes]
    m_coeff,b= np.polyfit(y_true,fst_x,1)
    
    return m_coeff, b, biased_pairwise, fst_x, y_true


def Fst_predict(vector_lib,m_coeff,b,n_comp= 5,pop_max= 8,Iter= 20,bias_range= [20,300],Eigen= False, Scale= False,Centre= True,ploidy= 1):
    ### Select pre and post processing measures. 
    
    length_haps= vector_lib.shape[1]
        
    print('length haps: {}, N iterations: {}, range pops: {}'.format(length_haps,Iter,pop_max))
    
    #### Predict
    predicted= []

    #def controled_fsts(vector_lib,Eigen,length_haps,Scale,Center,N_pops,n_comp,Iter,N_sims,MixL,MixP,Pairs):
    lengths_vector= []

    ### store distances between centroids
    biased_pairwise= []

    ### store PC projection:
    dist_PC_corrected= {x:[] for x in range(n_comp)}

    ### store fsts
    fst_store= []


    ### proceed.

    for rep in range(Iter):
        
        N_pops= np.random.choice(range(3,pop_max),1,replace= False)[0]
        
        ## Population Sizes and labels
        bias_scheme= np.random.choice(range(bias_range[0],bias_range[1]),N_pops,replace= False)
        
        bias_labels= np.repeat(np.array([x for x in range(N_pops)]),bias_scheme)
        
        ### triangular matrices extract.
        iu1= np.triu_indices(N_pops,1) # for centroid comparison

        iu_bias= np.triu_indices(sum(bias_scheme),1)

        iu_control= np.triu_indices(2,1)

        Pops= np.random.choice(vector_lib.shape[0],N_pops,replace= False)
        #print('Iter: {}, vectors selected: {}, hap length: {}'.format(rep,Pops,length_haps))
        ########## FST

        freqs_selected= vector_lib[Pops,:length_haps]
        Pairwise= Ste.return_fsts2(freqs_selected)

        #fsts_compare = scale(Pairwise.fst)
        fsts_compare= Pairwise.fst
        
        fst_store.extend(fsts_compare)

        ## lengths
        lengths_vector.extend([length_haps] * len(fsts_compare))
        
        #### generate data and perform PCA
        data= []

        for k in range(N_pops):

            probs= vector_lib[Pops[k],:]

            m= bias_scheme[k]
            Haps= [[np.random.choice([ploidy,0],p= [1-probs[x],probs[x]]) for x in range(length_haps)] for acc in range(m)]

            data.extend(Haps)

        data2= np.array(data)

        if Scale:
            data2= scale(data2)

        pca = PCA(n_components=n_comp, whiten=False,svd_solver='randomized').fit(data2)

        feat_bias= pca.transform(data2)

        if Eigen:
            feat_bias= feat_bias * pca.explained_variance_ratio_

        #### Centroid distances
        
        bias_centroids= [np.mean(feat_bias[[y for y in range(feat_bias.shape[0]) if bias_labels[y] == z],:],axis= 0) for z in range(N_pops)]
        bias_centroids= np.array(bias_centroids)

        bias_pair_dist= pairwise_distances(bias_centroids,metric= 'euclidean')
        bias_pair_dist= bias_pair_dist[iu1]
        #bias_pair_dist= scale(bias_pair_dist)
        
        fst_pred= [np.exp(m_coeff*np.log(x) + b) for x in bias_pair_dist]
        predicted.extend(fst_pred)
        

    
    fig= [go.Scatter(
        x= fst_store,
        y= predicted,
        mode= 'markers'
    )]

    layout = go.Layout(
        title= 'test of prediction',
        yaxis=dict(
            title='predicted Fst'),
        xaxis=dict(
            title='observed Fst')
    )

    fig= go.Figure(data=fig, layout=layout)
    iplot(fig)

