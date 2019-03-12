

import numpy as np
import pandas as pd
import itertools as it
import time

import scipy

from sklearn.neighbors import KernelDensity
from sklearn.decomposition import PCA
from sklearn.model_selection import GridSearchCV
from sklearn.cluster import estimate_bandwidth
from sklearn.cluster import MeanShift, estimate_bandwidth

from sklearn.metrics.pairwise import pairwise_distances
from sklearn.metrics.pairwise import euclidean_distances

import re
import matplotlib.pyplot as plt

from matplotlib.collections import BrokenBarHCollection

import collections

def recursively_default_dict():
        return collections.defaultdict(recursively_default_dict)


######################################################################################
### Load data 


def read_coords(coords_file):
    
    stock= recursively_default_dict()
    
    Input= open(coords_file,'r')
    front= []
    d= 0
    for line in Input:
        line= line.split()
        
        if d== 0:
            cols= line
            stock= {x:[] for x in line}
            d += 1
        else:
            for c in range(len(line)):
                stock[cols[c]].append(int(line[c]))
    
    Input.close()
    
    label_indexes= {z:[x for x in range(len(stock['label'])) if stock['label'][x] == z] for z in list(set(stock['label']))}
    order_by_label= {lab:[[stock['chrom'][x],stock['start'][x],stock['end'][x],stock['cluster'][x]] for x in label_indexes[lab]] for lab in label_indexes.keys()}
        
    return order_by_label


def read_books(books):
    library= []
    
    for shelf in books:
        card= shelf.split('/')
        
        cover= card[-1].split('_')
        Chr= int([re.findall(r'\d+',i)[0] for i in cover if re.search('CHR',i)][0])
        tag= cover[1]
        start= cover[2]
        library.append([shelf,tag,start,Chr])
    
    library= pd.DataFrame(library,columns= ['file','tag','start','Chr'])
    return library


def read_3D_profiles(Books):
    
    s0 = time.time()
    
    Blocks= recursively_default_dict()
    Profiles= recursively_default_dict()
    Out = recursively_default_dict()
    Names= []
    
    
    for Chr in Books.Chr.unique():
        
        Shelf= Books[(Books.Chr== Chr)]
        
        Blocks[Chr]= recursively_default_dict()
        Profiles[Chr]= recursively_default_dict()
        
        start= '3'
        
        for start in Shelf.start.unique():
            
            Series= Shelf[(Shelf.start== start)]
            
            if len(Series) != 2:
                print('{0} files were provided for analysis Chr: {1} - start: {2}. 2 files are expected.'.format(len(Series),Chr,start))
                break
            
            blocks_file= Series.loc[Series['tag'] == 'Request','file']
            blocks_file= blocks_file.iloc[0]
            
            profiles_file = Series.loc[Series['tag'] == 'profiles','file']
            profiles_file= profiles_file.iloc[0]
            ##### read blocks file
            Input= open(blocks_file,'r')
            d= 0
            
            for line in Input:
                line= line.split()
                
                if d== 0:
                    line=line[4:]
                    Names= line
                else:
                    Blocks[int(line[0])][int(line[1])][int(line[3])]= [float(x) for x in line[4:]]
                    if line[3] == '0':
                        Out[int(line[0])][int(line[1])] = int(line[2])
                d += 1
            
            Input.close()
            
            ##### read profiles file
            
            Input= open(profiles_file,'r')
            d= 0
            
            for line in Input:
                line= line.split()
                
                if d== 0:
                    line=line[3:]
                    Names= line
                else:
                    Profiles[int(line[0])][int(line[1])][int(line[2])]= [float(x) for x in line[3:]]
                
                d += 1
            
            Input.close()
    
    s1= time.time()
    
    print(str(s1-s0) + ' seconds elapsed.')
    print(str((s1-s0) / float(60)) + ' minutes elapsed.')
    return Blocks, Profiles, Names, Out


def Distance_profiles(vector,data_now,n_comp,label_select):
    pca = PCA(n_components=n_comp, whiten=False,svd_solver='randomized').fit(data_now)
    feats= pca.transform(data_now)
    
    N= 50
    bandwidth = estimate_bandwidth(feats, quantile=0.15)
    params = {'bandwidth': np.linspace(np.min(feats), np.max(feats),30)}
    grid = GridSearchCV(KernelDensity(algorithm = "ball_tree",breadth_first = False), params,verbose=0)
    
    ## perform MeanShift clustering.
    ms = MeanShift(bandwidth=bandwidth, bin_seeding=True, cluster_all=False, min_bin_freq=5)
    ms.fit(feats)
    labels1 = ms.labels_
    #label_select = {y:[x for x in range(len(labels1)) if labels1[x] == y] for y in sorted(list(set(labels1))) if y != -1}
    
    ## Extract the KDE of each cluster identified by MS.
    Proxy_data= []
    label_select_labels= [z for z in it.chain(*[[x] * len(label_select[x]) for x in label_select.keys()])]
    Center_store= {}
    Proxy_indexes= {}
    distance_vecs= []
    
    for lab in label_select.keys():
        if len(label_select[lab]) < 3:
            continue
            
        Quanted_set= feats[label_select[lab],:]
        grid.fit(Quanted_set)
        
        kde = grid.best_estimator_
        Extract= kde.sample(N)
        
        center= np.mean(Extract,axis= 0)
        Center_store[lab]= center
        Proxy_indexes[lab]= [x for x in range((len(Center_store) - 1) * N, len(Center_store) * N)]        
        Return= pca.inverse_transform(Extract)
        
        #Return= data_now[np.random.choice(label_select[lab],N),:]
        Proxy_data.extend(Return)
    
    ### get distances to other centers:
    Distances_vectors= []
    distances= euclidean_distances(vector.reshape(1,-1),feats)
    X_plot = np.linspace(-8, 2, 1000)
    
    kde = KernelDensity(kernel='gaussian', bandwidth=0.05).fit(np.array(distances).reshape(-1,1))
    log_dens = kde.score_samples(np.array(X_plot).reshape(-1,1))
       
    return log_dens



