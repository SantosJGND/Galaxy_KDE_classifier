# -*- coding: utf-8 -*-
"""
Created on Wed Oct 18 22:59:24 2017

@author: jgarcia
"""

"""
Created on Sun Mar 05 13:01:08 2017

@author: jgarcia
"""

import collections

import itertools as it

import numpy as np
from numpy import nan,mean

import time

from sklearn.neighbors import KernelDensity
from sklearn.decomposition import PCA
from sklearn.model_selection import GridSearchCV
from sklearn.cluster import MeanShift, estimate_bandwidth
from sklearn.ensemble import IsolationForest

import sys
sys.setrecursionlimit(10000)


import scipy
from scipy import stats
from scipy.signal import savgol_filter


from sklearn.cluster import DBSCAN

from random import randint

def recursively_default_dict():
        return collections.defaultdict(recursively_default_dict)


def OriginbySNMF(Geno_Q,t):
    """
    Classes individuals according to Group assignment by SNMF
    using user provided threshold (.8 advised). returns dict.
    """
    Geneo = open(Geno_Q,"r")
    Ind = 0
    Groups = recursively_default_dict()
    for line in Geneo:
        line = line.split()
        line = [float(x.strip("\n")) for x in line]
        if Ind == 0:
            Groups = {x:[] for x in range(1,len(line)+2)}
        bagged = 0
        for value in range(len(line)):
            if line[value] >= t:
                Groups[value+1] += [Ind]
                bagged += 1
        if bagged == 0:
            Groups[len(line)+1] += [Ind]
        Ind += 1
    return Groups




def Imiss_names(File):
    """
    Deprecated ID library function. used plink analysis output.
    """
    File = open(File,"r")
    Index = recursively_default_dict()
    d = 0
    for line in File:
        line = line.split()
        if line[0] == "FID":
            continue
        else:
            Index[d] = line[0]
            d += 1
    return Index





def Mass_merge(lmiss,Chr):
    """
    deprecated snp library funtion. used plink analysis output.
    """
    PosIndex = recursively_default_dict()
    lmiss = open(lmiss,"r")
    d = 1
    print("start")
    for line in lmiss:
        line = line.split()
        if line[0] == "CHR":
            continue
        if d == 1 and line[0] != "CHR":
            First = int(line[1])
        if int(line[0]) == int(Chr):
            PosIndex[d] = int(line[1]) - First
            d += 1
    lmiss.close()
    return PosIndex


def BIMread(bimFile):
    '''
    reads .bim file from plink genomic data.
    returns dictionary of {geno_snp_index: locus}
    '''
    File = open(bimFile,"r")
    #Nsnps = {x:recursively_default_dict() for x in range(1,13)}
    Nsnps = recursively_default_dict()
    Gindex = recursively_default_dict()
    for x in range(1,13):
        Nsnps[x] = recursively_default_dict()
    d = 0
    CHR = 0
    for line in File:
        line = line.split()
        if d == 0:
            CHR = int(line[0])
        if d != 0 and int(line[0]) != CHR:
            d = 0
            CHR = int(line[0])
        Nsnps[CHR][d] = [int(line[3]),line[4],line[5]]
        Gindex[CHR][int(line[3])] = [d]
        d += 1
    File.close()
    return Nsnps,Gindex



def FAMread(Famfile):
    '''
    reads plink .fam file for names of accesssions in geneo file.
    '''
    File = open(Famfile,"r")
    Inds = recursively_default_dict()
    d = 0
    for line in File:
        line = line.split()
        Inds[d] = line[0]
        Inds[line[0]] = d
        d += 1
    File.close()
    return Inds


def Geneo_CORE(FAM,Identifiers):
    '''
    Specific to CORE excel file where accessions to be kept are marked out.
    Returns dictionary of same structure as OriginbysNMF().
    '''
    Set = []
    Fam = {x:v for v,x in FAM.items()}
    IDs = open(Identifiers,'r')
    for line in IDs:
        line = line.split()
        if line[1] == 'c':
            Set.append(Fam[line[0]])
    return Set


def interpolate_nans(X):
    """Overwrite NaNs with column value interpolations."""
    col_means = np.array([int(np.nanmean(X[:,x])) for x in range(X.shape[1])])
    inds = np.where(np.isnan(X))
    X[inds] = np.take(col_means,inds[1])
    return X



def GenoQtoDict(Geno_Q):
    """
    Reads Q matrix output of sNMF analysis stores into dictionary.
    """
    GenoI = open(Geno_Q,"r")
    IndAx = {}
    Individual = 0
    
    for Id in GenoI:
         IndAx[Individual] = recursively_default_dict()
         Id = Id.split(" ")
         Id = [float(x.strip("\n")) for x in Id]
         for u in range(len(Id)):
             IndAx[Individual][u + 1] = Id[u]
         Individual += 1
    GenoI.close()
    return IndAx


def lognormalize(x):
    a = np.logaddexp.reduce(x)
    return np.exp(x - a)


def Progress(snp,Inter,kmax,partition,window):
    '''
    simple function to print out remaining time estimates.
    '''
    for step in range(partition + 1):
        if snp == int(step*kmax / float(partition)) + 1 and Inter:
            job = "[" + step * "-" + (partition - step) * " " + "]" + " estimated: "
            MeanT = str(float(mean(Inter)) * (kmax-snp) / (window * 3600)).split(".")
            Hours = str(int(MeanT[0]))
            Minutes = str(int(60 *float("0." + MeanT[1][0])))
            print(job + Hours + " hours and " + Minutes + " min remaining.")
        elif snp == int(step*kmax / float(partition)) + 1 and len(Inter) == 0:
            print("[" + step * "-" + (partition - step) * " " + "]")



### reads SIM conf file, or any file with those lines really.

def readSIM_conf(Fconf):
    '''
    Reads SIM conf file, or any file with those lines really.
    Only the number of markers and organization of 
    hybrid and parent populations in geno file are read.
    '''
    Orders = recursively_default_dict()
    Fconf = open(Fconf,"r")
    for line in Fconf:
        line = line.split("=")
        line = [x.strip() for x in line]
        if line[0] == "n.markers":
            Miss = {x:[x+1,"Irrelevant"] for x in range(int(line[1]))}
            Orders["Miss"] = Miss
        if line[0] == "mod.n.samplePerPop":
            Pops = filter(lambda ch: ch not in " c()", line[1]).split(",")
            print(Pops)
            Pops = [int(x) for x in Pops]
            Pops = {x+1:[y for y in range(sum(Pops[:x]),sum(Pops[:x+1]))] for x in range(len(Pops))}
            Orders["Geneo"] = Pops
    Fconf.close()
    return Orders


##
##
##

def read_NEWorder(OrderFile,Row):
    """
    Read slection file setup w/ JC.
    """
    Gps = []
    Ofile = open(OrderFile,"r")
    Row = Row-1
    d = 0
    for line in Ofile:
        if d == 0:
            d += 1
            continue
        line = line.split()
        Gps.append(int(line[Row]))
    Ofile.close()
    Geneo = {x:[] for x in list(set(Gps))}
    Ofile = open(OrderFile,"r")
    d = 0
    for line in Ofile:
        if d == 0:
            d += 1
            continue
        line = line.split()
        Geneo[int(line[Row])].append(Fam[line[0]])
    Ofile.close()
    return Geneo


########
########
########
########
########

##
## Second stage tools. 
## readBlocks function reads condensed output from Kernel_MarkIII.py
## Filter labels condense could be of more use.


def Gen_rand(Snp_lib,n,L):
    """
    create library of random genome windows.
    args: snp/index library (MissG)
    n: number of windows
    L: size of windows
    """
    from random import randint
    Seen = {x:recursively_default_dict() for x in range(1,13)}
    for i in range(n):
        CHR = randint(1,12)
        snp = randint(0,len(Snp_lib[CHR]) - L - 1)
        Seen[CHR][Snp_lib[CHR][snp][0]] = [Snp_lib[CHR][snp + L][0],'rand']
    return Seen


##### One shot Classifier 

def Filter_label(Parents,data,Geneo):
    """
    Filter label one shot.
    """
    Bandwidth_split = 20
    quantile_filter_1 = 20
    quantile_filter_2 = 0
    
    KDE_samples = 1000
    
    """
    Needs:
    
    Parents
    Aro
    Geneo
    
    """
    
    Likes = {x:[] for x in range(len(Parents))}
    Accurate = []
    # use grid search cross-validation to optimize the bandwidth
    params = {'bandwidth': np.linspace(np.min(data), np.max(data),Bandwidth_split)}
    grid = GridSearchCV(KernelDensity(algorithm = "ball_tree",breadth_first = False), params,verbose=0)
    #grid.fit(data[Aro.shape[0]:,:])
    #Likes = {x:[] for x in range(len(Parents))}
    
    for D in range(len(Parents)):
        Where = [sum([len(Geneo[x]) for x in Parents[0:D]])+Aro.shape[0],sum([len(Geneo[x]) for x in Parents[0:(D+1)]])+Aro.shape[0]]
        Where = [int(x) for x in Where]
        
        grid.fit(data[Where[0]:Where[1],:])
        #print("best bandwidth: {0}".format(grid.best_estimator_.bandwidth))
        
        # use the best estimator to compute the kernel density estimate
        kde = grid.best_estimator_
        
        #### For sampling based
        #Quanted_set = kde.sample(500)
        #P_dist = np.exp(kde.score_samples(Quanted_set))
        #y_pred = [int(x < stats.scoreatpercentile(P_dist,quantile_filter_1)) for x in P_dist]
        
        #### Outlier based
        Quanted_set = data[Where[0]:Where[1],:]
        #P_dist = np.exp(kde.score_samples(Quanted_set))
        bandwidth = estimate_bandwidth(data[Where[0]:Where[1],:], quantile=0.5)
        if bandwidth==0:
            bandwidth = 0.1
        ms = MeanShift(bandwidth = bandwidth,bin_seeding=None, cluster_all=False, min_bin_freq=int(.15*len(Quanted_set)))
        ms.fit(Quanted_set)
        labels = ms.labels_
        y_pred = [int(x == -1) for x in labels]
        #y_pred = scipy.stats.norm(mean(P_dist),np.std(P_dist)).cdf(P_dist)
        #Below = [x + Aro.shape[0] + sum([len(Geneo[x]) for x in Parents[0:D]]) for x in range(len(y_pred)) if y_pred[x] <= .2]
        Below = [x for x in range(len(y_pred)) if y_pred[x] == 1]
        print(Below)
        Indexes = range(Quanted_set.shape[0])
        if len(Below) < len(y_pred) and Below:
            Indexes = [x for x in Indexes if x not in Below]
            grid.fit(Quanted_set[Indexes,:])
            kde = grid.best_estimator_
            Quanted_set = kde.sample(1000)
            P_dist = kde.score_samples(Quanted_set)
        
        if len(Below) == 0:
            grid.fit(Quanted_set[Indexes,:])
            kde = grid.best_estimator_
            Quanted_set = kde.sample(1000)
            P_dist = kde.score_samples(Quanted_set)
        
        Fist = kde.score_samples(data)
        CDF = [len([y for y in P_dist if y <= x]) / float(len(P_dist)) for x in Fist]
        Fist = scipy.stats.norm(np.mean(P_dist),np.std(P_dist)).cdf(Fist)
        Dist = np.exp(kde.score_samples(data))
        Dist = Dist / max(Dist)
        Dist = [round(x,5) for x in Dist]
        Accurate.append(Fist)
        Likes[D].append(Fist)
        
    Test = [Likes[x][0] for x in range(len(Parents))]
    Test = [[[Test[x][y],0][int(Test[x][y]<1e-3)] for y in range(len(Test[x]))] for x in range(len(Test))]
    Test = np.array(Test).T
    Clear = [x for x in range(Test.shape[0]) if max(Test[x,:]) < 1e-03]
    labels2 = np.argmax(Test,axis = 1)
    
    #Method = "MeanShift"
    #bandwidth = estimate_bandwidth(Test, quantile=0.05)
    ###
    #if Method == "MeanShift":
    #    #### Mean Shift approach
    #    ## from sklearn.cluster import MeanShift, estimate_bandwidth
    #    #bandwidth = estimate_bandwidth(data[Aro.shape[0]:,:], quantile=0.2, n_samples=Daddy.shape[0])
    #    if bandwidth==0:
    #        bandwidth = 0.1
    #    ms = MeanShift(bandwidth=bandwidth, bin_seeding=True, cluster_all=True, min_bin_freq=25)
    #    ms.fit(Test)
    #    labels2 = ms.labels_
    
    maxim = labels2
    Consex = [x for x in it.combinations(range(len(Parents)),2)]
    
    for h in range(len(maxim)):
        CL = []
        for j in Consex:
            Diff = Test[h,j]
            if maxim[h] not in j or len([x for x in Diff if x == 0]) > 0:
                continue
            if max(Diff) <= 5e-2:
                Diff = 0
            else:
                Diff = abs(max(Diff)) / abs(min(Diff))
                Diff = int(Diff > 1.2)
            #Diff = len([x for x in Diff if x >= X_threshold])
            
            #print(Diff)
            if Diff == 0:
                CL.append(j)
        
        if len(CL) == 2:
            maxim[h] = 7
        if len(CL) == 1:
            maxim[h] = sum(CL[0]) + len(Parents)
    
    labels2[Clear] = 3
    #labels2 = [[0,1][int(Test[x,2] > .05)] for x in range(Test.shape[0])]
    return Test,labels2



def Filter_label2(Parents,data,Geneo,Outlier_method,KDE_tool,normalize,Outlier_threshold):
    """
    Kernel one shot processing. Receives DR coordinates.
    Updated to Kernel_MarkIII_clean.py 19-10-2017.
    """
    Bandwidth_split = 30
    quantile_filter_1 = 20
    X_threshold = Outlier_threshold
    Diff_threshold = 2
    KDE_samples = 1000
    
    """
    Needs:
    
    Parents
    Aro
    Geneo
    
    """
    Aro = data[:len([y for y in it.chain(*[Geneo[x] for x in Geneo.keys() if x not in Parents])]),:]
    
    Likes = {x:[] for x in range(len(Parents))}
    Accurate = []
    params = {'bandwidth': np.linspace(np.min(data), np.max(data),Bandwidth_split)}
    grid = GridSearchCV(KernelDensity(algorithm = "ball_tree",breadth_first = False), params,verbose=0)
    
    for D in range(len(Parents)):
        Where = [sum([len(Geneo[x]) for x in Parents[0:D]])+Aro.shape[0],sum([len(Geneo[x]) for x in Parents[0:(D+1)]])+Aro.shape[0]]
        Where = [int(x) for x in Where]
        
        ### FILTER OUT OUTLIERS
        ### Conservative methods are preferred since onluy extreme values 
        ### can reasonably be removed given an unreliable sampling method.
        ### In other words, no assumptions there.
        ##### DBSCAN
        if Outlier_method == 'DBSCAN':
            Quanted_set = data[Where[0]:Where[1],:]
            db = DBSCAN(eps=0.3, min_samples=5).fit(Quanted_set)
            core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
            core_samples_mask[db.core_sample_indices_] = True
            labels = db.labels_
            y_pred = [int(x == -1) for x in labels]
            Below = [x for x in range(len(y_pred)) if y_pred[x] == 1]
        ##### KDE percentile based
        if Outlier_method == 'Perc':
            Quanted_set = data[Where[0]:Where[1],:]
            kde = stats.gaussian_kde(Quanted_set.T)
            Quanted_set = kde.sample(500)
            P_dist = np.exp(kde.score_samples(Quanted_set))
            y_pred = [int(x < stats.scoreatpercentile(P_dist,quantile_filter_1)) for x in P_dist]
            Below =  [x for x in range(len(y_pred)) if y_pred[x] <= 5]
        ##### MEAN SHIFT by-product
        if Outlier_method == 'MS':
            Quanted_set = data[Where[0]:Where[1],:]
            Set_ref = np.vstack({tuple(row) for row in Quanted_set})
            if len(Set_ref) == 1:
                Below = []
            else:
                bandwidth = estimate_bandwidth(Quanted_set, quantile=0.5, n_samples=len(Quanted_set))
                if bandwidth==0:
                    bandwidth = 0.1
                ms = MeanShift(bandwidth=bandwidth,bin_seeding=True, cluster_all=False, min_bin_freq=5)
                ms.fit(Quanted_set)
                labels = ms.labels_
                y_pred = [int(x == -1) for x in labels]
                Below = [x for x in range(len(y_pred)) if y_pred[x] == 1]
        if Outlier_method == 'Z':
        ##### Normalized kde outliers (simpler is better)
            Quanted_set = data[Where[0]:Where[1],:]
            Set_ref = np.vstack({tuple(row) for row in Quanted_set})
            #if len(Set_ref) < 3:
            if 1 > 0:
                grid.fit(Quanted_set)
                kde = grid.best_estimator_
                #Quanted_set = kde.sample(KDE_samples)
                P_dist = kde.score_samples(Quanted_set)
            else:
                kde = stats.gaussian_kde(Quanted_set.T)
                Quanted_set = kde.resample(KDE_samples).T
                P_dist = np.log(kde(Quanted_set.T))
            
            y_pred = scipy.stats.norm(mean(P_dist),np.std(P_dist)).cdf(P_dist)
            Below = [x for x in range(len(y_pred)) if y_pred[x] <= 5e-3]
        #### ISOLATION FOREST ###########################
        if Outlier_method == 'Isolation_forest':
            rng = np.random.RandomState(42)
            Quanted_set = data[Where[0]:Where[1],:]
            clf = IsolationForest(max_samples=len(Quanted_set), random_state=rng)
            clf.fit(Quanted_set)
            y_pred = clf.predict(Quanted_set)
            Below = [x for x in range(len(y_pred)) if y_pred[x] == -1]
        else: 
            Quanted_set = data[Where[0]:Where[1],:]
            Below = []
        ###################            ###################
        ## KDE estimation on filtered parental data set ##
        ###################            ###################
        #print Below
        Indexes = range(Quanted_set.shape[0])
        Indexes = [x for x in Indexes if x not in Below]
        Set_ref = np.vstack({tuple(row) for row in Quanted_set[Indexes,:]})
        if KDE_tool == 'sklearn':
            grid.fit(Quanted_set[Indexes,:])
            kde = grid.best_estimator_
            #Quanted_set = kde.sample(KDE_samples)
            P_dist = kde.score_samples(Quanted_set)
            Fist = kde.score_samples(data)
        if KDE_tool == 'scipy':
            Craft = Quanted_set[Indexes,:]
            ### Resampling from estimate distribution is
            ### preferred for cdf estimates.
            kde = stats.gaussian_kde(Craft.T)
            Craft = kde.resample(KDE_samples).T
            P_dist = np.log(kde(Craft.T))
            Fist = np.log(kde(data.T))
        if sum(np.isnan(P_dist)) == len(P_dist):
            if len(Likes[D]) == 0:
                Likes[D].append([int(x in range(Where[0],Where[1])) for x in range(data.shape[0])])
                Accurate.append([int(x in range(Where[0],Where[1])) for x in range(data.shape[0])])
            else:
                Likes[D].append(Likes[D][-1])
                Accurate.append(Likes[D][-1])
            continue
        ### Neutrality tests of filtered reference pop KDE derived log-Likelihood.
        #Normality.append(scipy.stats.mstats.normaltest(P_dist)[1])
        ######
        CDF = [len([y for y in P_dist if y <= x]) / float(len(P_dist)) for x in Fist]
        #Fist = (Fist - np.mean(P_dist)) / np.std(P_dist)
        if normalize == 'CDF':
            Fist = scipy.stats.norm(np.mean(P_dist),np.std(P_dist)).cdf(Fist)
        else:
            Fist = kde(data.T)
            Fist = Fist / max(Fist)
            Fist = [round(x,5) for x in Fist]
        Accurate.append(Fist)
        Likes[D].append(Fist)
       
    
    Test = [Likes[x][0] for x in range(len(Parents))]
    Test = [[[Test[x][y],0][int(Test[x][y]<=X_threshold)] for y in range(len(Test[x]))] for x in range(len(Test))]
    Test = np.array(Test).T
    Clear = [x for x in range(Test.shape[0]) if max(Test[x,:]) < X_threshold]
    labels2 = np.argmax(Test,axis = 1)
    
    #Method = "MeanShift"
    #bandwidth = estimate_bandwidth(Test, quantile=0.05)
    ###
    #if Method == "MeanShift":
    #    #### Mean Shift approach
    #    ## from sklearn.cluster import MeanShift, estimate_bandwidth
    #    #bandwidth = estimate_bandwidth(data[Aro.shape[0]:,:], quantile=0.2, n_samples=Daddy.shape[0])
    #    if bandwidth==0:
    #        bandwidth = 0.1
    #    ms = MeanShift(bandwidth=bandwidth, bin_seeding=True, cluster_all=True, min_bin_freq=25)
    #    ms.fit(Test)
    #    labels2 = ms.labels_
    
    maxim = labels2
    Consex = [x for x in it.combinations(range(len(Parents)),2)]
    
    for h in range(len(maxim)):
        CL = []
        for j in Consex:
            Diff = Test[h,j]
            if maxim[h] not in j or len([x for x in Diff if x < X_threshold]) > 0:
                continue
            else:
                #Diff = int(len([x for x in Diff if x > X_threshold]) == 1)
                Diff = abs(max(Diff)) / abs(min(Diff))
                Diff = int(Diff > Diff_threshold)
            #Diff = len([x for x in Diff if x >= X_threshold])
            
            #print(Diff)
            if Diff == 0:
                CL.append(j)
        
        if len(CL) == 2:
            maxim[h] = 7
        if len(CL) == 1:
            maxim[h] = sum(CL[0]) + len(Parents)
    
    labels2[Clear] = 3
    #labels2 = [[0,1][int(Test[x,2] > .05)] for x in range(Test.shape[0])]
    return Test,labels2



def readBlocks(BLfile):
    """
    reads Blocks_ALL.txt. Stores by CHR/stat/ind - all integers. Returns dictionary.
    """
    Blocks = {x:recursively_default_dict() for x in range(1,13)}
    Sizes = {x:recursively_default_dict() for x in range(1,13)}
    BLfile = open(BLfile,"r")
    for line in BLfile:
        line = line.split()
        if line[0] == "CHR":
            continue
        else:
            Sizes[int(line[0])][int(line[1])] = int(line[2])
            #Blocks[int(line[0])][int(line[1])] = {Whom[c-3]:OriOrder.index(line[c]) + 1 for c in range(3,len(line)-1)}
            Blocks[int(line[0])][int(line[1])] = [int(x) for x in line[3:]]
    BLfile.close()
    Blocks = {x:Blocks[x] for x in Blocks.keys() if Blocks[x]}
    Sizes = {x:Sizes[x] for x in Sizes.keys() if Sizes[x]}
    return Blocks,Sizes


