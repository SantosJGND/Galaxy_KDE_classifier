# -*- coding: utf-8 -*-

"""
Created on Sun Mar 05 13:01:08 2017

@author: jgarcia
"""


import itertools as it

import numpy as np
from numpy import nan,mean

import time

from sklearn.neighbors import KernelDensity
from sklearn.decomposition import PCA
from sklearn.model_selection import GridSearchCV
from sklearn.cluster import MeanShift, estimate_bandwidth

import sys
sys.setrecursionlimit(10000)


import scipy
from scipy import stats
from scipy.signal import savgol_filter

from random import randint
import collections

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




def read_refs(index_file,Fam_lib):
    '''
    ref file indexes individuals to population code.
    '''
    indxs = recursively_default_dict()
    Absent= []
    
    Input = open(index_file,'r')
    for line in Input:
        line = line.split()
        if line[1] not in Fam_lib.keys():
            Absent.append(line[1])
            continue
        indxs[int(line[0])][Fam_lib[line[1]]] = []
    
    Input.close()
    
    indxs = {gop:[x for x in indxs[gop].keys()] for gop in indxs.keys()}
    
    return indxs, [x for x in sorted(indxs.keys())], Absent



def read_focus(index_file):
    indxs = []
    
    Input = open(index_file,'r')
    for line in Input:
        line = line.split()
        indxs.append(line[0])
    
    Input.close()
    
    return indxs



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
        Nsnps[CHR][d] = [float(line[3]),line[4],line[5]]
        if float(line[3]) in Gindex[CHR].keys():
            Gindex[CHR][float(line[3]) + .5] = [d]
        else:
            Gindex[CHR][float(line[3])] = [d]
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


#################################################
######### Other #################################

#### Local sampling correct
#### The use of this function was deprecated following the study of the impact of sampling on PCA projections.
#### see: https://github.com/SantosJGND/Genetic-data-analysis/blob/master/8.%20Controlling%20for%20size.ipynb

def local_sampling_correct(data):
    '''
    This function uses the MeanShift algorithm to identify clusters in the data and sample equally from each.
    The sampled observations are used to create a new PCA transformation that is applied to the original data.
    '''
    ncomp= data.shape[1]
    pca = PCA(n_components=ncomp, whiten=False,svd_solver='randomized')
    features= pca.fit_transform(data)
    
    N= 100
    bandwidth = estimate_bandwidth(features, quantile=0.2)
    if bandwidth <= 1e-3:
        bandwidth = 0.1
    
    params = {'bandwidth': np.linspace(np.min(features), np.max(features),20)}
    grid = GridSearchCV(KernelDensity(algorithm = "ball_tree",breadth_first = False), params,verbose=0)
    
    ## perform MeanShift clustering.
    ms = MeanShift(bandwidth=bandwidth, bin_seeding=True, cluster_all=False, min_bin_freq=20)
    ms.fit(features)
    labels1 = ms.labels_
    label_select = {y:[x for x in range(len(labels1)) if labels1[x] == y] for y in sorted(list(set(labels1)))}
    
    ## Extract the KDE of each cluster identified by MS.
    Proxy_data= []
    
    for lab in label_select.keys():
        if len(label_select[lab]) < 3:
            continue
        Quanted_set= features[label_select[lab],:]
        
        grid.fit(Quanted_set)
        
        kde = grid.best_estimator_
        Extract= kde.sample(N)
        Return= pca.inverse_transform(Extract)
        Proxy_data.extend(Return)
    
    Proxy_data= np.array(Proxy_data)
    
    pca = PCA(n_components=ncomp, whiten=False,svd_solver='randomized').fit(Proxy_data)
    New_features= pca.transform(data)
    
    return New_features


def lognormalize(x):
    a = np.logaddexp.reduce(x)
    return np.exp(x - a)


def interpolate_nans(X):
    """Overwrite NaNs with column value interpolations."""
    col_means = np.array([int(np.nanmean(X[:,x])) for x in range(X.shape[1])])
    inds = np.where(np.isnan(X))
    X[inds] = np.take(col_means,inds[1])
    return X


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



