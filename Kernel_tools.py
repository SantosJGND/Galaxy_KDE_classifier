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

import sys
sys.setrecursionlimit(10000)


import scipy
from scipy import stats
from scipy.signal import savgol_filter

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



