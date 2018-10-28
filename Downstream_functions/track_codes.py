# -*- coding: utf-8 -*-
"""
Created on Wed Apr 04 13:57:25 2018

@author: jgarcia
"""

from Kernel_tools import *
########## START HERE #############

import os
import itertools as it

import pandas as pd

import argparse
parser = argparse.ArgumentParser()

parser.add_argument("--ref",help = "reference accessions indexes in genofile.")

parser.add_argument('--cols',help= 'provide reference colors. optional')

args = parser.parse_args()



########## Complementary files.

def read_refs(index_file):
    indxs = recursively_default_dict()
    
    Input = open(index_file,'r')
    for line in Input:
        line = line.split()
        indxs[int(line[0])][int(line[1])] = []
    
    Input.close()
    
    indxs = {gop:[x for x in indxs[gop].keys()] for gop in indxs.keys()}
    
    return indxs, [x for x in sorted(indxs.keys())]


####

refs_lib, Parents = read_refs(args.ref)

color_list= ['red','yellow','blue','black','orange','purple','green','silver','red3','deepskyeblue','navy','chartreuse','darkorchid3','goldenrod2']



def codes(Parents):
    Store= []
    
    # Parent labels
    for x in range(len(Parents)):
        Store.append([x + 1,str(Parents[x]),color_list[x]])
    
    ## outlier label
    Store.append([len(Parents) + 1,','.join([str(x) for x in Parents]),color_list[len(Parents)]])
    
    ## Combination labels:
    combine= [x for x in it.combinations(Parents,2)]
    for i in range(len(combine)):
        Store.append([len(Parents) + 2 + i,','.join([str(x) for x in combine[i]]),color_list[len(Parents) + 1 + i]])
    
    if len(Parents) >=3:
        ## three way label    
        Store.append([len(Store) + 1,','.join([str(x) for x in Parents]),color_list[len(Store)]])
        
    Store= pd.DataFrame(Store,columns= ['code','target_labels','color'])
    return Store


cypher= codes(Parents)

cypher.to_csv('codes.txt',sep= '\t')

print('written color-code map to: codes.txt')
