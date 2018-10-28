# -*- coding: utf-8 -*-
"""
Created on Sun Feb 18 13:26:51 2018

@author: jgarcia
"""

from Kernel_tools import *
########## START HERE #############

import collections
import time
import itertools as it
import numpy as np
import re

import os
import argparse

import matplotlib
matplotlib.use('Agg')

from matplotlib import pyplot as plt
from sklearn.neighbors import KernelDensity
from sklearn.cluster import MeanShift, estimate_bandwidth

import Galaxy_summary_tools

parser = argparse.ArgumentParser()

### optional arguments
parser.add_argument("books",type=str,metavar= 'N',nargs= '+',
                    help = "Reference files to read. Any number can be given.")

parser.add_argument("--focus",type= str,help = "IDs of accessions to plot.")

parser.add_argument("--geneMerge",type= str,help = "Merged gene by block file.")

parser.add_argument("--id",type= str,help = "ID for this particular plot. if not given, name of last accession in focus.")

parser.add_argument("--target",type= str,default= '0,1,2',help = "target pops or combinations over which to sum. comma sep., will split digits between sep. for combinations.")

parser.add_argument("--outlier",type=float,default = 1e-3,help = "Outlier threshold")

parser.add_argument("--bornes",type= int,default= 1000,help = "bornes")

parser.add_argument("--individual",action='store_true',help= 'plot graphs by individual.')

parser.add_argument("--out",type= str,default= '',help = "output directory")

args = parser.parse_args()



Home= args.out

if len(Home) > 0:
    Home= Home + '/'

#################################################################################
##############################################################
####
#### Sum p-value proportions at every for every reference sample at every window.
####

def Target_pVals(Ref_profiles,focus_indicies,targets,Out,X_threshold):
    Blocks_genome = recursively_default_dict()
    print(targets)
    for CHR in Ref_profiles.keys():
        
        print(CHR)
        Points = sorted(Ref_profiles[CHR].keys())
        Likes = Ref_profiles[CHR]
        N_pops= len(Likes[[x for x in Likes.keys()][0]])
        
        print("number of reference populations: {0}".format(N_pops))
        Likes = {x:[Likes[bl][x] for bl in sorted(Likes.keys())] for x in range(N_pops)}
        Likes = {x:np.array(Likes[x]) for x in Likes.keys()}
        
        Topo = {ref:[] for ref in targets}
        
        #range_Parents = [x + Aro.shape[0] for x in range(Daddy.shape[0])]
        #range_Crossed = [x for x in range(Aro.shape[0])]
        
        
        
        for acc in focus_indicies:
            Guys = np.array([Likes[x][:,acc] for x in range(N_pops)])
            Guys = np.nan_to_num(Guys)
            Guys = [[[y,0][int(y<=X_threshold)] for y in x] for x in Guys]
            
            
            ##### Now grab the groups and lines that interest you. Stored in targets.
            for ref in targets:
                if len(ref) > 1:
                    background= [0]*len(Guys[0])
                    for ley in it.combinations(ref,2):
                        Test = [int(x <= X_threshold) for x in np.amax(np.array([Guys[x] for x in ley]),axis = 0)]
                        Inlier_indicies= [x for x in range(len(Test)) if Test[x] == 0]
                        Outlier_indicies= [x for x in range(len(Test)) if Test[x] == 1]
                        
                        
                        
                        combs= [1 - abs(Guys[ley[0]][x] - Guys[ley[1]][x]) / float(max([Guys[ley[0]][x],Guys[ley[1]][x]])) for x in Inlier_indicies]
                        for gimp in range(len(combs)):
                            background[Inlier_indicies[gimp]] +=  combs[gimp]
                    
                    combo= background
                    
                else:
                    combo= Guys[ref[0]]
                
                Topo[ref].append(combo)
            
            #####
        
        
        Topo = {ref:np.array(Topo[ref]) for ref in Topo.keys()}
        
        """
        if method == 'median':
            Topo= {ref:np.median(Topo[ref],axis= 0).reshape(1,-1) for ref in Topo.keys()}
        
        if args.coarse:
            Tail = {ref:savgol_filter(Topo[ref],BIN,args.sg_order,mode = "nearest") for ref in Topo.keys()}
        """
        
        Clove = {CHR:{Points[x]:{ref: Topo[ref][:,x] for ref in targets} for x in range(len(Points))}}
        
        Blocks_genome.update(Clove)
    
    return Blocks_genome



####
#### read block and profile files
####

print('To begin reading from: ')
print(args.books)

Ref_profiles, Names, Out = read_3D_profiles_list(args.books)


#### Read Block-gene merged file.
import pandas as pd

Gene_fit = pd.read_csv(args.geneMerge,sep= '\t')

####
#### some parameters
####

X_threshold= args.outlier


target= args.target
target= target.split(',')

for ref in range(len(target)):
    l = target[ref]
    if len(l) > 1:
        l= [int(x) for x in l]
        l= tuple(l)
        target[ref] = l
    else: 
        target[ref] = tuple([int(l)])


######
###### deciding who you're going to be looking at;
###### 

if args.focus:
    refs_lib, Parents = read_refs(args.focus)

    Focus= [z for z in it.chain(*[refs_lib[x] for x in refs_lib.keys()])]
    
    ### group to individual index
    norm_labels= {z:[x for x in range(len(Focus)) if Focus[x] in refs_lib[z]] for z in refs_lib.keys()}
    
else:
    Focus = Names


Absent= [x for x in Focus if x not in Names]

if len(Absent) > 0:
    print('The following individuals were not found in the files provided: {}'.format(Absent))
    print('Analysis will proceed without them.')
    Focus= [x for x in Focus if x in Names]


focus_indexes = [Names.index(x) for x in Focus]


### add-individual p-values in the desired combination, 
Merged_Pvalues= Target_pVals(Ref_profiles,focus_indexes,target,Out,args.outlier)



print("Number chromosomes selected: {0}".format(len(Merged_Pvalues)))

####
####

chromosomes= Merged_Pvalues.keys()

figsize= (18,7)

color_ref= ['red','yellow','blue','black','orange','purple','green','silver','red3','deepskyeblue','navy','chartreuse','darkorchid3','goldenrod2']

color_names= ['Indica','Aus','Japonica','outlier','Aus-Ind','Ind-Jap','Aus-Jap','Ind-Aus-Jap']

Genes_to_state= recursively_default_dict()

for Chr in chromosomes:
    
    CHR_genes= Gene_fit[(Gene_fit.chr == Chr)]
    
    for aim in target:
        
        ### First plot gene density.
        Genes= []
        for gen in CHR_genes.ID:
            
            lower_bound= [int(x) for x in CHR_genes[(CHR_genes.ID == gen)].start][0] - args.bornes
            upper_bound= [int(x) for x in CHR_genes[(CHR_genes.ID == gen)].end][0] + args.bornes
            bl= np.array([x for x in Merged_Pvalues[Chr].keys() if x <= upper_bound and Out[Chr][x] >= lower_bound])
            
            stat= []
            
            ## Get average estimate by reference gp
            for gp in norm_labels.keys():
                
                ## Get individual averages for windows overlapping with gen
                SummedP= [np.mean([Merged_Pvalues[Chr][z][aim][ind] for z in bl]) for ind in norm_labels[gp]]
                
                SummedP= np.mean(SummedP)
                stat.append(SummedP)
            
            ## average estimates across groups
            stat= np.mean(stat)
            
            Genes.append(stat)
            Genes_to_state[gen][aim]= stat
        
        X_plot= np.linspace(0,max(Genes),200)[:,np.newaxis]
        bandwidth = estimate_bandwidth(np.array(Genes).reshape(-1,1), quantile=0.15, n_samples=len(Genes))
        
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111)
        
        kde= KernelDensity(kernel= 'gaussian',bandwidth= bandwidth).fit(np.array(Genes).reshape(-1,1))
        log_dens= kde.score_samples(X_plot)
        ax.plot(X_plot[:,0],np.exp(log_dens),'-',label= 'label= {0}, band= {1}'.format(aim,bandwidth))
        ax.legend(loc='upper right')
        
        filename= Home + 'GeneValues_labs' + ''.join([str(x) for x in aim])+ '_CHR' + str(Chr).zfill(2) +'.png'
        os.makedirs(os.path.dirname(filename), exist_ok=True)
        plt.savefig(filename,bbox_inches = 'tight')
        
        #### Now plot window values
        
        Xplot= [x for x in sorted(Merged_Pvalues[Chr].keys())]
        X= [np.median(Merged_Pvalues[Chr][bl][aim]) for bl in sorted(Merged_Pvalues[Chr].keys())]
           
        
        #
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111)
        
        ax.plot(Xplot,X,'.',label= 'label= {0}, band= {1}'.format(aim,bandwidth))
        ax.legend(loc='upper left')
        
        
        ax.set_xticks([x for x in range(0,max(Merged_Pvalues[Chr].keys()),int(1e6))])
        plt.xticks(fontsize = 5,rotation = 90)
        ax.tick_params(axis = 'x',pad = 15)
        
        filename= Home + 'Compare_labs' + ''.join([str(x) for x in aim])+ '_CHR' + str(Chr).zfill(2) +'.png'
        os.makedirs(os.path.dirname(filename), exist_ok=True)
        plt.savefig(filename,bbox_inches = 'tight')
    
    ###############################################################
    #### combine
    ### First plot gene density.
    Genes= []
    for gen in CHR_genes.ID:
        lower_bound= [int(x) for x in CHR_genes[(CHR_genes.ID == gen)].start][0] - args.bornes
        upper_bound= [int(x) for x in CHR_genes[(CHR_genes.ID == gen)].end][0] + args.bornes
        bl= np.array([x for x in Merged_Pvalues[Chr].keys() if x <= upper_bound and Out[Chr][x] >= lower_bound])
        
        stat_total= 0
        
        for aim in target:
            stat= []
            ## Get average estimate by reference gp
            for gp in norm_labels.keys():
                
                SummedP= [np.mean([Merged_Pvalues[Chr][z][aim][ind] for z in bl]) for ind in norm_labels[gp]]
                
                SummedP= np.mean(SummedP)
                stat.append(SummedP)
            
            ## average estimates across groups
            stat= np.mean(stat)
            stat_total += stat
            
        Genes.append(stat)
        Genes_to_state[gen]['compound']= stat
    
    X_plot= np.linspace(0,max(Genes),200)[:,np.newaxis]
    bandwidth = estimate_bandwidth(np.array(Genes).reshape(-1,1), quantile=0.15, n_samples=len(Genes))
    
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)
    
    kde= KernelDensity(kernel= 'gaussian',bandwidth= bandwidth).fit(np.array(Genes).reshape(-1,1))
    log_dens= kde.score_samples(X_plot)
    ax.plot(X_plot[:,0],np.exp(log_dens),'-',label= 'label= {0}, band= {1}'.format(aim,bandwidth))
    ax.legend(loc='upper right')
    
    filename= Home + 'GeneValues_compound_' + args.target + 'CHR' + str(Chr).zfill(2)+'.png'
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    plt.savefig(filename,bbox_inches = 'tight')
    
    #####
    
    Xplot= [x for x in sorted(Merged_Pvalues[Chr].keys())]
    X= [sum([np.median(Merged_Pvalues[Chr][bl][aim]) for aim in target]) for bl in sorted(Merged_Pvalues[Chr].keys())]
    
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)
    
    ax.plot(Xplot,X,'.',label= 'label= {0}, band= {1}'.format(aim,bandwidth))
    ax.legend(loc='upper left')
    
    
    ax.set_xticks([x for x in range(0,max(Merged_Pvalues[Chr].keys()),int(1e6))])
    plt.xticks(fontsize = 5,rotation = 90)
    ax.tick_params(axis = 'x',pad = 15)
    
    filename= Home + 'Compound_labs' + args.target + '_CHR' + str(Chr).zfill(2) +'.png'
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    plt.savefig(filename,bbox_inches = 'tight')



if args.individual:
    
    for Chr in chromosomes:
        
        for aim in target:
            
            for inddy in range(len(focus_indexes)):
                Xplot= [x for x in sorted(Merged_Pvalues[Chr].keys())]
                X= [Merged_Pvalues[Chr][bl][aim][inddy] for bl in sorted(Merged_Pvalues[Chr].keys())]
                   
                
                #
                fig = plt.figure(figsize=figsize)
                ax = fig.add_subplot(111)
                
                ax.plot(Xplot,X,'.',label= 'label= {0}, band= {1}'.format(aim,bandwidth))
                ax.legend(loc='upper left')
                
                
                ax.set_xticks([x for x in range(0,max(Merged_Pvalues[Chr].keys()),int(1e6))])
                plt.xticks(fontsize = 5,rotation = 90)
                ax.tick_params(axis = 'x',pad = 15)
                
                filename= Home + Names[focus_indexes[inddy]]+'_labs' + ''.join([str(x) for x in aim])+ '_CHR' + str(Chr).zfill(2) +'.png'
                os.makedirs(os.path.dirname(filename), exist_ok=True)
                plt.savefig(filename,bbox_inches = 'tight')



##############################################################################
################ Output

filename= Home + 'Genes_labs' + args.target + '_cp.txt'
os.makedirs(os.path.dirname(filename), exist_ok=True)
Output= open(filename,'w')

Output.write('ID\t')
Output.write('\t'.join([''.join([str(x) for x in z]) for z in target]))
Output.write('\tcompound')
Output.write('\n')


for gen in Genes_to_state.keys():
    Output.write(gen + '\t')
    Output.write('\t'.join([str(round(Genes_to_state[gen][x],4)) for x in target]) + '\t')
    Output.write(str(Genes_to_state[gen]['compound']))
    Output.write('\n')

Output.close()

print('Done.')
