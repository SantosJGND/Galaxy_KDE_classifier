# -*- coding: utf-8 -*-
"""
Created on Wed Dec 06 18:53:36 2017

@author: jgarcia
"""

import collections
import itertools as it
import numpy as np
import pandas as pd
import scipy

from sklearn.neighbors import KernelDensity
from sklearn.decomposition import PCA
from sklearn.model_selection import GridSearchCV
from sklearn.cluster import estimate_bandwidth
from sklearn.cluster import MeanShift, estimate_bandwidth

from Kernel_tools import recursively_default_dict, read_refs, read_focus, BIMread, FAMread
from Galaxy_Ideogram_tools import chromosome_collections, plot_ideo2, compress_ideo, Merge_class
from Galaxy_summary_tools import read_books, read_3D_profiles

import time
import re

import os, sys
import argparse

###########################################
###########################################
import matplotlib
matplotlib.use('Agg')

from matplotlib import pyplot as plt
from matplotlib.collections import BrokenBarHCollection
import pandas as pd


parser = argparse.ArgumentParser()

parser.add_argument("books",type=str,metavar= 'N',nargs= '+',
                    help = "Reference files to read. Any number can be given.")

parser.add_argument("--focus",type= str,help = "reference accessions indexes in genofile.")

parser.add_argument("--id",type= str,help = "Name your analysis")

parser.add_argument("--target",type= str,help = "target population")

parser.add_argument("--Dr_var",type= str,default= 'target',help = "focus Dr on target ['target'], admix+target ['focus_inc'] or global variation ['all']")

parser.add_argument("--ref",type= str,help = "reference accessions indexes in genofile.")

parser.add_argument("--CHR",type= int,help = "chromosome to draw ideogram of.")

parser.add_argument("--random",type= int, default= 0, help = "if a number greater than 0 is provided, performs analysis on a random sample of targetted clusters of that size.")

parser.add_argument("--start",type= int,help = "where to begin, in markers. Only makes sense if --CHR is also used.")

parser.add_argument("--end",type= int,help = "where to end, in markers. Only makes sense if --CHR is also used.")

parser.add_argument("--fam",help = "accession name file. same order as in geno.")

parser.add_argument("--plot",action= "store_true",help = "if given reduces cluster points using KDE extraction per cluster identified.")

parser.add_argument("--info",type= str,help = "optional information file on accessions analysed. requires columns: ID and label, with IDs and reference populations.")

parser.add_argument("--app",action= "store_true",help = "if given writes dash application and accompanying files.")

parser.add_argument("--reduc",action= "store_true",help = "if given prints cluster stats.")

parser.add_argument("--liss_control",action= "store_true",help = "if given controls supervised classification used by one where data is smoothed. Used to ignore sporadic outliers.")

parser.add_argument("--coarse",action='store_false',help= 'to smooth or not to smooth.')

parser.add_argument("--bin",default = 5,type= int,help = "smoothing parameter, must be uneven [savgol filter]")

parser.add_argument("--sg_order",default = 3,type= int,help = "staviksy golay filter order")

parser.add_argument("--ms",type=float,default = .1,help = "cluster profile selection threshold")

parser.add_argument("--shared",type=float,default = .1,help = "Proportion of shared ancestry at given locus.")

parser.add_argument("--threshold",type = float,default = 2,help = "Intermediate classification threshold")

parser.add_argument("--outlier",type = float,default = 1e-3,help = "outlier threshold")

parser.add_argument("--chrom_height",type= float, default= 1, help= "height of ideograms")

parser.add_argument("--chrom_gap",type= float,default= 0,help= "gap between ideograms.")

parser.add_argument("--height",type= float, default= 10,help= "figure height, in inches.")

parser.add_argument("--width",type= float,default= 20,help= "figure width, in inches.")

parser.add_argument('--xticks',type= int,default= 1000000,help= 'xticks on final ideogram')


args = parser.parse_args()


Home= 'Analyses_' + args.id

########## Complementary files.

print('to begin reading from: ')
print(args.books)

Library= read_books(args.books)

print('library:')
print(Library.sort_values(by= ['Chr','start']))


Ref_profiles, Profiles, Names, Out= read_3D_profiles(Library)


######
###### deciding who you're going to obe looking at;
###### in the future, replace with input file only it should be easier.
Fam = FAMread(args.fam)


####
#### some parameters
####

Diff_threshold = args.threshold
X_threshold= args.outlier

refs_lib, Parents, Absent_refs = read_refs(args.ref,Fam)

if Absent_refs:
    print('The following IDs were not found in .fam file: ' + ', '.join([x for x in Absent_refs]))
    
    perc_missing= float(len(Absent_refs)) / (sum([len(x) for x in refs_lib.values()]) + len(Absent_refs))
    
    if perc_missing >= .7:
        print('{} % references missing from .fam file.'.format(round(perc_missing*100,2)))


if args.focus:
    Focus = read_focus(args.focus)
else:
    Focus = Names

compound_reference = [Names.index(x) for x in Focus]

focus_indexes= [x for x in range(len(Focus))]

Absent= [x for x in Focus if x not in Names]

if len(Absent) > 0:
    print('The following individuals were not found in the files provided: {}'.format(Absent))
    print('Analysis will proceed without them.')
    Focus= [x for x in Focus if x in Names]


Blocks, Npops = Merge_class(Ref_profiles,compound_reference,Out,args.threshold,args.bin,args.outlier,args.coarse,args.sg_order)

if args.liss_control:
    Blocks_liss, Npops = Merge_class(Ref_profiles,compound_reference,Out,args.threshold,args.bin,args.outlier,True,args.sg_order)

####
####

if args.CHR:
    if len([x for x in Blocks.keys() if x == args.CHR]) == 0:
        print('{} chromomes not found in input Blocks file.'.format([x for x in Blocks.keys() if x == args.CHR]))
    
    chromosomes= [args.CHR]
else:
    chromosomes= Blocks.keys()


if args.start:
    Blocks= {Chr:{x:Blocks[Chr][x] for x in Blocks[Chr].keys() if x >= args.start} for Chr in Blocks.keys()}
    if not args.CHR:
        print("start was selected with no CHR specification.")

if args.end:
    Blocks= {Chr:{x:Blocks[Chr][x] for x in Blocks[Chr].keys() if x <= args.end} for Chr in Blocks.keys()}
    if not args.CHR:
        print("end was selected with no CHR specification.")


####

threshold= args.shared
MS_threshold= args.ms


def codes(Parents):
    Store= recursively_default_dict()
    
    # Parent labels
    for x in range(len(Parents)):
        Store[x + 1]= [Parents[x]]
    
    ## outlier label
    Store[len(Parents) + 1]= [x for x in Parents]
    
    ## Combination labels:
    combine= [x for x in it.combinations(Parents,2)]
    for i in range(len(combine)):
        Store[len(Parents) + 2 + i]= [x for x in combine[i]]
    
    if len(Parents) >=3:
        ## three way label    
        Store[len(Store) + 1]= [x for x in Parents]
    
    return Store


code_back= codes(Parents)

if args.target:
    target= [int(x) for x in args.target.split(',')]
else:
    target= [x for x in code_back.keys()]

print('target: {0}'.format(target))

##### Selecting windows to focus on:

#### Chose_profiles: automatically chose clusters with at least one included 
#### focus accession of 'target' color. 
#### Cluster assignment is estimated as max pval above outlier threshold.

MS_inliers= recursively_default_dict()
Chose_profiles= recursively_default_dict()
Coordinates= []
Empty= []

for CHR in sorted(chromosomes):
    for bl in sorted(Blocks[CHR].keys()):
        if len([x for x in focus_indexes if Blocks[CHR][bl][x] in target]) / float(len(Focus)) < threshold:
            continue
        
        Bls= sorted(list(Profiles[CHR][bl].keys()))
        pVals= np.array([Profiles[CHR][bl][y] for y in Bls])
        
        max_vals= np.amax(pVals,axis= 0)
        max_indx= np.argmax(pVals,axis= 0)
        
        inlier= [x for x in compound_reference if max_vals[x] >= args.outlier and Blocks[CHR][bl][Focus.index(Names[x])] in target]
        
        if args.liss_control:
            inlier= [x for x in inlier if Blocks_liss[CHR][bl][Focus.index(Names[x])] in target]
        
        BL_select= list(set([max_indx[x] for x in inlier]))
        
        #print('clusters {} selected. {} %'.format(BL_select,len(BL_select)/float(len(Bls))))
        
        if not BL_select:
            Empty.append([CHR,bl])
            continue
        
        BL_select= {
            x: pVals[x] for x in BL_select
            }
        
        BLextract= list(BL_select.keys())
        
        Assignment= {
                Bls[b]: [Focus.index(Names[x]) for x in inlier if max_indx[x] == b] for b in BLextract
            }
        
        Chose_profiles[CHR][bl]= list(Assignment.keys())
        
        for bls in sorted(Assignment.keys()):
            Coordinates.append([CHR,bl,Out[CHR][bl],bls])
            MS_inliers[CHR][bl][bls]= Assignment[bls]



Coordinates= np.array(Coordinates)

Clover= [[[Profiles[CHR][bl][x] for x in Chose_profiles[CHR][bl]] for bl in sorted(Chose_profiles[CHR].keys())] for CHR in sorted(Chose_profiles.keys())]
Clover= [z for z in it.chain(*[y for y in it.chain(*Clover)])]
Clover= np.array(Clover)

print(Clover.shape)

if args.random > 0:
    Who= np.random.choice(range(Coordinates.shape[0]),args.random)
    
    Coordinates= Coordinates[Who,:]
    
    Clover= Clover[Who,:]
    

####
#### Pre processing and dimensionality reduction of matrix
#### of selected clusters.
####

from sklearn import preprocessing

Clover = np.nan_to_num(Clover)
preProc_Clover = Clover

print('Clover shape: ', Clover.shape)

Clover = preprocessing.scale(Clover,axis = 1)
#
print("Clover shape: ", Clover.shape)

reefer= [g for g in it.chain(*[refs_lib[y] for y in sorted(refs_lib.keys())])]

if not args.focus:
    mary= [Names.index(x) for x in Names if Fam[x] not in reefer]
else:
    mary= [Names.index(x) for x in Names if Fam[x] not in reefer and x in Focus]

reefer= [Names.index(Fam[x]) for x in reefer]

Subset= [x for x in it.chain(*[mary,reefer])]
Trend= np.repeat([0,1],[len(mary),len(reefer)])


## apply pca to reference accessions, transform the rest.

Dr_var= args.Dr_var
Dr_processes= ['target','focus_inc','all']

if Dr_var not in Dr_processes:
    print('Dr_process selected: {}, Dr_var processes available: {}'.format(Dr_var,Dr_processes))
    Dr_var= 'target'

print('focusing Dr on {}'.format(Dr_var))

if Dr_var== 'target':
    variation_focus= [Names.index(Fam[x]) for x in it.chain(*[refs_lib[z] for z in list(set(it.chain(*[code_back[y] for y in target])))])]

if Dr_var== 'focus_inc':
    variation_focus= [Names.index(x) for x in Focus]
    variation_focus.extend([Names.index(Fam[x]) for x in it.chain(*[refs_lib[z] for z in list(set(it.chain(*[code_back[y] for y in target])))])])

if Dr_var== 'all':
    variation_focus= Subset


### PCA
pca = PCA(n_components=10, whiten=False).fit(Clover[:,variation_focus].T)
X_se = pca.transform(Clover[:,Subset].T)
COMPS = pca.components_.T #*np.sqrt(pca.explained_variance_)


###############################################################################
########################### PAINTING SHIT!! ###################################
###############################################################################

## 
## CLUSTER EIGENVALUES
##

### Clustering on decomposition
    
#    bandwidth = estimate_bandwidth(COMPS, quantile=0.1)
#    if bandwidth==0:
#        bandwidth = 0.1
#    
#    ms = MeanShift(bandwidth=bandwidth, bin_seeding=False, cluster_all=True, min_bin_freq=1)
#    ms.fit(COMPS)
#    labels1 = ms.labels_
#    label_select = {y:[x for x in range(len(labels1)) if labels1[x] == y] for y in sorted(list(set(labels1)))}
   
### HDBSCAN
#    from sklearn.cluster import DBSCAN
#    db = DBSCAN(min_samples=35).fit(COMPS)
#    labels1= db.labels_
#    label_select = {y:[x for x in range(len(labels1)) if labels1[x] == y] for y in sorted(list(set(labels1))) if y != -1}

#    ### Ward

#    from sklearn.cluster import AgglomerativeClustering
#    clustering = AgglomerativeClustering(linkage='ward', n_clusters=10)
#    clustering.fit(COMPS)
#    labels1= clustering.labels_
#    label_select = {y:[x for x in range(len(labels1)) if labels1[x] == y] for y in sorted(list(set(labels1)))}
#    
#    #
#    ### K-means
from sklearn.cluster import KMeans
kmeans = KMeans(n_clusters=10, random_state=0).fit(COMPS[:,5])
labels1 = kmeans.labels_
label_select = {y:[x for x in range(len(labels1)) if labels1[x] == y] for y in sorted(list(set(labels1)))}



############################################################################
############################################################################
############ Bring out ideograms ##########################################

###
### DELETE Ref_profiles AND Profiles (trying to save space)
###
Ref_profiles= {}
Profiles= {}


###############################################################################
#### Average normalized likelihhod among clustered eigenvectors by haplotype #####
###############################################################################


Cameo = []

for cramp in sorted(label_select.keys()):
    Clamp = np.mean(preProc_Clover[label_select[cramp],:],axis = 0)
    Fry = [Clamp[x] for x in Subset]
    Cameo.append(Fry)

Cameo = np.array(Cameo).T


###########################################################################
### cosine of the clustered eigenvectors with haplotype coordinates ######## DEPRECATED
###########################################################################

#cos_threshold = .6
#
#
#from numpy import dot
#from numpy.linalg import norm
#
#SpaceY = recursively_default_dict()
#
#for g in label_select.keys():
#    Green = COMPS[label_select[g],:]
#    SpaceY[g] = [mean(dot(X_se[x],Green.T)/norm(Green,axis=1)/norm(X_se[x])) for x in range(X_se.shape[0])]
#
#Globe = np.array([SpaceY[x] for x in sorted(SpaceY.keys())]).T
#
#

######## Reducing the number of cluster profiles to print:


new_labs= labels1

if args.reduc:
    new_labs= []
    reduced_comp= []
    Size= 2000
    
    params = {'bandwidth': np.linspace(np.min(COMPS), np.max(COMPS),30)}
    grid = GridSearchCV(KernelDensity(algorithm = "ball_tree",breadth_first = False), params,verbose=0)
    
    for lab in label_select.keys():
        N_prop= round(len(label_select[lab]) * Size / float(sum([len(x) for x in label_select.values()])))
        
        if len(label_select[lab]) <= 3:
            reduced_comp.extend(COMPS[label_select[lab],:])
            new_labs.extend([lab]*N_prop)
            continue
        
        grid.fit(COMPS[label_select[lab],:])    
        
        kde = grid.best_estimator_
        
        new_data = kde.sample(N_prop, random_state=0)
        
        reduced_comp.extend(new_data)
        new_labs.extend([lab]*N_prop)
    
    
    COMPS= np.array(reduced_comp)



### print


CHR = [x for x in chromosomes][-1]
start= args.id

### Writting files for dash application;

if args.app:
    try:
        import cPickle as pickle
    except ImportError:  # python 3.x
        import pickle
    
    with open('app_files.p', 'rb') as fp:
        data = pickle.load(fp)
    
    elements= {
    'ID=': start,
    'Where=': CHR,
    'ref=': target[0],
    'pop_refs=': [x for x in it.chain(*[['unlabelled'],[str(x) for x in Parents]])]
    }
    
    for x in range(len(data['app.py'])):
        line= data['app.py'][x]
        for el in elements.keys():
            if re.search(el,line):
                line= line.split('=')
                if el == 'pop_refs=':
                    new_line= '{}= {}'.format(line[0],elements[el]) + '\n'
                else:
                    new_line= '{}= "{}"'.format(line[0],elements[el]) + '\n'
                data['app.py'][x]= new_line
                line= new_line
    
    for file_name in data:
        filename= ''.join([Home,'/',file_name])
        os.makedirs(os.path.dirname(filename), exist_ok=True)
        Output= open(filename,'w')
        for line in data[file_name]:
            Output.write(line)
        Output.close()
        
    if args.info:
        ordercore= list(open(args.info,'r'))
        filename= ''.join([Home,'/Order_core.txt'])
        os.makedirs(os.path.dirname(filename), exist_ok=True)
        Output= open(filename,'w')
        for line in ordercore:
            Output.write(line)
        Output.close()
    
    else:
        crafty= []
        ## set up an accession to group map.
        reference= recursively_default_dict()
        for gp in refs_lib.keys():
            for acc in refs_lib[gp]:
                reference[Fam[acc]]= Parents.index(gp) + 1
        for acc in Names:
            if acc not in reference.keys():
                reference[acc]= 0
                
        for race in Names:
            crafty.append([race,elements['pop_refs='][reference[race]],reference[race]])
        crafty= pd.DataFrame(crafty,columns= ['ID','label','code'])
        filename= ''.join([Home,'/Order_core.txt'])
        os.makedirs(os.path.dirname(filename), exist_ok=True)
        crafty.to_csv(filename,sep= '\t')




if args.plot == True:
    Ideo_home= Home + '/Ideos'
    
    Blancs= {aim:{Chr:{bl:[[-1,0][int(x in target)] for x in Blocks[Chr][bl]] for bl in Blocks[Chr].keys()} for Chr in Blocks.keys()} for aim in list(set(labels1))}
    
    target_block= {Chr:{bl:[[-1,1][int(x in target)] for x in Blocks[Chr][bl]] for bl in Blocks[Chr].keys()} for Chr in chromosomes}
    
    
    for Chr in chromosomes:
        plot_ideo2(target_block,[Chr],Focus,Out,'all',args.chrom_height,args.chrom_gap,args.height,args.width,Ideo_home,args.id,args.xticks)
        
    for n in range(len(Coordinates)):
        site= Coordinates[n]
        trigger= MS_inliers[site[0]][site[1]][site[3]]
        
        for v in trigger:
            Blancs[labels1[n]][site[0]][site[1]][v] = 1
    
    for aim in Blancs.keys():
        if len(Focus)== 1:
            plot_ideo2(Blancs[aim],Blocks.keys(),Focus,Out,aim,args.chrom_height,args.chrom_gap,args.height,args.width,Ideo_home,args.id,args.xticks)
        else:
            for Chr in chromosomes:
                plot_ideo2(Blancs[aim],[Chr],Focus,Out,aim,args.chrom_height,args.chrom_gap,args.height,args.width,Ideo_home,args.id,args.xticks)





filename= Home + "/ID_focus.txt"
os.makedirs(os.path.dirname(filename), exist_ok=True)

Output = open(filename,"w")

for ind in Focus:
    Output.write(ind)
    Output.write("\n")

Output.close()


filename= Home + "/Profile_" + str(target[0]) + "_CHR" + str(CHR) + "."+str(start) + ".txt"
os.makedirs(os.path.dirname(filename), exist_ok=True)

Output = open(filename,"w")

for title in sorted(label_select.keys()):
    Output.write("G" + str(title) + "\t")

Output.write("\n")

for Future in range(Cameo.shape[0]):
    for Machine in range(Cameo.shape[1]):
        Output.write(str(Cameo[Future,Machine]) + "\t")
    Output.write("\n")

Output.close()




filename= Home + "/Profile_coordinates_" + str(target[0]) + "_CHR" + str(CHR) + "."+str(start) + ".txt"
os.makedirs(os.path.dirname(filename), exist_ok=True)

Output = open(filename,"w")

Output.write('\t'.join(['chrom','start','end','cluster','members','label']))
Output.write('\n')

for axe in range(Coordinates.shape[0]):
    Output.write('\t'.join([str(x) for x in Coordinates[axe,:]]) + '\t')
    Output.write('.'.join([str(x) for x in MS_inliers[Coordinates[axe,0]][Coordinates[axe,1]][Coordinates[axe,3]]]) + '\t')
    
    Output.write(str(labels1[axe]) + '\n')

Output.close()



filename= Home + "/DIM_private_"+str(target[0])+"_request_CHR" + str(CHR) + "."+str(start)+".txt"
os.makedirs(os.path.dirname(filename), exist_ok=True)

Output = open(filename,"w")

for entry in range(X_se.shape[0]):
    Output.write(str(Trend[entry]) + "\t")
    
    Output.write(Names[Subset[entry]] + "\t")
    for scythe in range(X_se.shape[1]):
        Output.write(str(X_se[entry,scythe]) + "\t")
    Output.write("\n")

Output.close()



filename= Home + "/DIM_private_"+str(target[0])+"_comp_CHR" + str(CHR) + "."+str(start)+".txt"
os.makedirs(os.path.dirname(filename), exist_ok=True)

Output = open(filename,"w")
Output.write('0\t' + '\t'.join([str(x) for x in pca.explained_variance_ratio_]) + '\t\n')

for entry in range(COMPS.shape[0]):
    Output.write(str(new_labs[entry]) + "\t")
    for scythe in range(COMPS.shape[1]):
        Output.write(str(COMPS[entry,scythe]) + "\t")
    Output.write("\n")

Output.close()


print('Done.')