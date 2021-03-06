# -*- coding: utf-8 -*-
"""
Created on Sun Mar 05 13:01:08 2017

@author: jgarcia
"""

from Kernel_tools import *
########## START HERE #############

import os
import argparse
parser = argparse.ArgumentParser()

parser.add_argument("--CHR",type=int,help = "Chromosome number")

parser.add_argument("--geno",help = "Chromosome file to work on")

parser.add_argument("--fam",help = "accession name file. same order as in geno.")

parser.add_argument("--bim",help = "snp information bim format.")

parser.add_argument("--ref",help = "reference accessions indexes in genofile.")

parser.add_argument("--admx",help = "admixed accession indexes in geno file")

parser.add_argument("--proc",help = "number of processors requested")

parser.add_argument("--out",type= str,default= '',help = "output directory")

parser.add_argument("--MSprint",action= "store_true",help = "if given prints cluster stats.")
###
parser.add_argument("--PCvar",type= float,default= 0,help = "if above 0 retain dimensions as a function of PC variance")
###
parser.add_argument("--VARprint",action= "store_true",help = "if given prints PC explained variance per window. If PCA is not chosen just prints out 0's")
###
parser.add_argument("--id",type= str,default= '1',help = "Give your analysis an ID. default is set to integer 2")
###
parser.add_argument("-w",type = int,default = 200, help = "Window size - markers")
### 
parser.add_argument("--clustmethod",default = "MeanShift",choices = ["MeanShift","DBscan","HDBscan"],help = "Clustering method to extract reference specific clusters. MS, dbscan and hdbsan available. MS preferred")
###
parser.add_argument("--het",type = float,default = 5e-2, help = "Heterozygosity filter")
### 
parser.add_argument("--dr",default = 'NMF',help = "Dimensionality reduction. options: PCA, NMF")
###
parser.add_argument("--mono",action= "store_false", help= "if given includes monomorphic SNPs")
### 
parser.add_argument("--ncomp",type = int,default = 5,help = "Number of components kept in case of PCA reduction")
### 
parser.add_argument("--outmethod",default = "None",help = "Outlier filter of population refs method. options: DBSCAN, NMF, Perc, Z, Isolation_forest.")
###
parser.add_argument("--overlap",type= int,default = 100,help = "Overlap between windows, in snps")
###
parser.add_argument("--miss_cd",type= int,default = '9',help = "Missing data code")
###
parser.add_argument("--het_cd",type= int,default = '1',help = "Het allele code")
###
parser.add_argument("--hom_cd",type= int,default = '2',help = "Alt allele code")
###

args = parser.parse_args()


########## Complementary files.
#### read names (Fam file in same order as genotype data)
Fam = FAMread(args.fam)

#### read SNP info (.bim file, same order and number as in geno file)
MissG, Gindex = BIMread(args.bim)

GenoSUF = args.geno

#### read admix ind files.
admx_lib, Crossed, absent_admx = read_refs(args.admx,Fam)

if absent_admx:
    print(",".join([str(x) for x in absent_admx]) + ' absent from provided admx.')

if len(absent_admx) > (0.5 * len(Crossed)):
    print('over half the admx missing. Aborting run.')
    

#### read reference ind files.
refs_lib, Parents, absent_refs  = read_refs(args.ref,Fam)

if absent_refs:
    print(",".join([str(x) for x in absent_refs]) + ' absent from provided refs.')

if len(absent_admx) > (0.5 * len(Crossed)):
    print('over half the references missing. Aborting run.')
    


admx_lib.update(refs_lib)

Geneo = admx_lib

CHR = args.CHR

Home = args.out

if len(args.out) > 0:
    Home= args.out + '/'

print('Number of markers: {}'.format(len(MissG[CHR])))

print('Population labels: {}'.format(Geneo.keys()))

print('Population Sizes: {}'.format([len(x) for x in Geneo.values()]))
####
####


def Set_up(Chr0,Chr1,Sub,MissG):
    Settings= []
    for Chr in range(Chr0,Chr1):
        Size= len(MissG[Chr])
        Steps = [x for x in range(0,Size+1,int(Size/float(Sub)))]
        Steps[-1] = Size
        for i in range(len(Steps)-1):
            Settings.append([Chr,Steps[i],Steps[i+1]-1])
    return Settings




def Main_engine(Fam,MissG,Geneo,Parents,GenoSUF,CHR,start,end,args):
    '''
    Main Engine. This function performs calculations along a slidding window.
    
    Parameters used are provided by the user.
    '''
    
    #### Fam
    #### Whose
    #### Miss
    #### Geneo
    #### Parents
    #### GenoFile
    Miss= MissG[CHR]
    GenoFile= GenoSUF 
    ######################
    ######################
    
    ### Here define some things
    Window = args.w
    
    #### Het. filter
    
    Filter_Het = args.het
    
    
    ## Dimensionality reduction: PCA, NMF
    DIMr = args.dr
    n_comp = args.ncomp
    PC_var_lim = .01
    ##
    
    ## Chose clustering method
    Method = args.clustmethod
    
    ## Outlier filter method: DBSCAN, NMF, Perc, Z, Isolation_forest.
    Outlier_method = args.outmethod
    
    quantile_filter_1 = 20 ### quantile filter for Percentile based outlier.
    
    ## KDE estimation tool
    ## Attention: sklearn uses GridsearchCV. While this implies more parameters, it 
    ## nonetheless ensures the construction of a KDE when scipy KDE would break down.
    ## This is because scipy KDE will always estimate bandwidth itself, thus using the points given only,
    ## while sklearn KDE allows the bandwidth to be passed. This allows us in turn to set a minimum,
    ## and ensure it to never come to 0. 
    ##
    ## A population genetics note: At the local genomic level the bandwidth is our proxy of genetic differentiation.
    ## It is thus important to for us to control from which material is used to estimate it.
    
    ## sklearn, scipy
    
    KDE_tool = 'sklearn'
    Bandwidth_split = 30
    
    ### Output normalization: CDF, Max
    normalize = 'CDF'
    
    ### Resampling size for CDF estimates.
    KDE_samples = 1000
    
    ## Controlling for global likelihood.
    Control = False
    
    ##
    ##
    ### These things define themselves
    
    Crossed = [x for x in Geneo.keys() if x not in Parents]
        
    Whose = list(it.chain(*Geneo.values()))
    SequenceStore = {fy:[] for fy in Whose}
    
    Likes = {x:[] for x in range(len(Parents))}
    
    Accuracy = []
    
    Geno = open(GenoFile,"r")
    Points = []
    Points_out = []
    
    PC_var= []
    
    Win = 0
    Index = 0
    Intervals = []
    
    Construct = recursively_default_dict()
    
    Codes = [0,0,0,0,0,0,0,0,0,0]
    
    Codes[int(args.miss_cd)]= nan
    Codes[int(args.het_cd)]= nan
    Codes[int(args.hom_cd)]= int(args.hom_cd)
    
    for line in Geno:
        
        d = Miss[Index][0]
        
        step= 1
        
        if Index > end:
            break
        if len([x for x in Whose if line[x] == args.het_cd]) / float(len(Whose)) > Filter_Het:
            Index += 1
            continue
        if len(list(set([x for x in Whose if line[x] not in [args.het_cd,args.miss_cd]]))) == int(args.mono):
            Index += 1
            continue
        if Index >= start and Index <= end:
            for judas in SequenceStore.keys():
                SequenceStore[judas].append(Codes[int(line[judas])])
            Win += step
            if Win == Window:
                s1 = time.time()
                window_start = Index - Window + 1
                
                Seq = SequenceStore
                
                ### Nan to Numeric. 
                ### Important to pay attention to this when applying this method to other data sets.
                ### Allowed here for the homozygous nature of rice genomes and the high quality of the data.
                ### think of replacing by by more performant imputation procedures.

                Aro = np.array([Seq[x] for x in it.chain(*[Geneo[x] for x in Crossed])])
                Daddy = np.array([Seq[x] for x in it.chain(*[Geneo[x] for x in Parents])])
                b = np.concatenate((Aro,Daddy),axis = 0)

                ## missing / Na index 
                na_here= int(args.miss_cd)
                ind_thresh= .1

                nstat= np.array(b)
                nstat= nstat == na_here
                nstat= np.array(nstat,dtype= int)
                nstat= np.sum(nstat,axis= 1)
                nstat= nstat / b.shape[1]

                nstat= ind_thresh >= ind_thresh
                nstat= np.array(nstat,dtype= int)
                nidx= [x for x in range(len(nstat)) if nstat[x] == 1]
                ## nan to num.
                Aro= np.nan_to_num(Aro)
                Daddy= np.nan_to_num(Daddy)
                b= np.nan_to_num(b)

                ## add a 0 row. 
                data = np.zeros((b.shape[0],b.shape[1]+1))
                data[:,:-1] = b
                
                if DIMr == 'PCA':
                    pca = PCA(n_components= n_comp, whiten=False,svd_solver='randomized').fit(data)
                    local_pcvar= list(pca.explained_variance_ratio_)
                    data = pca.transform(data)
                    
                    PC_var.append(local_pcvar)
                
                if DIMr == 'NMF':
                    from sklearn.decomposition import NMF
                    data = NMF(n_components=n_comp, init='random', random_state=0).fit_transform(data)
                
                Accurate = []
                params = {'bandwidth': np.linspace(np.min(data), np.max(data),Bandwidth_split)}
                grid = GridSearchCV(KernelDensity(algorithm = "ball_tree",breadth_first = False), params,verbose=0)
                
                ######################################
                ####### TEST global Likelihood #######
                ######################################
                Focus_labels = range(data.shape[0])
                
                if Method == "MeanShift":
                    #### Mean Shift approach
                    ## from sklearn.cluster import MeanShift, estimate_bandwidth
                    
                    bandwidth = estimate_bandwidth(data[Focus_labels,:], quantile=0.2, n_samples=len(Focus_labels))
                    if bandwidth <= 1e-3:
                        bandwidth = 0.1
                    
                    ms = MeanShift(bandwidth=bandwidth, cluster_all=False, min_bin_freq=25)
                    ms.fit(data[Focus_labels,:])
                    labels = ms.labels_
                    
                    Tree = {x:[Focus_labels[y] for y in range(len(labels)) if labels[y] == x] for x in [g for g in list(set(labels)) if g != -1]}
                
                SpaceX = {x:data[Tree[x],:] for x in Tree.keys()}
                
                
                
                for hill in SpaceX.keys():
                    if len(Tree[hill]) <= 3:
                        continue
                    grid.fit(data[Tree[hill],:])
                    
                    # use the best estimator to compute the kernel density estimate
                    kde = grid.best_estimator_
                    
                    P_dist = kde.score_samples(data[Tree[hill],:])
                    Dist = kde.score_samples(data)
                    P_dist= np.nan_to_num(P_dist)
                    Dist= np.nan_to_num(Dist)
                    if np.std(P_dist) == 0:
                        Dist= np.array([int(Dist[x] in P_dist) for x in range(len(Dist))])
                    else:
                        Dist = scipy.stats.norm(mean(P_dist),np.std(P_dist)).cdf(Dist)
                    Dist= np.nan_to_num(Dist)
                    Construct[Miss[window_start][0]][hill] = Dist
                    
                
                
                    ######################################### 
                ############# Reference KDE treatment ################
                    #########################################
                
                
                for D in range(len(Parents)):
                    Where = [sum([len(Geneo[x]) for x in Parents[0:D]])+Aro.shape[0],sum([len(Geneo[x]) for x in Parents[0:(D+1)]])+Aro.shape[0]]
                    ## don't use accessions that don't pass the NA filter. 
                    Where = [int(x) for x in Where if x not in nidx]
                    ### FILTER OUT OUTLIERS
                    ### Conservative methods are preferred since only extreme values 
                    ### can reasonably be removed given an unreliable sampling method.
                    ### In other words, no assumptions there.
                    ##### DBSCAN
                    if Outlier_method == 'DBSCAN':
                        from sklearn.cluster import DBSCAN
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
                    ##### Normalized kde outliers (simpler is better)
                    if Outlier_method == 'Z':
                        Quanted_set = data[Where[0]:Where[1],:]
                        Set_ref = np.vstack({tuple(row) for row in Quanted_set})
                        
                        if KDE_tool == 'sklearn':
                            grid.fit(Quanted_set)
                            kde = grid.best_estimator_
                            
                            P_dist = kde.score_samples(Quanted_set)
                        else:
                            kde = stats.gaussian_kde(Quanted_set.T)
                            Quanted_set = kde.resample(KDE_samples).T
                            P_dist = np.log(kde(Quanted_set.T))
                        
                        y_pred = scipy.stats.norm(mean(P_dist),np.std(P_dist)).cdf(P_dist)
                        Below = [x for x in range(len(y_pred)) if y_pred[x] <= 5e-3]
                    
                    #### ISOLATION FOREST ###########################
                    if Outlier_method == 'Isolation_forest':
                        from sklearn.ensemble import IsolationForest
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
                    
                    Indexes = range(Quanted_set.shape[0])
                    
                    Indexes = [x for x in Indexes if x not in Below]
                    
                    Set_ref = np.vstack({tuple(row) for row in Quanted_set[Indexes,:]})
                    
                    if KDE_tool == 'sklearn':
                        grid.fit(Quanted_set[Indexes,:])
                        kde = grid.best_estimator_
                        P_dist = kde.score_samples(Quanted_set)
                        Fist = kde.score_samples(data)
                        P_dist = np.nan_to_num(P_dist)
                        Fist = np.nan_to_num(Fist)
                    
                    if KDE_tool == 'scipy':
                        Craft = Quanted_set[Indexes,:]
                        ### Resampling from estimate distribution is
                        ### preferred for cdf estimates.
                        kde = stats.gaussian_kde(Craft.T)
                        Craft = kde.resample(KDE_samples).T
                        P_dist = np.log(kde(Craft.T))
                        Fist = np.log(kde(data.T))
                        P_dist = np.nan_to_num(P_dist)
                        Fist = np.nan_to_num(Fist)
                    
                    ####
                    #### just add the copy the previous window if this one is completely NA. 
                    if sum(np.isnan(P_dist)) == len(P_dist):
                        if len(Likes[D]) == 0:
                            Likes[D].append([int(x in range(Where[0],Where[1])) for x in range(data.shape[0])])
                            Accurate.append([int(x in range(Where[0],Where[1])) for x in range(data.shape[0])])
                        else:
                            Likes[D].append(Likes[D][-1])
                            Accurate.append(Likes[D][-1])
                        continue
                    
                    ####
                    ### Neutrality tests of filtered reference pop KDE derived log-Likelihoods.
                    ## Normality.append(scipy.stats.mstats.normaltest(P_dist)[1])
                    ######
                    if normalize == 'CDF':
                        if np.std(P_dist)== 0:
                            Fist= np.array([int(Fist[x] in P_dist) for x in range(len(Fist))])
                        else:
                            Fist = scipy.stats.norm(np.mean(P_dist),np.std(P_dist)).pdf(Fist)
                    else:
                        Fist = kde(data.T)
                        Fist = Fist / max(Fist)
                        Fist = [round(x,5) for x in Fist]
                    if Control == True:
                        Fist = Fist * Mortal
                    
                    ###
                    ### set values of accessions with > individual NA threshold to 0.
                    for idx in nidx:
                        Fist[x]= 0

                    Accurate.append(Fist)
                    Likes[D].append(Fist)
                
                Points.append(Miss[window_start][0])
                Points_out.append(Miss[Index][0])
                s2 = time.time()
                elapsed = s2 - s1
                Intervals.append(elapsed)
                
                overlap= args.overlap
                
                SequenceStore = {fy:SequenceStore[fy][overlap:] for fy in Whose}
                Win = Window - overlap
                if end - Index < Window:
                    Window += end - Index
        
        Progress(Index,Intervals,end,50,Window)
        Index += 1
    
    Geno.close()
    
    if args.dr == 'PCA':
        PC_var= {CHR:{Points[star]:PC_var[star] for star in range(len(Points))}}
    else:
        PC_var= {CHR:{Points[star]:[0]*n_comp for star in range(len(Points))}}
    
    Out= {CHR:{Points[star]:Points_out[star] for star in range(len(Points))}}
    
    return {CHR:{start:Likes}},{CHR:Construct}, Out, PC_var


################  ###################  ###############################  ###################
###############    ###################  #############################  #####################

import multiprocessing as mp

nbProcs= int(args.proc)

parameters = Set_up(CHR,CHR + 1,nbProcs,MissG)

def what(job):
    try:
        rslt = Main_engine(*job)
    except Exception as e:
        print(e)
        rslt = 1
    finally:
        return rslt


print(parameters)

listJobs = []
for setting in parameters:
    listJobs.append([Fam,MissG,Geneo,Parents,GenoSUF,setting[0],setting[1],setting[2],args])

pool = mp.Pool(processes=nbProcs)
results = pool.map(what, listJobs)

### The current version is only applied to one genotype file at the time.
### It was original written to deal with multiple geno files.
### Different chromosomes still figure in the lines below, this is dealt 
### withn in the line 554: if Clover[repas]

Clover= {CHR: recursively_default_dict() for CHR in range(1,13)}
Construct= {CHR: recursively_default_dict() for CHR in range(1,13)}
Out= {CHR: recursively_default_dict() for CHR in range(1,13)}
PC_var= {CHR: recursively_default_dict() for CHR in range(1,13)}

for element in results:
    for k in element[0].keys():
        Clover[k].update(element[0][k])
        Out[k].update(element[2][k])
        Construct[k].update(element[1][k])
        PC_var[k].update(element[3][k])

Topo= {CHR:recursively_default_dict()}


for repas in Clover.keys():
    print(repas)
    if Clover[repas]:
        Topo[CHR] = {x:[] for x in range(len(Parents))}
        
        for block in sorted(Clover[repas].keys()):
            for reef in range(len(Parents)):
                Topo[CHR][reef].extend(Clover[repas][block][reef])


Points = sorted(Out[CHR].keys())

 ################## #################################################
############# WRITE ####################################################
 ############## ######################################################

start= args.id

print('writting analysis id:{0} to directory {1}'.format(args.id,Home))

Output = open(Home + "Blocks_Request_st"+str(start)+"_CHR" + str(CHR).zfill(2) + ".txt","w")

Output.write("CHR\tIn\tOut\tRef\t")

Crossed.extend(Parents)

for var in it.chain(*[Geneo[x] for x in Crossed]):
    Output.write(Fam[var] + "\t")

Output.write("\n")

for block in range(len(Topo[CHR][0])):
    for ref in Topo[CHR].keys():
        Output.write(str(CHR) + "\t")
        Output.write(str(int(Points[block])) + "\t")
        Output.write(str(int(Out[CHR][Points[block]])) + "\t")
        Output.write(str(ref) + '\t')
        for ass in range(len(Topo[CHR][ref][block])):
            Output.write(str(Topo[CHR][ref][block][ass]) + "\t")
        Output.write("\n")

Output.close()

#
#
#

if args.MSprint == True:
    Output= open(Home + 'Blocks_profiles_st'+str(start)+'_CHR'+ str(CHR).zfill(2)+ '.txt','w')
    
    Output.write('CHR\tIN\tcluster\t')
    
    for var in it.chain(*[Geneo[x] for x in Crossed]):
        Output.write(Fam[var] + "\t")
    
    Output.write("\n")
    
    for prf in Construct[CHR].keys():
        for cl in Construct[CHR][prf].keys():
            Output.write(str(CHR) + "\t")
            Output.write(str(int(prf)) + '\t')
            Output.write(str(cl) + '\t')
            Output.write('\t'.join([str(round(x,5)) for x in Construct[CHR][prf][cl]]))
            Output.write('\n')
    
    Output.close()


if args.VARprint == True:
    Output= open(Home + 'Blocks_ExVAR_st'+str(start)+'_CHR'+ str(CHR).zfill(2)+ '.txt','w')
    
    Output.write('CHR\tIN\t' + '\t'.join(['PC{0}'.format(x + 1) for x in range(args.ncomp)]))
    
    Output.write("\n")
    
    for prf in PC_var[CHR].keys():
        Output.write(str(CHR) + '\t')
        Output.write(str(int(prf)) + '\t')
        Output.write('\t'.join([str(round(x,5)) for x in PC_var[CHR][prf]]))
        Output.write('\n')
    
    Output.close()


print('Done.')
