# -*- coding: utf-8 -*-
"""
Created on Thu Feb 01 15:43:00 2018

@author: jgarcia
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Jan 23 14:23:15 2018

@author: jgarcia
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Dec 06 10:58:26 2017

@author: jgarcia
"""

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

parser.add_argument("--bin",default = 5,type= int,help = "smoothing parameter [savgol filter]")
###
parser.add_argument("--MSprint",action= "store_false",help = "if given prints cluster stats.")
###
parser.add_argument("-c",action = "store_true",help = "specific accession choice file")
###
parser.add_argument("-w",type = int,default = 200, help = "Window size - markers")
### 
parser.add_argument("--outlier",type=float,default = 1e-3,help = "Outlier threshold")
### 
parser.add_argument("--threshold",type = float,default = 1.6,help = "Intermediate classification threshold")
### 
parser.add_argument("--clustmethod",default = "MeanShift",choices = ["MeanShift","DBscan","HDBscan"],help = "Clustering method to extract reference specific clusters. MS, dbscan and hdbsan available. MS preferred")
###
parser.add_argument("--het",type = float,default = 5e-2, help = "Heterozygosity filter")
### 
parser.add_argument("--dr",default = 'NMF',help = "Dimensionality reduction. options: PCA, NMF")
### 
parser.add_argument("--ncomp",type = int,default = 4,help = "Number of components kept in case of PCA reduction")
### 
parser.add_argument("--outmethod",default = "None",help = "Outlier filter method. options: DBSCAN, NMF, Perc, Z, Isolation_forest.")
###
parser.add_argument("--overlap",type= int,default = 100,help = "Overlap between windows, in snps")
###

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
    
    return indxs, [x for x in indxs.keys()]


####

Fam = FAMread(args.fam)

MissG, Gindex = BIMread(args.bim)

GenoSUF = args.geno

admx_lib, Crossed = read_refs(args.admx)
refs_lib, Parents = read_refs(args.ref)

admx_lib.update(refs_lib)

Geneo = admx_lib

print(Geneo.keys())
print([len(x) for x in Geneo.values()])

CHR = args.CHR
BIN = args.bin

Home = ''

####
####


def Set_up(Chr0,Chr1,Sub,MissG):
    Settings= []
    for Chr in range(Chr0,Chr1):
        Size= len(MissG[Chr])
        Steps = range(0,Size+1,int(Size/float(Sub)))
        Steps[-1] = Size
        for i in range(len(Steps)-1):
            Settings.append([Chr,Steps[i],Steps[i+1]-1])
    return Settings



def Main_engine(Fam,MissG,Geneo,Parents,GenoSUF,CHR,start,end,args):
    
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
    ## Outlier threshold
    X_threshold = args.outlier
    #### Intermediate classification threshold
    Diff_threshold = args.threshold
    Filter_Het = args.het
    
    
    ## Dimensionality reduction: PCA, NMF
    DIMr = args.dr
    n_comp = args.ncomp
    PC_var_lim = .01
    ##
    ## savgol filter parameter
    BIN = args.bin
    
    ## Chose clustering method
    Method = args.clustmethod
    
    ## Outlier filter method: DBSCAN, NMF, Perc, Z, Isolation_forest.
    Outlier_method = args.outmethod
    
    
    quantile_filter_1 = 20 ### quantile filter for Percentile based outlier.
    
    ## KDE estimation tool
    # Attention: sklearn uses GridsearchCV. While this implies more parameters, it 
    ## nonetheless ensures the construction of a KDE when scipy KDE would break down.
    ## This is because scipy KDE will always estimate bandwidth itself, thus using the points given only,
    ## while sklearn KDE allows the bandwidth to be passed. This allows us in turn to set a minimum,
    ## and ensure it to never come to 0.
    ## sklearn, scipy
    
    KDE_tool = 'sklearn'
    Bandwidth_split = 30
    
    ### Output normalization: CDF, Max
    normalize = 'CDF'
    
    ### Resampling size for CDF estimates.
    KDE_samples = 1000
    
    # Controlling for global likelihood.
    Control = False
    
    
    ##
    ##
    ##
    ##
    ### These things define themselves
    
    
    Crossed = [x for x in Geneo.keys() if x not in Parents]
    
    t = [x for x in Miss.keys() if Miss[x][0] >= start and Miss[x][0] <= end]
    
    Whose = list(it.chain(*Geneo.values()))
    SequenceStore = {fy:[] for fy in Whose}
    
    Likes = {x:[] for x in range(len(Parents))}
    
    Accuracy = []
    
    Geno = open(GenoFile,"r")
    Points = []
    Points_out = []
    
    Win = 0
    Index = 0
    Intervals = []
    
    Construct = recursively_default_dict()
    
    for line in Geno:
        Codes = [0,1,2,0,0,0,0,0,0,nan]
        d = Miss[Index][0]
        if Index > end:
            break
        if len([x for x in Whose if line[x] =='1']) / float(len(Whose)) > Filter_Het:
            Index += 1
            continue
        if Index >= start and Index <= end:
            for judas in SequenceStore.keys():
                SequenceStore[judas].append(Codes[int(line[judas])])
            Win += 1
            if Win == Window and (end - Index - Window) > (Window/4):
                s1 = time.time()
                window_start = Index - Window + 1
                Seq = SequenceStore
                
                Aro = np.nan_to_num(np.array([Seq[x] for x in it.chain(*[Geneo[x] for x in Crossed])]))
                
                Daddy = np.nan_to_num(np.array([Seq[x] for x in it.chain(*[Geneo[x] for x in Parents])]))
                
                b = np.concatenate((Aro,Daddy),axis = 0)
                data = np.zeros((b.shape[0],b.shape[1]+1))
                data[:,:-1] = b
                
                if DIMr == 'PCA':
                    pca = PCA(n_components=n_comp, whiten=False,svd_solver='randomized')
                    data = pca.fit_transform(data)
                    data = data[:,[x for x in range(data.shape[1]) if pca.explained_variance_ratio_[x] >= PC_var_lim]]
                
                if DIMr == 'NMF':
                    from sklearn.decomposition import NMF
                    data = NMF(n_components=n_comp, init='random', random_state=0).fit_transform(data)
                
                Accurate = []
                params = {'bandwidth': np.linspace(np.min(data), np.max(data),Bandwidth_split)}
                grid = GridSearchCV(KernelDensity(algorithm = "ball_tree",breadth_first = False), params,verbose=0)
                
                #####################################
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
                    
                    #n_clusters_ = len(labels_unique)
                    Tree = {x:[Focus_labels[y] for y in range(len(labels)) if labels[y] == x] for x in [g for g in list(set(labels)) if g != -1]}
                
                SpaceX = {x:data[Tree[x],:] for x in Tree.keys()}
                
                
                
                for hill in SpaceX.keys():
                    #Crease = [sum([Apples[hill][y]==x for y in range(len(Apples[hill]))]) / float(len(Apples[hill])) for x in Geneo.keys()]
                    if len(Tree[hill]) <= 3:
                        continue
                    grid.fit(data[Tree[hill],:])
                    #print("best bandwidth: {0}".format(grid.best_estimator_.bandwidth))
                    
                    # use the best estimator to compute the kernel density estimate
                    kde = grid.best_estimator_
                    #Quanted_set = kde.sample(KDE_samples)
                    P_dist = kde.score_samples(data[Tree[hill],:])
                    Dist = kde.score_samples(data)
                    P_dist= np.nan_to_num(P_dist)
                    Dist= np.nan_to_num(Dist)
                    if np.std(P_dist) == 0:
                        Dist= [int(Dist[x] in P_dist) for x in range(len(Dist))]
                    else:
                        Dist = scipy.stats.norm(mean(P_dist),np.std(P_dist)).cdf(Dist)
                    Dist= np.nan_to_num(Dist)
                    Construct[Miss[window_start][0]][hill] = Dist
                    
                
                
                ######################################### 
                ############# TEST #####################
                #########################################
                
                for D in range(len(Parents)):
                    Where = [sum([len(Geneo[x]) for x in Parents[0:D]])+Aro.shape[0],sum([len(Geneo[x]) for x in Parents[0:(D+1)]])+Aro.shape[0]]
                    Where = [int(x) for x in Where]
                    
                    ### FILTER OUT OUTLIERS
                    ### Conservative methods are preferred since only extreme values 
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
                    if sum(np.isnan(P_dist)) == len(P_dist):
                        if len(Likes[D]) == 0:
                            Likes[D].append([int(x in range(Where[0],Where[1])) for x in range(data.shape[0])])
                            Accurate.append([int(x in range(Where[0],Where[1])) for x in range(data.shape[0])])
                        else:
                            Likes[D].append(Likes[D][-1])
                            Accurate.append(Likes[D][-1])
                        continue
                    
                    ### Neutrality tests of filtered reference pop KDE derived log-Likelihoods.
                    #Normality.append(scipy.stats.mstats.normaltest(P_dist)[1])
                    ######
                    CDF = [len([y for y in P_dist if y <= x]) / float(len(P_dist)) for x in Fist]
                    #Fist = (Fist - np.mean(P_dist)) / np.std(P_dist)
                    if normalize == 'CDF':
                        if np.std(P_dist)== 0:
                            Fist= np.zeros(len(Fist))
                        else:
                            Fist = scipy.stats.norm(np.mean(P_dist),np.std(P_dist)).cdf(Fist)
                    else:
                        Fist = kde(data.T)
                        Fist = Fist / max(Fist)
                        Fist = [round(x,5) for x in Fist]
                    if Control == True:
                        Fist = Fist * Mortal
                
                    Accurate.append(Fist)
                    Likes[D].append(Fist)
                   
                ### Not under use, but 'Accurate' is useful for local analysis/testing.
#                Accurate = np.vstack(Accurate).T
#                Accuracy = np.argmax(Accurate,axis = 1) 
#                Similar = [sum(Accurate[g,:] == Accurate[g,Accuracy[g]]) for g in range(Accurate.shape[0])]
#                Washed = {d:[c + Aro.shape[0] + sum([len(Geneo[x]) for x in Parents[0:d]]) for c in range(len(Geneo[Parents[d]])) \
#                #if Similar[c + sum([len(Geneo[x]) for x in Parents[0:d]]) + Aro.shape[0]] < 3 and Accuracy[c + sum([len(Geneo[x]) for x in Parents[0:d]]) + Aro.shape[0]] == d] for d in range(3)}
#                #if Accuracy[c + sum([len(Geneo[x]) for x in Parents[0:d]]) + Aro.shape[0]] == d \
#                if Accurate[c + sum([len(Geneo[x]) for x in Parents[0:d]]) + Aro.shape[0],d] >= 0.1] for d in range(len(Parents))}
                #Accurate[[x for x in range(Where[0],Where[1]) if x not in Washed[2]]]
                ######
                
                #JC_visual.append(np.hstack((Likings,pca.explained_variance_ratio_[:5])))
                
                Points.append(Miss[window_start][0])
                Points_out.append(Miss[Index][0])
                s2 = time.time()
                elapsed = s2 - s1
                Intervals.append(elapsed)
                SequenceStore = {fy:SequenceStore[fy][args.overlap:] for fy in Whose}
                Win = Window - args.overlap
        
        Progress(Index,Intervals,end,50,Window)
        Index += 1
    
    Geno.close()
    
    
    #rm /gs7k1/home/jgarcia/Likes/COMP/*
    #qsub -N Bubba -V -q normal.q -b y python -u Launch.py Launch_PO.py 5
    
    Out= {CHR:{Points[star]:Points_out[star] for star in range(len(Points))}}
    
    return {CHR:{start:Likes}},{CHR:Construct},Out


################  ###################  ###############################  ###################
###############    ###################  #############################  #####################

import multiprocessing as mp

nbProcs= int(args.proc)
#nbProcs= 4
parameters = Set_up(CHR,CHR + 1,nbProcs,MissG)


def what(job):
	try:
		rslt = Main_engine(*job)
	except Exception as e:
		print e
		rslt = 1
	finally:
		return rslt


#parameters= [[1,2000,4000],[2,200,2000],[2,2200,3000],[6,23000,26000]]


listJobs = []
for setting in parameters:
    listJobs.append([Fam,MissG,Geneo,Parents,GenoSUF,setting[0],setting[1],setting[2],args])

pool = mp.Pool(processes=nbProcs)
results = pool.map(what, listJobs)



###### Merge results

def Smooth_class(Likes,Indicies,Diff_threshold,BIN,X_threshold):
    Points = sorted(Indicies.keys())
    Likes = {x:np.array(Likes[x]) for x in Likes.keys()}
    
    Topo = []
    
    #range_Parents = [x + Aro.shape[0] for x in range(Daddy.shape[0])]
    #range_Crossed = [x for x in range(Aro.shape[0])]
    
    for acc in range(Likes[0].shape[1]):
        Guys = np.array([Likes[x][:,acc] for x in range(len(Parents))])
        Guys = np.nan_to_num(Guys)
        
        
        Test = [int(x < X_threshold) for x in np.amax(np.array(Guys),axis = 0)]
        Test = savgol_filter(Test,BIN,3,mode = "nearest")
        Test = [round(x) for x in Test]
        
        Guys = [[[y,0][int(y<=X_threshold)] for y in x] for x in Guys]
        Guys = [savgol_filter(x,BIN,3,mode = "nearest") for x in Guys]
        #    
        Guys = np.array(Guys).T
        
        maxim = np.argmax(Guys,axis = 1)
        where_X = [x for x in range(Guys.shape[0]) if Test[x] == 1]
        #where_X = [x for x in range(Guys.shape[0]) if len([c for c in Guys[x,:3] if c <= .0001]) == 3]
        Consex = [x for x in it.combinations(range(len(Parents)),2)]
        if Consex:
            for h in range(len(maxim)):
                CL = []
                for j in Consex:
                    Diff = Guys[h,j]
                    if maxim[h] not in j or len([x for x in Diff if x < X_threshold]) > 0:
                        continue
                    if max(Diff) <= X_threshold:
                        Diff = 0
                    else:
                        #Diff = int(len([x for x in Diff if x > 1e-5]) == 1)
                        Diff = abs(max(Diff)) / abs(min(Diff))
                        Diff = int(Diff > Diff_threshold)
                    
                    #print(Diff)
                    if Diff == 0:
                        CL.append(j)
                
                if len(CL) == 2:
                    maxim[h] = 7
                if len(CL) == 1:
                    maxim[h] = sum(CL[0]) + len(Parents)
            maxim[where_X] = len(Parents)
        
        if not Consex:
            for h in range(len(maxim)):
                maxim[h] = int(10*Guys[h,0])    
        
        #Similar = [sum(Guys[g,:] == np.amax(Guys,axis=1)[g]) for g in range(Guys.shape[0])]
        #Similar = [sum(Guys[g,:] == 0) for g in range(Guys.shape[0])]
        #Similar = [int(maxim[x] > 3 or sum(Guys[x,:] == 0) > 2 and maxim[x] != 3) + 1 for x in range(len(maxim))]
        Similar = [int(sum(Guys[x,:] <= 0) == Guys.shape[1] and maxim[x] != 3) + 1 for x in range(len(maxim))]
        
        peaks = [Points[x] for x in range(len(Points)) if Similar[x] == 1]
        where_peaks = [x for x in range(len(Points)) if Similar[x] == 1]
        
        d = 'none'
        if peaks:
            for l in range(len(Similar)):
                if Similar[l] == 1:
                    if d == 'none' and l > 0:
                        maxim[:l] = [maxim[l] for x in range(l)]
                    d = l
                if Similar[l] > 1:
                    if d != 'none':
                        #if max(Guys[l,:]) > 0:
                        #    Close = [x for x in range(Guys.shape[1]) if Guys[l,x] == max(Guys[l,:])]
                        #    maxim[l] = sum(Close) + 3
                        #else:
                        Distances = [abs(peaks[x] - Points[l]) for x in range(len(peaks))]
                        #print(maxim[l],maxim[where_peaks[Distances.index(min(Distances))]],Similar[l])
                        maxim[l] = maxim[where_peaks[Distances.index(min(Distances))]]                 
        
        ###
        ### PARK nber: 1
        ###    
        
        Topo.append(maxim + 1)
    
    
    Topo = np.array(Topo).T
    
    Clove = {Points[x]:Topo[x,] for x in range(len(Points))}
    
    return Clove



Clover= {CHR: recursively_default_dict() for CHR in range(1,13)}
Construct= {CHR: recursively_default_dict() for CHR in range(1,13)}
Out= {CHR: recursively_default_dict() for CHR in range(1,13)}


for element in results:
    for k in element[0].keys():
        Clover[k].update(element[0][k])
        Out[k].update(element[2][k])
        Construct[k].update(element[1][k])




Topo= {CHR:recursively_default_dict()}

for repas in Clover.keys():
    print(repas)
    if Clover[repas]:
        Topo[CHR] = {x:[] for x in range(len(Parents))}
        
        for block in sorted(Clover[repas].keys()):
            for reef in range(len(Parents)):
                Topo[CHR][reef].extend(Clover[repas][block][reef])
        
#        Topo[repas]= Smooth_class(Likes,Out[repas],args.threshold,args.bin,args.outlier)

Points = sorted(Out[CHR].keys())

 ################## #################################################
############# WRITE ####################################################
 ############## ######################################################

start= 1
print 'writting to directory ' + Home

if args.MSprint == True:
    Output = open(Home + "Blocks_Request_st"+str(start)+"_CHR" + str(CHR).zfill(2) + ".txt","w")
    
    Output.write("CHR\tIn\tOut\tRef\t")
    
    Crossed.extend(Parents)
    
    for var in it.chain(*[Geneo[x] for x in Crossed]):
        Output.write(Fam[var] + "\t")
    
    Output.write("\n")
    
    for block in range(len(Topo[CHR][0])):
        for ref in Topo[CHR].keys():
            Output.write(str(CHR) + "\t")
            Output.write(str(Points[block]) + "\t")
            Output.write(str(Out[CHR][Points[block]]) + "\t")
            Output.write(str(ref) + '\t')
            for ass in range(len(Topo[CHR][ref][block])):
                Output.write(str(Topo[CHR][ref][block][ass]) + "\t")
            Output.write("\n")
    
    Output.close()

#
#
#

Output= open(Home + 'Blocks_profiles_st'+str(start)+'_CHR'+ str(CHR).zfill(2)+ '.txt','w')

Output.write('CHR\tIN\tcluster\t')

for var in it.chain(*[Geneo[x] for x in Crossed]):
    Output.write(Fam[var] + "\t")

Output.write("\n")

for prf in Construct[CHR].keys():
    for cl in Construct[CHR][prf].keys():
        Output.write(str(CHR) + "\t")
        Output.write(str(prf) + '\t')
        Output.write(str(cl) + '\t')
        Output.write('\t'.join([str(round(x,5)) for x in Construct[CHR][prf][cl]]))
        Output.write('\n')

Output.close()



