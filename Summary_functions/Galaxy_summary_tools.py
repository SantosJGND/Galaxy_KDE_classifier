
import os
import argparse

import collections
import time
import itertools as it
import numpy as np
import re
import pandas as pd

##########################################
########## Load data.

def read_refs(index_file):
    indxs = recursively_default_dict()
    
    Input = open(index_file,'r')
    for line in Input:
        line = line.split()
        indxs[int(line[0])][line[1]] = []
    
    Input.close()
    
    indxs = {gop:[x for x in indxs[gop].keys()] for gop in indxs.keys()}
    
    return indxs, [x for x in sorted(indxs.keys())]



def read_focus(index_file):
    indxs = []
    
    Input = open(index_file,'r')
    for line in Input:
        line = line.split()
        indxs.append(line[0])
    
    Input.close()
    
    return indxs


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



def read_geno_books(books):
    '''
    reads files. required pattern: _chr(i)
    where i = chromosome number.
    Tag will be string preceding underscore.
    '''
    library= []
    
    for shelf in books:
        card= shelf.split('/')
        
        cover= card[-1].split('_')
        print(cover)
        Chr= int([re.findall(r'\d+',i)[0] for i in cover if re.search('chr',i)][0])
        print(Chr)
        tag= cover[0]
        
        library.append([shelf,tag,Chr])
    
    library= pd.DataFrame(library,columns= ['file','tag','Chr'])
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


def read_3D_profiles_list(File_list):
    
    s0 = time.time()
    
    Blocks= recursively_default_dict()
    Out = recursively_default_dict()
    Names= []        
    
    for File in File_list:
        
        ##### read blocks file
        Input= open(File,'r')
        d= 0
        
        for line in Input:
            line= line.split()
            
            if d== 0:
                line=line[4:]
                Names= line
            else:
                Blocks[int(line[0])][int(float(line[1]))][int(line[3])]= [float(x) for x in line[4:]]
                if line[3] == '0':
                    Out[int(line[0])][int(float(line[1]))] = int(float(line[2]))
            d += 1
        
        Input.close()
        #    
    s1= time.time()
    
    print(str(s1-s0) + ' seconds elapsed.')
    print(str((s1-s0) / float(60)) + ' minutes elapsed.')
    return Blocks, Names, Out



##########################################
######### Classification



def Merge_class(Ref_profiles,focus_indicies,Out,Diff_threshold,BIN,X_threshold,coarse,sg_order):
    Blocks_genome = recursively_default_dict()
    
    for CHR in Ref_profiles.keys():
        print(CHR)
        Points = sorted(Ref_profiles[CHR].keys())
        Likes = Ref_profiles[CHR]
        N_pops= len(Likes[[x for x in Likes.keys()][0]])
        
        print("number of reference populations: {0}".format(N_pops))
        Likes = {x:[Likes[bl][x] for bl in sorted(Likes.keys())] for x in range(N_pops)}
        Likes = {x:np.array(Likes[x]) for x in Likes.keys()}
        
        Topo = []
        
        #range_Parents = [x + Aro.shape[0] for x in range(Daddy.shape[0])]
        #range_Crossed = [x for x in range(Aro.shape[0])]
        
        for acc in focus_indicies:
            Guys = np.array([Likes[x][:,acc] for x in range(N_pops)])
            Guys = np.nan_to_num(Guys)
            Guys = [[[y,0][int(y<=X_threshold)] for y in x] for x in Guys]
            
            Test = [int(x <= X_threshold) for x in np.amax(np.array(Guys),axis = 0)]     
            
            if coarse:
                
                Guys = [savgol_filter(x,BIN,sg_order,mode = "nearest") for x in Guys]
                Test = savgol_filter(Test,BIN,sg_order,mode = "nearest")
                Test = [round(x) for x in Test]
            
            #
            Guys = np.array(Guys).T
            
            maxim = np.argmax(Guys,axis = 1)
            where_X = [x for x in range(Guys.shape[0]) if Test[x] == 1]
            
            #
            Consex = [x for x in it.combinations(range(N_pops),2)]
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
                            Diff = abs(max(Diff)) / abs(min(Diff))
                            Diff = int(Diff > Diff_threshold)
                        
                        if Diff == 0:
                            CL.append(j)
                    
                    if len(CL) == 2:
                        maxim[h] = 7
                    if len(CL) == 1:
                        maxim[h] = sum(CL[0]) + N_pops
            
            maxim[where_X] = N_pops
            
            if not Consex:
                for h in range(len(maxim)):
                    maxim[h] = int(10*Guys[h,0])    
            
            
            Topo.append(maxim + 1)
            
        Topo = np.array(Topo).T
        
        Clove = {CHR:{Points[x]:Topo[x,] for x in range(len(Points))}}
        
        Blocks_genome.update(Clove)
    
    return Blocks_genome


##########################################
######### Classification processing


def compress_ideo(df,Out,chromosome_list):
    '''
    Compress classification matrix by individual.
    
    '''
    new_set = []
    
    for CHR in range(len(chromosome_list)):
        
        Chr = int(re.search('chr(.+?)_',chromosome_list[CHR]).group(1))
        sub = df[df.chrom == chromosome_list[CHR]]
        Id= '_'.join(chromosome_list[CHR].split('_')[1:])
        Coordinates = sorted(sub.start)
        Size = sub.shape[0]
        start = min(sub.start)
        First = sub.gieStain.iloc[0]
        for index in range(len(Coordinates)):
            row = sub[sub.start == Coordinates[index]]
            if index == 0:
                continue
            if index == (Size - 1):
                if row.gieStain.iloc[0] == First:
                    new_set.append([Id,Chr,start,Out[Chr][max(sub.start)],First])
                else:
                    new_set.append([Id,Chr,start,Out[Chr][max(sub.start)],First])
                    First = row.gieStain.iloc[0]
                    start = row.start.iloc[0]
                    new_set.append([Id,Chr,start,Out[Chr][max(sub.start)],First])
            else:
                if row.gieStain.iloc[0] == First:
                    continue
                else:
                    new_set.append([Id,Chr,start,row.start.iloc[0]-1,First])
                    First = row.gieStain.iloc[0]
                    start = row.start.iloc[0]
        
    new_set = pd.DataFrame(new_set,columns = ['ID','chrom', 'start', 'end', 'gieStain'])
    return new_set



##########################################
######### Complementary

def recursively_default_dict():
        return collections.defaultdict(recursively_default_dict)

