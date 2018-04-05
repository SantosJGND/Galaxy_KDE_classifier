# -*- coding: utf-8 -*-
"""
Created on Tue Apr 03 09:44:25 2018

@author: jgarcia
"""


from Kernel_tools import *
import collections
import time
import re

import os
import argparse
import pandas as pd


parser = argparse.ArgumentParser()

### optional arguments
parser.add_argument("books",type=str,metavar= 'N',nargs= '+',
                    help = "Reference files to read from genome crawl. Any number can be given.")

parser.add_argument("--gff",type= str,help = "gff file ")

parser.add_argument("--extract",type= str,default= 'gene',help = 'what do you want to merge with the widow crawl? have to be elements in gff file (gene, miRNA). Comma delimited')

parser.add_argument("--interest",type= str,default= 'ID,Name,Note',help = 'information elements to retain on targeted extractions. Comma delimited.')

parser.add_argument("-o",type= str,default= '',help = 'Output directory')

parser.add_argument("--id",type= str,default= 'MSU7_blockMerge',help = "Name your analysis. Output directory.")

args = parser.parse_args()


def recursively_default_dict():
        return collections.defaultdict(recursively_default_dict)


Home= args.o

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
        
        
        for start in Shelf.start.unique():
            
            Series= Shelf[(Shelf.start== start)]
            
            blocks_file= Series.loc[Series['tag'] == 'Request','file']
            blocks_file= blocks_file.iloc[0]
            
            ##### read blocks file
            Input= open(blocks_file,'r')
            d= 0
            
            for line in Input:
                line= line.split()
                
                if d== 0:
                    line=line[4:]
                    Names= line
                else:
                    if line[3] == '0':
                        Out[int(line[0])][int(line[1])] = int(line[2])
                d += 1
            
            Input.close()
            
            ##### read profiles file
            
    
    s1= time.time()
    
    print(str(s1-s0) + ' seconds elapsed.')
    print(str((s1-s0) / float(60)) + ' minutes elapsed.')
    return Names, Out



### Genic information - connection with gff3 file.

def extract_genes(Out,gffile,extract,interest):
    Genes = []
    Input = open(gffile,"r")
    for line in Input:
        if line[0] == '#':
            continue
        line = line.split('\t')
        info = []
        if line[2] in extract:
            vault = {x:[] for x in interest}
            
            Chr= int(line[0].replace('chr',''))
            
            start= int(line[3])
            end= int(line[4])
            
            block_tag= [bl for bl in Out[Chr].keys() if bl <= start and Out[Chr][bl] >= start]
            
            
            loot = line[8].split(';')
            for box in loot:
                box = box.split('=')
                if box[0] in interest:
                    vault[box[0]].append(box[1])
            if block_tag:
                bill = [block_tag[0],Chr,start,end]
                
                for v in interest:
                    if len(vault[v]) == 0:
                        vault[v].append('-')
                    if len(vault[v]) > 1:
                        vault[v] = ['_'.join(vault[v])]
                    bill.extend(vault[v])
                
                Genes.append(bill)
    
    Input.close()
    colnames = ['tag','chr','start','end']
    colnames.extend(interest)
    
    Genes = pd.DataFrame(Genes,columns = colnames)
    return Genes



### extract windows surveyed during genome crawl
print('to begin reading from: ')
print(args.books)

Library= read_books(args.books)

print('library:')
print(Library.sort_values(by= ['Chr','start']))

filename= ''.join([Home,'/',args.id,'.txt'])

print('eventually writting to: {}'.format(filename))
os.makedirs(os.path.dirname(filename), exist_ok=True)


## prepare extraction:
extract= args.extract
extract= extract.split(',')

interest= args.interest
interest= interest.split(',')

print('extracting columns {} for target {}'.format(interest,extract))


Names, Out = read_3D_profiles(Library)


print([len(x) for x in Out.values()])

## Merge files
Gene_fit = extract_genes(Out,args.gff,extract,interest)

print(Gene_fit.head())
## printout

Output = open(filename,'w')

Output.write('\t'.join(Gene_fit.columns))
Output.write('\n')

for gene in range(len(Gene_fit)):
    Output.write('\t'.join([str(c) for c in Gene_fit.iloc[gene,:]]))
    Output.write('\n')

Output.close()

