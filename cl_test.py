# -*- coding: utf-8 -*-
"""
Created on Wed Feb 07 17:45:43 2018

@author: jgarcia
"""

from Kernel_tools import *
import collections
import time
import itertools as it
import numpy as np
import re

import sys
import argparse
parser = argparse.ArgumentParser()

parser.add_argument("books",type=str,metavar= 'N',nargs= '+',
                    help = "Fam file, ind id's, same as for kernel analysis")

parser.add_argument("--number",type=int,help = "int test")

args = parser.parse_args()


print(args.books)
print(args.number)
