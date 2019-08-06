
import os
import sys
from glob import glob 
import subprocess 

#Request a minimum version of snakemake 

from snakemake.utils import min_version 
min_version("5.4.0")

cwd= os.getcwd()
print(cwd)

########## --- 1) Get the genome information 
# x = -1 
# if x < 0:
# 	raise Exception("Sorry no numbers below zero")


# x = "hello"

# if not type(x) is int:
#   raise TypeError("Only integers are allowed")

print("module load STAR/2.5.2a-foss-2016a".split(" "))





