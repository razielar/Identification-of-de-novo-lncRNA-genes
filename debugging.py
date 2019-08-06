
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
genome='/users/rg/projects/references/Genome/D.melanogaster/dm6/dm6.fa'
genome_base= os.path.basename(genome)
print(genome_base)






