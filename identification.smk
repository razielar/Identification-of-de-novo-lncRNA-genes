#!/usr/bin/env python3 

"""
Snakemake workflow to identify the novo lncRNAs 

22/07/2019

"""
#Authorship and License information: 
__author__ = "Raziel Amador Rios"
__copyright__ = "Copyright 2019"
__license__ = "GNU General Public License"
__maintainer__ = "Raziel Amador Rios"
__email__ = "raziel.amador@crg.eu"
__status__ = "Development"
__version__ = "1.0"

### Import modules: 
import os
import sys
from glob import glob
import subprocess 

#Request a minimum version of snakemake 
from snakemake.utils import min_version 
min_version("5.4.0")

cwd= os.getcwd()

########## --- 1) Get the genome information 
# get genome:
genome=os.path.abspath(config["genome"])
if not os.path.exists(genome):
	print("ERROR: The genome file cannot be accessed: {0}".format(genome))
	sys.exit()

# check if genome is gzipped and exit if true: 
genome_base=os.path.basename(genome)
gzipped= None
gzipped= True if genome.lower().endswith(".gz") else False
if not gzipped:
	print("The genome file is not gzipped: {0}".format(genome))
	sys.exit()



























