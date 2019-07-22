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




















