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
from decimal import Decimal  

#Request a minimum version of snakemake 
from snakemake.utils import min_version 
min_version("5.4.0")

cwd= os.getcwd()

########## --- 1) Get the genome information --- ########## 
# get genome:
genome=os.path.abspath(config["genome"])
if not os.path.exists(genome):
	print("ERROR: The genome file cannot be accessed: {0}".format(genome))
	sys.exit()

# check if genome is gzipped and exit if true: 
genome_base=os.path.basename(genome)
gzipped= None
gzipped= True if genome.lower().endswith(".gz") else False
if gzipped:
	print("ERROR: {0} file should not be gzipped, unzip it and then return".format(genome))
	sys.exit()

# check the genome size:
cmd= "awk 'BEGIN {total=0} {if($0 !~ /^>/) {total+=length}} END{print total}' " + genome
p= subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, universal_newlines=True)
(result,error)=p.communicate()
exit_code=p.returncode
if exit_code:
	raise subprocess.CalledProcessError(exit_code, cmd)
result= int(result) #convert str to int
print("The worked out genome length is: {:.2E}".format(Decimal(result)))

##### ---- 2) Load the tools needed to run the pipeline --- ##########

# STAR_version= None
# load= config["load"]

# for tool, source in load.items():
# 	if tool in ("STAR"):
# 		source_dump, STAR_version= tool.split("")

##### ---- 3) Define the output folder --- ##########

OUTPUT= os.path.abspath(config["output"])

##### ---- 4) Get sample information --- ##########

#work out if samples are gzipped or not 
rna_seq_samples = []
read_dir= os.path.join(OUTPUT, "data", "reads")

all_samples= config["samples"]
gzipped= None
for sample_name, reads in all_samples.items():
	rna_seq_samples.append(sample_name)
	reads_dict= {}
	reads_dict["R1"]= reads[0]
	reads_dict["R2"]= reads[1]

for x,y in reads_dict.items():
	print("keys:{0} -- values:{1}".format(x,y))






















