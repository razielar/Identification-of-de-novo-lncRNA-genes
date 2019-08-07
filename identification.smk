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

########## --- 2) Get the indexed-genome information --- ##########

#get indexed-genome:
index_genome= os.path.abspath(config["index"])
if not os.path.exists(index_genome):
	print("ERROR: The indexed-genome file cannot be accessed: '{0}'".format(index_genome))
	sys.exit()

index_end= None
index_end= True if index_genome.lower().endswith(".fa.fai") else False
if not index_end:
	print("ERROR: indexed-genome: '{0}' should end: '.fa.fai' ".format(index_genome)) 

########## --- 3) Get the GTF file information --- ##########
gtf_path= os.path.abspath(config["gtf"])
if not os.path.exists(gtf_path):
	print("ERROR: The GTF file cannot be accessed: '{0}'".format(gtf_path))
	sys.exit()





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
for sample_name, reads in all_samples.items(): #variables={"sample_name": "Sample16", "reads": ["reads[0]","reads[1]"]}
	rna_seq_samples.append(sample_name)
	reads_dict= {}
	reads_dict["R1"]= reads[0]
	reads_dict["R2"]= reads[1]
	for pair, path in reads_dict.items():
		r_path= os.path.abspath(path)
		gzipped= True if r_path.lower().endswith(".gz") else False
		if not os.path.exists(r_path):
			print("ERROR: The pair '{0}' of read '{1}' for sample '{2}' cannot be accessed".format(pair,r_path,sample_name))
			sys.exit()
		if not os.path.exists(read_dir):
			os.makedirs(read_dir)
		new_read_name_list= (sample_name,pair,"fastq.gz") if gzipped else (sample_name,pair,"fastq")
		new_read_name= ".".join(new_read_name_list)
		cmd= "cd " + read_dir + " && ln -s "+ r_path + " " + new_read_name
		if not os.path.exists(os.path.join(read_dir, new_read_name)):
			os.system(cmd)


#######################
# RULES STARTS HERE
#######################

#######################
# ALIGNMENT & ASSEMBLY WORKFLOW
#######################


##### ---- Strandness information --- ##########  

# RSeQC: infer_experiment.py: http://rseqc.sourceforge.net/#infer-experiment-py

##### ---- 4) create logs folder --- ##########
cluster_logs_dir= os.path.join(cwd,"logs","cluster")
if not os.path.exists(cluster_logs_dir):
	os.makedirs(cluster_logs_dir)








