
import os
import sys
from glob import glob 
import subprocess 

#Request a minimum version of snakemake 

from snakemake.utils import min_version 
min_version("5.4.0")

cwd= os.getcwd()
# print(cwd)

########## --- 1) Get the genome information 
# x = -1 
# if x < 0:
# 	raise Exception("Sorry no numbers below zero")


# x = "hello"

# if not type(x) is int:
#   raise TypeError("Only integers are allowed")

# print("module load STAR/2.5.2a-foss-2016a".split(" "))

thisdict =	{
  "brand": "Ford",
  "model": "Mustang",
  "year": 1964
}

print(thisdict)

x = thisdict["model"]
print(x)


############################################# ---- APPENDIX from identification.smk ---- #####################################################

##### ---- Load the tools needed to run the pipeline --- ##########

# STAR_version= None
# load= config["load"]

# for tool, source in load.items():
# 	if tool in ("STAR"):
# 		source_dump, STAR_version= tool.split("")

##### ---- Generate symbolik links of gtf and genome.fa inside genome_index_folder --- ##########

# #GTF file
# cmd= "cd " + genome_index_folder + " && ln -s " + gtf_path + " " + "gtf"
# if not os.path.exists(os.path.join(genome_index_folder,"gtf")):
# 	os.system(cmd)
# #Genome.fa
# cmd= "cd " + genome_index_folder + " && ln -s "+ genome + " " + "genome.fa"
# if not os.path.exists(os.path.join(genome_index_folder, "genome.fa")):
# 	os.system(cmd)









