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

####################################### Import modules: 
import os
import sys
from glob import glob
import subprocess
from decimal import Decimal 
from snakemake.utils import min_version 
#######################################

min_version("5.4.0") #Request a minimum version of snakemake
cwd= os.getcwd()

############################################################################# --- 1) Get Genome information --- ########## 
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

### 1.1) Check Genome size:
cmd= "awk 'BEGIN {total=0} {if($0 !~ /^>/) {total+=length}} END{print total}' " + genome
p= subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, universal_newlines=True)
(result,error)=p.communicate()
exit_code=p.returncode
if exit_code:
	raise subprocess.CalledProcessError(exit_code, cmd)
result= int(result) #convert str to int
print("The worked out genome length is: {:.2E}".format(Decimal(result)))

############################################################################# --- 2) Get indexed-genome information --- ##########
#get indexed-genome:
index_genome= os.path.abspath(config["index"])
if not os.path.exists(index_genome):
	print("ERROR: The indexed-genome file cannot be accessed: '{0}'".format(index_genome))
	sys.exit()

index_end= None
index_end= True if index_genome.lower().endswith(".fa.fai") else False
if not index_end:
	print("ERROR: indexed-genome: '{0}' should end: '.fa.fai' ".format(index_genome))
	sys.exit()

############################################################################# --- 3) Get the GTF file information --- ##########
#get gtf:
gtf_path= os.path.abspath(config["gtf"])
if not os.path.exists(gtf_path):
	print("ERROR: The GTF file cannot be accessed: '{0}'".format(gtf_path))
	sys.exit()

gtf_end= None
gtf_end= True if gtf_path.lower().endswith(".gtf") else False
if not gtf_end:
	print("ERROR: GTF file: '{0}' should end: '.gtf'".format(gtf_path))
	sys.exit()

### 3.1) Check Number of Genes:
cmd= "awk -F \"\\t\" '$3 ~ \"gene\" {print $0}'  " + gtf_path + " | wc -l" 
p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, universal_newlines=True)
(result,error) = p.communicate()
exit_code = p.returncode
if exit_code:
    raise subprocess.CalledProcessError(exit_code,cmd)
result= int(result)
print("The worked out GTF file contains: {} genes".format(result))

### 3.2) Check Number of Transcripts:
cmd= "awk -F \"\\t\" '$3 ~ \"transcript\" {print $0}'  " + gtf_path + " | wc -l"
p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, universal_newlines=True)
(result,error) = p.communicate()
exit_code = p.returncode
if exit_code:
    raise subprocess.CalledProcessError(exit_code,cmd)
result= int(result)
print("The worked out GTF file contains: {} transcripts".format(result))

############################################################################# ---- 4) Define the output folder --- ##########
#output folder:
OUTPUT= os.path.abspath(config["output"])

############################################################################# ---- 5) Get sample information --- ##########
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

############################################################################# ---- 6) Create Cluster_logs folder --- ##########
cluster_logs_dir= os.path.join(cwd,"logs","cluster")
if not os.path.exists(cluster_logs_dir):
	os.makedirs(cluster_logs_dir)

#######################
# RULES STARTS HERE
#######################

rule all: 
	input:
		expand(os.path.join(OUTPUT, "genome_index")),
		expand(os.path.join(OUTPUT, "{sample}.test.txt"), sample= rna_seq_samples) 

#######################
# ALIGNMENT & ASSEMBLY WORKFLOW
#######################

# run STAR_index
# ------------------
rule STAR_index:
	input:
		genome= genome,
		gtf= gtf_path
	output:
		directory(expand(os.path.join(OUTPUT, "genome_index")))
	log:
		expand(os.path.join(OUTPUT,"genome_index", "STAR_index.log"))
	threads: 4
	shell:
		"(mkdir {output} && "
		"STAR --runThreadN {threads} "
		"--runMode genomeGenerate "
		"--genomeDir {output} "
		"--genomeFastaFiles {input.genome} "
		"--sjdbGTFfile {input.gtf} "
		"--sjdbOverhang 100) 2> {log}"

# run STAR_mapping
# ------------------
def get_r1(wildcards):
	"""Code that returns a list of fastq files for read 1 based on *wildcards.sample* e.g."""
	return glob(os.path.join(read_dir,wildcards.sample + "*R1.fastq*"))

def get_r2(wildcards):
	"""Code that returns a list of fastq files for read 1 based on *wildcards.sample* e.g."""
	return glob(os.path.join(read_dir,wildcards.sample + "*R2.fastq*"))

rule STAR_mapping:
	input:
		lambda wildcards: config["samples"][wildcards.sample], 
		r1= get_r1,
		r2= get_r2,
		idx= os.path.join(OUTPUT, "genome_index")
	output:
		os.path.join(OUTPUT, "{sample}.test.txt")
	threads: 8
	shell:
		"echo \"Hello\" > {output}"


		# "STAR --runThreadN 8 --genomeDir {input.idx} " 
		# "--readFilesIn {input.r1} {input.r2} --outSAMunmapped Within "
		# "--outFilterType BySJout --outSAMattributes NH HI AS NM MD "
		# "--outFilterMultimapNmax 10 --outFilterMismatchNmax 999 "
		# "--outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 "
		# "--alignIntronMax 1000000 --alignMatesGapMax 1000000 "
		# "--alignSJoverhangMin 8 --alignSJDBoverhangMin 1 "
		# "--readFilesCommand pigz -p8 -dc --outSAMtype BAM SortedByCoordinate "
		# "--quantMode TranscriptomeSAM "
		# "--limitBAMsortRAM 33285996544"

#Mismatches: 4 and multimapping: 10
#explanation of SAM tags: https://samtools.github.io/hts-specs/SAMtags.pdf
#pigz: Parallel gzip, -p=threads
#@RG information: https://samtools.github.io/hts-specs/SAMv1.pdf

######## Customize: between --quantMode & --limitBAMsortRAM 
# --outSAMattrRGline ID:fEWP.2 PU:fEWP.2 SM:fEWP.2 PL:ILLUMINA CN:CRG --limitBAMsortRAM 33285996544
# ID: Read group identifier, PU: Platform Unit , SM: sample , PL: Platform, CN: Name of the Sequencing center 


##### ---- Strandness information --- ##########  

# RSeQC: infer_experiment.py: http://rseqc.sourceforge.net/#infer-experiment-py









