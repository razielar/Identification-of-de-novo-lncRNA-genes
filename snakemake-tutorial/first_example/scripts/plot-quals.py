import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt 
from pysam import VariantFile

quals=[i.qual for i in VariantFile(snakemake.input[0])]

#Histogram:
plt.hist(quals)

plt.savefig(snakemake.output[0])

