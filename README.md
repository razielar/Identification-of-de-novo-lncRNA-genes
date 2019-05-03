# Identification of *de novo* lncRNA genes

## 1) Genome-guided transcriptome assembly

[CLASS (2013)](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-14-S5-S14): (Contraint-based Local Assembly and Selection of Splice variants) a transcript selection scheme that takes into account **contiguity constrains** from read pairs and spliced reads and, where available, knowledge about gene structure (cDNA sequence databases). Do not estimate transcript abundance can be passed to RSEM, etc.

[StringTie (2015)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4643835/): applies a **network flow algorithm**, together with optional *de novo* transcriptome assembly. The reference-based uses alignments of reads to identify clusters of reads that represent potential transcripts. If paired-end reads, they improve the ability of the assembler to link together exons belonging to the same transcript.

[Strawberry (2017)](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005851): consists of two modules: assembly and quantification. The novelty is that the two modules use different optimization frameworks but utilize the same **data graph substructure**. The assembly module parses aligned reads into splicing graphs, and uses network flow algorithms. The quantification module corrects for sequencing bias through an EM algorithm.

[Scallop (2017)](https://www.nature.com/articles/nbt.4020): is a reference-based transcript assembler that improves reconstruction of of **multiexon** and **lowly expressed transcripts**. Scallop minimizes the read coverage deviation and minimizes the number of expressed transcripts by iteratively decomposing vertices of the splice graph.

## 2) *De novo* transcriptome assembly

[Trinity (2011)](https://www.nature.com/articles/nbt.1883): efficiently constructs and analyses sets of **de Brujin graphs**. Trinity fully reconstruction a large fraction of transcripts, transcripts from recently duplicated genes with a sensitivity similar to methods that rely on genome alignments.
