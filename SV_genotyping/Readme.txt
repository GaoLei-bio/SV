####################################################################
SV calling scripts
By Lei Gao (lg397@cornell.edu; leigao@wbgcas.cn)
Fei lab (http://bioinfo.bti.cornell.edu/)
Version 1.1
Date: 08/29/2020

This pipeline is designed to genotype a known list of structural variants (SVs) in one given sample based on its resequencing data.
The SVs are detected between high-quality genomes (Reference and Query genomes) of two closely related species/varieties with coordinates on both genomes.
For SV detection, please refer our SV_calling pipeline (https://github.com/GaoLei-bio/SV/tree/master/SV_calling). 

####################################################################
Dependencies:
   1. samtools (v1.5 or higher, http://www.htslib.org/)
   2. Python 2.7

####################################################################
How to run:
The resequencing reads should be aligned to Reference and Query genomes, respectively, by any read mapping program (e.g. bwa) and then sort the bam file.
To run this pipeline, we assume you already have the read alignments files, i.e. two sorted bam files for one sample.

Parameters:
INPUT=Example_SV.tsv          # Path to your input SV file. A example of input SV file is provided.
sample=Sample_name            # Sample name, prefix of outputs
Reference=Reference_genome    # Path to Reference genome, fasta format, indexed by samtools faidx
Query=Query_genome            # Path to Query genome , fasta format, indexed by samtools faidx
Ref_bam=Reference_base_bam    # Sorted bam file with reads aligned on Reference genome
Qry_bam=Query_base_bam        # Sorted bam file with reads aligned on Query genome
Mismath=3                     # Allowed mismath percentage of aligned read


Command:

python SV_genotyping.py $INPUT $sample $Reference $Query $Ref_bam $Qry_bam $Mismath

Outputs:
    Sample_name.GT.txt        # SV genotyping result

	  

