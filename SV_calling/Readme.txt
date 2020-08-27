####################################################################
SV calling scripts
By Lei Gao (lg397@cornell.edu or leigao@wbgcas.cn)
Fei lab (http://bioinfo.bti.cornell.edu/)
Version 1.1
Date: 07/04/2020
This pipeline is designed to detect structural variants (SVs) between high-quality genomes of two closely related species/varieties.
One is called Reference genome, the other is Query genome.

SV_calling.sh:
Step 1. Genome alignment by minimap2	
Step 2. Unique anchor/alignment detection
Step 3. SV calling based on genome comparison, and then filtering by Illumina read mapping 

SV_PacBio.sh (Optional): 
Step 4: Combining indels called by PacBio read mapping

####################################################################
How to run:
Please (1) copy all scripts to path_to_SV_calling_script (absolute path);
       (2) put all inputs in working directory;
       (3) install all dependencies

Step 1-3:
  SV_calling.sh
  Command:      
      bash path_to_SV_calling_script/SV_calling.sh \
           path_to_SV_calling_script \
           <Reference_genome_file>   \
           <Query_genome_file>  \
           <Prefix_for_outputs> \
           <number of threads>  \
           <min_SV_size>   \
           <max_SV_size>  \
           <assemblytics_path>
  Inputs:
     Reference_genome_file   Fasta format, see "Genome file format requirement" for detail
     Query_genome_file       Fasta format, see "Genome file format requirement" for detail
     QryRead2Ref.bam         Sorted bam file, Query Illumina read alignments to Reference genome
     RefRead2Qry.bam         Sorted bam file, Reference Illumina read alignments to Query genome
     Ref_self.bam            Sorted bam file, Reference Illumina read alignments to Reference genome
     Qry_self.bam            Sorted bam file, Query Illumina read alignments to Query genome
     Note: the four bam files should be named to the ones exactly shown (QryRead2Ref.bam, RefRead2Qry.bam, Ref_self.bam, Qry_self.bam) and these bam files and the corresponding index files should be put in working directory.

  For example:
      bash path_to_SV_calling_script/SV_calling.sh \
           path_to_SV_calling_script \
           SL4.0.genome.fasta \
           Pimp_v1.4.fasta \
           SP2SL 24 10 1000000 \
           path_to_assemblytics_scripts

  Final result for this example (Same SVs in 2 formats):
      SP2SL.Genome_comp_SV.tsv
      SP2SL.NR.bed (Assemblytics output format. This one is the input file for next step)

======================
Step 4:
  SV_PacBio.sh
  Command:
      bash path_to_SV_calling_script/SV_PacBio.sh \
           path_to_SV_calling_script \
           <Reference_genome_file>   \
           <Query_genome_file>  \
           <Prefix_for_outputs> \
           <number of threads> \
           <Ref_base_pbsv_vcf> \
           <Qry_base_pbsv_vcf>

  Inputs:
     Prefix_for_outputs.NR.bed      Result of SV_calling.sh
     Reference_genome_file          Fasta format
     Query_genome_file              Fasta format
     Ref_base_pbsv_vcf              vcf file. SV calling based on Query PacBio read alignments to Reference genome by pbsv (https://github.com/PacificBiosciences/pbsv)
     Qry_base_pbsv_vcf              vcf file. SV calling based on Reference PacBio read alignments to Query genome by pbsv (https://github.com/PacificBiosciences/pbsv)
 
For example:
      bash path_to_SV_calling_script/SV_PacBio.sh \
           path_to_SV_calling_script \
           SL4.0.genome.fasta  \
           Pimp_v1.4.fasta   \
           SP2SL  24  \
           PimpReads2SL4.0.var.vcf  \
           HeinzReads2Pimp.var.vcf

Final result for this example:
      SP2SL.Master_list.tsv

####################################################################
Dependencies:
   1. minimap2 (v2.11 or higher, https://github.com/lh3/minimap2)
   2. sam2delta.py (A python script from RaGOO package, https://github.com/malonge/RaGOO/raw/master/sam2delta.py)
   3. Assemblytics (Assemblytics, https://github.com/MariaNattestad/Assemblytics)
   4. samtools (v1.5 or higher, http://www.htslib.org/)
   5. blast+ 
   6. Python 2.7

####################################################################
Genome file format requirement:
    1. fasta
    2. Chromosome name should be look like xxxch01, xxxch02 et al. 
    3. The contigs not anchored into pseudochromosome can be merged as xxxch00.
       For example, SL4.0ch00, SL4.0ch01, SL4.0ch02...SL4.0ch12 in Heinz 1706 genome;
                    PIMPch00, SPIMPch01, SPIMPch02...SPIMPch12 in LA2093 genome.
    4. SV calling will run on normal chromosomes first and then ch00.

####################################################################
Note:
  1. SV_calling.sh
     1.1 This pipeline calls SVs based on the anchors identified by Assemblytics (https://pubmed.ncbi.nlm.nih.gov/27318204/).
     1.2 Because the genome assembly may be incomplete, some unassembled regions might be identified as deletion. To solve this issue, this pipeline uses Illumina reads to filter detected SVs.
         To do this, please align the Illumina reads to the two reference genomes, respectively, by any read mapping program (e.g. bwa), and named the bam file as instructed in “How to run” section.

   2. SV_PacBio.sh
     If you have PacBio reads for the two genomes, you may call SVs using these reads, and combined the resulted SVs into the above genome comparison SVs. 
     Please put the vcf files based on the two genomes under working directory. A sample vcf file is provided (Example.vcf).


