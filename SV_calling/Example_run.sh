# By Lei Gao
# genome
ln -s path_to_genome/SL4.0.genome.fasta .
ln -s path_to_genome/Pimp_v1.4.fasta .
# illumina bam
ln -s path_to_bam/LA2093_to_Heinz.bam      QryRead2Ref.bam           # Align query reads to reference genome
ln -s path_to_bam/LA2093_to_Heinz.bam.bai  QryRead2Ref.bam.bai
ln -s path_to_bam/Heinz_to_Heinz.bam       Ref_self.bam
ln -s path_to_bam/Heinz_to_Heinz.bam.bai   Ref_self.bam.bai

ln -s path_to_bam/Heinz_to_LA2093.bam      RefRead2Qry.bam               # Align reference reads to query genome
ln -s path_to_bam/Heinz_to_LA2093.bam.bai  RefRead2Qry.bam.bai
ln -s path_to_bam/LA2093_to_LA2093.bam     Qry_self.bam
ln -s path_to_bam/LA2093_to_LA2093.bam.bai Qry_self.bam.bai
# pbsv vcf
ln path_to_pbsv/PimpReads2SL4.0.var.vcf .
ln path_to_pbsv/HeinzReads2Pimp.var.vcf .

# path_to_SV_calling_script
# please change it to the real path 
path_to_SV_calling_script=xxx    

# SV calling by genome comparison
      bash path_to_SV_calling_script/SV_calling.sh path_to_SV_calling_script SL4.0.genome.fasta Pimp_v1.4.fasta   SP2SL 24 10 1000000 /data/lgVB/miniconda3/envs/py27/libexec/assemblytics
# If no PacBio SVs, just stop here.

# Combine PacBio SVs
      bash path_to_SV_calling_script/SV_PacBio.sh  path_to_SV_calling_script SL4.0.genome.fasta  Pimp_v1.4.fasta   SP2SL 24  PimpReads2SL4.0.var.vcf  HeinzReads2Pimp.var.vcf






