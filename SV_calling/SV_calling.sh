# By Lei Gao
path_to_SV_calling_script=$1
reference=$2
query=$3
prefix=$4
cpuN=$5
min_SV_size=$6
max_SV_size=$7
assemblytics_path=$8
Anchor=10000     # Anchor length. Please see detail of this parameter in Assemblytics paper (http://www.ncbi.nlm.nih.gov/pubmed/27318204)

wget https://github.com/malonge/RaGOO/raw/master/sam2delta.py

####################################################################
# Step 0. Genome files
date 
echo "Now, prepare genome files..."

# samtools index
samtools faidx $query
samtools faidx $reference


####################################################################
# Step 1. Genome alignment
echo "Step 1, genome alignment"
# Align query genome to reference genome
date
echo -e "Reference genome: " $reference
echo -e "Query genome: " $query
minimap2 -ax asm5 --cs -t$cpuN  $reference $query > Qry2Ref.sam
samtools view -bS Qry2Ref.sam | samtools sort --threads $cpuN -o Qry2Ref.bam /dev/stdin
samtools index Qry2Ref.bam
# Align reference genome to query genome
date
echo -e "Reference genome: " $query
echo -e "Query genome: " $reference
minimap2 -ax asm5 --cs -t$cpuN  $query $reference > Ref2Qry.sam
samtools view -bS Ref2Qry.sam | samtools sort --threads $cpuN -o Ref2Qry.bam /dev/stdin
samtools index Ref2Qry.bam

####################################################################
# Step 2. Unique anchor detection
# Step 2.1 Convert sam file to delta file
#          sam2delta.py is a python script from RaGOO
date 
echo "Step 2, get unique alignments"
echo "Step 2.1, convert sam file to delta file"
python sam2delta.py  Qry2Ref.sam   # Output Qry2Ref.sam.delta

#Assemblytics Qry2Ref.sam.delta $prefix $Anchor $min_SV_size $max_SV_size

# Step 2.2 Get unique alignments as anchors for SV calling
echo "Step 2.1, Step 2.2 Get unique alignments as anchors for SV calling"
python $assemblytics_path/Assemblytics_uniq_anchor.py --delta Qry2Ref.sam.delta \
                            --unique-length $Anchor \
							--out $prefix \
							--keep-small-uniques 
# Outputs:
#       $prefix.coords.csv    All anchors               My script for between_anchors_SVs will base on this list; I will get the Synteny blocks on same chromsomes
#  		$prefix.coords.tab    Unique anchors only       My script for within_anchor_SVs will base on this list; I will also get back the REPETITIVE anchors within good Synteny blocks


####################################################################
# Step 3. Call between_anchors_SVs based on my Synteny blocks

# Step 3.1 Sort anchor list 
python $path_to_SV_calling_script/Sort_anchor.py $prefix.coords.csv > $prefix.Sorted_anchor.csv

# Step 3.2 Good anchors in Synteny blocks; INV, BND and DUP calling
python $path_to_SV_calling_script/Synteny_block.py --Input_coords $prefix.Sorted_anchor.csv \
                          --Prefix  $prefix  
# Output:
#     Ox2KoV7.a10000.Unique_anchor.sam_chr.Synteny_check.megaSV_flanking       SV anchors on same chromosomes
#     Ox2KoV7.a10000.My_SV.tsv                                                 INV, DUP and BND 



# Step 3.2 Call SVs between good anchors within and between Synteny blocks
#          1. Keep synteny block	
#          2. If some anchors locate between the processed anchors/synteny blocks, get them back	
#             1) if only one, just get back
#             2) if more than one, just merge them for between SV calling;

# My between anchor SVs
Start_ID=`cat $prefix.Sorted_anchor.csv | wc -l`

python $path_to_SV_calling_script/Call_SV_between_anchors.py --Input_coords   $prefix.Sorted_anchor.csv  \
                                  --Input_bam      Qry2Ref.bam \
								  --Reads_bam      QryRead2Ref.bam \
								  --Rev_bam        Ref2Qry.bam \
								  --Rev_Reads_bam  RefRead2Qry.bam\
								  --Ref_self_bam   Ref_self.bam   \
								  --Query_self_bam Qry_self.bam  \
                                  --Input_block    $prefix.Unique_anchor.sam_chr.Synteny_check.megaSV_flanking \
								  --Max_anchor_distance  $max_SV_size \
								  --Start_ID       $Start_ID \
								  --ref_genome     $reference \
								  --query_genome   $query  \
								  --Prefix  $prefix > $prefix.Call_SV_between_anchors.log


# Rescue REPETITIVE anchors within good Synteny blocks for assemblitics SV calling between anchors
cp $prefix.coords.tab  $prefix.rskRep.coords.tab
python $path_to_SV_calling_script/Rescue_repetitive_anchor.py --Initial_block  $prefix.Unique_anchor.sam_chr.Synteny_check.megaSV_flanking \
                                  --Final_block  $prefix.sam_chr.Synteny.final \
								  --Raw_Unique   $prefix.coords.tab >> $prefix.rskRep.coords.tab

echo -e "#reference\tref_start\tref_stop\tID\tsize\tstrand\ttype\tref_gap_size\tquery_gap_size\tquery_coordinates\tmethod" > $prefix.variants_between_alignments.bed
perl $assemblytics_path/Assemblytics_between_alignments.pl $prefix.rskRep.coords.tab 1 $max_SV_size all-chromosomes exclude-longrange bed >> $prefix.variants_between_alignments.bed 

# Example output:
#      Ox2KoV7.a10000.my_same_chr_between_anchor.bed         My between anchor SVs
#      Ox2KoV7.a10000.sam_chr.Synteny.final                  My Synteny anchors on same chromosomes


# Get anchors for SV calling within anchors
python $path_to_SV_calling_script/Get_anchor_for_within_SV.py --Initial_block  $prefix.Unique_anchor.sam_chr.Synteny_check.megaSV_flanking \
                                   --Final_block  $prefix.sam_chr.Synteny.final \
                                   --Raw_Unique   $prefix.coords.tab \
								   --Prefix       $prefix

# run my indel calling script
# SVs on my uniue achoors on same chromosomes
# Output $prefix.my_anchor.variants_within_alignments.bed
python $path_to_SV_calling_script/Call_Indel_within_alignment.py --sam_file    Qry2Ref.sam \
                              --ref_genome   $reference \
							  --query_genome $query \
							  --uniqe_anchor $prefix.my_anchor.tab \
							  --min_size     10  \
							  --max_size     100000 \
							  --prefix       $prefix.my_anchor  >  $prefix.my_anchor.SV_within_anchor.log  

# SVs on the uniqe anchors only in assemblitics output
# Output $prefix.raw_unique_only.variants_within_alignments.bed
python $path_to_SV_calling_script/Call_Indel_within_alignment.py --sam_file    Qry2Ref.sam \
                              --ref_genome   $reference \
							  --query_genome $query \
							  --uniqe_anchor $prefix.raw_unique_only.tab \
							  --min_size     10  \
							  --max_size     100000 \
							  --prefix       $prefix.raw_unique_only  >  $prefix.raw_unique_only.SV_within_anchor.log
							  
							  

#############################################
# Merge different list
# strict overlaping check between My_between and Raw_between

# Merge SVs 
#   1. $prefix.my_same_chr_between_anchor.bed, keep it
#   2. $prefix.variants_between_alignments.bed:
#        remove redundancy;
#        quality check
#   3. $prefix.my_anchor.variants_within_alignments.bed, remove redundancy
#   4. $prefix.raw_unique_only.variants_within_alignments.bed, remove redundancy

python $path_to_SV_calling_script/Merge_SVs_strict.py  --Within       $prefix.my_anchor.variants_within_alignments.bed  \
                     --Bad_within          $prefix.raw_unique_only.variants_within_alignments.bed \
                     --My_between          $prefix.my_same_chr_between_anchor.bed \
					 --Raw_between         $prefix.variants_between_alignments.bed \
					 --Raw_unique_anchor   $prefix.coords.tab  \
					 --My_same_chr_anchor  $prefix.sam_chr.Synteny.final \
					 --ref_genome          $reference \
					 --query_genome        $query  \
					 --Input_bam           Qry2Ref.bam \
					 --Rev_bam             Ref2Qry.bam \
					 --Reads_bam           QryRead2Ref.bam \
					 --Rev_Reads_bam       RefRead2Qry.bam \
					 --Ref_self_bam        Ref_self.bam  \
					 --Query_self_bam      Qry_self.bam  \
					 --Output_SV           $prefix.NR.bed  > $prefix.NR.log 

python $path_to_SV_calling_script/Formating_rawSVs.py $prefix.NR.bed > $prefix.Genome_comp_SV.tsv





