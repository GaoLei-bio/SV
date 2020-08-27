# By Lei Gao
# Final genome comparison SV list
#  $prefix.NR.bed
path_to_SV_calling_script=$1
reference=$2
query=$3
prefix=$4
cpuN=$5
Ref_base_pbsv_vcf=$6
Qry_base_pbsv_vcf=$7

input_SV=$prefix.NR.bed
my_within_sv=$prefix.my_anchor.variants_within_alignments.bed
anchor_file=$prefix.coords.tab

ref_pbsv=Ref_base.Indel
query_pbsv=Qry_base.Indel

#################################
# Prepare PacBio indel lists based on reference and query genomes, respectively

python $path_to_SV_calling_script/Extract_pbsv_indel.py    --vcf_file        $Qry_base_pbsv_vcf \
                                --ref_genome      $query     \
                		    	--query_genome    $reference  \
							    --Output          Qry_base >  Qry_base.log  &
# Use Heinz as ref genome
python $path_to_SV_calling_script/Extract_pbsv_indel.py    --vcf_file        $Ref_base_pbsv_vcf \
                                --ref_genome      $reference     \
                			    --query_genome    $query  \
							    --Output          Ref_base >  Ref_base.log  & 


wait


#################################
# Make blast database for each chr 
for chr in `cut -f 1 "$reference".fai`
do
	samtools faidx $reference $chr > $chr.fasta 
	samtools faidx $chr.fasta
	makeblastdb -dbtype nucl -in $chr.fasta
done

for chr in `cut -f 1 "$query".fai`
do
	samtools faidx $query $chr > $chr.fasta 
	samtools faidx $chr.fasta
	makeblastdb -dbtype nucl -in $chr.fasta
done


####################################################################
# Step 5. Combined pbsv
# The input pbsv indels have already filtered gap-related ones
# SV_combine_pbsv: support Substitution
echo "#####################################"
date 
echo "Step 5. Combined pbsv"


head -1  $input_SV > ${prefix}.temp
#tail -n +2 $input_SV | sort >> ${prefix}.temp
grep within_alignment   $input_SV >> ${prefix}.temp
grep between_alignments $input_SV | grep Forward >> ${prefix}.temp
grep between_alignments $input_SV | grep Reverse >> ${prefix}.temp
grep between_alignments $input_SV | grep -v Forward | grep -v Reverse |awk '{print $0"\tForward"}' >>  ${prefix}.temp

python $path_to_SV_calling_script/Split_file.py --Infile ${prefix}.temp --Split_Number $cpuN --Suffix_length 4 --Out_prefix Piece

for file in `ls Piece_????`
do

python $path_to_SV_calling_script/SV_combine_pbsv.py --SV_file      $file \
                          --anchor_file  $anchor_file \
                          --ref_genome   $reference \
						  --query_genome $query \
						  --ref_pbsv     $ref_pbsv \
						  --query_pbsv   $query_pbsv \
						  --gap_flank    50 \
					  	  --blast_flank  5000 \
						  --temp_dir     ${file}_comb      \
						  --prefix       ${file}  > $file.SV_combine_pbsv.log &

done

wait 

mv $ref_pbsv.bad_pbsv.Piece_0000   $ref_pbsv.bad_pbsv
mv $query_pbsv.bad_pbsv.Piece_0000 $query_pbsv.bad_pbsv

cat Piece_????.added_pbsv.tab > ${prefix}.filter.pass.added_pbsv
cat Piece_????.formated.tsv > ${prefix}.formated.tsv

rm Piece_????* -rf
rm *.bad_pbsv.*
# Keep all within SVs
tail -n +2 $my_within_sv  | cut -f 4 > Good_within.list 
python $path_to_SV_calling_script/SV_combine_summary.py  --bed_input ${prefix}.filter.pass.added_pbsv --tsv_out ${prefix}.formated.tsv  --Good_within  Good_within.list  --prefix ${prefix} 

# outputs:
#    ${prefix}.formated.NR.tsv
#    ${prefix}.conbine_pbsv.summary

####################################################################
# Step 6. Get the uncombined between SVs back 
echo "#####################################"
date 
echo "Step 6. Get the uncombined between SVs back"

python $path_to_SV_calling_script/Get_uncombined_between_SV_back.py --Raw_SVs      ${prefix}.filter.pass.added_pbsv  \
                                         --Combined_NR  ${prefix}.formated.NR.tsv \
										 --Out_prefix   ${prefix} 
# Output:  
#    SP2SL.Master_list.tsv     Final master list 
#    SP2SL.Check_comb.tsv      Check conbination results
#    SP2SL.Within_only.tsv        SVs Within alignments only
#    SP2SL.Redefined.tsv          Successifully redefined SVs between alignments
#    SP2SL.Raw_for_redefined.tsv  The rawSV for those successifully redefined
#    SP2SL.Non_redefined.tsv      The rawSV for those not redefined












