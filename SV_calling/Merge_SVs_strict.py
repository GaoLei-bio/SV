"""
By Lei Gao
usage: Merge_SVs_strict.py [-h] --Within WITHIN --Bad_within BAD_WITHIN
                           --My_between MY_BETWEEN --Raw_between RAW_BETWEEN
                           --Raw_unique_anchor RAW_UNIQUE_ANCHOR
                           --My_same_chr_anchor MY_SAME_CHR_ANCHOR
                           --ref_genome REF_GENOME --query_genome QUERY_GENOME
                           --Input_bam INPUT_BAM --Rev_bam REV_BAM --Reads_bam
                           READS_BAM --Rev_Reads_bam REV_READS_BAM
                           --Ref_self_bam REF_SELF_BAM --Query_self_bam
                           QUERY_SELF_BAM --Output_SV OUTPUT_SV

optional arguments:
  -h, --help            show this help message and exit
  --Within WITHIN       SVs within alignments. $prefix.a$Anchor.b$Between.vari
                        ants_within_alignments.bed
  --Bad_within BAD_WITHIN
                        SVs within the anchors only unique in Assemblytics
                        outputs. $prefix.a$Anchor.raw_uniqye_only.variants_wit
                        hin_alignments.bed
  --My_between MY_BETWEEN
                        SVs between alignments by my scripts. $prefix.a$Anchor
                        .b$Between.my_same_chr_between_anchor.bed
  --Raw_between RAW_BETWEEN
                        SVs between alignments by Assemblytics. $prefix.a$Anch
                        or.b$Between.variants_between_alignments.bed
  --Raw_unique_anchor RAW_UNIQUE_ANCHOR
                        Unique anchors found by Assemblytics.
                        $prefix.a$Anchor.coords.tab
  --My_same_chr_anchor MY_SAME_CHR_ANCHOR
                        The anchors used for my between SV calling.
  --ref_genome REF_GENOME
                        Reference genome.
  --query_genome QUERY_GENOME
                        Query genome.
  --Input_bam INPUT_BAM
                        Bam file. Align Query to Reference genome.
  --Rev_bam REV_BAM     Bam file. Align Reference to Query genome.
  --Reads_bam READS_BAM
                        Bam file. Align Query reads to Reference genome.
  --Rev_Reads_bam REV_READS_BAM
                        Bam file. Align Reference reads to Query genome.
  --Ref_self_bam REF_SELF_BAM
                        Bam file. Align Reference reads to Reference genome.
  --Query_self_bam QUERY_SELF_BAM
                        Bam file. Align Query reads to Query genome.
  --Output_SV OUTPUT_SV
                        Merged non-redundant SV list.

"""
import argparse
from pathlib import Path
from operator import itemgetter
import re
import os
import subprocess
parser = argparse.ArgumentParser()
parser.add_argument("--Within", type=str, help="SVs within alignments. $prefix.a$Anchor.b$Between.variants_within_alignments.bed", required=True, default="")
parser.add_argument("--Bad_within", type=str, help="SVs within the anchors only unique in Assemblytics outputs. $prefix.a$Anchor.raw_uniqye_only.variants_within_alignments.bed", required=True, default="")
parser.add_argument("--My_between", type=str, help="SVs between alignments by my scripts. $prefix.a$Anchor.b$Between.my_same_chr_between_anchor.bed", required=True, default="")
parser.add_argument("--Raw_between", type=str, help="SVs between alignments by Assemblytics. $prefix.a$Anchor.b$Between.variants_between_alignments.bed", required=True, default="")
parser.add_argument("--Raw_unique_anchor", type=str, help="Unique anchors found by Assemblytics. $prefix.a$Anchor.coords.tab", required=True, default="")
parser.add_argument("--My_same_chr_anchor", type=str, help="The anchors used for my between SV calling.", required=True, default="")
parser.add_argument("--ref_genome", type=str, help="Reference genome.", required=True, default="")
parser.add_argument("--query_genome", type=str, help="Query genome.", required=True, default="")
parser.add_argument("--Input_bam", type=str, help="Bam file. Align Query to Reference genome.", required=True, default="")
parser.add_argument("--Rev_bam", type=str, help="Bam file. Align Reference to Query  genome.", required=True, default="")
parser.add_argument("--Reads_bam", type=str, help="Bam file. Align Query reads to Reference genome.", required=True, default="")
parser.add_argument("--Rev_Reads_bam", type=str, help="Bam file. Align Reference reads to Query  genome.", required=True, default="")
parser.add_argument("--Ref_self_bam", type=str, help="Bam file. Align Reference reads to Reference  genome.", required=True, default="")
parser.add_argument("--Query_self_bam", type=str, help="Bam file. Align Query reads to Query  genome.", required=True, default="")


#parser.add_argument("--PBSV_indel_2_blast_dir", type=str, help="Directory including blastouts for flanking region check for PBSV_indel_2.", required=True, default="")

parser.add_argument("--Output_SV", type=str, help="Merged non-redundant SV list.", required=True, default="")

args = parser.parse_args()

Within = args.Within
Bad_within = args.Bad_within
My_between = args.My_between
Raw_between = args.Raw_between
Raw_unique_anchor = args.Raw_unique_anchor
My_same_chr_anchor = args.My_same_chr_anchor
ref_genome = args.ref_genome
query_genome = args.query_genome
ref_genome = args.ref_genome
query_genome = args.query_genome
Input_bam = args.Input_bam
Rev_bam = args.Rev_bam
Reads_bam = args.Reads_bam
Rev_Reads_bam = args.Rev_Reads_bam
Ref_self_bam = args.Ref_self_bam
Query_self_bam = args.Query_self_bam
#PBSV_indel_2_blast_dir = args.PBSV_indel_2_blast_dir
Output_SV = args.Output_SV

''' Function '''
def convert_int(fields):
    for i in range(0,len(fields)):
        if fields[i].isdigit() or (fields[i].startswith("-") and fields[i][1:].isdigit()):
            fields[i] = int(fields[i])
    return fields

def get_qry_loc(cells):
    Start = min(cells[5],cells[6])
    End   = max(cells[5],cells[6])
    Size = End - Start + 1
    return Start, End, Size

def get_ref_loc(cells):
    Start = min(cells[1],cells[2])
    End   = max(cells[1],cells[2])
    Size = End - Start + 1
    return Start, End, Size

def add_to_dict(Key,Value,Dict):
    if Key not in Dict:
        Dict[Key] = []
        Dict[Key].append(Value)
    else:
        Dict[Key].append(Value)

def keep_NR_sv(NR_by_ref_coord,NR_by_qry_coord,cells):
    ref_key = "_".join(cells[:3])
    qry_key = cells[9][:-2]
    NR_by_ref_coord[ref_key] = cells
    NR_by_qry_coord[qry_key] = cells

def keep_NR_on_chr(NR_by_ref_coord,NR_by_qry_coord,cells):
    ref_chr = cells[0]
    qry_chr = cells[9].split(":")[0]
    add_to_dict(ref_chr,cells,NR_by_ref_coord)
    add_to_dict(qry_chr,cells,NR_by_qry_coord)

def gap_count(Chr,start,stop,genome):
    region = Chr + ":" + str(min(start,stop)) + "-" + str(max(start,stop))
    #Gap = "".join(subprocess.check_output("samtools faidx " + genome + " " +  region, shell=True).split("\n")[1:-1]).upper().count("N")
    gap_seq = "".join(subprocess.check_output("samtools faidx " + genome + " " +  region, shell=True).split("\n")[1:-1]).upper()
    Gap = gap_seq.count("N")
    gap_to_sta = "nan"
    gap_to_end = "nan"
    if Gap > 0:
        # check gap start and end
        Sta = min(start,stop) - 1
        End = max(start,stop) + 1
        for x in gap_seq:
            Sta += 1
            if x == "N":
                gap_sta = Sta
                break
        for x in gap_seq[::-1]:
            End -= 1
            if x == "N":
                gap_end = End
                break
        gap_to_sta = gap_sta - min(start,stop)
        gap_to_end = max(start,stop) - gap_end
    return Gap, gap_to_sta, gap_to_end

def cov_count(Chr,Sta,End,Bam,reads_bam,Gap_size):
    Size = (End - Sta) - Gap_size
    region = Chr + ":" + str(Sta+1) + "-" + str(End-1)
    if Size > 0:
        Cov1 = float(subprocess.check_output("samtools depth " + Bam + " -r " +  region + "| awk '$3>0' | wc -l", shell=True)[:-1])/Size
        Cov2 = float(subprocess.check_output("samtools depth " + reads_bam + " -r " +  region + "| awk '$3>4' | wc -l", shell=True)[:-1])/Size
        return max(Cov1,Cov2)
    else:
        return 1

def self_cov_count(Chr,Sta,End,reads_bam,Gap_size):
    Size = (End - Sta) - Gap_size
    region = Chr + ":" + str(Sta+1) + "-" + str(End-1)
    if Size > 0:
        Cov2 = float(subprocess.check_output("samtools depth " + reads_bam + " -r " +  region + "| awk '$3>4' | wc -l", shell=True)[:-1])/Size
        return Cov2
    else:
        return 1

def overlap_check(Key,NR_by_chr,Genome):
    coord = Key.replace("_","\t").replace(":","\t").replace("-","\t").split("\t")
    Chr = coord[0]
    Over = False
    if Chr in NR_by_chr:
        Sta = min(int(coord[1]), int(coord[2]))
        End = max(int(coord[1]), int(coord[2]))

        if Genome == "Ref":
            for cells in NR_by_chr[Chr]:
                A_sta = int(cells[1])
                A_end = int(cells[2])
                if not (A_end < Sta or A_sta > End):
                    Over = True
                    break
        else:
            for cells in NR_by_chr[Chr]:
                q = cells[9].replace(":","\t").replace("-","\t").split("\t")
                A_sta = int(q[1])
                A_end = int(q[2])
                if not (A_end < Sta or A_sta > End):
                    Over = True
                    break

    return Over

''' Step 0.0 Get Raw uniqe anchors '''
raw_anchor_dict = {}
with open(Raw_unique_anchor) as infile:
    for line in infile:
        cells = line.strip().split("\t")
        anchor_key = "_".join([cells[6],cells[0],cells[1],cells[7],cells[2],cells[3]])
        raw_anchor_dict[anchor_key] = cells

''' Step 0.1 Get My uniqe anchors '''
my_anchor_dict = {}
i = -1
with open(My_same_chr_anchor) as infile:
    for line in infile:
        i += 1
        if i > 0:
            cells = line.strip().split("\t")
            anchor_key = "_".join([cells[2],cells[3],cells[4],cells[6],cells[7],cells[8]])
            my_anchor_dict[anchor_key] = cells


''' Step 1.0 Add my SVs between alignments '''
# NR within for next step
NR_by_ref_coord = {}
NR_by_qry_coord = {}
NR_by_ref_chr = {}
NR_by_qry_chr = {}

out_head = open(Within).readline().strip()

i = -1
with open(My_between) as infile, open(Output_SV, "w") as outfile:
    for line in infile:
        i += 1
        if i == 0:
            outfile.write(out_head + "\n")
        else:
            cells = line.strip().split("\t")
            ref_key = "_".join(cells[:3])
            qry_key = cells[9][:-2]

            NR_by_ref_coord[ref_key] = cells
            NR_by_qry_coord[qry_key] = cells
            outfile.write("\t".join(cells) + "\n")
            keep_NR_on_chr(NR_by_ref_chr,NR_by_qry_chr,cells)

''' Step 3.0 Add raw SVs between alignments after (1) size and (2) gap checks'''
i = -1
keep_by_ref_coord = {}
keep_by_qry_coord = {}
with open(Raw_between) as infile:
    for line in infile:
        i += 1
        if i > 0:
            cells = line.strip().split("\t")
            ref_key = "_".join(cells[:3])
            qry_key = cells[9][:-2]
            # not in previous list
            # not covered by previous list
            if ref_key in NR_by_ref_coord and qry_key in NR_by_qry_coord:
                print "both_cov_by_NR\t" + line.strip()
            elif ref_key in NR_by_ref_coord:
                print "Ref_cov_by_NR\t" + line.strip()
            elif qry_key in NR_by_qry_coord:
                print "Qry_cov_by_NR\t" + line.strip()
            elif overlap_check(ref_key,NR_by_ref_chr,"Ref") and overlap_check(qry_key,NR_by_qry_chr,"Qry"):
                print "both_Over_by_NR\t" + line.strip()
            elif overlap_check(ref_key,NR_by_ref_chr,"Ref"):
                print "Ref_Over_by_NR\t" + line.strip()
            elif overlap_check(qry_key,NR_by_qry_chr,"Qry"):
                print "Qry_Over_by_NR\t" + line.strip()
            else:
                # size check
                ref_gap_size = int(cells[7])
                query_gap_size = int(cells[8])
                sv_size = query_gap_size - ref_gap_size

                if ref_gap_size > 0 and ref_gap_size - 9 > query_gap_size and (-3 < query_gap_size < 3 or ref_gap_size / float(abs(query_gap_size)) >= 1000):
                    type = "Deletion"
                elif query_gap_size > 0 and ref_gap_size < query_gap_size - 9 and (-3 < ref_gap_size < 3 or query_gap_size / float(abs(ref_gap_size)) >= 1000):
                    type = "Insertion"
                elif ref_gap_size > query_gap_size > 0 and ref_gap_size > 9:
                    type = "Repeat_contraction"
                    if abs(query_gap_size - ref_gap_size) < 10 or abs(query_gap_size - ref_gap_size) / float(min(abs(query_gap_size),abs(query_gap_size))) < 0.05:
                        type = "Substitution"
                elif 0 < ref_gap_size < query_gap_size and query_gap_size > 9:
                    type = "Repeat_expansion"
                    if abs(query_gap_size - ref_gap_size) < 10 or abs(query_gap_size - ref_gap_size) / float(min(abs(query_gap_size),abs(query_gap_size))) < 0.05:
                        type = "Substitution"
                elif query_gap_size < 0 and query_gap_size + 9 < ref_gap_size:
                    type = "Tandem_contraction"
                    if abs(query_gap_size - ref_gap_size) / float(abs(query_gap_size)) < 0.25:
                        type = "error"
                elif ref_gap_size < 0 and ref_gap_size + 9 < query_gap_size:
                    type = "Tandem_expansion"
                    if abs(query_gap_size - ref_gap_size) / float(abs(ref_gap_size)) < 0.25:
                        type = "error"
                elif ref_gap_size == query_gap_size > 9:
                    type = "Substitution"
                #elif ref_gap_size == query_gap_size < 0:
                #    type = "Duplication"
                else:
                    type = "error"

                if type == "error":
                    print "Bad_Between_SV_calling\t" + line.strip()
                else:
                    # gap check
                    ref_Chr = cells[0]
                    ref_start = int(cells[1])
                    ref_stop = int(cells[2])
                    qry_coord = cells[9].replace("-",":").split(":")
                    qry_Chr = qry_coord[0]
                    qry_start = int(qry_coord[1])
                    qry_stop = int(qry_coord[2])

                    Ref_Gap = 0
                    Qry_Gap = 0
                    Ref_gap_to_Sta = 'nan'
                    Ref_gap_to_End = 'nan'
                    Qry_gap_to_Sta = 'nan'
                    Qry_gap_to_End = 'nan'

                    gap_check = "No"
                    Ref_Gap,Ref_gap_to_Sta,Ref_gap_to_End = gap_count(ref_Chr,ref_start,ref_stop,ref_genome)
                    if (abs(ref_gap_size) > 0 and Ref_Gap / float(abs(ref_gap_size)) > 0.3) or ("ch00" in ref_Chr and Ref_Gap > 0):
                        gap_check = "Yes"
                        print "Ref_Gap\t" + line.strip() + "\t" + str(Ref_Gap)

                    if gap_check == "No":
                        Qry_Gap,Qry_gap_to_Sta,Qry_gap_to_End = gap_count(qry_Chr,qry_start,qry_stop,query_genome)
                        if (abs(query_gap_size) > 0 and Qry_Gap / float(abs(query_gap_size)) > 0.3) or ("ch00" in qry_Chr and Qry_Gap > 0):
                            gap_check = "Yes"
                            print "Qry_Gap\t" + line.strip() + "\t" + str(Qry_Gap)

                    if gap_check == "No":
                        cells[6] = type
                        Ref_cov = 0
                        Qry_cov = 0

                        if ref_gap_size > 1:
                            Ref_cov = cov_count(ref_Chr,ref_start,ref_stop,Input_bam,Reads_bam,Ref_Gap)
                        if query_gap_size > 1:
                            Qry_cov = cov_count(qry_Chr,qry_start,qry_stop,Rev_bam,Rev_Reads_bam,Qry_Gap)

                        if ((type == "Deletion" or type == "Tandem_contraction") and Ref_cov > 0.75) or ((type == "Insertion" or type == "Tandem_expansion") and Qry_cov > 0.75) or \
                           ((type == "Repeat_contraction" or type == "Repeat_expansion" or type == "Substitution") and Ref_cov > 0.75 and Qry_cov > 0.75):
                            # bad SVs
                            Outlist = cells + [Ref_cov, Qry_cov, "bad"]
                            #outfile.write("\t".join(map(str,Outlist)) + "\n")
                            print "Raw_Between_SV_Covered_by_other_anchors\t" + "\t".join(map(str,Outlist))
                        else:
                            Ref_self_cov = 0
                            Qry_self_cov = 0
                            if ref_gap_size > 1:
                                Ref_self_cov = self_cov_count(ref_Chr,ref_start,ref_stop,Ref_self_bam,Ref_Gap)
                            if query_gap_size > 1:
                                Qry_self_cov = self_cov_count(qry_Chr,qry_start,qry_stop,Query_self_bam,Qry_Gap)

                            if ((type == "Deletion" or type == "Tandem_contraction") and (Ref_self_cov - Ref_cov < 0.3 or (Ref_cov > 0 and Ref_self_cov / Ref_cov < 1.5))) \
                            or ((type == "Insertion" or type == "Tandem_expansion") and (Qry_self_cov - Qry_cov < 0.3 or (Qry_cov > 0 and Qry_self_cov / Qry_cov < 1.5))) or \
                               ((type == "Repeat_contraction" or type == "Repeat_expansion" or type == "Substitution") and (Ref_self_cov - Ref_cov < 0.3 or (Ref_cov > 0 and Ref_self_cov / Ref_cov < 1.5)) and (Qry_self_cov - Qry_cov < 0.3 or (Qry_cov > 0 and Qry_self_cov / Qry_cov < 1.5))):
                                # bad SVs
                                Outlist = cells + [Ref_cov, Qry_cov, "bad", Ref_self_cov, Qry_self_cov]
                                #outfile.write("\t".join(map(str,Outlist)) + "\n")
                                print "No_sig_diff\t" + "\t".join(map(str,Outlist))
                            else:
                                add_to_dict(ref_key,cells,keep_by_ref_coord)
                                add_to_dict(qry_key,cells,keep_by_qry_coord)



# keep the SVs with unique ref_key and qry_key
used_qry_keys = set()
used_ref_keys = set()
with open(Output_SV, "a+") as outfile:
    # From ref
    for ref_key in sorted(keep_by_ref_coord.keys()):
        used_ref_keys.add(ref_key)
        for cells in keep_by_ref_coord[ref_key]:
            qry_key = cells[9][:-2]
            if qry_key not in used_qry_keys:
                used_qry_keys.add(qry_key)
                #outfile.write("\t".join(cells) + "\n")

                NR_by_ref_coord[ref_key] = cells
                NR_by_qry_coord[qry_key] = cells
                outfile.write("\t".join(cells) + "\n")
                keep_NR_on_chr(NR_by_ref_chr,NR_by_qry_chr,cells)
                break




''' Step 2.0 Filter SVs Within alignments '''
# keep by keys
keep_by_ref_coord = {}
keep_by_qry_coord = {}
i = -1
with open(Within) as infile, open(Output_SV, "a+") as outfile:
    for line in infile:
        i += 1
        if i > 0:
            cells = line.strip().split("\t")
            ref_key = "_".join(cells[:3])
            qry_key = cells[9][:-2]

            if ref_key in NR_by_ref_coord and qry_key in NR_by_qry_coord:
                print "both_cov_by_MySV\t" + line.strip()
            elif ref_key in NR_by_ref_coord:
                print "Ref_cov_by_MySV\t" + line.strip()
            elif qry_key in NR_by_qry_coord:
                print "Qry_cov_by_MySV\t" + line.strip()
            else:
                NR_by_ref_coord[ref_key] = cells
                NR_by_qry_coord[qry_key] = cells
                outfile.write("\t".join(cells) + "\n")
                keep_NR_on_chr(NR_by_ref_chr,NR_by_qry_chr,cells)

my_within_count = i + 1

''' Step 2.1 Filter SVs Within bad alignments '''
# keep by keys
keep_by_ref_coord = {}
keep_by_qry_coord = {}
i = -1
with open(Bad_within) as infile, open(Output_SV, "a+") as outfile:
    for line in infile:
        i += 1
        if i > 0:
            cells = line.strip().split("\t")
            ref_key = "_".join(cells[:3])
            qry_key = cells[9][:-2]
            cells[3] = "IND_w_" + str(int(cells[3].replace("IND_w_","")) + my_within_count + 1000)
            if ref_key in NR_by_ref_coord and qry_key in NR_by_qry_coord:
                print "both_cov_by_MySV\t" + line.strip()
            elif ref_key in NR_by_ref_coord:
                print "Ref_cov_by_MySV\t" + line.strip()
            elif qry_key in NR_by_qry_coord:
                print "Qry_cov_by_MySV\t" + line.strip()
            elif overlap_check(ref_key,NR_by_ref_chr,"Ref") and overlap_check(qry_key,NR_by_qry_chr,"Qry"):
                print "both_Over_by_MySV\t" + line.strip()
            elif overlap_check(ref_key,NR_by_ref_chr,"Ref"):
                print "Ref_Over_by_MySV\t" + line.strip()
            elif overlap_check(qry_key,NR_by_qry_chr,"Qry"):
                print "Qry_Over_by_MySV\t" + line.strip()
            else:
                add_to_dict(ref_key,cells,keep_by_ref_coord)

    for ref_key in keep_by_ref_coord.keys():
        if len(keep_by_ref_coord[ref_key]) == 1:
            cells = keep_by_ref_coord[ref_key][0]
            qry_key = cells[9][:-2]
            add_to_dict(qry_key,cells,keep_by_qry_coord)
        else:
            found_in_my_unique = "No"
            for cells in keep_by_ref_coord[ref_key]:
                anchor_key = "_".join([cells[18],cells[12],cells[13],cells[19],cells[14],cells[15]])
                if anchor_key in my_anchor_dict:
                    found_in_my_unique = "Yes"
                    qry_key = cells[9][:-2]
                    add_to_dict(qry_key,cells,keep_by_qry_coord)
                    break
            if found_in_my_unique == "No":
                # not in my anchor, just randomly keep one
                cells = keep_by_ref_coord[ref_key][0]
                qry_key = cells[9][:-2]
                add_to_dict(qry_key,cells,keep_by_qry_coord)

    # NR within for next step
    for qry_key in sorted(keep_by_qry_coord.keys()):
        if len(keep_by_qry_coord[qry_key]) == 1:
            outfile.write("\t".join(keep_by_qry_coord[qry_key][0]) + "\n")
            keep_NR_sv(NR_by_ref_coord,NR_by_qry_coord,keep_by_qry_coord[qry_key][0])
            keep_NR_on_chr(NR_by_ref_chr,NR_by_qry_chr,keep_by_qry_coord[qry_key][0])
        else:
            found_in_my_unique = "No"
            for cells in keep_by_qry_coord[qry_key]:
                anchor_key = "_".join([cells[18],cells[12],cells[13],cells[19],cells[14],cells[15]])
                if anchor_key in my_anchor_dict:
                    found_in_my_unique = "Yes"
                    outfile.write("\t".join(cells) + "\n")
                    keep_NR_sv(NR_by_ref_coord,NR_by_qry_coord,cells)
                    keep_NR_on_chr(NR_by_ref_chr,NR_by_qry_chr,cells)
                    break
            if found_in_my_unique == "No":
                # not in my anchor, just randomly keep one
                outfile.write("\t".join(keep_by_qry_coord[qry_key][0]) + "\n")
                keep_NR_sv(NR_by_ref_coord,NR_by_qry_coord,keep_by_qry_coord[qry_key][0])
                keep_NR_on_chr(NR_by_ref_chr,NR_by_qry_chr,keep_by_qry_coord[qry_key][0])



#srf
