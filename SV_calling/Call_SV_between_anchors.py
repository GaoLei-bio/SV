"""
By Lei Gao
usage: Call_SV_between_anchors.py [-h] --Input_coords INPUT_COORDS --Input_bam
                                  INPUT_BAM --Reads_bam READS_BAM --Rev_bam
                                  REV_BAM --Rev_Reads_bam REV_READS_BAM
                                  --Ref_self_bam REF_SELF_BAM --Query_self_bam
                                  QUERY_SELF_BAM --Input_block INPUT_BLOCK
                                  --Max_anchor_distance MAX_ANCHOR_DISTANCE
                                  --Start_ID START_ID --ref_genome REF_GENOME
                                  --query_genome QUERY_GENOME --Prefix PREFIX

optional arguments:
  -h, --help            show this help message and exit
  --Input_coords INPUT_COORDS
                        Oriented coords from Assemblytics
  --Input_bam INPUT_BAM
                        Bam file. Align Query genome to Reference genome.
  --Reads_bam READS_BAM
                        Bam file. Align Query reads to Reference genome.
  --Rev_bam REV_BAM     Bam file. Align Reference genome to Query genome.
  --Rev_Reads_bam REV_READS_BAM
                        Bam file. Align Reference reads to Query genome.
  --Ref_self_bam REF_SELF_BAM
                        Bam file. Align Reference reads to Reference genome.
  --Query_self_bam QUERY_SELF_BAM
                        Bam file. Align Query reads to Query genome.
  --Input_block INPUT_BLOCK
                        Synteny blocks on same chromosomes
  --Max_anchor_distance MAX_ANCHOR_DISTANCE
                        Max distance between 2 anchors allowed for SV calling
  --Start_ID START_ID   The start ID for my SVs
  --ref_genome REF_GENOME
                        Reference genome.
  --query_genome QUERY_GENOME
                        Query genome.
  --Prefix PREFIX       Prefix for output files


"""
import argparse
from pathlib import Path
from operator import itemgetter
import re
import os
import subprocess
parser = argparse.ArgumentParser()
parser.add_argument("--Input_coords", type=str, help="Oriented coords from Assemblytics", required=True, default="")
parser.add_argument("--Input_bam", type=str, help="Bam file. Align Query genome to Reference genome.", required=True, default="")
parser.add_argument("--Reads_bam", type=str, help="Bam file. Align Query reads to Reference genome.", required=True, default="")
parser.add_argument("--Rev_bam", type=str, help="Bam file. Align Reference genome to Query  genome.", required=True, default="")
parser.add_argument("--Rev_Reads_bam", type=str, help="Bam file. Align Reference reads to Query  genome.", required=True, default="")
parser.add_argument("--Ref_self_bam", type=str, help="Bam file. Align Reference reads to Reference  genome.", required=True, default="")
parser.add_argument("--Query_self_bam", type=str, help="Bam file. Align Query reads to Query  genome.", required=True, default="")
parser.add_argument("--Input_block", type=str, help="Synteny blocks on same chromosomes", required=True, default="")
parser.add_argument("--Max_anchor_distance", type=str, help="Max distance between 2 anchors allowed for SV calling", required=True, default="")
parser.add_argument("--Start_ID", type=str, help="The start ID for my SVs", required=True, default="")
parser.add_argument("--ref_genome", type=str, help="Reference genome.", required=True, default="")
parser.add_argument("--query_genome", type=str, help="Query genome.", required=True, default="")

#parser.add_argument("--PBSV_indel_2_blast_dir", type=str, help="Directory including blastouts for flanking region check for PBSV_indel_2.", required=True, default="")

parser.add_argument("--Prefix", type=str, help="Prefix for output files.", required=True, default="")

args = parser.parse_args()

Input_coords = args.Input_coords
Input_bam = args.Input_bam
Reads_bam = args.Reads_bam
Rev_bam = args.Rev_bam
Rev_Reads_bam = args.Rev_Reads_bam
Ref_self_bam = args.Ref_self_bam
Query_self_bam = args.Query_self_bam
Input_block = args.Input_block
Max_anchor_distance = int(args.Max_anchor_distance)
ID = int(args.Start_ID)
ref_genome = args.ref_genome
query_genome = args.query_genome
#PBSV_indel_2_blast_dir = args.PBSV_indel_2_blast_dir
Prefix = args.Prefix

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

def best_candidates(Prev_index,This_index,raw_tab_dict,Direction):
    candidate_index = []
    if Direction == "Forward":
        Left = raw_tab_dict[Prev_index]
        Right = raw_tab_dict[This_index]
    else:
        Left = raw_tab_dict[This_index]
        Right = raw_tab_dict[Prev_index]
    Left_Start, Left_End, Left_Size = get_qry_loc(Left)
    Right_Start, Right_End, Right_Size = get_qry_loc(Right)

    Left_R_Start, Left_R_End, Left_R_Size = get_ref_loc(Left)
    Right_R_Start, Right_R_End, Right_R_Size = get_ref_loc(Right)


    for index in range(Prev_index+1,This_index):
        This_anchor = raw_tab_dict[index]
        This_Start, This_End, This_Size = get_qry_loc(This_anchor)
        This_R_Start, This_R_End, This_R_Size = get_ref_loc(This_anchor)
        if Left_End <= This_Start < This_End <= Right_Start and Left_R_End <= This_R_Start < This_R_End <= Right_R_Start:
            candidate_index.append(index)
    best_index = get_best_one(candidate_index,raw_tab_dict)
    return best_index

def get_best_one(candidate_index,raw_tab_dict):
    add_index = ""
    if len(candidate_index) == 0:
        add_index = ""
    elif len(candidate_index) == 1:
        add_index = candidate_index[0]
    elif len(candidate_index) >= 1:
        # Get best one
        minGap = 100
        for index in candidate_index:
            Q_Start, Q_end, Q_size = get_qry_loc(raw_tab_dict[index])
            R_Start, R_end, R_size = get_ref_loc(raw_tab_dict[index])
            gap = abs(Q_size - R_size) * 100.0 / min(Q_size,R_size)
            if gap < minGap:
                add_index = index
    return add_index

def ok_candidates(Prev_index,This_index,raw_tab_dict,Direction):
    candidate_index = []
    if Direction == "Forward":
        Left = raw_tab_dict[Prev_index]
        Right = raw_tab_dict[This_index]
    else:
        Left = raw_tab_dict[This_index]
        Right = raw_tab_dict[Prev_index]
    Left_Start, Left_End, Left_Size = get_qry_loc(Left)
    Right_Start, Right_End, Right_Size = get_qry_loc(Right)
    L = Left_End - Left_Size * 0.25
    R = Right_Start + Right_Size * 0.25

    Left_R_Start, Left_R_End, Left_R_Size = get_ref_loc(Left)
    Right_R_Start, Right_R_End, Right_R_Size = get_ref_loc(Right)
    Lref = Left_R_End - Left_R_Size * 0.25
    Rref = Right_R_Start + Right_R_Size * 0.25

    for index in range(Prev_index+1,This_index):
        This_anchor = raw_tab_dict[index]
        This_Start, This_End, This_Size = get_qry_loc(This_anchor)
        This_R_Start, This_R_End, This_R_Size = get_ref_loc(This_anchor)
        if (L < This_Start < This_End < R and Left_End < This_End and This_Start < Right_Start) and (Lref < This_R_Start < This_R_End < Rref and Left_R_End < This_R_End and This_R_Start < Right_R_Start):
            if This_Start < Left_End and This_R_Start < Left_R_End:
                # if left overlap
                if Left_End - This_Start < This_Size * 0.25 and Left_R_End - This_R_Start < This_R_Size * 0.25:
                    candidate_index.append(index)
            elif This_End > Right_Start and This_R_End > Right_R_Start:
                # if right overlap
                if This_End - Right_Start < This_Size * 0.25 and This_R_End - Right_R_Start < This_R_Size * 0.25:
                    candidate_index.append(index)

    best_index = get_best_one(candidate_index,raw_tab_dict)
    return best_index



def try_rescue_anchor(clu_table):
    new_table = []
    first = clu_table[0]
    new_table.append(first)
    Direction = first[12]
    Prev_index = raw_key_dict["_".join(map(str, first[2:10]))]
    for cells in clu_table[1:]:
        This_index = raw_key_dict["_".join(map(str, cells[2:10]))]
        if This_index - Prev_index > 1:
            add_index = best_candidates(Prev_index,This_index,raw_tab_dict,Direction)
            if add_index != "":
                add_anchor = [cells[0],'best'] + raw_tab_dict[add_index] + ["NAN","NAN"] + cells[12:]
                new_table.append(add_anchor)
            else:
                add_index = ok_candidates(Prev_index,This_index,raw_tab_dict,Direction)
                if add_index != "":
                    add_anchor = [cells[0],'okay'] + raw_tab_dict[add_index] + ["NAN","NAN"] + cells[12:]
                    new_table.append(add_anchor)
        # keep anchor
        # next round
        new_table.append(cells)
        Prev_index = raw_key_dict["_".join(map(str, cells[2:10]))]
    return new_table

def get_anchor_coord(cells):
    # ref: ref_End always > ref_Sta
    ref_Sta = cells[1]
    ref_End = cells[2]
    # qry
    qry_Sta = cells[5]
    qry_End = cells[6]
    if qry_End > qry_Sta:
        Str = "plus"
    else:
        Str = "minus"

    return ref_Sta,ref_End,qry_Sta,qry_End,Str


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

def call_between_SV(SV_id,Left,Right,Direction,outfile):
    #["#reference", "ref_start", "ref_stop", "ID", "size", "strand", "type", "ref_gap_size", "query_gap_size", "query_coordinates", "method"]
    # anchor infor
    ref_Chr = Left[0]
    qry_Chr = Left[4]
    L_ref_Sta, L_ref_End, L_qry_Sta, L_qry_End, L_Str = get_anchor_coord(Left)
    R_ref_Sta, R_ref_End, R_qry_Sta, R_qry_End, R_Str = get_anchor_coord(Right)
    if not (L_Str == R_Str and (L_qry_Sta == R_qry_Sta or L_qry_End == R_qry_End)):
        # Call SV
        if Direction == "Forward":
            ref_start = L_ref_End
            ref_stop = R_ref_Sta
            qry_start = max(L_qry_Sta, L_qry_End)
            qry_stop = min(R_qry_Sta, R_qry_End)
        else:
            ref_start = L_ref_End
            ref_stop = R_ref_Sta
            qry_start = max(R_qry_Sta, R_qry_End)
            qry_stop = min(L_qry_Sta, L_qry_End)

        ref_gap_size = ref_stop - ref_start
        query_gap_size = qry_stop - qry_start
        sv_size = query_gap_size - ref_gap_size

        if query_gap_size < 0:
            query_coordinates = qry_Chr + ":" + str(qry_stop) + "-" + str(qry_start) + ":" + "-"
        else:
            query_coordinates = qry_Chr + ":" + str(qry_start) + "-" + str(qry_stop) + ":" + "+"

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
            Outlist = [ref_Chr, min(ref_start,ref_stop), max(ref_start,ref_stop), SV_id, abs(sv_size), "+", type, ref_gap_size, query_gap_size, query_coordinates, "between_alignments",Direction]
            print "\t".join(map(str,Outlist))
        else:
            Ref_Gap = 0
            Qry_Gap = 0
            Ref_gap_to_Sta = 'nan'
            Ref_gap_to_End = 'nan'
            Qry_gap_to_Sta = 'nan'
            Qry_gap_to_End = 'nan'
            gap_check = "No"
            if ref_gap_size > 0:
                Ref_Gap,Ref_gap_to_Sta,Ref_gap_to_End = gap_count(ref_Chr,ref_start,ref_stop,ref_genome)
                if Ref_Gap / float(ref_gap_size) > 0.3:
                    gap_check = "Yes"
            if gap_check == "No" and query_gap_size > 0:
                Qry_Gap,Qry_gap_to_Sta,Qry_gap_to_End = gap_count(qry_Chr,qry_start,qry_stop,query_genome)
                if Qry_Gap / float(query_gap_size) > 0.3:
                    gap_check = "Yes"
            if gap_check == "No":
                Ref_cov = 0
                Qry_cov = 0

                if ref_gap_size > 1:
                    Ref_cov = cov_count(ref_Chr,ref_start,ref_stop,Input_bam,Reads_bam,Ref_Gap)
                if query_gap_size > 1:
                    Qry_cov = cov_count(qry_Chr,qry_start,qry_stop,Rev_bam,Rev_Reads_bam,Qry_Gap)

                '''Outlist = [ref_Chr, min(ref_start,ref_stop), max(ref_start,ref_stop), SV_id, abs(sv_size), "+", type, ref_gap_size, query_gap_size, query_coordinates, "between_alignments", Direction, Ref_cov, Qry_cov]
                outfile.write("\t".join(map(str,Outlist)) + "\n")'''
                if ((type == "Deletion" or type == "Tandem_contraction") and Ref_cov > 0.75) or ((type == "Insertion" or type == "Tandem_expansion") and Qry_cov > 0.75) or \
                   ((type == "Repeat_contraction" or type == "Repeat_expansion" or type == "Substitution") and Ref_cov > 0.75 and Qry_cov > 0.75):
                    # bad SVs
                    Outlist = [ref_Chr, min(ref_start,ref_stop), max(ref_start,ref_stop), SV_id, abs(sv_size), "+", type, ref_gap_size, query_gap_size, query_coordinates, "between_alignments", Direction, Ref_cov, Qry_cov, "bad"]
                    #outfile.write("\t".join(map(str,Outlist)) + "\n")
                    print "Covered_by_other_anchors\t" + "\t".join(map(str,Outlist))
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
                        Outlist = [ref_Chr, min(ref_start,ref_stop), max(ref_start,ref_stop), SV_id, abs(sv_size), "+", type, ref_gap_size, query_gap_size, query_coordinates, "between_alignments", Direction, Ref_cov, Qry_cov, "bad", Ref_self_cov, Qry_self_cov]
                        #outfile.write("\t".join(map(str,Outlist)) + "\n")
                        print "No_sig_diff\t" + "\t".join(map(str,Outlist))
                    else:
                        Outlist = [ref_Chr, min(ref_start,ref_stop), max(ref_start,ref_stop), SV_id, abs(sv_size), "+", type, ref_gap_size, query_gap_size, query_coordinates, "between_alignments", Direction]
                        outfile.write("\t".join(map(str,Outlist)) + "\n")


def call_chrEnd_SV(SV_id,outfile,ref_Chr,ref_start,ref_stop,qry_Chr,qry_start,qry_stop):
    ref_gap_size = ref_stop - ref_start
    query_gap_size = qry_stop - qry_start
    sv_size = query_gap_size - ref_gap_size

    query_coordinates = qry_Chr + ":" + str(qry_start) + "-" + str(qry_stop) + ":" + "+"

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
    elif ref_gap_size == query_gap_size > 9:
        type = "Substitution"
    #elif ref_gap_size == query_gap_size < 0:
    #    type = "Duplication"
    else:
        type = "error"

    if type == "error":
        Outlist = [ref_Chr, min(ref_start,ref_stop), max(ref_start,ref_stop), SV_id, abs(sv_size), "+", type, ref_gap_size, query_gap_size, query_coordinates, "chrEnd", "Forward"]
        print "\t".join(map(str,Outlist))
    else:
        Ref_Gap = 0
        Qry_Gap = 0
        Ref_gap_to_Sta = 'nan'
        Ref_gap_to_End = 'nan'
        Qry_gap_to_Sta = 'nan'
        Qry_gap_to_End = 'nan'
        gap_check = "No"
        if ref_gap_size > 0:
            Ref_Gap,Ref_gap_to_Sta,Ref_gap_to_End = gap_count(ref_Chr,ref_start,ref_stop,ref_genome)
            if Ref_Gap / float(ref_gap_size) > 0.3:
                gap_check = "Yes"
        if gap_check == "No" and query_gap_size > 0:
            Qry_Gap,Qry_gap_to_Sta,Qry_gap_to_End = gap_count(qry_Chr,qry_start,qry_stop,query_genome)
            if Qry_Gap / float(query_gap_size) > 0.3:
                gap_check = "Yes"
        if gap_check == "No":
            Ref_cov = 0
            Qry_cov = 0

            if ref_gap_size > 1:
                Ref_cov = cov_count(ref_Chr,ref_start,ref_stop,Input_bam,Reads_bam,Ref_Gap)
            if query_gap_size > 1:
                Qry_cov = cov_count(qry_Chr,qry_start,qry_stop,Rev_bam,Rev_Reads_bam,Qry_Gap)

            '''Outlist = [ref_Chr, min(ref_start,ref_stop), max(ref_start,ref_stop), SV_id, abs(sv_size), "+", type, ref_gap_size, query_gap_size, query_coordinates, "between_alignments", Direction, Ref_cov, Qry_cov]
            outfile.write("\t".join(map(str,Outlist)) + "\n")'''
            if ((type == "Deletion" or type == "Tandem_contraction") and Ref_cov > 0.75) or ((type == "Insertion" or type == "Tandem_expansion") and Qry_cov > 0.75) or \
               ((type == "Repeat_contraction" or type == "Repeat_expansion" or type == "Substitution") and Ref_cov > 0.75 and Qry_cov > 0.75):
                # bad SVs
                Outlist = [ref_Chr, min(ref_start,ref_stop), max(ref_start,ref_stop), SV_id, abs(sv_size), "+", type, ref_gap_size, query_gap_size, query_coordinates, "chrEnd", "Forward", Ref_cov, Qry_cov, "bad"]
                #outfile.write("\t".join(map(str,Outlist)) + "\n")
                print "Covered_by_reads\t" + "\t".join(map(str,Outlist))
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
                    Outlist = [ref_Chr, min(ref_start,ref_stop), max(ref_start,ref_stop), SV_id, abs(sv_size), "+", type, ref_gap_size, query_gap_size, query_coordinates, "between_alignments", "Forward", Ref_cov, Qry_cov, "bad", Ref_self_cov, Qry_self_cov]
                    #outfile.write("\t".join(map(str,Outlist)) + "\n")
                    print "No_sig_diff\t" + "\t".join(map(str,Outlist))
                else:
                    Outlist = [ref_Chr, min(ref_start,ref_stop), max(ref_start,ref_stop), SV_id, abs(sv_size), "+", type, ref_gap_size, query_gap_size, query_coordinates, "between_alignments", "Forward"]
                    outfile.write("\t".join(map(str,Outlist)) + "\n")

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


def merge_cluster(clu_table):
    # [Ref_chr, Ref_start, Ref_end, Ref_size, Qry_chr, Qry_start, Qry_end, Qry_size]
    cells = clu_table[0]
    Ref_chr = cells[2]
    Ref_start = min([cells[3] , cells[4]])
    Ref_end = max([cells[3] , cells[4]])
    #Ref_size = Ref_end - Ref_start + 1
    Qry_chr = cells[6]
    Qry_start = min([cells[7] , cells[8]])
    Qry_end = max([cells[7] , cells[8]])
    #Qry_size = Qry_end - Qry_start + 1

    for cells in clu_table:
        if Ref_start > min([cells[3] , cells[4]]):
            Ref_start = min([cells[3] , cells[4]])
        if Ref_end < max([cells[3] , cells[4]]):
            Ref_end = max([cells[3] , cells[4]])
        if Qry_start > min([cells[7] , cells[8]]):
            Qry_start = min([cells[7] , cells[8]])
        if Qry_end < max([cells[7] , cells[8]]):
            Qry_end = max([cells[7] , cells[8]])

    Ref_size = Ref_end - Ref_start + 1
    Qry_size = Qry_end - Qry_start + 1

    return [Ref_chr, Ref_start, Ref_end, Ref_size, Qry_chr, Qry_start, Qry_end, Qry_size]

def add_to_dict(Key,Value,Dict):
    if Key not in Dict:
        Dict[Key] = []
        Dict[Key].append(Value)
    else:
        Dict[Key].append(Value)

def check_chr_cov(infor_list,Chr_coved_dict):
    # infor_list = [chr, sta, end]
    Chr = infor_list[0]
    if Chr not in Chr_coved_dict:
        Chr_coved_dict[Chr] = infor_list
    else:
        if infor_list[1] < Chr_coved_dict[Chr][1]:
            Chr_coved_dict[Chr][1] = infor_list[1]
        if infor_list[2] > Chr_coved_dict[Chr][2]:
            Chr_coved_dict[Chr][2] = infor_list[2]



''' Step 0.0 Get sorted anchors '''
raw_key_dict = {}
raw_tab_dict = {}


index = -1
with open(Input_coords) as input_file:
    for line in input_file:
        index += 1
        if index > 0:
            cells = line.strip().split("\t")
            out_cells = [cells[6],int(cells[0]),int(cells[1]),int(cells[1]) - int(cells[0]) + 1,cells[7],int(cells[2]),int(cells[3]), int(cells[3]) - int(cells[2]) + 1]
            key = "_".join(map(str, out_cells))
            raw_key_dict[key] = index
            raw_tab_dict[index] = out_cells


''' Step 1.0 Get Input_block '''

block_dict = {}
i = -1
n = 0    # count not synteny block
c = 0

with open(Input_block) as infile, open(Prefix + ".sam_chr.Synteny.clean", "w") as outfile:
    for line in infile:
        cells = convert_int(line.strip().split("\t"))
        i += 1
        if i == 0:
            outfile.write("Cluster\tReason\t" + line)
            Synteny_head = "Cluster\tReason\t" + line
        else:
            if cells[0][-4:] != cells[4][-4:]:
                print "Diff_chr\t" + line.strip()
            #elif cells[10] == "No":
                # not synteny region
            #    n += 1
            else:
                Chr = cells[0][-4:]
                B_index = cells[12]
                if Chr not in block_dict:
                    block_dict[Chr] = {}
                if B_index not in block_dict[Chr]:
                    block_dict[Chr][B_index] = []
                    block_dict[Chr][B_index].append(cells)
                else:
                    block_dict[Chr][B_index].append(cells)
    ''' Step 2.0 Get Check Input_block '''
    Cluster_dict = {}
    RM_anchors = {}
    for Chr in sorted(block_dict.keys()):
        Cluster_dict[Chr] = {}
        for B_index in sorted(block_dict[Chr].keys()):
            raw_table = block_dict[Chr][B_index]
            if block_dict[Chr][B_index][0][10] != "No":
                c += 1
                Cluster_dict[Chr][c] = []
                # add first anchor
                Cluster_dict[Chr][c].append([c,"Begin"] + raw_table[0])
                outfile.write("\t".join(map(str,[c,"Begin"] + raw_table[0])) + "\n")
            else:
                RM_anchors[c] = []
                RM_anchors[c].append([c,"Begin_RM"] + raw_table[0])
                outfile.write("\t".join(map(str,[c,"Begin_RM"] + raw_table[0])) + "\n")


            # keep first anchor coords
            Direct = raw_table[0][10]
            Prev_start, Prev_end, Prev_size = get_qry_loc(raw_table[0])

            for cells in raw_table[1:]:
                reason = ""
                This_start, This_end, This_size = get_qry_loc(cells)
                if Direct == "Forward":
                    # This anchor should be after the previous one
                    if This_start < Prev_start:
                        # break type 1; latter is former
                        c += 1
                        reason = "Transcloation"
                    elif This_start - Prev_end > Max_anchor_distance:
                        c += 1
                        reason = This_start - Prev_end
                elif  Direct ==  "Reverse":
                    # This anchor should be before the previous one
                    if This_start > Prev_start:
                        # break type 1; latter is former
                        c += 1
                        reason = "Transcloation"
                    elif Prev_start - This_end > Max_anchor_distance:
                        c += 1
                        reason = Prev_start - This_end
                else:
                    reason = "RM"
                if reason != "RM":
                    if c not in Cluster_dict[Chr]:
                        Cluster_dict[Chr][c] = []
                        Cluster_dict[Chr][c].append([c,reason] + cells)
                    else:
                        Cluster_dict[Chr][c].append([c,reason] + cells)
                else:
                    RM_anchors[c].append([c,reason] + cells)

                outfile.write("\t".join(map(str,[c,reason] + cells)) + "\n")
                # next round
                Prev_start, Prev_end, Prev_size = get_qry_loc(cells)

''' Step 3.0 Call SV in cluster '''
#ID = Start_ID
with open(Prefix + ".my_same_chr_between_anchor.bed", "w") as outfile, open(Prefix + ".sam_chr.Synteny.final", "w") as anchor_file:
    #outhead = ["#reference", "ref_start", "ref_stop", "ID", "size", "strand", "type", "ref_gap_size", "query_gap_size", "query_coordinates", "method", "Ref_Gap_check", "Qry_Gap_check", "Ref_gap_to_Sta", "Ref_gap_to_End", "Qry_gap_to_Sta", "Qry_gap_to_End"]
    outhead = ["#reference", "ref_start", "ref_stop", "ID", "size", "strand", "type", "ref_gap_size", "query_gap_size", "query_coordinates", "method", "Anchor_direction"]
    outfile.write("\t".join(outhead) + "\n")
    anchor_file.write(Synteny_head)
    for Chr in sorted(Cluster_dict.keys()):
        '''SV calling in cluster'''
        for c in sorted(Cluster_dict[Chr].keys()):
            clu_table = Cluster_dict[Chr][c]
            if len(clu_table) > 1:
                new_table = try_rescue_anchor(clu_table)
                while len(new_table) > len(clu_table):
                    clu_table = new_table
                    new_table = try_rescue_anchor(clu_table)
                Cluster_dict[Chr][c] = new_table
            for cells in Cluster_dict[Chr][c]:
                anchor_file.write("\t".join(map(str,cells)) + "\n")
            # SV calling in cluster
            clu_table = Cluster_dict[Chr][c]
            if len(clu_table) > 1:
                first = clu_table[0]
                new_table.append(first)
                Direction = first[12]
                Prev_index = raw_key_dict["_".join(map(str, first[2:10]))]
                for cells in clu_table[1:-1]:
                    This_index = raw_key_dict["_".join(map(str, cells[2:10]))]
                    ID += 1
                    call_between_SV("Assemblytics_b_" + str(ID),raw_tab_dict[Prev_index],raw_tab_dict[This_index],Direction,outfile)
                    #print "Assemblytics_b_" + str(ID) + "\tLeft\t" + "\t".join(map(str, raw_tab_dict[Prev_index])) + "\tRight\t" + "\t".join(map(str, raw_tab_dict[This_index]))
                    # Next round
                    Prev_index = raw_key_dict["_".join(map(str, cells[2:10]))]
    for Chr in sorted(Cluster_dict.keys()):
        '''SV calling between cluster'''
        c_list = Cluster_dict[Chr].keys()
        first_c = c_list[0]
        Prev_cluster = merge_cluster(Cluster_dict[Chr][first_c])

        for c in c_list[1:]:
            This_cluster = merge_cluster(Cluster_dict[Chr][c])
            # check direction and distance
            #print "Left\t" + "\t".join(map(str, Prev_cluster)) + "\tRight\t" + "\t".join(map(str, This_cluster))
            if c-1 not in RM_anchors:
                # no removed anchor between these 2 clusters
                if Prev_cluster[1] < This_cluster[1] and Prev_cluster[2] < This_cluster[2]  and Prev_cluster[5] < This_cluster[5] and Prev_cluster[6] < This_cluster[6]:
                    if -100 < This_cluster[5] - Prev_cluster[6] < Max_anchor_distance and -100 < This_cluster[1] - Prev_cluster[2] < Max_anchor_distance:
                        ID += 1
                        call_between_SV("Assemblytics_b_" + str(ID),Prev_cluster,This_cluster,"Forward",outfile)
                        print "Assemblytics_b_" + str(ID) + "\t" + str(c - 1) + "\t" + str(c) + "\tLeft\t" + "\t".join(map(str, Prev_cluster)) + "\tRight\t" + "\t".join(map(str, This_cluster))
            else:
                # have removed anchor between these 2 clusters
                # check them
                get_back_anchors = []
                for cells in RM_anchors[c-1]:
                    Ref_chr = cells[2]
                    Ref_start = min([cells[3] , cells[4]])
                    Ref_end = max([cells[3] , cells[4]])
                    Qry_chr = cells[6]
                    Qry_start = min([cells[7] , cells[8]])
                    Qry_end = max([cells[7] , cells[8]])
                    # return [Ref_chr, Ref_start, Ref_end, Ref_size, Qry_chr, Qry_start, Qry_end, Qry_size]
                    if Prev_cluster[1] < Ref_start < This_cluster[1] and Prev_cluster[2] < Ref_end < This_cluster[2]  and Prev_cluster[5] < Qry_start < This_cluster[5] and Prev_cluster[6] < Qry_end < This_cluster[6]:
                        get_back_anchors.append(cells)
                if len(get_back_anchors) > 0:
                    Mid_cluster = merge_cluster(get_back_anchors)
                    # first half
                    if Prev_cluster[1] < Mid_cluster[1] and Prev_cluster[2] < Mid_cluster[2]  and Prev_cluster[5] < Mid_cluster[5] and Prev_cluster[6] < Mid_cluster[6]:
                        if -100 < Mid_cluster[5] - Prev_cluster[6] < Max_anchor_distance and -100 < Mid_cluster[1] - Prev_cluster[2] < Max_anchor_distance:
                            ID += 1
                            call_between_SV("Assemblytics_b_" + str(ID),Prev_cluster,Mid_cluster,"Forward",outfile)
                            print "Assemblytics_b_" + str(ID) + "\t" + str(c - 1) + "\t" + str(c) + "\tLeft\t" + "\t".join(map(str, Prev_cluster)) + "\tRight\t" + "\t".join(map(str, Mid_cluster))
                    # 2nd half
                    if Mid_cluster[1] < This_cluster[1] and Mid_cluster[2] < This_cluster[2]  and Mid_cluster[5] < This_cluster[5] and Mid_cluster[6] < This_cluster[6]:
                        if -100 < This_cluster[5] - Mid_cluster[6] < Max_anchor_distance and -100 < This_cluster[1] - Mid_cluster[2] < Max_anchor_distance:
                            ID += 1
                            call_between_SV("Assemblytics_b_" + str(ID),Mid_cluster,This_cluster,"Forward",outfile)
                            print "Assemblytics_b_" + str(ID) + "\t" + str(c - 1) + "\t" + str(c) + "\tLeft\t" + "\t".join(map(str, Mid_cluster)) + "\tRight\t" + "\t".join(map(str, This_cluster))
                else:
                    if Prev_cluster[1] < This_cluster[1] and Prev_cluster[2] < This_cluster[2]  and Prev_cluster[5] < This_cluster[5] and Prev_cluster[6] < This_cluster[6]:
                        if -100 < This_cluster[5] - Prev_cluster[6] < Max_anchor_distance and -100 < This_cluster[1] - Prev_cluster[2] < Max_anchor_distance:
                            ID += 1
                            call_between_SV("Assemblytics_b_" + str(ID),Prev_cluster,This_cluster,"Forward",outfile)
                            print "Assemblytics_b_" + str(ID) + "\t" + str(c - 1) + "\t" + str(c) + "\tLeft\t" + "\t".join(map(str, Prev_cluster)) + "\tRight\t" + "\t".join(map(str, This_cluster))

            # next round
            Prev_cluster = merge_cluster(Cluster_dict[Chr][c])


''' Step 4.0 SV on chromosome ends '''
Ref_size_dict = {}
Qry_size_dict = {}
Chr_coved_dict = {}

index = -1
with open(Input_coords) as input_file, open(Prefix + ".my_same_chr_between_anchor.bed", "a+") as outfile:
    for line in input_file:
        index += 1
        if index > 0:
            cells = convert_int(line.strip().split("\t"))
            ref_list = [cells[6],min(cells[0],cells[1]),max(cells[0],cells[1])]
            qry_list = [cells[7],min(cells[2],cells[3]),max(cells[2],cells[3])]
            check_chr_cov(ref_list,Chr_coved_dict)
            check_chr_cov(qry_list,Chr_coved_dict)
            if cells[6] not in Ref_size_dict:
                Ref_size_dict[cells[6]] = cells[4]
            if cells[7] not in Qry_size_dict:
                Qry_size_dict[cells[7]] = cells[5]
    for ref_chr in sorted(Ref_size_dict.keys()):
        if ref_chr[-4:] != "ch00":
            for qry_chr in sorted(Qry_size_dict.keys()):
                if ref_chr[-4:] == qry_chr[-4:]:
                    # infor_list = [chr, sta, end]
                    ref_coved = Chr_coved_dict[ref_chr]
                    qry_coved = Chr_coved_dict[qry_chr]

                    if not (ref_coved[1] < 100 and qry_coved[1] < 100):
                        # chr start
                        ref_start = 1
                        ref_stop = ref_coved[1]
                        qry_start = 1
                        qry_stop = qry_coved[1]
                        ID += 1
                        call_chrEnd_SV("Assemblytics_b_" + str(ID),outfile,ref_chr,ref_start,ref_stop,qry_chr,qry_start,qry_stop)
                    if not(Ref_size_dict[ref_chr] - ref_coved[2] < 100 and Qry_size_dict[qry_chr] - ref_coved[2] < 100):
                        # chr end
                        ref_start = ref_coved[2]
                        ref_stop = Ref_size_dict[ref_chr]
                        qry_start = qry_coved[2]
                        qry_stop = Qry_size_dict[qry_chr]
                        ID += 1
                        call_chrEnd_SV("Assemblytics_b_" + str(ID),outfile,ref_chr,ref_start,ref_stop,qry_chr,qry_start,qry_stop)




#srf
