"""
By Lei Gao
usage: SV_combine_pbsv.py [-h] --SV_file SV_FILE --anchor_file ANCHOR_FILE
                          --ref_genome REF_GENOME --query_genome QUERY_GENOME
                          --ref_pbsv REF_PBSV --query_pbsv QUERY_PBSV
                          --gap_flank GAP_FLANK --blast_flank BLAST_FLANK
                          --temp_dir TEMP_DIR --prefix PREFIX

optional arguments:
  -h, --help            show this help message and exit
  --SV_file SV_FILE     Filtered Assemblytics output SV file, bed format
  --anchor_file ANCHOR_FILE
                        Assemblytics sorted unique anchor file, tab format.
                        e.g.
  --ref_genome REF_GENOME
                        Reference genome file, fasta format. Indexed by
                        samtools
  --query_genome QUERY_GENOME
                        Query species genome file, fasta format. Indexed by
                        samtools
  --ref_pbsv REF_PBSV   pbsv calling based on reference genome. bed file.
  --query_pbsv QUERY_PBSV
                        pbsv calling based on query genome. bed file.
  --gap_flank GAP_FLANK
                        Require SVs are not very close to Gap. This is the
                        size of flank regions not allowed with gaps. For
                        within_alignment indels, this is the minimal distance
                        to anchor end.
  --blast_flank BLAST_FLANK
                        The size for flanking regions for blast check.
  --temp_dir TEMP_DIR   Temporary directory for flanking fasta seqs.
  --prefix PREFIX       Prefix for output files.


"""
import argparse
import subprocess
import re
from operator import itemgetter

parser = argparse.ArgumentParser()
parser.add_argument("--SV_file", type=str, help="Filtered Assemblytics output SV file, bed format", required=True, default="")
parser.add_argument("--anchor_file", type=str, help="Assemblytics sorted unique anchor file, tab format. e.g. ", required=True, default="")
parser.add_argument("--ref_genome", type=str, help="Reference genome file, fasta format. Indexed by samtools", required=True, default="")
parser.add_argument("--query_genome", type=str, help="Query species genome file, fasta format. Indexed by samtools", required=True, default="")
parser.add_argument("--ref_pbsv", type=str, help="pbsv calling based on reference genome. bed file.", required=True, default="")
parser.add_argument("--query_pbsv", type=str, help="pbsv calling based on query genome. bed file.", required=True, default="")
parser.add_argument("--gap_flank", type=str, help="Require SVs are not very close to Gap. This is the size of flank regions not allowed with gaps. For within_alignment indels, this is the minimal distance to anchor end.", required=True, default="")
parser.add_argument("--blast_flank", type=str, help="The size for flanking regions for blast check.", required=True, default="")
parser.add_argument("--temp_dir", type=str, help="Temporary directory for flanking fasta seqs.", required=True, default="")
parser.add_argument("--prefix", type=str, help="Prefix for output files.", required=True, default="")

args = parser.parse_args()

SV_file = args.SV_file
anchor_file = args.anchor_file
ref_genome = args.ref_genome
query_genome = args.query_genome
ref_pbsv = args.ref_pbsv
query_pbsv = args.query_pbsv
gap_flank = int(args.gap_flank)
blast_flank = int(args.blast_flank)
temp_dir = args.temp_dir
prefix = args.prefix

blast_parameter = "  -perc_identity 90 -dust no -num_threads 1 -outfmt '7 qseqid sseqid qlen slen length qstart qend sstart send sstrand pident nident mismatch gapopen gaps qcovhsp qcovs score bitscore evalue'  -evalue 1e-5 -gapopen 5 -gapextend 3 "
letters = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'

# if temp dir does not exist, creat it
file_list = set(subprocess.check_output("ls", shell=True).decode().split("\n"))
if temp_dir not in file_list:
    subprocess.check_output("mkdir " + temp_dir, shell=True)

''' Functions'''
def convert_int(fields):
    for i in range(0,len(fields)):
        if fields[i].isdigit() or (fields[i].startswith("-") and fields[i][1:].isdigit()):
            fields[i] = int(fields[i])
    return fields

def add_to_dict(key,value,dict):
    if key in dict:
        dict[key].append(value)
    else:
        dict[key] = []
        dict[key].append(value)

def gap_count(coord,genome):
    ext_seq = "".join(subprocess.check_output("samtools faidx " + genome + " " +  coord, shell=True).decode().split("\n")[1:-1]).upper()
    N_num = ext_seq.count("N")
    return N_num

def SV_2_AnchorEnd(ref_sta,ref_end,qry_sta,qry_end,anchor_info):
    anchor_ref_sta = anchor_info[0]
    anchor_ref_end = anchor_info[1]
    anchor_qry_sta = min(anchor_info[2],anchor_info[3])
    anchor_qry_end = max(anchor_info[2],anchor_info[3])
    distance2end = min( ref_sta - anchor_ref_sta, anchor_ref_end - ref_end, min(qry_sta,qry_end) - anchor_qry_sta, anchor_qry_end - max(qry_sta,qry_end))
    return distance2end

def get_pbsv_list(pbsv_in,pbsv_by_chr,pbsv_by_ID):
    with open(pbsv_in) as infile, open(pbsv_in + ".bad_pbsv." + SV_file, "w") as outfile:
        for line in infile:
            if line.startswith("SV_id"):
                outfile.write(line)
            else:
                cells = convert_int(line.strip().split("\t"))
                sv_id = cells[0]
                sv_type = cells[1]
                if sv_type not in pbsv_by_chr:
                    pbsv_by_chr[sv_type] = {}
                chrID = cells[3]
                ref_reads = int(cells[2].split(":")[1].split(",")[0])
                alt_reads = int(cells[2].split(":")[1].split(",")[1])
                if alt_reads > 3 and alt_reads * 100.0 / (ref_reads + alt_reads) >= 20:
                    add_to_dict(chrID,cells,pbsv_by_chr[sv_type])
                    pbsv_by_ID[sv_id] = cells
                else:
                    outfile.write(line)

def within_pbsv_check(cells,ref_pbsv_by_chr,query_pbsv_by_chr):
    SV_type = cells[6]
    SV_size = cells[4]
    Ref_chr = cells[0]
    Ref_start = cells[1]
    Ref_end = cells[2]
    Ref_Size = cells[7]
    Qry_chr = cells[9].split(":")[0]
    Qry_start = int(cells[9].split(":")[1].split("-")[0])
    Qry_end = int(cells[9].split(":")[1].split("-")[1])
    Qry_Size = cells[8]

    if SV_type.startswith("Deletion"):
        # DEL in query relative to ref
        # Ref_Size: DEL_size
        # Qry_Size: INS_Size
        Ref_pbsv_IDs = overlap_w_DEL(Ref_chr,Ref_start,Ref_end,Ref_Size,Qry_Size,ref_pbsv_by_chr)
        Qry_pbsv_IDs = overlap_w_INS(Qry_chr,Qry_start,Qry_end,Ref_Size,Qry_Size,query_pbsv_by_chr)
    else:
        # DEL in ref relative to query
        # Qry_Size: DEL_size
        # Ref_Size: INS_Size
        Ref_pbsv_IDs = overlap_w_INS(Ref_chr,Ref_start,Ref_end,Qry_Size,Ref_Size,ref_pbsv_by_chr)
        Qry_pbsv_IDs = overlap_w_DEL(Qry_chr,Qry_start,Qry_end,Qry_Size,Ref_Size,query_pbsv_by_chr)
    return Ref_pbsv_IDs,Qry_pbsv_IDs


def overlap_w_DEL(chrID,DEL_start,DEL_end,DEL_size,INS_Size,pbsv_by_chr):
    # pbsv information:
    # SV_id   Type    GT      Chromosome      Start   End     Ref_seq Query_seq       Ref_size        Query_size
    hit_pbsv_IDs = []
    DEL_size -= INS_Size
    for cells in pbsv_by_chr["DEL"][chrID]:
        if not (max(DEL_start,DEL_end) < min(cells[4],cells[5]) or min(DEL_start,DEL_end) > max(cells[4],cells[5])):
            overlap_size = len(set(range(DEL_start,DEL_end)).intersection(range(cells[4],cells[5])))
            pbsv_size = len(cells[6])
            if overlap_size * 100.0 / max(pbsv_size,DEL_size) > 50:
                # overlap >50%
                hit_pbsv_IDs.append(cells[0])
    if len(hit_pbsv_IDs) == 0:
        hit_pbsv_IDs = "None"
    else:
        hit_pbsv_IDs = ",".join(hit_pbsv_IDs)
    return hit_pbsv_IDs

def overlap_w_INS(chrID,INS_start,INS_end,DEL_size,INS_Size,pbsv_by_chr):
    # pbsv information:
    # SV_id   Type    GT      Chromosome      Start   End     Ref_seq Query_seq       Ref_size        Query_size
    hit_pbsv_IDs = []
    INS_Size = DEL_size - INS_Size
    for cells in pbsv_by_chr["INS"][chrID]:
        pbsv_pos = cells[4]
        if abs(pbsv_pos - INS_start) <= 2 or abs(pbsv_pos - INS_end) <= 2 :
            # +- 2 bp
            pbsv_size = len(cells[7])
            if abs(pbsv_size - INS_Size) * 100.0 / min(INS_Size,pbsv_size) <= 20:
                # size diff <= 20
                hit_pbsv_IDs.append(cells[0])
    if len(hit_pbsv_IDs) == 0:
        hit_pbsv_IDs = "None"
    else:
        hit_pbsv_IDs = ",".join(hit_pbsv_IDs)
    return hit_pbsv_IDs


def between_pbsv_check(cells,ref_pbsv_by_chr,query_pbsv_by_chr):
    SV_type = cells[6]
    SV_size = cells[4]
    Ref_chr = cells[0]
    Ref_start = cells[1]
    Ref_end = cells[2]
    Ref_Size = cells[7]
    Qry_chr = cells[9].split(":")[0]
    Qry_start = int(cells[9].split(":")[1].split("-")[0])
    Qry_end = int(cells[9].split(":")[1].split("-")[1])
    Qry_Size = cells[8]

    if SV_type == "Deletion":
        # DEL in query relative to ref
        # Ref_Size: DEL_size
        # Qry_Size: INS_Size
        Ref_pbsv_IDs = overlap_w_DEL(Ref_chr,Ref_start,Ref_end,Ref_Size,Qry_Size,ref_pbsv_by_chr)
        Qry_pbsv_IDs = overlap_b_INS(Qry_chr,Qry_start,Qry_end,Ref_Size,Qry_Size,query_pbsv_by_chr)
    elif SV_type == "Insertion":
        # DEL in ref relative to query
        # Qry_Size: DEL_size
        # Ref_Size: INS_Size
        Ref_pbsv_IDs = overlap_b_INS(Ref_chr,Ref_start,Ref_end,Qry_Size,Ref_Size,ref_pbsv_by_chr)
        Qry_pbsv_IDs = overlap_w_DEL(Qry_chr,Qry_start,Qry_end,Qry_Size,Ref_Size,query_pbsv_by_chr)
    elif SV_type == "Repeat_contraction":
        # DEL in query relative to ref
        # Ref_Size: DEL_size
        # Qry_Size: INS_Size
        Ref_pbsv_IDs = Repeat_sv_DEL(Ref_chr,Ref_start,Ref_end,Ref_Size,Qry_Size,ref_pbsv_by_chr)
        Qry_pbsv_IDs = Repeat_sv_INS(Qry_chr,Qry_start,Qry_end,Ref_Size,Qry_Size,query_pbsv_by_chr)
    elif SV_type == "Repeat_expansion":
        # DEL in ref relative to query
        # Qry_Size: DEL_size
        # Ref_Size: INS_Size
        Ref_pbsv_IDs = Repeat_sv_INS(Ref_chr,Ref_start,Ref_end,Qry_Size,Ref_Size,ref_pbsv_by_chr)
        Qry_pbsv_IDs = Repeat_sv_DEL(Qry_chr,Qry_start,Qry_end,Qry_Size,Ref_Size,query_pbsv_by_chr)
    elif SV_type == "Tandem_contraction":
        # DEL in query relative to ref
        # Ref_Size: DEL_size
        # Qry_Size: INS_Size
        Ref_pbsv_IDs = Repeat_sv_DEL(Ref_chr,Ref_start,Ref_end,Ref_Size,Qry_Size,ref_pbsv_by_chr)
        Qry_pbsv_IDs = Repeat_sv_INS(Qry_chr,Qry_start,Qry_end,Ref_Size,Qry_Size,query_pbsv_by_chr)
    elif SV_type == "Tandem_expansion":
        # DEL in ref relative to query
        # Qry_Size: DEL_size
        # Ref_Size: INS_Size
        Ref_pbsv_IDs = Repeat_sv_INS(Ref_chr,Ref_start,Ref_end,Qry_Size,Ref_Size,ref_pbsv_by_chr)
        Qry_pbsv_IDs = Repeat_sv_DEL(Qry_chr,Qry_start,Qry_end,Qry_Size,Ref_Size,query_pbsv_by_chr)
    else:
        # Substitution
        Ref_pbsv_IDs = Repeat_sv_DEL(Ref_chr,Ref_start,Ref_end,Ref_Size,Qry_Size,ref_pbsv_by_chr)
        Qry_pbsv_IDs = Repeat_sv_INS(Qry_chr,Qry_start,Qry_end,Ref_Size,Qry_Size,query_pbsv_by_chr)

    return Ref_pbsv_IDs,Qry_pbsv_IDs

def overlap_b_INS(chrID,INS_start,INS_end,DEL_size,INS_Size,pbsv_by_chr):
    # pbsv information:
    # SV_id   Type    GT      Chromosome      Start   End     Ref_seq Query_seq       Ref_size        Query_size
    hit_pbsv_IDs = []
    INS_Size = DEL_size - INS_Size
    if abs(DEL_size) < 10:
        # simple indel
        max_distance = 3
    else:
        max_distance = int(abs(DEL_size)/3)
    exp_ins_sites = set(range(INS_start-max_distance,INS_end+max_distance))

    for cells in pbsv_by_chr["INS"][chrID]:
        pbsv_pos = cells[4]
        if pbsv_pos in exp_ins_sites:
            pbsv_size = len(cells[7])
            if abs(pbsv_size - INS_Size) * 100.0 / min(INS_Size,pbsv_size) <= 20:
                # size diff <= 20
                hit_pbsv_IDs.append(cells[0])
    if len(hit_pbsv_IDs) == 0:
        hit_pbsv_IDs = "None"
    else:
        hit_pbsv_IDs = ",".join(hit_pbsv_IDs)
    return hit_pbsv_IDs

def Repeat_sv_DEL(chrID,DEL_start,DEL_end,DEL_size,INS_Size,pbsv_by_chr):
    # pbsv information:
    # SV_id   Type    GT      Chromosome      Start   End     Ref_seq Query_seq       Ref_size        Query_size
    hit_pbsv_IDs = []
    DEL_size -= INS_Size
    for cells in pbsv_by_chr["DEL"][chrID]:
        if not (max(DEL_start,DEL_end) < min(cells[4],cells[5]) or min(DEL_start,DEL_end) > max(cells[4],cells[5])):
            overlap_size = len(set(range(DEL_start,DEL_end)).intersection(range(cells[4],cells[5])))
            pbsv_size = len(cells[6])
            if overlap_size * 100.0 / pbsv_size > 50:
                # overlap >50% of pbsv_size only
                hit_pbsv_IDs.append(cells[0])
    if len(hit_pbsv_IDs) == 0:
        hit_pbsv_IDs = "None"
    else:
        hit_pbsv_IDs = ",".join(hit_pbsv_IDs)
    return hit_pbsv_IDs

def Repeat_sv_INS(chrID,INS_start,INS_end,DEL_size,INS_Size,pbsv_by_chr):
    # pbsv information:
    # SV_id   Type    GT      Chromosome      Start   End     Ref_seq Query_seq       Ref_size        Query_size
    hit_pbsv_IDs = []
    INS_Size = DEL_size - INS_Size       # insetion size
    max_distance = int(abs(INS_Size)/3)  # All DEL_size >= 50
    exp_ins_sites = set(range(INS_start-max_distance,INS_end+max_distance))

    for cells in pbsv_by_chr["INS"][chrID]:
        pbsv_pos = cells[4]
        if pbsv_pos in exp_ins_sites:
            pbsv_size = len(cells[7])
            if pbsv_size * 100.0 / INS_Size < 120:
                # not too large
                hit_pbsv_IDs.append(cells[0])
    if len(hit_pbsv_IDs) == 0:
        hit_pbsv_IDs = "None"
    else:
        hit_pbsv_IDs = ",".join(hit_pbsv_IDs)
    return hit_pbsv_IDs



def redefine_between_SV(cells,Ref_pbsv_IDs,Qry_pbsv_IDs,ref_pbsv_by_ID,query_pbsv_by_ID,formated_out,raw_prefix,new_prefix):
    SV_type = cells[6]
    SV_size = cells[4]
    Ref_chr = cells[0]
    Ref_start = cells[1]
    Ref_end = cells[2]
    Ref_Size = cells[7]
    Qry_chr = cells[9].split(":")[0]
    Qry_start = int(cells[9].split(":")[1].split("-")[0])
    Qry_end = int(cells[9].split(":")[1].split("-")[1])
    Qry_Size = cells[8]
    Assemblytics_ID = cells[3]
    Anchor_direction = cells[11]
    between_order = 0
    # simple indels directly output
    if (SV_type == "Deletion" and abs(Qry_Size) <= 3) or (SV_type == "Insertion" and abs(Ref_Size) <= 3):
        Final_ID = Assemblytics_ID.replace(raw_prefix,new_prefix)
        formated_list = [SV_type,Ref_chr,Ref_start,Ref_end,Ref_Size,Qry_chr,Qry_start,Qry_end,Qry_Size,Assemblytics_ID,Ref_pbsv_IDs,Qry_pbsv_IDs,Final_ID]
        formated_out.write("\t".join(map(str,formated_list)) + "\n")
    else:
        Ref_redefined = []
        Qry_redefined = []
        if Ref_pbsv_IDs != "None":
            # redefine by ref_base pbsv
            Ref_redefined = redefine_by_blast_flank(Ref_pbsv_IDs,ref_pbsv_by_ID,Anchor_direction,Ref_chr,Ref_start,Ref_end,Ref_Size,Qry_chr,Qry_start,Qry_end,Qry_Size)
            '''
            # check redefine results
            if Ref_redefined != []:
                for Tab in Ref_redefined:
                    print ("\t".join(map(str,cells)) + "\tRef_pbsv\t" + "\t".join(map(str,Tab)))
            '''
        if Qry_pbsv_IDs != "None":
            # redefine by qry_base pbsv
            Qry_redefined = redefine_by_blast_flank(Qry_pbsv_IDs,query_pbsv_by_ID,Anchor_direction,Qry_chr,Qry_start,Qry_end,Qry_Size,Ref_chr,Ref_start,Ref_end,Ref_Size)
            '''
            if Qry_redefined != []:
                for Tab in Qry_redefined:
                    print ("\t".join(map(str,cells)) + "\tQry_pbsv\t" + "\t".join(map(str,Tab)))
            '''

        if Ref_redefined != [] and Qry_redefined != []:
            # conbine redefine results
            combined_output_pair = {}
            conbined_qry_pbsv_ids = []
            for ref_tab in Ref_redefined:
                for qry_tab in Qry_redefined:
                    ref_base_Ref_sta = ref_tab[2]
                    ref_base_Ref_end = ref_tab[3]
                    ref_base_Ref_size = ref_base_Ref_end - ref_base_Ref_sta
                    ref_base_Qry_sta = ref_tab[5]
                    ref_base_Qry_end = ref_tab[6]
                    ref_base_Qry_size = ref_base_Qry_end - ref_base_Qry_sta

                    qry_base_Ref_sta = qry_tab[5]
                    qry_base_Ref_end = qry_tab[6]
                    qry_base_Ref_size = qry_base_Ref_end - qry_base_Ref_sta
                    qry_base_Qry_sta = qry_tab[2]
                    qry_base_Qry_end = qry_tab[3]
                    qry_base_Qry_size = qry_base_Qry_end - qry_base_Qry_sta
                    if ref_base_Ref_size > ref_base_Qry_size and \
                     not (max(qry_base_Ref_sta,qry_base_Ref_end) < min(ref_base_Ref_sta,ref_base_Ref_end) or min(qry_base_Ref_sta,qry_base_Ref_end) > max(ref_base_Ref_sta,ref_base_Ref_end) )  and  \
                     len(set(range(min(ref_base_Ref_sta,ref_base_Ref_end),max(ref_base_Ref_sta,ref_base_Ref_end)+1)).intersection(range(min(qry_base_Ref_sta,qry_base_Ref_end),max(qry_base_Ref_sta,qry_base_Ref_end)+1))) * 100.0 / min(ref_base_Ref_size,qry_base_Ref_size) > 50:
                        # Del, check Ref overlap
                        add_to_dict(ref_tab[0],qry_tab[0],combined_output_pair)
                        conbined_qry_pbsv_ids.append(qry_tab[0])
                    elif ref_base_Ref_size < ref_base_Qry_size and   \
                      not (max(qry_base_Qry_sta,qry_base_Qry_end) < min(ref_base_Qry_sta,ref_base_Qry_end) or min(qry_base_Qry_sta,qry_base_Qry_end) > max(ref_base_Qry_sta,ref_base_Qry_end) )  and  \
                      len(set(range(min(ref_base_Qry_sta,ref_base_Qry_end),max(ref_base_Qry_sta,ref_base_Qry_end)+1)).intersection(range(min(qry_base_Qry_sta,qry_base_Qry_end),max(qry_base_Qry_sta,qry_base_Qry_end)+1))) * 100.0 / min(ref_base_Qry_size,qry_base_Qry_size) > 50:
                        # ins, check Qry overlap
                        add_to_dict(ref_tab[0],qry_tab[0],combined_output_pair)
                        conbined_qry_pbsv_ids.append(qry_tab[0])
            if len(combined_output_pair) == 0:
                # can not combine
                # output each respectively
                between_order = output_Ref_redefined(Ref_redefined,cells,formated_out,between_order,raw_prefix,new_prefix)
                between_order = output_Qry_redefined(Qry_redefined,cells,formated_out,between_order,raw_prefix,new_prefix)
            else:
                # at least 1 can be combined
                left_Ref_redefined = []
                left_Qry_redefined = []
                for ref_tab in Ref_redefined:
                    ref_pbsv_id = ref_tab[0]
                    if ref_pbsv_id in combined_output_pair:
                        qry_pbsv_id = ",".join(combined_output_pair[ref_pbsv_id])
                        # [pbsv_id,Ref_chr,real_ref_sta,real_ref_end,Qry_chr,real_qry_sta,real_qry_end]
                        #Final_ID = Assemblytics_ID.replace("Assemblytics_","IND_")
                        #formated_list = [SV_type,Ref_chr,Ref_start,Ref_end,Ref_Size,Qry_chr,Qry_start,Qry_end,Qry_Size,Assemblytics_ID,Ref_pbsv_IDs,Qry_pbsv_IDs,Final_ID]
                        #formated_out.write("\t".join(map(str,formated_list)) + "\n")
                        real_ref_sta = ref_tab[2]
                        real_ref_end = ref_tab[3]
                        real_ref_size = real_ref_end - real_ref_sta
                        real_qry_sta = ref_tab[5]
                        real_qry_end = ref_tab[6]
                        real_qry_size = real_qry_end - real_qry_sta
                        # creat ID
                        if between_order < 52:
                            Final_ID = Assemblytics_ID.replace(raw_prefix,new_prefix) + letters[between_order]
                        else:
                            round = int(between_order / 52)
                            Final_ID = Assemblytics_ID.replace(raw_prefix,new_prefix) + letters[between_order] + str(round)
                        between_order += 1
                        # check type
                        SV_type = cells[6] + "/" + re.sub('\d*','',ref_pbsv_id.split(",")[0])
                        formated_list = [SV_type,Ref_chr,real_ref_sta,real_ref_end,real_ref_size,Qry_chr,real_qry_sta,real_qry_end,real_qry_size,Assemblytics_ID,ref_pbsv_id,qry_pbsv_id,Final_ID]
                        formated_out.write("\t".join(map(str,formated_list)) + "\n")
                    else:
                        left_Ref_redefined.append(ref_tab)

                for qry_tab in Qry_redefined:
                    if qry_tab[0] not in set(conbined_qry_pbsv_ids):
                        left_Qry_redefined.append(qry_tab)

                if left_Ref_redefined != []:
                    between_order = output_Ref_redefined(left_Ref_redefined,cells,formated_out,between_order,raw_prefix,new_prefix)
                if left_Qry_redefined != []:
                    between_order = output_Qry_redefined(left_Qry_redefined,cells,formated_out,between_order,raw_prefix,new_prefix)

        elif Ref_redefined != []:
            # redefine by ref_base pbsv is final results
            between_order = output_Ref_redefined(Ref_redefined,cells,formated_out,between_order,raw_prefix,new_prefix)
        elif Qry_redefined != []:
            # redefine by qry_base pbsv is final results
            between_order = output_Qry_redefined(Qry_redefined,cells,formated_out,between_order,raw_prefix,new_prefix)
        else:
            if SV_type == "Deletion" or SV_type == "Insertion":
                # Even if fail to redefine simple, output it
                Final_ID = Assemblytics_ID.replace(raw_prefix,new_prefix)
                formated_list = [SV_type,Ref_chr,Ref_start,Ref_end,Ref_Size,Qry_chr,Qry_start,Qry_end,Qry_Size,Assemblytics_ID,Ref_pbsv_IDs,Qry_pbsv_IDs,Final_ID]
                formated_out.write("\t".join(map(str,formated_list)) + "\n")

def output_Ref_redefined(Ref_redefined,cells,formated_out,between_order,raw_prefix,new_prefix):
    Ref_chr = cells[0]
    Qry_chr = cells[9].split(":")[0]
    Assemblytics_ID = cells[3]

    for ref_tab in Ref_redefined:
        SV_type = cells[6]
        qry_pbsv_id = "None"
        ref_pbsv_id = ref_tab[0]
        real_ref_sta = ref_tab[2]
        real_ref_end = ref_tab[3]
        real_ref_size = real_ref_end - real_ref_sta
        real_qry_sta = ref_tab[5]
        real_qry_end = ref_tab[6]
        real_qry_size = real_qry_end - real_qry_sta

        # creat ID
        if between_order < 52:
            Final_ID = Assemblytics_ID.replace(raw_prefix,new_prefix) + letters[between_order]
        else:
            round = int(between_order / 52)
            Final_ID = Assemblytics_ID.replace(raw_prefix,new_prefix) + letters[between_order] + str(round)
        between_order += 1

        # check type
        SV_type += "/" + re.sub('\d*','',ref_pbsv_id.split(",")[0])
        formated_list = [SV_type,Ref_chr,real_ref_sta,real_ref_end,real_ref_size,Qry_chr,real_qry_sta,real_qry_end,real_qry_size,Assemblytics_ID,ref_pbsv_id,qry_pbsv_id,Final_ID]
        formated_out.write("\t".join(map(str,formated_list)) + "\n")

    return between_order

def output_Qry_redefined(Qry_redefined,cells,formated_out,between_order,raw_prefix,new_prefix):
    Ref_chr = cells[0]
    Qry_chr = cells[9].split(":")[0]
    Assemblytics_ID = cells[3]

    for qry_tab in Qry_redefined:
        SV_type = cells[6]
        qry_pbsv_id = qry_tab[0]
        ref_pbsv_id = "None"

        real_ref_sta = qry_tab[5]
        real_ref_end = qry_tab[6]
        real_ref_size = real_ref_end - real_ref_sta
        real_qry_sta = qry_tab[2]
        real_qry_end = qry_tab[3]
        real_qry_size = real_qry_end - real_qry_sta

        # creat ID
        if between_order < 52:
            Final_ID = Assemblytics_ID.replace(raw_prefix,new_prefix) + letters[between_order]
        else:
            round = int(between_order / 52)
            Final_ID = Assemblytics_ID.replace(raw_prefix,new_prefix) + letters[between_order] + str(round)
        between_order += 1

        # check type
        if re.sub('\d*','',qry_pbsv_id.split(",")[0]) == "INS":
            SV_type += "/DEL"
        else:
            SV_type += "/INS"
        formated_list = [SV_type,Ref_chr,real_ref_sta,real_ref_end,real_ref_size,Qry_chr,real_qry_sta,real_qry_end,real_qry_size,Assemblytics_ID,ref_pbsv_id,qry_pbsv_id,Final_ID]
        formated_out.write("\t".join(map(str,formated_list)) + "\n")

    return between_order



def redefine_by_blast_flank(Ref_pbsv_IDs,ref_pbsv_by_ID,Anchor_direction,Ref_chr,Ref_start,Ref_end,Ref_Size,Qry_chr,Qry_start,Qry_end,Qry_Size):
    # redefine between sv by pbsv indels
    # pbsv information:
    # SV_id   Type    GT      Chromosome      Start   End     Ref_seq Query_seq       Ref_size        Query_size
    Ref_redefined = []
    for pbsv_id in Ref_pbsv_IDs.split(","):
        fields = ref_pbsv_by_ID[pbsv_id]
        pbsv_type = fields[1]
        Chromosome = fields[3]
        Start = fields[4]
        End = fields[5]
        exp_gap_size = len(fields[7])
        pbsv_size = max(len(fields[6]),len(fields[7]))
        if Start < Ref_start:
            extra_left = Ref_start - Start
        else:
            extra_left = 0
        if End > Ref_end:
            extra_right = End - Ref_end
        else:
            extra_right = 0
        exp_sites = set(range(min(Qry_start,Qry_end) - extra_left - 10,max(Qry_start,Qry_end) + extra_right + 11))
        " check left flank "
        if Start >= blast_flank:
            left_region = Chromosome + ":" + str(Start - (blast_flank-1)) + "-" + str(Start)
        elif Start > 1:
            left_region = Chromosome + ":" + "1" + "-" + str(Start)
        if Start > 1:
            left_seq = re.sub('.*N', '', "".join(subprocess.check_output("samtools faidx " + Ref_chr + ".fasta " +  left_region, shell=True).decode().split("\n")[1:-1]))
            if len(left_seq) >= gap_flank:
                with open(temp_dir + "/" + pbsv_id + ".left.fa", 'w') as fa_file:
                    fa_file.write(">" + pbsv_id + "_L_" + left_region + "\n")
                    fa_file.write(left_seq + "\n")
                # blast
                left_blast_output = subprocess.check_output("blastn -query " + temp_dir + "/" + pbsv_id + ".left.fa " + " -db " +  Qry_chr + ".fasta " + blast_parameter, shell=True).decode().split("\n")[:-1]
                left_blast = get_blast_hits(left_blast_output,Anchor_direction,Qry_chr)
                ''' if left blast success, check right flank '''
                if len(left_blast) > 0:
                    right_region = Chromosome + ":" + str(End+1) + "-" + str(End+blast_flank)
                    right_seq = re.sub('N.*', '', "".join(subprocess.check_output("samtools faidx " + Ref_chr + ".fasta " +  right_region, shell=True).decode().split("\n")[1:-1]))
                    if len(right_seq) >= gap_flank:
                        with open(temp_dir + "/" + pbsv_id + ".right.fa", 'w') as fa_file:
                            fa_file.write(">" + pbsv_id + "_R_" + right_region +  "\n")
                            fa_file.write(right_seq + "\n")
                        right_blast_output = subprocess.check_output("blastn -query " + temp_dir + "/" + pbsv_id + ".right.fa " + " -db " +  Qry_chr + ".fasta " + blast_parameter, shell=True).decode().split("\n")[:-1]
                        right_blast = get_blast_hits(right_blast_output,Anchor_direction,Qry_chr)
                        ''' get blast pair '''
                        size_diff = 100000000000000
                        best_pair = []
                        if len(right_blast) > 0:
                            for L in left_blast:
                                if L[2] - L[6] < 100 and L[8] in exp_sites:
                                    for R in right_blast:
                                        if R[5] < 100 and R[7] in exp_sites:
                                            if Anchor_direction == "Forward":
                                                obs_gap_size = R[7] - L[8] - (R[5] - 1) - (L[2] - L[6])
                                                real_ref_sta = Start - (L[2] - L[6])
                                                real_ref_end = End + (R[5] - 1)
                                                real_qry_sta = L[8]
                                                real_qry_end = R[7] - 1
                                            else:
                                                obs_gap_size = L[8] - R[7] + (R[5] - 1) + (L[2] - L[6])
                                                real_ref_sta = Start - (L[2] - L[6])
                                                real_ref_end = End + (R[5] - 1)
                                                real_qry_sta = R[7] + 1
                                                real_qry_end = L[8]
                                            if (pbsv_type == "INS" and obs_gap_size >= 10 and abs(obs_gap_size - exp_gap_size)*100.0 / min(obs_gap_size, exp_gap_size) <=20) or \
                                               (pbsv_type == "DEL" and (abs(obs_gap_size) <= 3 or abs(obs_gap_size - exp_gap_size)*100.0 / pbsv_size <=20)):
                                                if abs(obs_gap_size - exp_gap_size) < size_diff:
                                                    size_diff = abs(obs_gap_size - exp_gap_size)
                                                    best_pair = [pbsv_id,Ref_chr,real_ref_sta,real_ref_end,Qry_chr,real_qry_sta,real_qry_end]
                                                    align_len = (L[6] - L[5]) + (R[6] - R[5])
                                                elif abs(obs_gap_size - exp_gap_size) == size_diff:
                                                    # same size diff, select the one with larger alignments length
                                                    if (L[6] - L[5]) + (R[6] - R[5]) > align_len:
                                                        size_diff = abs(obs_gap_size - exp_gap_size)
                                                        best_pair = [pbsv_id,Ref_chr,real_ref_sta,real_ref_end,Qry_chr,real_qry_sta,real_qry_end]
                                                        align_len = (L[6] - L[5]) + (R[6] - R[5])
                        if best_pair != []:
                            # found best pair
                            Ref_redefined.append(best_pair)
    ''' Remove redundancy '''
    if len(Ref_redefined) < 2:
        # no or only one
        return Ref_redefined
    else:
        if pbsv_type == "DEL":
            sort_table = sorted(Ref_redefined, key = itemgetter(2,3))
        else:
            sort_table = sorted(Ref_redefined, key = itemgetter(5,6))
        pbsv_id = sort_table[0][0]
        real_ref_sta = sort_table[0][2]
        real_ref_end = sort_table[0][3]
        real_qry_sta = sort_table[0][5]
        real_qry_end = sort_table[0][6]
        merged_table = []
        for fields in sort_table[1:]:
            previous_ref_size = real_ref_end - real_ref_sta
            previous_qry_size = real_qry_end - real_qry_sta
            this_id = fields[0]
            this_ref_sta = fields[2]
            this_ref_end = fields[3]
            this_qry_sta = fields[5]
            this_qry_end = fields[6]
            this_ref_size = this_ref_end - this_ref_sta
            this_qry_size = this_qry_end - this_qry_sta
            if pbsv_type == "DEL":
                # del, check ref
                if len(set(range(this_ref_sta,this_ref_end+1)).intersection(range(real_ref_sta,real_ref_end+1))) * 100.0 / min(previous_ref_size,this_ref_size) > 50:
                    # overlap > 50%
                    # take large one
                    pbsv_id = pbsv_id + ',' + this_id
                    if this_ref_size > previous_ref_size:
                        real_ref_sta = this_ref_sta
                        real_ref_end = this_ref_end
                        real_ref_sta = this_ref_sta
                        real_ref_end = this_ref_end
                else:
                    # no enough overlap
                    merged_table.append([pbsv_id,Ref_chr,real_ref_sta,real_ref_end,Qry_chr,real_qry_sta,real_qry_end])
                    pbsv_id = fields[0]
                    real_ref_sta = fields[2]
                    real_ref_end = fields[3]
                    real_qry_sta = fields[5]
                    real_qry_end = fields[6]

            else:
                # ins, check qry
                if len(set(range(this_qry_sta,this_qry_end+1)).intersection(range(real_qry_sta,real_qry_end+1))) * 100.0 / min(previous_qry_size,this_qry_size) > 50:
                    # overlap > 50%
                    # take large one
                    pbsv_id = pbsv_id + ',' + this_id
                    if this_qry_size > previous_qry_size:
                        real_ref_sta = this_ref_sta
                        real_ref_end = this_ref_end
                        real_qry_sta = this_qry_sta
                        real_qry_end = this_qry_end
                else:
                    # no enough overlap
                    merged_table.append([pbsv_id,Ref_chr,real_ref_sta,real_ref_end,Qry_chr,real_qry_sta,real_qry_end])
                    pbsv_id = fields[0]
                    real_ref_sta = fields[2]
                    real_ref_end = fields[3]
                    real_qry_sta = fields[5]
                    real_qry_end = fields[6]
        merged_table.append([pbsv_id,Ref_chr,real_ref_sta,real_ref_end,Qry_chr,real_qry_sta,real_qry_end])
        return merged_table



def get_blast_hits(blast_output,Anchor_direction,hit_chr):
    blast_table = []
    if Anchor_direction == "Forward":
        hit_str = "plus"
    else:
        hit_str = "minus"
    for line in blast_output:
        line = line.strip()
        if not line.startswith("#"):
            cells = convert_int(line.split("\t"))
            if cells[9] == hit_str and cells[1] == hit_chr:
                blast_table.append(cells)
    return blast_table

def add_merged_pbsv(pbsv_IDs,Merged_dict,cells):
    if pbsv_IDs != "None":
        for pbsv_id in pbsv_IDs.split(","):
            add_to_dict(pbsv_id,cells,Merged_dict)

'''
def between_anchor_other_pbsv_check(cells,ref_pbsv_by_chr,query_pbsv_by_chr,ref_base_type,Merged_Ref_pbsv,Merged_Qry_pbsv):
    SV_type = cells[6]
    SV_size = cells[4]
    Ref_chr = cells[0]
    Ref_start = cells[1]
    Ref_end = cells[2]
    Ref_Size = cells[7]
    Qry_chr = cells[9].split(":")[0]
    Qry_start = int(cells[9].split(":")[1].split("-")[0])
    Qry_end = int(cells[9].split(":")[1].split("-")[1])
    Qry_Size = cells[8]
    if ref_base_type == "DEL":
        Ref_pbsv_IDs = Repeat_overlap_DEL(Ref_chr,Ref_start,Ref_end,Ref_Size,Qry_Size,ref_pbsv_by_chr,Merged_Ref_pbsv)
        Qry_pbsv_IDs = Repeat_overlap_INS(Qry_chr,Qry_start,Qry_end,Ref_Size,Qry_Size,query_pbsv_by_chr,Merged_Qry_pbsv)
    else:
        Ref_pbsv_IDs = Repeat_overlap_INS(Ref_chr,Ref_start,Ref_end,Qry_Size,Ref_Size,ref_pbsv_by_chr,Merged_Ref_pbsv)
        Qry_pbsv_IDs = Repeat_overlap_DEL(Qry_chr,Qry_start,Qry_end,Qry_Size,Ref_Size,query_pbsv_by_chr,Merged_Qry_pbsv)
    return Ref_pbsv_IDs,Qry_pbsv_IDs

def Repeat_overlap_DEL(chrID,DEL_start,DEL_end,DEL_size,INS_Size,pbsv_by_chr,already_merged):
    # pbsv information:
    # SV_id   Type    GT      Chromosome      Start   End     Ref_seq Query_seq       Ref_size        Query_size
    hit_pbsv_IDs = []
    for cells in pbsv_by_chr["DEL"][chrID]:
        if cells[0] not in already_merged:
            if not (max(DEL_start,DEL_end) < min(cells[4],cells[5]) or min(DEL_start,DEL_end) > max(cells[4],cells[5])):
                overlap_size = len(set(range(DEL_start,DEL_end)).intersection(range(cells[4],cells[5])))
                pbsv_size = len(cells[6])
                if overlap_size * 100.0 / pbsv_size > 50:
                    # overlap >50% of pbsv_size only
                    hit_pbsv_IDs.append(cells[0])

    if len(hit_pbsv_IDs) == 0:
        hit_pbsv_IDs = "None"
    else:
        hit_pbsv_IDs = ",".join(hit_pbsv_IDs)
    return hit_pbsv_IDs

def Repeat_overlap_INS(chrID,INS_start,INS_end,DEL_size,INS_Size,pbsv_by_chr,already_merged):
    # pbsv information:
    # SV_id   Type    GT      Chromosome      Start   End     Ref_seq Query_seq       Ref_size        Query_size
    hit_pbsv_IDs = []
    exp_ins_sites = set(range(INS_start-3,INS_end+4))

    for cells in pbsv_by_chr["INS"][chrID]:
        pbsv_pos = cells[4]
        if cells[0] not in already_merged and pbsv_pos in exp_ins_sites:
            hit_pbsv_IDs.append(cells[0])

    if len(hit_pbsv_IDs) == 0:
        hit_pbsv_IDs = "None"
    else:
        hit_pbsv_IDs = ",".join(hit_pbsv_IDs)
    return hit_pbsv_IDs
'''

''' Get anchors '''
anchor_dict = {}
anchor_table = []

with open(anchor_file) as infile:
    for line in infile:
        line = line.strip()
        cells = convert_int(line.split("\t"))
        ref_chr = cells[6]
        qry_chr = cells[7]
        add_to_dict(ref_chr + "-" + qry_chr,cells,anchor_dict)
        anchor_table.append(cells)


''' Prepare pbsv input '''
ref_pbsv_by_chr = {}
ref_pbsv_by_ID = {}
get_pbsv_list(ref_pbsv,ref_pbsv_by_chr,ref_pbsv_by_ID)

query_pbsv_by_chr = {}
query_pbsv_by_ID = {}
get_pbsv_list(query_pbsv,query_pbsv_by_chr,query_pbsv_by_ID)


''' Main Program'''
Merged_Ref_pbsv = {}
Merged_Qry_pbsv = {}


"""
Step 1. Find overlapping pbsv and Assemblytics SVs
    1.1. keep within_alignment Assemblytics SVs
    1.2. redefine between_alignment Assemblytics SVs
"""
with open(SV_file) as infile, open(prefix + ".added_pbsv.tab", "w") as outfile, open(prefix + ".formated.tsv", "w") as formated_out:
    for line in infile:
        line = line.strip()
        cells = convert_int(line.split("\t"))
        if line.startswith("#"):
            outfile.write(line + "\t" + ref_pbsv + "\t" + query_pbsv + "\n")
            formated_list = ["SV_type","Ref_chr","Ref_start","Ref_end","Ref_Size","Qry_chr","Qry_start","Qry_end","Qry_Size","Assemblytics",ref_pbsv,query_pbsv,"SV_ID"]
            formated_out.write("\t".join(formated_list) + "\n")
        else:
            SV_type = cells[6]
            SV_size = cells[4]
            Ref_chr = cells[0]
            Ref_start = cells[1]
            Ref_end = cells[2]
            Ref_Size = cells[7]
            Qry_chr = cells[9].split(":")[0]
            Qry_start = int(cells[9].split(":")[1].split("-")[0])
            Qry_end = int(cells[9].split(":")[1].split("-")[1])
            Qry_Size = cells[8]
            Final_ID = cells[3]
            Assemblytics_ID = cells[3].replace("IND_","Assemblytics_")
            Anchor_direction = cells[11]
            Anchor_detail = cells[12:-1]
            if cells[10] == "within_alignment":
                # Find pbsvs overlaping with within_alignment indels
                Ref_pbsv_IDs,Qry_pbsv_IDs = within_pbsv_check(cells,ref_pbsv_by_chr,query_pbsv_by_chr)
                outfile.write(line + "\t" + Ref_pbsv_IDs + "\t" + Qry_pbsv_IDs + "\n")
                formated_list = [SV_type,Ref_chr,Ref_start,Ref_end,Ref_Size,Qry_chr,Qry_start,Qry_end,Qry_Size,Assemblytics_ID,Ref_pbsv_IDs,Qry_pbsv_IDs,Final_ID]
                formated_out.write("\t".join(map(str,formated_list)) + "\n")
                add_merged_pbsv(Ref_pbsv_IDs,Merged_Ref_pbsv,cells)
                add_merged_pbsv(Qry_pbsv_IDs,Merged_Qry_pbsv,cells)
            else:
                # Find pbsvs overlaping with between_alignment indels
                Ref_pbsv_IDs,Qry_pbsv_IDs = between_pbsv_check(cells,ref_pbsv_by_chr,query_pbsv_by_chr)
                outfile.write(line + "\t" + Ref_pbsv_IDs + "\t" + Qry_pbsv_IDs + "\n")
                add_merged_pbsv(Ref_pbsv_IDs,Merged_Ref_pbsv,cells)
                add_merged_pbsv(Qry_pbsv_IDs,Merged_Qry_pbsv,cells)
                redefine_between_SV(cells,Ref_pbsv_IDs,Qry_pbsv_IDs,ref_pbsv_by_ID,query_pbsv_by_ID,formated_out,"Assemblytics_b_","IND_b_")
    """
        Step 2. For merged pbsv, try to find those overlap with between_alignment Assemblytics SVs, then define it in the other genome based on Assemblytics SVs region.

    for line in infile:
        line = line.strip()
        cells = convert_int(line.split("\t"))
        SV_type = cells[6]
        cells[3] = cells[3].replace("Assemblytics_b_","Anchor_b_")
        if cells[10] == "between_alignments" and SV_type != "Insertion" and SV_type != "Deletion":
            # ref based DEL and qry besed INS
            Ref_pbsv_IDs,Qry_pbsv_IDs = between_anchor_other_pbsv_check(cells,ref_pbsv_by_chr,query_pbsv_by_chr,"DEL",Merged_Ref_pbsv,Merged_Qry_pbsv)
            if Ref_pbsv_IDs != "None" or Qry_pbsv_IDs != "None":
                outfile.write(line + "\t" + Ref_pbsv_IDs + "\t" + Qry_pbsv_IDs + "\n")
                add_merged_pbsv(Ref_pbsv_IDs,Merged_Ref_pbsv,cells)
                add_merged_pbsv(Qry_pbsv_IDs,Merged_Qry_pbsv,cells)
                redefine_between_SV(cells,Ref_pbsv_IDs,Qry_pbsv_IDs,ref_pbsv_by_ID,query_pbsv_by_ID,formated_out,"Anchor_b_","IND_AD_")

            # ref based INS and qry besed DEL
            Ref_pbsv_IDs,Qry_pbsv_IDs = between_anchor_other_pbsv_check(cells,ref_pbsv_by_chr,query_pbsv_by_chr,"INS",Merged_Ref_pbsv,Merged_Qry_pbsv)
            if Ref_pbsv_IDs != "None" or Qry_pbsv_IDs != "None":
                outfile.write(line + "\t" + Ref_pbsv_IDs + "\t" + Qry_pbsv_IDs + "\n")
                add_merged_pbsv(Ref_pbsv_IDs,Merged_Ref_pbsv,cells)
                add_merged_pbsv(Qry_pbsv_IDs,Merged_Qry_pbsv,cells)
                redefine_between_SV(cells,Ref_pbsv_IDs,Qry_pbsv_IDs,ref_pbsv_by_ID,query_pbsv_by_ID,formated_out,"Anchor_b_","IND_AI_")
    """



















'''
print ("Good_Ref_pbsv:\t" + str(len(ref_pbsv_by_ID)))
print ("Merged_Ref_pbsv:\t" + str(len(Merged_Ref_pbsv)))
print ("Good_Qry_pbsv:\t" + str(len(query_pbsv_by_ID)))
print ("Merged_Qry_pbsv:\t" + str(len(Merged_Qry_pbsv)))
'''























#srf
