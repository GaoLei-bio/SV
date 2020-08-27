"""
By Lei Gao
usage: Call_Indel_within_alignment.py [-h] --sam_file SAM_FILE --ref_genome
                                      REF_GENOME --query_genome QUERY_GENOME
                                      --uniqe_anchor UNIQE_ANCHOR --min_size
                                      MIN_SIZE --max_size MAX_SIZE --prefix
                                      PREFIX

optional arguments:
  -h, --help            show this help message and exit
  --sam_file SAM_FILE   Input sam file. Produced by minmap2.
  --ref_genome REF_GENOME
                        Reference genome file, fasta format. Indexed by
                        samtools
  --query_genome QUERY_GENOME
                        Query species genome file, fasta format. Indexed by
                        samtools
  --uniqe_anchor UNIQE_ANCHOR
                        Unique anchors selected by RaGOO
                        (Assemblytics_uniq_anchor.py). We can only call Indel
                        based on unique anchors.
  --min_size MIN_SIZE   Minimal indel size
  --max_size MAX_SIZE   Maximal indel size
  --prefix PREFIX       Prefix for output files.



"""
import argparse
import subprocess
import re

parser = argparse.ArgumentParser()
parser.add_argument("--sam_file", type=str, help="Input sam file. Produced by minmap2.", required=True, default="")
parser.add_argument("--ref_genome", type=str, help="Reference genome file, fasta format. Indexed by samtools", required=True, default="")
parser.add_argument("--query_genome", type=str, help="Query species genome file, fasta format. Indexed by samtools", required=True, default="")
parser.add_argument("--uniqe_anchor", type=str, help="Unique anchors selected by RaGOO (Assemblytics_uniq_anchor.py). We can only call Indel based on unique anchors.", required=True, default="")
parser.add_argument("--min_size", type=str, help="Minimal indel size", required=True, default="")
parser.add_argument("--max_size", type=str, help="Maximal indel size", required=True, default="")
parser.add_argument("--prefix", type=str, help="Prefix for output files.", required=True, default="")

args = parser.parse_args()

sam_file = args.sam_file
ref_genome = args.ref_genome
query_genome = args.query_genome
UNIQE_ANCHOR = args.uniqe_anchor
min_size = int(args.min_size)
max_size = int(args.max_size)
prefix = args.prefix

''' Functions'''

def get_chr_size(chr_size,genome):
    with open(genome + ".fai") as infile:
        for line in infile:
            cells = line.strip().split("\t")
            chr_size[cells[0]] = int(cells[1])

def get_forward_anchor(cells):
    ref_chr = cells[2]
    ref_chr_len = ref_chr_size_dict[ref_chr]
    ref_sta = int(cells[3])
    ref_end = ref_sta - 1       # initial value

    qry_chr = cells[0]
    qry_chr_len = qry_chr_size_dict[qry_chr]

    CIGAR_STRs = re.findall(r"\d+\w",cells[5])
    # 1st
    if "H" in CIGAR_STRs[0] or "S" in CIGAR_STRs[0]:
        qry_sta = int(CIGAR_STRs[0][:-1]) + 1
        qry_end = qry_sta - 1  # initial value
    else:
        # "M"
        ref_end += int(CIGAR_STRs[0][:-1])
        qry_sta = 1
        qry_end = int(CIGAR_STRs[0][:-1])
    # 2nd to last 2nd
    for cigar in CIGAR_STRs[1:-1]:
        if "M" in cigar:
            ref_end += int(cigar[:-1])
            qry_end += int(cigar[:-1])
        elif "H" in cigar or "S" in cigar or "I" in cigar:
            qry_end += int(cigar[:-1])
        elif "D" in cigar:
            ref_end += int(cigar[:-1])
    # last
    if "M" in CIGAR_STRs[-1]:
        ref_end += int(CIGAR_STRs[-1][:-1])
        qry_end += int(CIGAR_STRs[-1][:-1])

    anchor_infor = [ref_sta,ref_end,qry_sta,qry_end,ref_chr_len,qry_chr_len,ref_chr,qry_chr]
    return anchor_infor

def get_reverse_anchor(cells):
    ref_chr = cells[2]
    ref_chr_len = ref_chr_size_dict[ref_chr]
    ref_sta = int(cells[3])
    ref_end = ref_sta - 1       # initial value

    qry_chr = cells[0]
    qry_chr_len = qry_chr_size_dict[qry_chr]

    CIGAR_STRs = re.findall(r"\d+\w",cells[5])

    # get ref end
    # 1st
    if "M" in CIGAR_STRs[0]:
        ref_end += int(CIGAR_STRs[0][:-1])
    # 2nd to last 2nd
    for cigar in CIGAR_STRs[1:-1]:
        if "M" in cigar or "D" in cigar:
            ref_end += int(cigar[:-1])
    # last
    if "M" in CIGAR_STRs[-1]:
        ref_end += int(CIGAR_STRs[-1][:-1])

    # get qry end
    # last
    if "H" in CIGAR_STRs[-1] or "S" in CIGAR_STRs[-1]:
        qry_sta = int(CIGAR_STRs[-1][:-1]) + 1
        qry_end = qry_sta - 1  # initial value
    else:
        # "M"
        qry_sta = 1
        qry_end = int(CIGAR_STRs[-1][:-1])
    # 2nd to last 2nd
    for cigar in CIGAR_STRs[1:-1]:
        if "M" in cigar or "H" in cigar or "S" in cigar or "I" in cigar:
            qry_end += int(cigar[:-1])
    # first
    if "M" in CIGAR_STRs[0]:
        qry_end += int(CIGAR_STRs[0][:-1])

    anchor_infor = [ref_sta,ref_end,qry_end,qry_sta,ref_chr_len,qry_chr_len,ref_chr,qry_chr]
    return anchor_infor


def SV_calling(cells,sv_order,anchor_infor,Anchor_direction,sv_bed_file):
    if Anchor_direction == "Forward":
        sv_order = Forward_SV_calling(cells,sv_order,anchor_infor,Anchor_direction,sv_bed_file)
    else:
        sv_order = Reverse_SV_calling(cells,sv_order,anchor_infor,Anchor_direction,sv_bed_file)
    return sv_order

def Forward_SV_calling(cells,sv_order,anchor_infor,Anchor_direction,sv_bed_file):
    CIGAR_STRs = re.findall(r"\d+\w",cells[5])
    ref_pos = anchor_infor[0]
    qry_pos = anchor_infor[2]
    ref_chr = anchor_infor[-2]
    qry_chr = anchor_infor[-1]
    # first one
    if "M" in CIGAR_STRs[0]:
        ref_pos += int(CIGAR_STRs[0][:-1])
        qry_pos += int(CIGAR_STRs[0][:-1])
    # from 2nd
    i = 1
    checked_cigars = set()
    for cigar in CIGAR_STRs[1:-1]:
        if "M" in cigar:
            ref_pos += int(cigar[:-1])
            qry_pos += int(cigar[:-1])
        elif "D" in cigar:
            if max_size >= int(cigar[:-1]) >= min_size and i not in checked_cigars:
                """Deletion"""
                ref_gap_size = int(cigar[:-1])
                qry_gap_size = 0
                sv_type = "Deletion"
                """ref coord"""
                ref_sv_sta = ref_pos
                ref_sv_end = ref_pos + int(cigar[:-1])
                """qry coord"""
                qry_sv_sta = qry_pos
                qry_sv_end = qry_pos
                if "I" in CIGAR_STRs[i+1]:
                    # if next is "I"
                    qry_sv_end += int(CIGAR_STRs[i+1][:-1])
                    qry_gap_size += int(CIGAR_STRs[i+1][:-1])
                    checked_cigars.add(i+1)
                if "I" in CIGAR_STRs[i-1] and int(CIGAR_STRs[i-1][:-1]) < min_size:
                    # previous is samll I
                    qry_sv_sta -= int(CIGAR_STRs[i-1][:-1])
                    qry_gap_size += int(CIGAR_STRs[i-1][:-1])
                query_coordinates = qry_chr + ":" + str(qry_sv_sta) + "-" + str(qry_sv_end) + ":+"
                """output"""
                indel_size = ref_gap_size + qry_gap_size
                if ref_gap_size < qry_gap_size:
                    sv_type = "Insertion"
                if not (100.0 * abs(ref_gap_size - qry_gap_size) / max(ref_gap_size, qry_gap_size) < 10 and max(ref_gap_size, qry_gap_size) >= 25):
                    # minmap2 may make mistake. It put a "I" immediately following "D" with same or nealy same size. I found if the size is larger than ~25 bp, acturaly it insert and delete same seqs. So we should discard this "SV".
                    # For small "I" and "D" pair, they are real Substitution
                    if min(ref_gap_size, qry_gap_size) > 0:
                        sv_type += "/Substitution"
                    sv_infor = [ref_chr, ref_sv_sta, ref_sv_end, "IND_w_" + str(sv_order), indel_size, "+", sv_type, ref_gap_size, qry_gap_size, query_coordinates, "within_alignment",Anchor_direction]
                    out_list = sv_infor + anchor_infor
                    sv_bed_file.write("\t".join(map(str,out_list)) + "\n")
                    sv_order += 1
            # get new ref_pos
            ref_pos += int(cigar[:-1])
        elif "I" in cigar:
            if max_size >= int(cigar[:-1]) >= min_size and i not in checked_cigars:
                """Insertion"""
                ref_gap_size = 0
                qry_gap_size = int(cigar[:-1])
                sv_type = "Insertion"
                """ref coord"""
                ref_sv_sta = ref_pos
                ref_sv_end = ref_pos
                if "D" in CIGAR_STRs[i+1]:
                    # if next is "D"
                    ref_sv_end += int(CIGAR_STRs[i+1][:-1])
                    ref_gap_size += int(CIGAR_STRs[i+1][:-1])
                    checked_cigars.add(i+1)
                if "D" in CIGAR_STRs[i-1] and int(CIGAR_STRs[i-1][:-1]) < min_size:
                    # if previous is a small D
                    ref_sv_sta -= int(CIGAR_STRs[i-1][:-1])
                    ref_gap_size += int(CIGAR_STRs[i-1][:-1])
                """qry coord"""
                qry_sv_sta = qry_pos
                qry_sv_end = qry_pos + int(cigar[:-1])
                query_coordinates = qry_chr + ":" + str(qry_sv_sta) + "-" + str(qry_sv_end) + ":+"
                """output"""
                indel_size = ref_gap_size + qry_gap_size
                if ref_gap_size > qry_gap_size:
                    sv_type = "Deletion"
                if not (100.0 * abs(ref_gap_size - qry_gap_size) / max(ref_gap_size, qry_gap_size) < 10 and max(ref_gap_size, qry_gap_size) >= 25):
                    if min(ref_gap_size, qry_gap_size) > 0:
                        sv_type += "/Substitution"
                    sv_infor = [ref_chr, ref_sv_sta, ref_sv_end, "IND_w_" + str(sv_order), indel_size, "+", sv_type, ref_gap_size, qry_gap_size, query_coordinates, "within_alignment",Anchor_direction]
                    out_list = sv_infor + anchor_infor
                    sv_bed_file.write("\t".join(map(str,out_list)) + "\n")
                    sv_order += 1
            # get new qry_pos
            qry_pos += int(cigar[:-1])
        # keep checked cigars
        checked_cigars.add(i)
        i += 1
    return sv_order



def Reverse_SV_calling(cells,sv_order,anchor_infor,Anchor_direction,sv_bed_file):
    CIGAR_STRs = re.findall(r"\d+\w",cells[5])
    ref_pos = anchor_infor[0]
    qry_pos = anchor_infor[2]
    ref_chr = anchor_infor[-2]
    qry_chr = anchor_infor[-1]
    # first one
    if "M" in CIGAR_STRs[0]:
        ref_pos += int(CIGAR_STRs[0][:-1])
        qry_pos -= int(CIGAR_STRs[0][:-1])
    # from 2nd
    i = 1
    checked_cigars = set()
    for cigar in CIGAR_STRs[1:-1]:
        if "M" in cigar:
            ref_pos += int(cigar[:-1])
            qry_pos -= int(cigar[:-1])
        elif "D" in cigar:
            if max_size >= int(cigar[:-1]) >= min_size and i not in checked_cigars:
                """Deletion"""
                ref_gap_size = int(cigar[:-1])
                qry_gap_size = 0
                sv_type = "Deletion"
                """ref coord"""
                ref_sv_sta = ref_pos
                ref_sv_end = ref_pos + int(cigar[:-1])
                """qry coord"""
                qry_sv_sta = qry_pos
                qry_sv_end = qry_pos
                if "I" in CIGAR_STRs[i+1]:
                    # if next is "I"
                    qry_sv_end -= int(CIGAR_STRs[i+1][:-1])
                    qry_gap_size += int(CIGAR_STRs[i+1][:-1])
                    checked_cigars.add(i+1)
                if "I" in CIGAR_STRs[i-1] and int(CIGAR_STRs[i-1][:-1]) < min_size:
                    # previous is samll I
                    qry_sv_sta += int(CIGAR_STRs[i-1][:-1])
                    qry_gap_size += int(CIGAR_STRs[i-1][:-1])
                query_coordinates = qry_chr + ":" + str(qry_sv_end) + "-" + str(qry_sv_sta) + ":-"
                """output"""
                indel_size = ref_gap_size + qry_gap_size
                if ref_gap_size < qry_gap_size:
                    sv_type = "Insertion"
                if not (100.0 * abs(ref_gap_size - qry_gap_size) / max(ref_gap_size, qry_gap_size) < 10 and max(ref_gap_size, qry_gap_size) >= 25):
                    if min(ref_gap_size, qry_gap_size) > 0:
                        sv_type += "/Substitution"
                    sv_infor = [ref_chr, ref_sv_sta, ref_sv_end, "IND_w_" + str(sv_order), indel_size, "-", sv_type, ref_gap_size, qry_gap_size, query_coordinates, "within_alignment",Anchor_direction]
                    out_list = sv_infor + anchor_infor
                    sv_bed_file.write("\t".join(map(str,out_list)) + "\n")
                    sv_order += 1
            # get new ref_pos
            ref_pos += int(cigar[:-1])
        elif "I" in cigar:
            if max_size >= int(cigar[:-1]) >= min_size and i not in checked_cigars:
                """Insertion"""
                ref_gap_size = 0
                qry_gap_size = int(cigar[:-1])
                sv_type = "Insertion"
                """ref coord"""
                ref_sv_sta = ref_pos
                ref_sv_end = ref_pos
                if "D" in CIGAR_STRs[i+1]:
                    # if next is "D"
                    ref_sv_end += int(CIGAR_STRs[i+1][:-1])
                    ref_gap_size += int(CIGAR_STRs[i+1][:-1])
                    checked_cigars.add(i+1)
                if "D" in CIGAR_STRs[i-1] and int(CIGAR_STRs[i-1][:-1]) < min_size:
                    # if previous is a small D
                    ref_sv_sta -= int(CIGAR_STRs[i-1][:-1])
                    ref_gap_size += int(CIGAR_STRs[i-1][:-1])
                """qry coord"""
                qry_sv_sta = qry_pos
                qry_sv_end = qry_pos - int(cigar[:-1])
                query_coordinates = qry_chr + ":" + str(qry_sv_end) + "-" + str(qry_sv_sta) + ":-"
                """output"""
                indel_size = ref_gap_size + qry_gap_size
                if ref_gap_size > qry_gap_size:
                    sv_type = "Deletion"
                if not (100.0 * abs(ref_gap_size - qry_gap_size) / max(ref_gap_size, qry_gap_size) < 10 and max(ref_gap_size, qry_gap_size) >= 25):
                    if min(ref_gap_size, qry_gap_size) > 0:
                        sv_type += "/Substitution"
                    sv_infor = [ref_chr, ref_sv_sta, ref_sv_end, "IND_w_" + str(sv_order), indel_size, "-", sv_type, ref_gap_size, qry_gap_size, query_coordinates, "within_alignment",Anchor_direction]
                    out_list = sv_infor + anchor_infor
                    sv_bed_file.write("\t".join(map(str,out_list)) + "\n")
                    sv_order += 1
            # get new qry_pos
            qry_pos -= int(cigar[:-1])
        # keep checked cigars
        checked_cigars.add(i)
        i += 1
    return sv_order


''' Get chromosome size '''
ref_chr_size_dict = {}
qry_chr_size_dict = {}

get_chr_size(ref_chr_size_dict,ref_genome)
get_chr_size(qry_chr_size_dict,query_genome)

''' Get uniqe_anchor '''
uniqe_anchor_set = set(open(UNIQE_ANCHOR).read().split("\n"))

''' Get alignments '''
forward_alignments = subprocess.check_output("samtools view -F 16  " + sam_file + "| cut -f 1-6" , shell=True).decode().split("\n")
reverse_alignments = subprocess.check_output("samtools view -f 16  " + sam_file + "| cut -f 1-6" , shell=True).decode().split("\n")

Uniq_For = 0
Uniq_Rev = 0

sv_order = 1

with open(prefix + ".variants_within_alignments.bed", "w") as sv_bed_file:
    bed_list = ["#reference", "ref_start", "ref_stop", "ID", "size", "strand", "type", "ref_gap_size", "query_gap_size", "query_coordinates", "method","Anchor_direction","ref_start","ref_end","query_start","query_end","ref_length","query_length","ref","query"]
    sv_bed_file.write("\t".join(bed_list) + "\n")
    # according to Assemblytics format
    # add anchor information
    for line in forward_alignments:
        line = line.strip()
        if line != "":
            cells = line.split("\t")
            anchor_infor = get_forward_anchor(cells)
            if "\t".join(map(str,anchor_infor)) in uniqe_anchor_set:
                Uniq_For += 1
                sv_order = SV_calling(cells,sv_order,anchor_infor,"Forward",sv_bed_file)

    SV_on_Forward = sv_order - 1

    for line in reverse_alignments:
        line = line.strip()
        if line != "":
            cells = line.split("\t")
            anchor_infor = get_reverse_anchor(cells)
            if "\t".join(map(str,anchor_infor)) in uniqe_anchor_set:
                Uniq_Rev += 1
                sv_order = SV_calling(cells,sv_order,anchor_infor,"Reverse",sv_bed_file)


# log
# log
print ("#" * 30)
print ("Input sam file:\t" + sam_file)
print ("Reference genome file:\t" + ref_genome)
print ("Query genome file:\t" + query_genome)
print ("Unique anchors for Indel calling:\t" + UNIQE_ANCHOR)
print ("#" * 20)
print ("Evidence Anchor check:")
print ("Forward alignment:\t" + str(len(forward_alignments)))
print ("Unique Forward alignment:\t" + str(Uniq_For))
print ("Repetitive Forward alignment:\t" + str(len(forward_alignments) - Uniq_For))
print ("Reverse alignment:\t" + str(len(reverse_alignments)))
print ("Uniqe Reverse alignment:\t" + str(Uniq_Rev))
print ("Repetitive Reverse alignment:\t" + str(len(reverse_alignments) - Uniq_Rev))
print ("Total unique anchors:\t" + str(Uniq_For + Uniq_Rev))
print ("#" * 20)
print ("SV calling:")
print ("Detected indels:\t" + str(sv_order - 1))
print ("On Forward alignment:\t" + str(SV_on_Forward))
print ("On Reverse alignment:\t" + str(sv_order - 1 - SV_on_Forward))


























#srf
