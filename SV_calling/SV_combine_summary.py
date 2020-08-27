"""
By Lei Gao
usage: SV_combine_summary.py [-h] --bed_input BED_INPUT --tsv_out TSV_OUT
                             --Good_within GOOD_WITHIN --prefix PREFIX

optional arguments:
  -h, --help            show this help message and exit
  --bed_input BED_INPUT
                        Filtered Assemblytics output SV file, bed format
  --tsv_out TSV_OUT     Combined pbsv indels. For complex SVs, only keep those
                        can be redefined as indels
  --Good_within GOOD_WITHIN
                        My within SV list
  --prefix PREFIX       Prefix for output files.

"""
import argparse
import subprocess
import re
from operator import itemgetter

parser = argparse.ArgumentParser()
parser.add_argument("--bed_input", type=str, help="Filtered Assemblytics output SV file, bed format", required=True, default="")
parser.add_argument("--tsv_out", type=str, help="Combined pbsv indels. For complex SVs, only keep those can be redefined as indels", required=True, default="")
parser.add_argument("--Good_within", type=str, help="My within SV list", required=True, default="")
parser.add_argument("--prefix", type=str, help="Prefix for output files.", required=True, default="")

args = parser.parse_args()

bed_input = args.bed_input
tsv_out = args.tsv_out
prefix = args.prefix
Good_within = set(open(args.Good_within).read().split())

input_sv_number = len(open(bed_input).read().split("\n")) - 1


''' Function '''
def convert_int(fields):
    for i in range(0,len(fields)):
        if fields[i].isdigit() or (fields[i].startswith("-") and fields[i][1:].isdigit()):
            fields[i] = int(fields[i])
    return fields

def check_overlap(temp_pass_table):
    overlaped_ids = []
    sort_table = sorted(temp_pass_table, key = itemgetter(0,1,2))
    sv_chr = sort_table[0][0]
    sv_sta = sort_table[0][1]
    sv_end = sort_table[0][2]
    sv_id = sort_table[0][3]
    sv_size = max(sv_sta,sv_end) - min(sv_sta,sv_end)
    for sv in sort_table[1:]:
        this_chr = sv[0]
        this_sta = sv[1]
        this_end = sv[2]
        this_id = sv[3]
        this_size = max(this_sta,this_end) - min(this_sta,this_end)
        if this_chr == sv_chr and sv_sta <= this_sta <= sv_end:
            # overlaping
            overlap_size = len(set(range(min(sv_sta,sv_end),max(sv_sta,sv_end))).intersection(range(min(this_sta,this_end),max(this_sta,this_end))))
            if overlap_size * 100.0 / sv_size > 50:
                overlaped_ids.append(sv_id)
            if overlap_size * 100.0 / this_size > 50:
                overlaped_ids.append(this_id)
        sv_chr = this_chr
        sv_sta = this_sta
        sv_end = this_end
        sv_id = this_id
        sv_size = this_size
    return overlaped_ids

def add_to_dict(Key,Value,Dict):
    if Key not in Dict:
        Dict[Key] = []
        Dict[Key].append(Value)
    else:
        Dict[Key].append(Value)


''' Main program '''
temp_pass_ID_list = []
temp_pass_by_ID = {}

temp_pass_by_Ref = []
temp_pass_by_Qry = []

both_pbsv = 0
ref_pbsv = 0
qry_pbsv = 0
no_pbsv = 0
Assemblytics_ids = []
ref_pbsv_ids = []
qry_pbsv_ids = []

with open(tsv_out) as infile, open(prefix + ".formated.NR.tsv", "w") as nr_out, open(prefix + ".formated.redundancy.tsv", "w") as red_out, open(prefix + ".conbine_pbsv.summary", "w") as summary_file:
    for line in infile:
        if line.startswith("SV_type"):
            nr_out.write(line)
            red_out.write(line)
        else:
            fields = convert_int(line.strip().split("\t"))
            sv_id = fields[-1]
            temp_pass_ID_list.append(sv_id)
            temp_pass_by_ID[sv_id] = line
            ref_size = fields[4]
            qry_size = fields[8]
            if ref_size > qry_size:
                temp_pass_by_Ref.append(fields[1:4] + [sv_id])
            else:
                temp_pass_by_Qry.append(fields[5:8] + [sv_id])

    ''' check overlaping '''
    overlaped_by_Ref = check_overlap(temp_pass_by_Ref)
    overlaped_by_Qry = check_overlap(temp_pass_by_Qry)
    overlaped_IDs = set(overlaped_by_Ref + overlaped_by_Qry)

    ''' get NR list '''
    keep_by_ref_coord = {}
    keep_by_qry_coord = {}
    out_IDs = set()

    for sv_id in temp_pass_ID_list:
        line = temp_pass_by_ID[sv_id]
        if sv_id in overlaped_IDs and sv_id not in Good_within:
            red_out.write(line)
            out_IDs.add(sv_id)
        else:
            cells = line.strip().split("\t")
            ref_key = "_".join(cells[1:4])
            qry_key = "_".join(cells[5:8])
            add_to_dict(ref_key,cells,keep_by_ref_coord)

    for ref_key in keep_by_ref_coord.keys():
        if len(keep_by_ref_coord[ref_key]) == 1:
            cells = keep_by_ref_coord[ref_key][0]
            qry_key = "_".join(cells[5:8])
            add_to_dict(qry_key,cells,keep_by_qry_coord)
        else:
            found_IND_w = "No"
            for cells in keep_by_ref_coord[ref_key]:
                if cells[-1] in Good_within:
                    found_IND_w = "Yes"
                    qry_key = "_".join(cells[5:8])
                    add_to_dict(qry_key,cells,keep_by_qry_coord)
                    break
            if found_IND_w == "No":
                best_support = 0
                for fields in keep_by_ref_coord[ref_key]:
                    if fields[10] != "None" and fields[11] != "None":
                        pbsv = 2
                    elif fields[10] != "None" or fields[11] != "None":
                        pbsv = 1
                    else:
                        pbsv = 0
                    if pbsv >= best_support:
                        best_support = pbsv
                        best_sv = fields
                cells = best_sv
                qry_key = "_".join(cells[5:8])
                add_to_dict(qry_key,cells,keep_by_qry_coord)
    NR_count = 0
    for qry_key in keep_by_qry_coord.keys():
        if len(keep_by_qry_coord[qry_key]) == 1:
            fields = keep_by_qry_coord[qry_key][0]

            NR_count += 1
            nr_out.write("\t".join(fields) + "\n")
            out_IDs.add(fields[-1])
            if fields[10] != "None" and fields[11] != "None":
                both_pbsv += 1
            elif fields[10] != "None":
                ref_pbsv += 1
            elif fields[11] != "None":
                qry_pbsv += 1
            else:
                no_pbsv += 1

            Assemblytics_ids.append(fields[9])
            if fields[10] != "None":
                ref_pbsv_ids += fields[10].split(",")
            if fields[11] != "None":
                qry_pbsv_ids += fields[11].split(",")
        else:
            found_IND_w = "No"
            for fields in keep_by_qry_coord[qry_key]:
                if fields[-1] in Good_within and found_IND_w == "No":
                    found_IND_w = "Yes"
                    NR_count += 1
                    nr_out.write("\t".join(fields) + "\n")
                    out_IDs.add(fields[-1])
                    if fields[10] != "None" and fields[11] != "None":
                        both_pbsv += 1
                    elif fields[10] != "None":
                        ref_pbsv += 1
                    elif fields[11] != "None":
                        qry_pbsv += 1
                    else:
                        no_pbsv += 1

                    Assemblytics_ids.append(fields[9])
                    if fields[10] != "None":
                        ref_pbsv_ids += fields[10].split(",")
                    if fields[11] != "None":
                        qry_pbsv_ids += fields[11].split(",")

                    break

            if found_IND_w == "No":
                best_support = 0
                for fields in keep_by_qry_coord[qry_key]:
                    if fields[10] != "None" and fields[11] != "None":
                        pbsv = 2
                    elif fields[10] != "None" or fields[11] != "None":
                        pbsv = 1
                    else:
                        pbsv = 0
                    if pbsv >= best_support:
                        best_support = pbsv
                        best_sv = fields
                fields = best_sv
                NR_count += 1
                nr_out.write("\t".join(fields) + "\n")
                out_IDs.add(fields[-1])
                if fields[10] != "None" and fields[11] != "None":
                    both_pbsv += 1
                elif fields[10] != "None":
                    ref_pbsv += 1
                elif fields[11] != "None":
                    qry_pbsv += 1
                else:
                    no_pbsv += 1

                Assemblytics_ids.append(fields[9])
                if fields[10] != "None":
                    ref_pbsv_ids += fields[10].split(",")
                if fields[11] != "None":
                    qry_pbsv_ids += fields[11].split(",")


    for sv_id in sorted(list(set(temp_pass_ID_list).difference(out_IDs))):
        line = temp_pass_by_ID[sv_id]
        red_out.write(line)

    summary_file.write("Input SV number:\t" + str(input_sv_number) + "\n")
    summary_file.write("Raw and redefined indel number:\t" + str(len(temp_pass_ID_list)) + "\n")
    summary_file.write("Amonng them, " + str(len(overlaped_IDs)) + " with >50% being covered by other SVs\n")
    #summary_file.write("Kept non-redundant (NR) indel number:\t" + str(len(temp_pass_ID_list) - len(overlaped_IDs)) + "\n")
    summary_file.write("Kept non-redundant (NR) indel number:\t" + str(NR_count) + "\n")
    summary_file.write("#############################:\n")
    summary_file.write("In these final NR indel list:\n")
    summary_file.write("\t" + str(both_pbsv) + " have both reference- and query-based pbsv supports\n")
    summary_file.write("\t" + str(ref_pbsv) + " have only reference-based pbsv support\n")
    summary_file.write("\t" + str(qry_pbsv) + " have only query-based pbsv support\n")
    summary_file.write("\t" + str(no_pbsv) + " have no pbsv supports\n")
    summary_file.write("\t A total of " + str(len(set(Assemblytics_ids))) + " Assemblytics SVs have been combined \n")
    summary_file.write("\t A total of " + str(len(set(ref_pbsv_ids))) + " reference-based pbsv indels have been combined \n")
    summary_file.write("\t A total of " + str(len(set(qry_pbsv_ids))) + " query-based pbsv indels have been combined \n")

























#srf
