"""
By Lei Gao
usage: Get_uncombined_between_SV_back.py [-h] --Raw_SVs RAW_SVS --Combined_NR
                                         COMBINED_NR --Out_prefix OUT_PREFIX

optional arguments:
  -h, --help            show this help message and exit
  --Raw_SVs RAW_SVS     Raw SV list. Assemblytics format
  --Combined_NR COMBINED_NR
                        After combination of psbv. Formated for genotyping.
  --Out_prefix OUT_PREFIX
                        Prefix for output file.

Keep:
Within_only         SVs within alignments
Redefined           The SVs between alignments and redefined by pbsv
Raw_for_redefined   Raw SVs for those redefined
Non_redefined       Raw SVs for those not redefined
"""
import argparse
from pathlib import Path
from operator import itemgetter
import re
import os
import subprocess
parser = argparse.ArgumentParser()
parser.add_argument("--Raw_SVs", type=str, help="Raw SV list. Assemblytics format", required=True, default="")
parser.add_argument("--Combined_NR", type=str, help="After combination of psbv. Formated for genotyping.", required=True, default="")
parser.add_argument("--Out_prefix", type=str, help="Prefix for output file.", required=True, default="")

args = parser.parse_args()

Raw_SVs = args.Raw_SVs
Combined_NR = args.Combined_NR
Out_prefix = args.Out_prefix

letters = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'

''' Function '''
def convert_int(fields):
    for i in range(0,len(fields)):
        if fields[i].isdigit() or (fields[i].startswith("-") and fields[i][1:].isdigit()):
            fields[i] = int(fields[i])
    return fields

def adjust_coord(cells):
    Ref_start = cells[2]
    Ref_end = cells[3]
    cells[2] = min(Ref_start,Ref_end)
    cells[3] = max(Ref_start,Ref_end)
    Qry_start = cells[6]
    Qry_end = cells[7]
    cells[6] = min(Qry_start,Qry_end)
    cells[7] = max(Qry_start,Qry_end)
    return cells


def Cover_check(Assemblytics_ID,cells):
    # Raw SV information
    fields = raw_SVs[Assemblytics_ID]
    Raw_ref = fields[0:3]
    Raw_qry = convert_int(fields[9].replace("-",":").split(":"))
    Raw_size = fields[4]
    Raw_type = fields[6]
    ref_gap_size = fields[7]
    query_gap_size = fields[8]
    # Formated SVs
    Fmt_ref = cells[1:4]
    Fmt_qry = cells[5:8]
    # Ref
    Ref_Loc, Ref_cov = cov_count(Raw_ref,Fmt_ref)
    # Qry
    Qry_Loc, Qry_cov = cov_count(Raw_qry,Fmt_qry)

    return [Ref_Loc, Ref_cov, Qry_Loc, Qry_cov]

def cov_count(Raw_coord,Fmt_coord):
    Raw_sta = min(Raw_coord[1],Raw_coord[2])
    Raw_end = max(Raw_coord[1],Raw_coord[2])
    Fmt_sta = min(Fmt_coord[1],Fmt_coord[2])
    Fmt_end = max(Fmt_coord[1],Fmt_coord[2])
    # if no overlap, give distance
    # else, give overlapped size
    coords = sorted([Raw_sta,Raw_end,Fmt_sta,Fmt_end])
    X =  coords[2] - coords[1]
    if Raw_end < Fmt_sta or Fmt_end < Raw_sta:
        Loc = "Outside"
    else:
        Loc = "Overlap"
    return Loc, X

def overlap_cov(gap_size,Overlapped):
    if gap_size == 0:
        Over = 1
    else:
        Over = float(Overlapped) / abs(gap_size)
    return Over


''' Get raw between SV information '''
raw_SVs = {}

with open(Raw_SVs) as infile:
    for line in infile:
        cells = convert_int(line.strip().split("\t"))
        if cells[10] == "between_alignments":
            raw_SVs[cells[3]] = cells

''' Get raw between SV information '''
Formated_SVs = set()
Good_defined = set()
with open(Combined_NR) as infile, open(Out_prefix + ".Check_comb.tsv", "w") as outfile, \
     open(Out_prefix + ".Master_list.tsv", "w") as final_out, \
     open(Out_prefix + ".Within_only.tsv", "w") as Within_only, \
     open(Out_prefix + ".Redefined.tsv", "w") as Redefined, \
     open(Out_prefix + ".Raw_for_redefined.tsv", "w") as Raw_for_redefined, \
     open(Out_prefix + ".Non_redefined.tsv", "w") as Non_redefined:
    i = -1
    for line in infile:
        i += 1
        line = line.strip()
        if i == 0:
            outfile.write(line + "\tNote\n")
            final_out.write(line + "\n")
            Within_only.write(line + "\n")
            Redefined.write(line + "\n")
            Raw_for_redefined.write(line + "\n")
            Non_redefined.write(line + "\n")
        else:
            line = line.strip().replace("IND_","SV_")
            cells = adjust_coord(convert_int(line.strip().split("\t")))
            line = "\t".join(map(str,cells))
            if "_w_" in cells[9]:
                # SV within alignment
                Note = "within_alignment"
                outfile.write(line + "\t" + Note + "\n")
                final_out.write(line + "\n")
                Within_only.write(line + "\n")
            else:
                # between alignments
                Assemblytics_ID = cells[9]
                fields = raw_SVs[Assemblytics_ID]
                Raw_type = fields[6]
                Raw_size = fields[4]
                ref_gap_size = fields[7]
                query_gap_size = fields[8]
                if cells[10] == "None" and cells[11] == "None":
                    # No pbsv evidence
                    Note = "No_Pbsv"
                    outfile.write(line + "\t" + Note + "\n")
                    #Covered_out = ["NAN"] * 4
                elif Raw_type == "Substitution":
                    Note = "Substitution"
                    outfile.write(line + "\t" + Note + "\n")
                    #Covered_out = ["NAN"] * 4
                else:
                    # Note = Cover_check(Assemblytics_ID,cells)
                    # outfile.write(line + "\t" + Note + "\n")
                    # covered size
                    Covered_out = Cover_check(Assemblytics_ID,cells)
                    # [Ref_Loc, Ref_cov, Qry_Loc, Qry_cov]
                    if Covered_out[0] == Covered_out[2] == "Overlap":
                        Note = "Both_Over"
                        #Ref_cov = overlap_cov(ref_gap_size,Covered_out[1])
                        #Qry_cov = overlap_cov(query_gap_size,Covered_out[3])
                        if Raw_type == "Insertion":
                            # check query
                            Qry_cov = overlap_cov(Raw_size,Covered_out[3])
                            if Qry_cov > 0.5:
                                Note = "Good_INS\t" + str(Qry_cov)
                            else:
                                Note = "Bad_INS\t" + str(Qry_cov)
                        elif Raw_type == "Deletion":
                            Ref_cov = overlap_cov(Raw_size,Covered_out[1])
                            if Ref_cov > 0.5:
                                Note = "Good_DEL\t" + str(Ref_cov)
                            else:
                                Note = "Bad_DEL\t" + str(Ref_cov)
                        else:
                            Ref_cov = overlap_cov(Raw_size,Covered_out[1])
                            Qry_cov = overlap_cov(Raw_size,Covered_out[3])
                            if max(Ref_cov,Qry_cov) > 0.5:
                                Note = "Good_REP\t" + str(Ref_cov) + "\t" + str(Qry_cov)
                            else:
                                Note = "Bad_REP\t" + str(Ref_cov) + "\t" + str(Qry_cov)
                    elif Covered_out[0] == Covered_out[2] == "Outside":
                        # bad SVs
                        Note = "Both_Outside"
                    else:
                        Note = "One_Overlap"
                        if Covered_out[2] == "Overlap":
                            Cov = overlap_cov(Raw_size,Covered_out[3])
                            distance = overlap_cov(Raw_size,Covered_out[1])
                        else:
                            Cov = overlap_cov(Raw_size,Covered_out[1])
                            distance = overlap_cov(Raw_size,Covered_out[3])

                        if Covered_out[2] == "Overlap":
                            if Cov > 0.75 and distance < 0.25:
                                Note = "Good_red_INS_\t" + str(Cov) + "\t_Dis\t" + str(distance)
                            else:
                                Note = "Bad_red_INS\t" + str(Cov) + "\t_Dis\t" + str(distance)
                        elif Covered_out[0] == "Overlap":
                            if Cov > 0.75 and distance < 0.25:
                                Note = "Good_red_DEL\t" + str(Cov) + "\t_Dis\t" + str(distance)
                            else:
                                Note = "Bad_red_DEL\t" + str(Cov) + "\t_Dis\t" + str(distance)

                    outlist = cells + Covered_out + [Note] + fields
                    #outlist = cells + [Note] + Covered_out + [Ref_cov,Qry_cov]
                    outfile.write("\t".join(map(str,outlist)) + "\n")

                if Note.startswith("Good") or Note == "No_Pbsv":
                    # Formated
                    Formated_SVs.add(Assemblytics_ID)
                    final_out.write(line + "\n")
                    if Note.startswith("Good"):
                        Redefined.write(line + "\n")
                        Good_defined.add(Assemblytics_ID)
                    else:
                        # No_Pbsv
                        Non_redefined.write(line + "\n")


with open(Out_prefix + ".Master_list.tsv", "a+") as final_out, open(Out_prefix + ".Non_redefined.tsv", "a+") as Non_redefined, open(Out_prefix + ".Raw_for_redefined.tsv", "a+") as Raw_for_redefined:
    for Assemblytics_ID in sorted(raw_SVs.keys()):
        cells = raw_SVs[Assemblytics_ID]
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
        Final_ID = cells[3].replace("Assemblytics_","SV_")
        Ref_pbsv_IDs = cells[-2]
        Qry_pbsv_IDs = cells[-1]
        formated_list = [SV_type,Ref_chr,Ref_start,Ref_end,Ref_Size,Qry_chr,Qry_start,Qry_end,Qry_Size,Assemblytics_ID,Ref_pbsv_IDs,Qry_pbsv_IDs,Final_ID]
        if Assemblytics_ID not in Formated_SVs:
            final_out.write("\t".join(map(str,formated_list)) + "\n")
            Non_redefined.write("\t".join(map(str,formated_list)) + "\n")
        if Assemblytics_ID in Good_defined:
            Raw_for_redefined.write("\t".join(map(str,formated_list)) + "\n")
