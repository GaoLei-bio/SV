'''
By Lei Gao
usage: python Formating_rawSVs.py prefix.NR.bed > prefix.Genome_comp_SV.tsv
'''

import sys

Input = sys.argv[1]

i = -1
with open(Input) as infile:
    for line in infile:
        i += 1
        line = line.strip()
        cells = line.split("\t")
        if i == 0:
            formated_list = ["SV_type","Ref_chr","Ref_start","Ref_end","Ref_Size","Qry_chr","Qry_start","Qry_end","Qry_Size","SV_ID", "Method"]
            print "\t".join(formated_list)
        else:
            SV_type = cells[6].split("/")[0]
            SV_size = cells[4]
            Ref_chr = cells[0]
            Ref_start = cells[1]
            Ref_end = cells[2]
            Ref_Size = cells[7]
            Qry_chr = cells[9].split(":")[0]
            Qry_start = int(cells[9].split(":")[1].split("-")[0])
            Qry_end = int(cells[9].split(":")[1].split("-")[1])
            Qry_Size = cells[8]
            Assemblytics_ID = cells[3].replace("Assemblytics_","IND_")
            if "_b_" in Assemblytics_ID:
                formated_list = [SV_type,Ref_chr,Ref_start,Ref_end,Ref_Size,Qry_chr,Qry_start,Qry_end,Qry_Size,Assemblytics_ID,"between_alignments"]
            else:
                formated_list = [SV_type,Ref_chr,Ref_start,Ref_end,Ref_Size,Qry_chr,Qry_start,Qry_end,Qry_Size,Assemblytics_ID,"within_alignment"]
            print "\t".join(map(str,formated_list))


#
