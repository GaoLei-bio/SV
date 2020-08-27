"""
By Lei Gao
usage: Synteny_block.py [-h] --Input_coords INPUT_COORDS --Prefix PREFIX

optional arguments:
  -h, --help            show this help message and exit
  --Input_coords INPUT_COORDS
                        Oriented coords from Assemblytics
  --Prefix PREFIX       My SV list. (1) Get filtered unique anchors from
                        Input_coords; (2) Filter Input_SV by unique anchors;
                        (3) Derive INV, BND and DUP from unique anchors.


"""
import argparse
from pathlib import Path
from operator import itemgetter
import re
import os

parser = argparse.ArgumentParser()
parser.add_argument("--Input_coords", type=str, help="Oriented coords from Assemblytics", required=True, default="")
parser.add_argument("--Prefix", type=str, help="My SV list. (1) Get filtered unique anchors from Input_coords; (2) Filter Input_SV by unique anchors; (3) Derive INV, BND and DUP from unique anchors.", required=True, default="")

args = parser.parse_args()

Input_coords = args.Input_coords
Prefix = args.Prefix

def put_best_first(NXGL_pbsv_Evi):
    if NXGL_pbsv_Evi[3] == "Self":
        NXGL_pbsv_Evi_hits = NXGL_pbsv_Evi[0]
    elif NXGL_pbsv_Evi[3] == "None":
        NXGL_pbsv_Evi_hits = NXGL_pbsv_Evi[3]
    else:
        best_hit = NXGL_pbsv_Evi[0]
        NXGL_pbsv_Evi_list = NXGL_pbsv_Evi[3].split("|")
        sorted_hits = [best_hit]
        for hit in NXGL_pbsv_Evi_list:
            if hit != best_hit:
                sorted_hits.append(hit)
        NXGL_pbsv_Evi_hits = "|".join(sorted_hits)
    return NXGL_pbsv_Evi_hits


def output_list_to_file(input_list,output_file):
    with open(output_file, "w") as OutFile:
        for id in input_list:
            OutFile.write(id + "\n")


def count_item(item,item_count):
    if item in item_count:
        item_count[item] += 1
    else:
        item_count[item] = 1

def add_to_dict(key,value,dict):
    if key in dict:
        dict[key].append(value)
    else:
        dict[key] = []
        dict[key].append(value)

def check_covered_by_neighbour(out_cells,i,step):
    conclusion = "None"
    dup_pair = []
    j = i + step
    next_cells = sam_chr_table[j]
    if out_cells[0] == next_cells[0]:
        # check in same chr
        this_sta = out_cells[1]
        this_end = out_cells[2]
        next_sta = next_cells[1]
        next_end = next_cells[2]
        if this_end <= next_sta or next_end <= this_sta:
            # not overlaping
            conclusion = "None"
        elif this_sta >= next_sta and this_end <= next_end:
            # This anchor is in next anchor
            conclusion = [i]
        elif this_sta <= next_sta and this_end >= next_end:
            # Next anchor is in this anchor
            conclusion = [j]
        else:
            # Partially overlaping
            this_len = this_end - this_sta
            next_len = next_end - next_sta
            anchor_coords = sorted([this_sta,this_end,next_sta,next_end])
            overlaping_size = anchor_coords[2] - anchor_coords[1]
            if float(overlaping_size) / this_len > 0.5 and float(overlaping_size) / next_len > 0.5:
                conclusion = [i,j]
                if float(overlaping_size) / this_len > 0.9 and float(overlaping_size) / next_len > 0.9:
                    dup_pair = [i,j]
            elif float(overlaping_size) / this_len > 0.5:
                conclusion = [i]
            elif float(overlaping_size) / next_len > 0.5:
                conclusion = [j]
    return conclusion,dup_pair


def check_cover_by_good_anchor(redundant_anchor,prev_good_anchor_Index,next_good_anchor_Index):
    conclusion = "None"
    check_index_list = sorted(list(set(range(prev_good_anchor_Index - 3, prev_good_anchor_Index + 3) +  range(next_good_anchor_Index - 3, next_good_anchor_Index + 3))))
    out_cells = redundant_anchor
    this_chr = out_cells[0]
    this_sta = out_cells[1]
    this_end = out_cells[2]
    this_sites = range(this_sta,this_end + 1)
    for index in check_index_list:
        next_cells = clean_sam_chr_table_With_orders_Syn_Block[index]
        next_chr = next_cells[0]
        next_sta = next_cells[1]
        next_end = next_cells[2]
        next_sites = range(next_sta,next_end + 1)
        if this_chr == next_chr:
            overlaping_size = len(set(this_sites).intersection(next_sites))
            if float(overlaping_size) / len(this_sites) > 0.5 or float(overlaping_size) / len(next_sites) > 0.5:
                conclusion = "Overlapped"
    return conclusion




def synteny_check(out_cells,ref_index,qry_index):
    # clean_sam_chr_table             Using Ref index
    # clean_sam_chr_table_By_qry      Using Qry index
    synteny = ''
    prev_ref_index = ref_index - 1
    next_ref_index = ref_index + 1

    if ref_index == 0 or out_cells[0] != clean_sam_chr_table[prev_ref_index][0]:
        # 1st anchor for each chr
        # Only need check next anchor
        next_cells = clean_sam_chr_table[next_ref_index]
        key = "_".join(map(str, next_cells))
        next_qry_index = Index_by_Qry[key]
        synteny = check_qry_index(qry_index,next_qry_index)

    elif ref_index == len(clean_sam_chr_table) - 1 or out_cells[0] != clean_sam_chr_table[next_ref_index][0]:
        # last anchor for each chr
        # Only need check previous anchor
        prev_cells = clean_sam_chr_table[prev_ref_index]
        key = "_".join(map(str, prev_cells))
        prev_qry_index = Index_by_Qry[key]
        synteny = check_qry_index(prev_qry_index,qry_index)

    else:
        # In middle
        # Compare both
        prev_cells = clean_sam_chr_table[prev_ref_index]
        key = "_".join(map(str, prev_cells))
        prev_qry_index = Index_by_Qry[key]
        Prev_synteny = check_qry_index(prev_qry_index,qry_index)

        next_cells = clean_sam_chr_table[next_ref_index]
        key = "_".join(map(str, next_cells))
        next_qry_index = Index_by_Qry[key]
        next_synteny = check_qry_index(qry_index,next_qry_index)

        if Prev_synteny != "No" or next_synteny != "No":
            if Prev_synteny != "No":
                synteny = Prev_synteny
            else:
                synteny = next_synteny
        else:
            synteny = "No"

    return synteny


def check_qry_index(qry_index,next_qry_index):
    if next_qry_index - qry_index == 1:
        synteny = "Forward"
    elif next_qry_index - qry_index == -1:
        synteny = "Reverse"
    else:
        synteny = "No"
    return synteny


def covered_check(covered,sta,end,sta2,end2):
    covered = "No"
    dup_check = ''
    Sta = min(int(sta),int(end))
    End = max(int(sta),int(end))
    Sta2 = min(int(sta2),int(end2))
    End2 = max(int(sta2),int(end2))
    if  Sta <= Sta2 <= End or Sta <= End2 <= End or Sta2 <= Sta <= End2 or Sta2 <= End <= End2:
        # overlaping
        this_sit  = range(sta, end  + 1)
        other_sit = range(sta2,end2 + 1)
        overlaping_size = len(set(this_sit).intersection(other_sit))
        if float(overlaping_size) / len(this_sit) > 0.5:
            covered = "Yes"
            if float(overlaping_size) / len(this_sit) > 0.9 and float(overlaping_size) / len(other_sit) > 0.9:
                dup_check = "Yes"
    return covered,dup_check

def creat_BND_id(bnd_order,up_good_Q_index,this_anchor,down_good_Q_index):
    if up_good_Q_index < this_anchor[9] < down_good_Q_index:
        # on site sv
        sv_id = "cBND_" + str(bnd_order)
    else:
        # simple BND
        sv_id = "BND_" + str(bnd_order)
    return sv_id


def get_megaSV_flanking(up_ref_index,ref_chr):
    if up_ref_index in sam_chr_SV_R_Index_By_chr[ref_chr]:
        # Not beginning of chr
        up_anchor = sam_chr_SV_By_Ref_index[up_ref_index]
        left_R_flanking = up_anchor[0] + ":" + str(up_anchor[1]) + "-" + str(up_anchor[2])
        left_Q_flanking = up_anchor[4] + ":" + str(up_anchor[5]) + "-" + str(up_anchor[6])
    else:
        # beginning of chr
        left_R_flanking = "None"
        left_Q_flanking = "None"
        up_ref_index = "None"
    return up_ref_index,left_R_flanking,left_Q_flanking



''' Step 1.0 Get tables '''
# count anchor
ref_count = {}
qry_count = {}
raw_table = []
raw_table_index = {}

diff_chr_count = {}


i = -1
with open(Input_coords) as input_file:
    for line in input_file:
        if not line.startswith("ref_start"):
            i += 1
            cells = line.strip().split("\t")
            out_cells = [cells[6],int(cells[0]),int(cells[1]),int(cells[1]) - int(cells[0]) + 1,cells[7],int(cells[2]),int(cells[3]), int(cells[3]) - int(cells[2]) + 1]
            raw_table.append(out_cells)
            key = "_".join(map(str, out_cells))
            raw_table_index[key] = i
            ref = "_".join([cells[6],cells[0],cells[1]])
            qry = "_".join([cells[7],cells[2],cells[3]])
            if cells[6][-4:] != "ch00" and cells[7][-4:] != "ch00":
                if cells[6][-4:] == cells[7][-4:]:
                    # count anchor on same chr
                    count_item(ref,ref_count)
                    count_item(qry,qry_count)
                else:
                    # count anchor on diff chr
                    count_item(ref,diff_chr_count)




""" Only keep unique anchors """

whole_table = []
sam_chr_table = []
dif_chr_table = []
redundancy_dict = {}


with open(Input_coords) as input_file:
    for line in input_file:
        if not line.startswith("ref_start"):
            cells = line.strip().split("\t")
            ref = "_".join([cells[6],cells[0],cells[1]])
            qry = "_".join([cells[7],cells[2],cells[3]])
            out_cells = [cells[6],int(cells[0]),int(cells[1]),int(cells[1]) - int(cells[0]) + 1,cells[7],int(cells[2]),int(cells[3]), int(cells[3]) - int(cells[2]) + 1]
            key = "_".join(map(str, out_cells))
            whole_table.append(out_cells)

            if cells[6][-4:] != "ch00" and cells[7][-4:] != "ch00":
                if cells[6][-4:] == cells[7][-4:]:
                    # same chr
                    if ref_count[ref] == 1:
                        # unique
                        sam_chr_table.append(out_cells)
                    else:
                        # not unique on same chr
                        add_to_dict(ref,out_cells,redundancy_dict)
                else:
                    if ref not in ref_count and diff_chr_count[ref] == 1:
                        # do not have hit on same chr
                        # only 1 hit one diff chr
                        dif_chr_table.append(out_cells)
                        #dif_chr.write(output)


""" check the anchor covered by neighbour """
i = -1
covered_by_neighbour_index = []
almost_dup_indexes = []

for out_cells in sam_chr_table:
    i += 1
    for step in range(1,10):
        if i < len(sam_chr_table) - step:
            # if one anchor covered by its neighbour, keep its index
            covered_by_neighbour,dup_pair = check_covered_by_neighbour(out_cells,i,step)
            if covered_by_neighbour != "None":
                covered_by_neighbour_index += covered_by_neighbour
                if dup_pair != []:
                    almost_dup_indexes += dup_pair

almost_dup_table = []
almost_dup_indexes = sorted(list(set(almost_dup_indexes)))
print "#" * 30
print "number of almost_dup_pair anchors from sam_chr_table:" + str(len(almost_dup_indexes))
for index in almost_dup_indexes:
    almost_dup_table.append(sam_chr_table[index])

"""
Output sam_chr table:
"""
clean_sam_chr_table = []

covered_by_neighbour_index = set(covered_by_neighbour_index)
i = -1
for out_cells in sam_chr_table:
    i += 1
    #output = "\t".join(map(str, out_cells)) + "\n"
    if i not in covered_by_neighbour_index:
        clean_sam_chr_table.append(out_cells)
        #sam_chr.write(output)



# sorted by qry coords
clean_sam_chr_table_By_qry = sorted(clean_sam_chr_table, key = itemgetter(4,5,6))
Index_by_Qry = {}

i = -1
for out_cells in clean_sam_chr_table_By_qry:
    i += 1
    key = "_".join(map(str, out_cells))
    Index_by_Qry[key] = i


clean_sam_chr_table_With_orders = []

i = -1
for out_cells in clean_sam_chr_table:
    i += 1
    key = "_".join(map(str, out_cells))
    ref_index = i
    qry_index = Index_by_Qry[key]
    synteny = synteny_check(out_cells,ref_index,qry_index)
    output_list = out_cells + [ref_index,qry_index,synteny]
    clean_sam_chr_table_With_orders.append(output_list)

""" Synteny_check """
clean_sam_chr_table_With_orders_Syn = []

for output_list in clean_sam_chr_table_With_orders:
    if output_list[7] < 0:
        strand = "minus"
    else:
        strand = "plus"

    ref_index = output_list[8]
    qry_index = output_list[9]
    prev_ref_index = ref_index - 1
    next_ref_index = ref_index + 1

    if output_list[10] == "No":
        prev_ref_index = ref_index - 1
        next_ref_index = ref_index + 1
        if prev_ref_index in clean_sam_chr_table_With_orders and next_ref_index in clean_sam_chr_table_With_orders:
            prev_qry_index = clean_sam_chr_table_With_orders[prev_ref_index][9]
            next_qry_index = clean_sam_chr_table_With_orders[next_ref_index][9]

            if not (ref_index == 0 or output_list[0] != clean_sam_chr_table[prev_ref_index][0]) and not (ref_index == len(clean_sam_chr_table) - 1 or output_list[0] != clean_sam_chr_table[next_ref_index][0]):
                # not 1st & not last
                prev_synteny = clean_sam_chr_table_With_orders[prev_ref_index][10]
                next_synteny = clean_sam_chr_table_With_orders[next_ref_index][10]
                if prev_synteny != "No" and next_synteny != "No":
                    if prev_qry_index < qry_index < next_qry_index or prev_qry_index > qry_index > next_qry_index:
                        if prev_synteny == next_synteny:
                            output_list[10] = prev_synteny


    Note = output_list[10] + "_" + strand
    output_list.append(Note)

    clean_sam_chr_table_With_orders_Syn.append(output_list)


""" block check """
clean_sam_chr_table_With_orders_Syn_Block = []

Synteny_Block = 0
Note_block = 0
Synteny_Block_count = {}
Note_block_count = {}

for output_list in clean_sam_chr_table_With_orders_Syn:
    ref_index = output_list[8]
    qry_index = output_list[9]
    prev_ref_index = ref_index - 1
    if ref_index == 0 or output_list[0] != clean_sam_chr_table_With_orders_Syn[prev_ref_index][0]:
        # 1st anchor for each chr
        # start a new block at the beginning of each chr
        Synteny_Block += 1
        Note_block += 1
    else:
        if output_list[10] != clean_sam_chr_table_With_orders_Syn[prev_ref_index][10]:
            Synteny_Block += 1
        if output_list[11] != clean_sam_chr_table_With_orders_Syn[prev_ref_index][11]:
            Note_block += 1

    count_item(Synteny_Block,Synteny_Block_count)
    count_item(Note_block,Note_block_count)
    output_list += [Synteny_Block,Note_block]
    clean_sam_chr_table_With_orders_Syn_Block.append(output_list)


""" Keep index for each blocked anchor """
blocked_anchor_index_dict = {}

i = -1
for output_list in clean_sam_chr_table_With_orders_Syn_Block:
    i += 1
    Synteny_Block = output_list[12]
    Note_block = output_list[13]
    output_list += [Synteny_Block_count[Synteny_Block], Note_block_count[Note_block]]
    key = "_".join(map(str, output_list[0:8]))
    blocked_anchor_index_dict[key] = i


''' Get back redundant anchors on within good block '''
redundant_ref_Keys = sorted(redundancy_dict.keys())
rescued_anchors = {}
rescued_anchors_count = {}

for ref_Key in redundant_ref_Keys:
    redundant_anchors = redundancy_dict[ref_Key]
    for redundant_anchor in redundant_anchors:
        redundant_key = "_".join(map(str, redundant_anchor))
        # get index in raw tables
        index_in_raw_table = raw_table_index[redundant_key]
        # get its neighbour good anchors
        prev_index_in_raw_table = index_in_raw_table - 1
        next_index_in_raw_table = index_in_raw_table + 1
        if prev_index_in_raw_table in raw_table and next_index_in_raw_table in raw_table:
            while "_".join(map(str, raw_table[prev_index_in_raw_table])) not in blocked_anchor_index_dict or \
                   clean_sam_chr_table_With_orders_Syn_Block[blocked_anchor_index_dict["_".join(map(str, raw_table[prev_index_in_raw_table]))]][10] == "No":
                prev_index_in_raw_table -= 1

            while "_".join(map(str, raw_table[next_index_in_raw_table])) not in blocked_anchor_index_dict or \
                   clean_sam_chr_table_With_orders_Syn_Block[blocked_anchor_index_dict["_".join(map(str, raw_table[next_index_in_raw_table]))]][10] == "No":
                next_index_in_raw_table += 1

            prev_good_anchor_Index = blocked_anchor_index_dict["_".join(map(str, raw_table[prev_index_in_raw_table]))]
            next_good_anchor_Index = blocked_anchor_index_dict["_".join(map(str, raw_table[next_index_in_raw_table]))]

            cover_by_good_anchor = check_cover_by_good_anchor(redundant_anchor,prev_good_anchor_Index,next_good_anchor_Index)

            if cover_by_good_anchor == "None":
                # Not covered by good anchor
                prev_good_anchor = clean_sam_chr_table_With_orders_Syn_Block[prev_good_anchor_Index]
                next_good_anchor = clean_sam_chr_table_With_orders_Syn_Block[next_good_anchor_Index]
                if prev_good_anchor[11] == next_good_anchor[11] and prev_good_anchor[0] == next_good_anchor[0]:
                    # same blocktype and same chr
                    # check within good block
                    if (prev_good_anchor[11] == "Forward_plus"  and prev_good_anchor[5] < redundant_anchor[5] < next_good_anchor[5]) or \
                       (prev_good_anchor[11] == "Reverse_minus" and prev_good_anchor[5] > redundant_anchor[5] > next_good_anchor[5]):
                        add_to_dict(prev_good_anchor_Index,redundant_anchor,rescued_anchors)
                        count_item(ref_Key,rescued_anchors_count)



""" Get final good anchors """
Final_good_keys = set()            # Good uniqe same chr anchor

kept_anchor_by_Ref_chr = {}        # Good uniqe same chr anchor & Dupliccation. For filter anchor related ch00
kept_anchor_by_Qry_chr = {}        # Good uniqe same chr anchor & Dupliccation. For filter anchor related ch00

Duplication_table = []

with open(Prefix + ".Duplication.Exact.txt", "w") as duplication:
    output_list = ['Ref_chr','Ref_start','Ref_end','Ref_size','Qry_chr','Qry_start','Qry_end','Qry_size']
    output = "\t".join(map(str, output_list)) + "\n"
    duplication.write(output)
    i = -1
    for output_list in clean_sam_chr_table_With_orders_Syn_Block:
        i += 1
        key = "_".join(map(str, output_list[0:8]))
        Final_good_keys.add(key)
        add_to_dict(output_list[0],output_list[0:8],kept_anchor_by_Ref_chr)
        add_to_dict(output_list[4],output_list[0:8],kept_anchor_by_Qry_chr)
        if i in rescued_anchors:
            for redundant_anchor in rescued_anchors[i]:
                ref_Key = "_".join(map(str, redundant_anchor[0:3]))
                add_to_dict(redundant_anchor[0],redundant_anchor,kept_anchor_by_Ref_chr)
                add_to_dict(redundant_anchor[4],redundant_anchor,kept_anchor_by_Qry_chr)
                if ref_Key in rescued_anchors_count and rescued_anchors_count[ref_Key] == 1:
                    key = "_".join(map(str, redundant_anchor))
                    Final_good_keys.add(key)
                elif ref_Key in rescued_anchors_count and rescued_anchors_count[ref_Key] > 1:
                    output = "\t".join(map(str, redundant_anchor)) + "\n"
                    duplication.write(output)
                    Duplication_table.append(redundant_anchor)


"""
Step 1.1 Get clean_sam_chr_table
Output clean sam_chr table
    ------- Again -------
 """
clean_sam_chr_table = []

with open(Prefix + ".Unique_anchor.sam_chr", "w") as sam_chr:
    output_list = ['Ref_chr','Ref_start','Ref_end','Ref_size','Qry_chr','Qry_start','Qry_end','Qry_size']
    output = "\t".join(map(str, output_list)) + "\n"
    sam_chr.write(output)
    for out_cells in raw_table:
        key = "_".join(map(str, out_cells))
        if key in Final_good_keys:
            clean_sam_chr_table.append(out_cells)
            output = "\t".join(map(str, out_cells)) + "\n"
            sam_chr.write(output)

"""sorted by qry coords"""
clean_sam_chr_table_By_qry = sorted(clean_sam_chr_table, key = itemgetter(4,5,6))
Index_by_Qry = {}

i = -1
for out_cells in clean_sam_chr_table_By_qry:
    i += 1
    key = "_".join(map(str, out_cells))
    Index_by_Qry[key] = i


clean_sam_chr_table_With_orders = []

i = -1
for out_cells in clean_sam_chr_table:
    i += 1
    key = "_".join(map(str, out_cells))
    ref_index = i
    qry_index = Index_by_Qry[key]
    synteny = synteny_check(out_cells,ref_index,qry_index)
    output_list = out_cells + [ref_index,qry_index,synteny]
    clean_sam_chr_table_With_orders.append(output_list)

""" Synteny_check """
clean_sam_chr_table_With_orders_Syn = []

for output_list in clean_sam_chr_table_With_orders:
    if output_list[7] < 0:
        strand = "minus"
    else:
        strand = "plus"
    ref_index = output_list[8]
    qry_index = output_list[9]
    prev_ref_index = ref_index - 1
    next_ref_index = ref_index + 1

    if output_list[10] == "No":
        prev_ref_index = ref_index - 1
        next_ref_index = ref_index + 1
        if prev_ref_index in clean_sam_chr_table_With_orders and next_ref_index in clean_sam_chr_table_With_orders:
            prev_qry_index = clean_sam_chr_table_With_orders[prev_ref_index][9]
            next_qry_index = clean_sam_chr_table_With_orders[next_ref_index][9]
            if prev_qry_index in clean_sam_chr_table and next_qry_index in clean_sam_chr_table:
                if not (ref_index == 0 or output_list[0] != clean_sam_chr_table[prev_ref_index][0]) and not (ref_index == len(clean_sam_chr_table) - 1 or output_list[0] != clean_sam_chr_table[next_ref_index][0]):
                    # not 1st & not last
                    prev_synteny = clean_sam_chr_table_With_orders[prev_ref_index][10]
                    next_synteny = clean_sam_chr_table_With_orders[next_ref_index][10]
                    if prev_synteny != "No" and next_synteny != "No":
                        if prev_qry_index < qry_index < next_qry_index or prev_qry_index > qry_index > next_qry_index:
                            if prev_synteny == next_synteny:
                                output_list[10] = prev_synteny


    Note = output_list[10] + "_" + strand
    output_list.append(Note)
    clean_sam_chr_table_With_orders_Syn.append(output_list)



""" block check """
clean_sam_chr_table_With_orders_Syn_Block = []

Synteny_Block = 0
Note_block = 0
Synteny_Block_count = {}
Note_block_count = {}

for output_list in clean_sam_chr_table_With_orders_Syn:
    ref_index = output_list[8]
    qry_index = output_list[9]
    prev_ref_index = ref_index - 1
    if ref_index == 0 or output_list[0] != clean_sam_chr_table_With_orders_Syn[prev_ref_index][0]:
        # 1st anchor for each chr
        # start a new block at the beginning of each chr
        Synteny_Block += 1
        Note_block += 1
    else:
        if output_list[10] != clean_sam_chr_table_With_orders_Syn[prev_ref_index][10]:
            Synteny_Block += 1
        if output_list[11] != clean_sam_chr_table_With_orders_Syn[prev_ref_index][11]:
            Note_block += 1

    count_item(Synteny_Block,Synteny_Block_count)
    count_item(Note_block,Note_block_count)
    output_list += [Synteny_Block,Note_block]
    clean_sam_chr_table_With_orders_Syn_Block.append(output_list)


"""
Step 1.1.1
      clean_sam_chr_table_With_orders_Syn_Block
         count each blocked anchor
"""
blocked_anchor_index_dict = {}

i = -1
for output_list in clean_sam_chr_table_With_orders_Syn_Block:
    i += 1
    Synteny_Block = output_list[12]
    Note_block = output_list[13]
    output_list += [Synteny_Block_count[Synteny_Block], Note_block_count[Note_block]]
    #output = "\t".join(map(str, output_list)) + "\n"
    #output_file.write(output)
    key = "_".join(map(str, output_list[0:8]))
    blocked_anchor_index_dict[key] = i


"""
Step 1.2  Clean dif_chr_table

"""
dif_chr_table_keys = {}
i = -1
for out_cells in dif_chr_table:
    i += 1
    key = "_".join(map(str, out_cells))
    dif_chr_table_keys[key] = i

clean_dif_chr_table = []
almost_dif_chr_dup_key_pair = []
covered_by_better_Evi_keys = set()

with open(Prefix + ".Unique_anchor.dif_chr", "w") as output_file:
    output_list = ['Ref_chr','Ref_start','Ref_end','Ref_size','Qry_chr','Qry_start','Qry_end','Qry_size']
    output = "\t".join(map(str, output_list)) + "\n"
    output_file.write(output)
    for out_cells in dif_chr_table:
        key = "_".join(map(str, out_cells))
        covered = "No"
        overlaped = ""
        ref_chr = out_cells[0]
        ref_sta = out_cells[1]
        ref_end = out_cells[2]
        qry_chr = out_cells[4]
        qry_sta = min(out_cells[5],out_cells[6])
        qry_end = max(out_cells[5],out_cells[6])
        check_ahchor_keys = set()

        for anchor in dif_chr_table + kept_anchor_by_Ref_chr[ref_chr] + kept_anchor_by_Qry_chr[qry_chr] + almost_dup_table:
            #ahchor_key = "_".join(map(str, anchor[0-8]))
            ahchor_key = "_".join(map(str, anchor))
            if key != ahchor_key and ahchor_key not in check_ahchor_keys:
                check_ahchor_keys.add(ahchor_key)
                ref_chr2 = anchor[0]
                ref_sta2 = anchor[1]
                ref_end2 = anchor[2]
                qry_chr2 = anchor[4]
                qry_sta2 = min(anchor[5],anchor[6])
                qry_end2 = max(anchor[5],anchor[6])
                if ref_chr == ref_chr2:
                    covered,dup_check = covered_check(covered,ref_sta,ref_end,ref_sta2,ref_end2)
                    if covered == "Yes":
                        overlaped = "Yes"
                        if ahchor_key not in dif_chr_table_keys:
                            # covered by other good Evidence
                            covered_by_better_Evi_keys.add(key)
                        elif dup_check == "Yes" and key in dif_chr_table_keys and ahchor_key in dif_chr_table_keys:
                            almost_dif_chr_dup_key_pair.append([key,ahchor_key])
                        break
                elif covered == "No" and qry_chr == qry_chr2:
                    covered,dup_check = covered_check(covered,qry_sta,qry_end,qry_sta2,qry_end2)
                    if covered == "Yes":
                        overlaped = "Yes"
                        if ahchor_key not in dif_chr_table_keys:
                            # covered by other good Evidence
                            covered_by_better_Evi_keys.add(key)
                        elif dup_check == "Yes" and key in dif_chr_table_keys and ahchor_key in dif_chr_table_keys:
                            almost_dif_chr_dup_key_pair.append([key,ahchor_key])
                        break
                if overlaped == "Yes":
                    break
        if overlaped == "":
            clean_dif_chr_table.append(out_cells)
            output = "\t".join(map(str, out_cells)) + "\n"
            output_file.write(output)


almost_dif_chr_dup_index = set()
good = 0
for keys in almost_dif_chr_dup_key_pair:
    if keys[0] not in covered_by_better_Evi_keys and keys[1] not in covered_by_better_Evi_keys:
        good += 1
        almost_dif_chr_dup_index.add(dif_chr_table_keys[keys[0]])
        almost_dif_chr_dup_index.add(dif_chr_table_keys[keys[1]])


"""
Step 1.2.1  almost_dup_table
From same chr OR diff chr
No ch00
"""

print "almost_dup anchors number in dif_chr_table: "  + str(len(almost_dif_chr_dup_index))

for index in sorted(list(almost_dif_chr_dup_index)):
    almost_dup_table.append(dif_chr_table[index])

#print "Check almost_dup_table"
almost_dup_table = sorted(almost_dup_table, key = itemgetter(0,1,2,3))

with open(Prefix + ".Duplication.cov90.txt", "w") as output:
    for anchor in almost_dup_table:
        output.write("\t".join(map(str, anchor)) + "\n")
        # add to kept_anchor_by_Ref_chr
        #        kept_anchor_by_Qry_chr
        # for filter ch00 anchor
        add_to_dict(anchor[0],anchor,kept_anchor_by_Ref_chr)
        add_to_dict(anchor[4],anchor,kept_anchor_by_Qry_chr)


""" Get ch00 related anchors
    filter by:
       kept_anchor_by_Ref_chr
       kept_anchor_by_Qry_chr
"""

ch00_anchors = []

for out_cells in raw_table:
    if "ch00" in out_cells[0] and "ch00" in out_cells[4]:
        # both on ch00
        ch00_anchors.append(out_cells)
    elif "ch00" in out_cells[4]:
        # qry is ch00
        # check whether ref is covered by non-ch00 anchors
        ref_chr = out_cells[0]
        ref_sta = out_cells[1]
        ref_end = out_cells[2]
        overlapped = "No"
        for good_anchor in kept_anchor_by_Ref_chr[ref_chr]:
            good_sta = good_anchor[1]
            good_end = good_anchor[2]
            if  good_sta <= ref_sta <= good_end or good_sta <= ref_end <= good_end or ref_sta <= good_sta <= ref_end or ref_sta <= good_end <= ref_end:
                # overlapping
                ref_sit = range(ref_sta,ref_end+1)
                ref_len = len(ref_sit)
                good_sit = range(good_sta,good_end + 1)
                overlaping_size = len(set(ref_sit).intersection(good_sit))
                if float(overlaping_size) / ref_len > 0.3:
                    overlapped = "Yes"
                    break
        if overlapped == "No":
            ch00_anchors.append(out_cells)
    elif "ch00" in out_cells[0]:
        # ref in ch00
        # check whether qry is covered by non-ch00 anchors
        qry_chr = out_cells[4]
        qry_sta = min(out_cells[5],out_cells[6])
        qry_end = max(out_cells[5],out_cells[6])
        overlapped = "No"
        for good_anchor in kept_anchor_by_Qry_chr[qry_chr]:
            good_sta = min(good_anchor[5],good_anchor[6])
            good_end = max(good_anchor[5],good_anchor[6])
            if  good_sta <= qry_sta <= good_end or good_sta <= qry_end <= good_end or qry_sta <= good_sta <= qry_end or qry_sta <= good_end <= qry_end:
                # overlapping
                qry_sit = range(qry_sta,qry_end+1)
                qry_len = len(qry_sit)
                good_sit = range(good_sta,good_end + 1)
                overlaping_size = len(set(qry_sit).intersection(good_sit))
                if float(overlaping_size) / qry_len > 0.3:
                    overlapped = "Yes"
                    break
        if overlapped == "No":
            ch00_anchors.append(out_cells)

'''
 Step 1.3 Get clean_ch00_anchors
   check covered by other ch00_anchors
   if any one coveded > 0.3 by others, removed it
'''
clean_ch00_anchors = []

with open(Prefix + ".Unique_anchor.ch00", "w") as output_file:
    output_list = ['Ref_chr','Ref_start','Ref_end','Ref_size','Qry_chr','Qry_start','Qry_end','Qry_size']
    output = "\t".join(map(str, output_list)) + "\n"
    output_file.write(output)
    for out_cells in ch00_anchors:
        key = "_".join(map(str, out_cells))
        covered = "No"
        overlaped = ""
        ref_chr = out_cells[0]
        ref_sta = out_cells[1]
        ref_end = out_cells[2]
        qry_chr = out_cells[4]
        qry_sta = min(out_cells[5],out_cells[6])
        qry_end = max(out_cells[5],out_cells[6])
        for anchor in ch00_anchors + clean_dif_chr_table:
            if key != "_".join(map(str, anchor)):
                ref_chr2 = anchor[0]
                ref_sta2 = anchor[1]
                ref_end2 = anchor[2]
                qry_chr2 = anchor[4]
                qry_sta2 = min(anchor[5],anchor[6])
                qry_end2 = max(anchor[5],anchor[6])
                if ref_chr == ref_chr2:
                    covered,dup_check = covered_check(covered,ref_sta,ref_end,ref_sta2,ref_end2)
                    if covered == "Yes":
                        overlaped = "Yes"
                        break
                elif covered == "No" and qry_chr == qry_chr2:
                    covered,dup_check = covered_check(covered,qry_sta,qry_end,qry_sta2,qry_end2)
                    if covered == "Yes":
                        overlaped = "Yes"
                        break
                if overlaped == "Yes":
                    break
        if  overlaped == "":
            clean_ch00_anchors.append(out_cells)
            output = "\t".join(map(str, out_cells)) + "\n"
            output_file.write(output)



"""
#########################################################
All reliable anchors for SV check

      clean_sam_chr_table     Step 1.1
      clean_dif_chr_table     Step 1.2
      clean_ch00_anchors      Step 1.3


print len(clean_sam_chr_table)
print clean_sam_chr_table[0]
print len(clean_dif_chr_table)
print clean_dif_chr_table[0]
print len(clean_ch00_anchors)
print clean_ch00_anchors[0]

"""


'''
  Step 2.0
  output INV and BND based on good unique anchor on same chr

  Evidence:
      BND       clean_dif_chr_table
                clean_sam_chr_table_With_orders_Syn_Block      Step 1.1.1

      INV       clean_sam_chr_table_With_orders_Syn_Block      Step 1.1.1


print len(clean_sam_chr_table_With_orders_Syn_Block)
print clean_sam_chr_table_With_orders_Syn_Block[0]
print clean_sam_chr_table_With_orders_Syn_Block[-1]

17484
['NXGLv5ch01', 4860, 9044, 4185, 'BOv5ch01', 19649, 23846, 4198, 0, 0, 'Forward', 'Forward_plus', 1, 1, 159, 159]
['NXGLv5ch09', 66422010, 66423117, 1108, 'BOv5ch09', 64417967, 64419074, 1108, 17483, 17483, 'Forward', 'Forward_plus', 371, 455, 74, 74]

'''


# split clean_sam_chr_table_With_orders_Syn_Block into each chr
clean_sam_chr_table_Block_by_R_chr = {}
clean_sam_chr_table_Block_by_R_index = {}
clean_sam_chr_table_Block_by_Q_chr = {}
clean_sam_chr_table_Block_by_Q_index = {}


for output_list in clean_sam_chr_table_With_orders_Syn_Block:
    ref_chr = output_list[0]
    qry_chr = output_list[4]
    add_to_dict(ref_chr,output_list,clean_sam_chr_table_Block_by_R_chr)
    add_to_dict(qry_chr,output_list,clean_sam_chr_table_Block_by_Q_chr)
    ref_index = output_list[8]
    qry_index = output_list[9]
    clean_sam_chr_table_Block_by_R_index[ref_index] = output_list
    clean_sam_chr_table_Block_by_Q_index[qry_index] = output_list


sv_order = 0
checked_sv_indexes = set()
inv_order = 0
bnd_order = 0
complex_order = 0

clean_sam_chr_SV = []


for output_list in clean_sam_chr_table_With_orders_Syn_Block:
    ref_chr = output_list[0]
    the_Q_index = output_list[9]
    this_R_index = output_list[8]
    next_Q_index = the_Q_index + 1
    while next_Q_index in clean_sam_chr_table_Block_by_Q_index and clean_sam_chr_table_Block_by_Q_index[next_Q_index][11] !=  "Forward_plus":
        next_Q_index += 1
    if this_R_index not in checked_sv_indexes and output_list[11] == "Forward_plus":
        checked_sv_indexes.add(this_R_index)
        if next_Q_index in clean_sam_chr_table_Block_by_Q_index and clean_sam_chr_table_Block_by_Q_index[next_Q_index][0] == ref_chr:
            next_Qry_R_index = clean_sam_chr_table_Block_by_Q_index[next_Q_index][8]
        else:
            next_Qry_R_index = 1000000000
        distance = next_Qry_R_index - this_R_index
        #output = "\t".join(map(str, output_list)) + "\t" + str(distance) + "\n"
        #output_file.write(output)
        clean_sam_chr_SV.append(output_list + [distance])
        if distance > 1 and next_Qry_R_index != 1000000000:
            if this_R_index < len(clean_sam_chr_table_With_orders_Syn_Block) - 1:
                neighbour_anchor = clean_sam_chr_table_With_orders_Syn_Block[this_R_index + 1]
                if ref_chr == neighbour_anchor[0] and neighbour_anchor[11] != "Forward_plus":
                    # found SV region
                    # 1. something is found between this and next qry region
                    # 2. not last anchor in chr
                    # 3. Not in a continual syytenic region

                    sv_anchor_indexes = range(this_R_index + 1,next_Qry_R_index)
                    up_good_R_index = sv_anchor_indexes[0] - 1
                    up_good_Q_index = clean_sam_chr_table_Block_by_R_index[up_good_R_index][9]
                    down_good_R_index = sv_anchor_indexes[-1] + 1
                    down_good_Q_index = clean_sam_chr_table_Block_by_R_index[down_good_R_index][9]

                    block_types = set()
                    sv_id_dict = {}
                    for index in sv_anchor_indexes:
                        checked_sv_indexes.add(index)
                        block_types.add(clean_sam_chr_table_With_orders_Syn_Block[index][11])
                    if len(block_types) > 1:
                        complex_order += 1
                        complex_id = "ComplexSV_" + str(complex_order)
                        take_out_index = set()
                        for index in sv_anchor_indexes:
                            if index not in take_out_index:
                                prev_anchor = clean_sam_chr_table_With_orders_Syn_Block[index-1]
                                this_anchor = clean_sam_chr_table_With_orders_Syn_Block[index]
                                next_anchor = clean_sam_chr_table_With_orders_Syn_Block[index+1]
                                if this_anchor[11] == "Forward_minus":
                                    if this_anchor[15] == 1:
                                        if (this_anchor[9] - prev_anchor[9] == 1  and index - 1 not in sv_anchor_indexes) or (this_anchor[9] - next_anchor[9] == -1 and index + 1 not in sv_anchor_indexes):
                                            inv_order += 1
                                            sv_id = "INV_" + str(inv_order)
                                            take_out_index.add(index)
                                        else:
                                            bnd_order += 1
                                            sv_id = creat_BND_id(bnd_order,up_good_Q_index,this_anchor,down_good_Q_index)
                                            sv_id_dict[index] = sv_id
                                            take_out_index.add(index)
                                    else:
                                        for index in range(index,index+this_anchor[15]):
                                            bnd_order += 1
                                            sv_id = creat_BND_id(bnd_order,up_good_Q_index,this_anchor,down_good_Q_index)
                                            sv_id_dict[index] = sv_id
                                            take_out_index.add(index)
                                elif this_anchor[10] == "No" or this_anchor[11] == "Reverse_plus":
                                    bnd_order += 1
                                    sv_id = creat_BND_id(bnd_order,up_good_Q_index,this_anchor,down_good_Q_index)
                                    sv_id_dict[index] = sv_id
                                    take_out_index.add(index)
                                elif this_anchor[11] == "Forward_plus":
                                    bnd_order += 1
                                    sv_id = creat_BND_id(bnd_order,up_good_Q_index,this_anchor,down_good_Q_index)
                                    for index in range(index,index+this_anchor[15]):
                                        sv_id_dict[index] = sv_id
                                        take_out_index.add(index)
                                else:
                                    sv_id = complex_id
                                sv_id_dict[index] = sv_id

                        if len(take_out_index) > 0:
                            # have removed head or tail Forward_minus
                            new_block_types = set()
                            new_sv_anchor_indexes = []
                            for index in sv_anchor_indexes:
                                if index not in take_out_index:
                                    new_sv_anchor_indexes.append(index)
                            for index in new_sv_anchor_indexes:
                                new_block_types.add(clean_sam_chr_table_With_orders_Syn_Block[index][11])
                            if len(new_block_types) == 1:
                                if list(new_block_types)[0] == "Reverse_minus":
                                    # merge all on site Inv together
                                    confirmed_inv_indexs = set()
                                    for index in new_sv_anchor_indexes:
                                        if index not in confirmed_inv_indexs:
                                            check_anchor = clean_sam_chr_table_With_orders_Syn_Block[index]
                                            if up_good_Q_index <  check_anchor[9] < down_good_Q_index:
                                                inv_order += 1
                                                sv_id = "cINV_" + str(inv_order)
                                                for index in range(index,index+check_anchor[15]):
                                                    sv_id_dict[index] = sv_id
                                                    confirmed_inv_indexs.add(index)
                                            else:
                                                bnd_order += 1
                                                sv_id = "BND_" + str(bnd_order)
                                                for index in range(index,index+check_anchor[15]):
                                                    sv_id_dict[index] = sv_id
                                                    confirmed_inv_indexs.add(index)
                                else:
                                    bnd_order += 1
                                    sv_id = "Problem_" + str(bnd_order)
                                    for index in new_sv_anchor_indexes:
                                        sv_id_dict[index] = sv_id
                    elif len(block_types) == 1:
                        if list(block_types)[0] == "Reverse_minus":
                            confirmed_inv_indexs = set()
                            for index in sv_anchor_indexes:
                                if index not in confirmed_inv_indexs:
                                    check_anchor = clean_sam_chr_table_With_orders_Syn_Block[index]
                                    if up_good_Q_index <  check_anchor[9] < down_good_Q_index:
                                        inv_order += 1
                                        sv_id = "INV_" + str(inv_order)
                                        for index in range(index,index+check_anchor[15]):
                                            sv_id_dict[index] = sv_id
                                            confirmed_inv_indexs.add(index)
                                    else:
                                        bnd_order += 1
                                        sv_id = "BND_" + str(bnd_order)
                                        for index in range(index,index+check_anchor[15]):
                                            sv_id_dict[index] = sv_id
                                            confirmed_inv_indexs.add(index)
                        elif list(block_types)[0] == "Reverse_plus":
                            for index in sv_anchor_indexes:
                                bnd_order += 1
                                sv_id = "cBND_" + str(bnd_order)
                                sv_id_dict[index] = sv_id
                        else:
                            bnd_order += 1
                            sv_id = "BND_" + str(bnd_order)
                            for index in sv_anchor_indexes:
                                sv_id_dict[index] = sv_id

                    for index in sv_anchor_indexes:
                        sv_id = sv_id_dict[index]
                        #output = "\t".join(map(str, clean_sam_chr_table_With_orders_Syn_Block[index])) + "\tNot_check\t" + sv_id + "\n"
                        #output_file.write(output)
                        clean_sam_chr_SV.append(clean_sam_chr_table_With_orders_Syn_Block[index] + ["Not_check",sv_id])
    elif this_R_index not in checked_sv_indexes:
        block_size = output_list[15]
        if output_list[11] == "Reverse_minus":
            inv_order += 1
            sv_id = "INV_" + str(inv_order)
        elif output_list[10] == "No":
            bnd_order += 1
            sv_id = "BND_" + str(bnd_order)
        elif output_list[11] == "Forward_minus" and  output_list[15] == 1:
            inv_order += 1
            sv_id = "INV_" + str(inv_order)
        else:
            sv_id = "NotSure"
        sv_anchor_indexes = range(this_R_index,this_R_index + block_size)
        for index in sv_anchor_indexes:
            checked_sv_indexes.add(index)
            #output = "\t".join(map(str, clean_sam_chr_table_With_orders_Syn_Block[index])) + "\tNot_check\t" + sv_id + "\n"
            #output_file.write(output)
            clean_sam_chr_SV.append(clean_sam_chr_table_With_orders_Syn_Block[index] + ["Not_check",sv_id])



"""
Step 2.1
     Get Check clean_sam_chr_SV
Manual correct simple inv "Forward_minus
"""


for output_list in clean_sam_chr_SV:
    ref_chr = output_list[0]
    this_R_index = output_list[8]
    this_Q_index = output_list[9]
    if 0 < this_R_index < len(clean_sam_chr_SV) - 1 and   output_list[11] == "Forward_minus" and output_list[15] == 1:
        prev_anchor = clean_sam_chr_table_With_orders_Syn_Block[this_R_index-1]
        next_anchor = clean_sam_chr_table_With_orders_Syn_Block[this_R_index+1]
        if  prev_anchor[0] == ref_chr == next_anchor[0] and prev_anchor[11] == next_anchor[11] == "Forward_plus" and (prev_anchor[9] + 1 == this_Q_index or next_anchor[9] - 1 ==this_Q_index):
            inv_order += 1
            sv_id = "INV_" + str(inv_order)
            output_list[17] = sv_id




"""
Step 2.2
    Output same_chr_SVs based on:
                             clean_sam_chr_SV
     1. Split into each chr
     2. Split into each SV
     3. Split into each Mega SV regions, for output flanking region
"""
sam_chr_SV_By_chr = {}              # value is a table
sam_chr_SV_R_Index_By_chr = {}      # value is a list of index

sam_chr_SV_By_Ref_index = {}        # value is anchor
sam_chr_SV_By_Qry_index = {}        # value is anchor

sam_chr_SV_Key_in_chr = {}          # value is list of key in each chr, used for check flanking anchors for DUP



# split table
for cells in clean_sam_chr_SV:
    ref_chr = cells[0]
    ref_index = cells[8]
    qry_index = cells[9]
    key = "_".join(map(str, cells[0:8]))
    add_to_dict(ref_chr,cells,sam_chr_SV_By_chr)
    add_to_dict(ref_chr,ref_index,sam_chr_SV_R_Index_By_chr)
    add_to_dict(ref_chr,key,sam_chr_SV_R_Index_By_chr)
    sam_chr_SV_By_Ref_index[ref_index] = cells
    sam_chr_SV_By_Qry_index[qry_index] = cells


# get megaSV
sorted_chr_list = sorted(sam_chr_SV_By_chr.keys())

sam_chr_SV_By_SV = {}               # only SV lines, value is a table
sam_chr_SV_By_megaSV = {}           # only SV lines, value is a table

SV_list = []
SV_vs_megaSV = {}

mega_sv_order = 0


for chr in sorted_chr_list:
    SV_line_check = "No"
    for cells in sam_chr_SV_By_chr[chr]:
        ref_chr = cells[0]
        ref_index = cells[8]
        qry_index = cells[9]
        if len(cells) == 18:
            if SV_line_check == "No":
                # new SV region
                SV_line_check = "Yes"
                mega_sv_order += 1
                mega_sv_id = "Mega_" + str(mega_sv_order)
                sv_id = cells[17]
                cells.append(mega_sv_id)
                add_to_dict(sv_id,cells,sam_chr_SV_By_SV)
                add_to_dict(mega_sv_id,cells,sam_chr_SV_By_megaSV)
                if sv_id not in set(SV_list):
                    SV_list.append(sv_id)
                    SV_vs_megaSV[sv_id] = mega_sv_id
            else:
                # continue SV region
                sv_id = cells[17]
                cells.append(mega_sv_id)
                add_to_dict(sv_id,cells,sam_chr_SV_By_SV)
                add_to_dict(mega_sv_id,cells,sam_chr_SV_By_megaSV)
                if sv_id not in set(SV_list):
                    SV_list.append(sv_id)
                    SV_vs_megaSV[sv_id] = mega_sv_id
        elif len(cells) == 17 and SV_line_check == "Yes":
            SV_line_check = "No"


""" Get flanking region for megaSV """
megaSV_flanking = {}
megaSV_flanking_indexs = {}
for megaSV in sam_chr_SV_By_megaSV.keys():
    sv_lines = sam_chr_SV_By_megaSV[megaSV]
    # get left anchor
    first_anchor = sv_lines[0]
    ref_chr = first_anchor[0]
    qry_chr = first_anchor[4]

    first_ref_index = first_anchor[8]
    up_ref_index = first_ref_index - 1
    up_ref_index,left_R_flanking,left_Q_flanking, = get_megaSV_flanking(up_ref_index,ref_chr)

    last_anchor = sv_lines[-1]
    last_ref_index  = last_anchor[8]
    down_ref_index = last_ref_index + 1
    down_ref_index,right_R_flanking,right_Q_flanking = get_megaSV_flanking(down_ref_index,ref_chr)

    megaSV_flanking[megaSV] = [up_ref_index,down_ref_index,left_R_flanking,left_Q_flanking,right_R_flanking,right_Q_flanking]


'''
Step 2.3
   Get sam_chr_SV_By_chr

Check INV in megaSV region
If it should be transcloation, add "t"


 '''
revised_SV_list = []
revised_SV_flanking = {}
sam_chr_SV_By_SV = {}
SV_vs_megaSV = {}
sam_chr_SV_By_Key = {}   # for DUP flanking

for chr in sorted_chr_list:
    for cells in sam_chr_SV_By_chr[chr]:
        key = "_".join(map(str, cells[0:8]))

        if len(cells) > 17:
            qry_index = cells[9]
            flanking = megaSV_flanking[cells[18]]
            up_ref_index = flanking[0]
            down_ref_index = flanking[1]
            if cells[17].startswith("INV") and ((up_ref_index != "None" and qry_index < sam_chr_SV_By_Ref_index[up_ref_index][9]) or (down_ref_index != "None" and qry_index > sam_chr_SV_By_Ref_index[down_ref_index][9])):
                cells[17] = "t" + cells[17]
            elif cells[17].startswith("BND") and ((up_ref_index != "None" and qry_index > sam_chr_SV_By_Ref_index[up_ref_index][9]) and (down_ref_index != "None" and qry_index < sam_chr_SV_By_Ref_index[down_ref_index][9])):
                cells[17] = "c" + cells[17]
            add_to_dict(cells[17],cells,sam_chr_SV_By_SV)
            if cells[17] not in SV_vs_megaSV:
                revised_SV_flanking[cells[17]] = flanking
                revised_SV_list.append(cells[17])
                SV_vs_megaSV[cells[17]] = cells[18]

        sam_chr_SV_By_Key[key] = cells


'''
Check
print len(megaSV_flanking)
print len(SV_list)

print len(revised_SV_flanking)
print len(sam_chr_SV_By_SV)
print len(SV_vs_megaSV)
print len(revised_SV_list)
'''

"""
########################################################################################
########################################################################################
Step 2.4
    Output my SV list:
     1. BND and INV         Step 2.3
            revised_SV_list = []
            revised_SV_flanking = {}
            sam_chr_SV_By_SV = {}
            sam_chr_SV_By_chr[chr]      Will reivse here ---> Final version   ---> Step 2.5 Output for checking

     2. BND based on        Step 1.2
            clean_dif_chr_table

     3. DUP
            Duplication_table        just before Step 1.1

     4. Assemblytics SVs
            raw_SV                   Input
          Fillter by:
            sam_chr_SV_By_chr[chr]   Here        Both within_alignment and between_alignments
            clean_dif_chr_table      Step 1.2    Both within_alignment and between_alignments
            clean_ch00_anchors       Step 1.3    within_alignment SV only
########################################################################################
########################################################################################
"""
final_SV_list = []
My_SV_table = []

dup_ref_key_count = {}
dup_ref_key_VS_id = {}
dup_order = 0

combined_table = []
combined_dict_By_ID = {}

bad_Assemblytics_id_list = []

letters = "0ABCDEFGHIJKLMNOPQRSTUVWXYZ"

with open(Prefix + ".Other_SV.tsv" , "w") as output_file, open(Prefix + ".Other_detail.txt", "w") as MY_SV_detail:
    output_list = ["SV_type","Ref_chr","Ref_start","Ref_end","Ref_Size","Qry_chr","Qry_start","Qry_end","Qry_Size","Left_Ref_Flank","Right_Ref_Flank","Left_Qry_Flank","Right_Qry_Flank","Assemblytics","Anchor_Number","Ref_aligned_Size/Anchor_strand","Qry_aligned_Size"]
    output = "\t".join(map(str, output_list)) + "\n"
    output_file.write(output)

    """
    (1) Output my BND and INV on same chr
        Based on anchors on same chrs
    """

    for sv_id in revised_SV_list:
        if "BND" in sv_id:
            SV_type = "BND"
        else:
            SV_type = "INV"
        sv_anchors = sam_chr_SV_By_SV[sv_id]
        Ref_chr   = sv_anchors[0][0]
        Ref_start = sv_anchors[0][1]
        Ref_end  = sv_anchors[-1][2]
        Ref_Size = Ref_end - Ref_start + 1
        Qry_chr   = sv_anchors[0][4]
        Qry_start = sv_anchors[0][5]
        Qry_end  = sv_anchors[-1][6]
        if Qry_end > Qry_start:
            Qry_Size = Qry_end - Qry_start + 1
        else:
            Qry_Size = Qry_end - Qry_start - 1
        conserved_flanking = revised_SV_flanking[sv_id][2:]

        ref_anchor_size = 0
        qry_anchor_size = 0
        for cells in sv_anchors:
            ref_anchor_size += cells[3]
            qry_anchor_size += cells[7]

        if len(sv_anchors) == 2 and float(ref_anchor_size) / Ref_Size < 0.20:
            # only 2 anchor
            # anchor size is too small
            for cells in sv_anchors:
                if SV_type == "BND":
                    bnd_order += 1
                    new_id = re.sub("_.*","_",sv_id) + str(bnd_order)
                else:
                    inv_order += 1
                    new_id = re.sub("_.*","_",sv_id) + str(inv_order)

                final_SV_list.append(new_id)
                revised_SV_flanking[new_id] = revised_SV_flanking[sv_id]
                sam_chr_SV_By_SV[new_id] = [cells]
                # revise sam_chr_SV_By_chr
                for anchor in sam_chr_SV_By_chr[Ref_chr]:
                    if anchor[8] == cells[8]:
                        anchor[17] = new_id
                        break
                # Output
                output_list = [SV_type]
                output_list += cells[0:8]
                output_list += conserved_flanking
                output_list.append(new_id)
                output_list += [1,cells[3],cells[7]]
                output = "\t".join(map(str, output_list)) + "\n"
                output_file.write(output)
                # keep my SV table
                My_SV_table.append(output_list)
                for cells in sv_anchors:
                    output = "\t".join(map(str, output_list)) + "\t" + "\t".join(map(str, cells)) + "\n"
                    MY_SV_detail.write(output)

        else:
            final_SV_list.append(sv_id)
            output_list = [SV_type, Ref_chr, Ref_start, Ref_end, Ref_Size, Qry_chr, Qry_start, Qry_end, Qry_Size]
            output_list += conserved_flanking
            output_list.append(sv_id)
            output_list.append(len(sv_anchors))
            output_list.append(ref_anchor_size)
            output_list.append(qry_anchor_size)
            output = "\t".join(map(str, output_list)) + "\n"
            output_file.write(output)
            # keep my SV table
            My_SV_table.append(output_list)
            for cells in sv_anchors:
                output = "\t".join(map(str, output_list)) + "\t" + "\t".join(map(str, cells)) + "\n"
                MY_SV_detail.write(output)
    """
    (2) Output BND based on Step 1.2
           clean_dif_chr_table
           Simply output, no more anchors need
    """
    for dif_cells in clean_dif_chr_table:
        bnd_order += 1
        bnd_order += 1
        sv_id = "BND_" + str(bnd_order)
        output_list = ["BND"] + dif_cells + ["None"] * 4 + [sv_id] + [1,dif_cells[3],dif_cells[7]]
        final_SV_list.append(sv_id)

        output = "\t".join(map(str, output_list)) + "\n"
        output_file.write(output)
        MY_SV_detail.write(output)        # This is already detail
        # keep my SV table
        My_SV_table.append(output_list)
        # add to SV dictionary
        sam_chr_SV_By_SV[sv_id]  = output_list

    """
    (3) Output Exact DUP based on:
            Duplication_table        DUP anchors                                    just before Step 1.1
            raw_table = []           Raw table including all oriented coords        Step 1.0
            raw_table_index = {}     Key & Index for raw_table                      Step 1.0
            sam_chr_SV_By_Key = {}   Same chr SV calling table (Note: No new_id)    Step 2.3
    """
    exact_dup_count = 0
    for dup_cells in Duplication_table:
        key = "_".join(map(str, dup_cells))
        ref_key = "_".join(map(str, dup_cells[0:3]))
        index_in_raw_table = raw_table_index[key]
        ref_chr = dup_cells[0]

        # Get conserved flanking
        prev_index_in_raw_table = index_in_raw_table - 1
        while "_".join(map(str, raw_table[prev_index_in_raw_table])) not in sam_chr_SV_By_Key or \
              len(sam_chr_SV_By_Key["_".join(map(str, raw_table[prev_index_in_raw_table]))]) > 17:
            prev_index_in_raw_table -= 1

        up_ref_index = sam_chr_SV_By_Key["_".join(map(str, raw_table[prev_index_in_raw_table]))][8]
        up_ref_index,left_R_flanking,left_Q_flanking, = get_megaSV_flanking(up_ref_index,ref_chr)

        next_index_in_raw_table = index_in_raw_table + 1
        while "_".join(map(str, raw_table[next_index_in_raw_table])) not in sam_chr_SV_By_Key or \
              len(sam_chr_SV_By_Key["_".join(map(str, raw_table[next_index_in_raw_table]))]) > 17:
            next_index_in_raw_table += 1

        down_ref_index = sam_chr_SV_By_Key["_".join(map(str, raw_table[next_index_in_raw_table]))][8]
        down_ref_index,right_R_flanking,right_Q_flanking = get_megaSV_flanking(down_ref_index,ref_chr)

        # Creat id
        if ref_key not in dup_ref_key_count:
            dup_order += 1
            exact_dup_count += 1
            base_id = "DUP_" + str(dup_order)
            final_SV_list.append(base_id)
            dup_ref_key_count[ref_key] = 1
            dup_ref_key_VS_id[ref_key] = base_id
            copy = 1
            dup_id = base_id + letters[copy]
        else:
            dup_ref_key_count[ref_key] += 1
            base_id = dup_ref_key_VS_id[ref_key]
            copy = dup_ref_key_count[ref_key]
            dup_id = base_id + letters[copy]


        # Get output
        output_list = ["DUP"] + dup_cells
        output_list += [left_R_flanking,left_Q_flanking,right_R_flanking,right_Q_flanking]
        output_list.append(dup_id)
        output_list.append(1)
        output_list.append(dup_cells[3])
        output_list.append(dup_cells[7])

        output = "\t".join(map(str, output_list)) + "\n"
        output_file.write(output)
        MY_SV_detail.write(output)        # This is already detail
        # keep my SV table
        My_SV_table.append(output_list)
        # add to SV dictionary
        dup_infor = dup_cells + [up_ref_index,down_ref_index]
        sam_chr_SV_By_SV[dup_id]  = dup_infor
    """
    (4) Output 90% DUP based on:
            almost_dup_table         DUP anchors                                    Step 1.2.1
            raw_table = []           Raw table including all oriented coords        Step 1.0
            raw_table_index = {}     Key & Index for raw_table                      Step 1.0
            sam_chr_SV_By_Key = {}   Same chr SV calling table (Note: No new_id)    Step 2.3
    """
    cov90_dup_count = 0
    sam_chr_cov90_anchors = 0
    sam_chr_cov90_dup = 0
    dup_id_count = {}
    '''# first line'''
    dup_cells = almost_dup_table[0]
    key = "_".join(map(str, dup_cells))
    index_in_raw_table = raw_table_index[key]
    ref_chr = dup_cells[0]
    # Get conserved flanking
    if dup_cells[0][-4:] == dup_cells[4][-4:]:
        sam_chr_cov90_anchors += 1
        prev_index_in_raw_table = index_in_raw_table - 1
        while "_".join(map(str, raw_table[prev_index_in_raw_table])) not in sam_chr_SV_By_Key or \
              len(sam_chr_SV_By_Key["_".join(map(str, raw_table[prev_index_in_raw_table]))]) > 17:
            prev_index_in_raw_table -= 1

        up_ref_index = sam_chr_SV_By_Key["_".join(map(str, raw_table[prev_index_in_raw_table]))][8]
        up_ref_index,left_R_flanking,left_Q_flanking, = get_megaSV_flanking(up_ref_index,ref_chr)

        next_index_in_raw_table = index_in_raw_table + 1
        while "_".join(map(str, raw_table[next_index_in_raw_table])) not in sam_chr_SV_By_Key or \
              len(sam_chr_SV_By_Key["_".join(map(str, raw_table[next_index_in_raw_table]))]) > 17:
            next_index_in_raw_table += 1

        down_ref_index = sam_chr_SV_By_Key["_".join(map(str, raw_table[next_index_in_raw_table]))][8]
        down_ref_index,right_R_flanking,right_Q_flanking = get_megaSV_flanking(down_ref_index,ref_chr)
    else:
        left_R_flanking = "Not_Check"
        left_Q_flanking = "Not_Check"
        right_R_flanking = "Not_Check"
        right_Q_flanking = "Not_Check"
    # creat id
    dup_order += 1
    cov90_dup_count += 1
    if dup_cells[0][-4:] == dup_cells[4][-4:]:
        base_id = "DUP_" + str(dup_order)
        sam_chr_cov90_dup += 1
    else:
        base_id = "dup_" + str(dup_order)
    final_SV_list.append(base_id)
    dup_id_count[base_id] = 1
    copy = dup_id_count[base_id]
    dup_id = base_id + letters[dup_id_count[base_id]]

    # Get output
    output_list = ["DUP"] + dup_cells
    output_list += [left_R_flanking,left_Q_flanking,right_R_flanking,right_Q_flanking]
    output_list.append(dup_id)
    output_list.append(1)
    output_list.append(dup_cells[3])
    output_list.append(dup_cells[7])

    output = "\t".join(map(str, output_list)) + "\n"
    output_file.write(output)
    MY_SV_detail.write(output)        # This is already detail
    # keep my SV table
    My_SV_table.append(output_list)
    # add to SV dictionary
    dup_infor = dup_cells + [up_ref_index,down_ref_index]
    sam_chr_SV_By_SV[dup_id]  = dup_infor

    '''# 2nd and after'''
    i = 0
    for dup_cells in almost_dup_table[1:]:
        i += 1
        key = "_".join(map(str, dup_cells))
        index_in_raw_table = raw_table_index[key]
        ref_chr = dup_cells[0]
        # Get conserved flanking
        if dup_cells[0][-4:] == dup_cells[4][-4:]:
            sam_chr_cov90_anchors += 1
            prev_index_in_raw_table = index_in_raw_table - 1
            while "_".join(map(str, raw_table[prev_index_in_raw_table])) not in sam_chr_SV_By_Key or \
                  len(sam_chr_SV_By_Key["_".join(map(str, raw_table[prev_index_in_raw_table]))]) > 17:
                prev_index_in_raw_table -= 1

            up_ref_index = sam_chr_SV_By_Key["_".join(map(str, raw_table[prev_index_in_raw_table]))][8]
            up_ref_index,left_R_flanking,left_Q_flanking, = get_megaSV_flanking(up_ref_index,ref_chr)

            next_index_in_raw_table = index_in_raw_table + 1
            while "_".join(map(str, raw_table[next_index_in_raw_table])) not in sam_chr_SV_By_Key or \
                  len(sam_chr_SV_By_Key["_".join(map(str, raw_table[next_index_in_raw_table]))]) > 17:
                next_index_in_raw_table += 1

            down_ref_index = sam_chr_SV_By_Key["_".join(map(str, raw_table[next_index_in_raw_table]))][8]
            down_ref_index,right_R_flanking,right_Q_flanking = get_megaSV_flanking(down_ref_index,ref_chr)
        else:
            left_R_flanking = "Not_Check"
            left_Q_flanking = "Not_Check"
            right_R_flanking = "Not_Check"
            right_Q_flanking = "Not_Check"

        # Creat id
        # check whether overlap with precious anchor
        prev_index = i - 1
        if almost_dup_table[prev_index][0] == dup_cells[0]:
            # chr match
            # check overlap
            this_sta = dup_cells[1]
            this_end = dup_cells[2]
            prev_sta = almost_dup_table[prev_index][1]
            prev_end = almost_dup_table[prev_index][2]
            covered,dup_check = covered_check(covered,this_sta,this_end,prev_sta,prev_end)
        else:
            # chr not match
            # don't check covered
            covered = "No"
        if covered == "No":
            # new dup
            dup_order += 1
            cov90_dup_count += 1
            if dup_cells[0][-4:] == dup_cells[4][-4:]:
                base_id = "DUP_" + str(dup_order)
                sam_chr_cov90_dup += 1
            else:
                base_id = "dup_" + str(dup_order)
            final_SV_list.append(base_id)
            dup_id_count[base_id] = 1
            copy = dup_id_count[base_id]
            dup_id = base_id + letters[dup_id_count[base_id]]
        else:
            # next copy
            dup_id_count[base_id] += 1
            copy = dup_id_count[base_id]
            dup_id = base_id + letters[dup_id_count[base_id]]

        # Get output
        output_list = ["DUP"] + dup_cells
        output_list += [left_R_flanking,left_Q_flanking,right_R_flanking,right_Q_flanking]
        output_list.append(dup_id)
        output_list.append(1)
        output_list.append(dup_cells[3])
        output_list.append(dup_cells[7])

        output = "\t".join(map(str, output_list)) + "\n"
        output_file.write(output)
        MY_SV_detail.write(output)        # This is already detail
        # keep my SV table
        My_SV_table.append(output_list)
        # add to SV dictionary
        dup_infor = dup_cells + [up_ref_index,down_ref_index]
        sam_chr_SV_By_SV[dup_id]  = dup_infor






"""
##############################
           Summary
##############################
"""
bnd_count = 0
inv_count = 0

for sv_id in final_SV_list:
    if "BND" in sv_id:
        bnd_count += 1
    elif 'INV' in sv_id:
        inv_count += 1



print "#" * 30
print "Summary of Assemblytics results processing:"
print "\t1. Assemblytics reported:"
print "\t\tUnique anchor: " + str(len(raw_table)) + " (File name: " + Input_coords + ")"
#print "\t\tIndels: " + str(raw_SV_count) + " (File name: " + Input_SV  + ")"

print "\t2. Clean unique anchors after my filtering:"
print "\t\t" + str(len(sam_chr_SV_By_Ref_index)) + " anchors on same chromosome (not ch00) in both genome. File: Unique_anchor.sam_chr.Synteny_check.megaSV_flanking"
print "\t\t" + str(len(clean_dif_chr_table)) + " anchors on different chromosomes (not ch00) in the 2 genome. File: Unique_anchor.dif_chr"
print "\t\t" + str(len(clean_ch00_anchors)) + " anchors are related to ch00. File name: Unique_anchor.ch00"
print "\t\t" + str(len(Duplication_table)) + " anchors for  " + str(exact_dup_count) + "  exact full-length duplication events (prefix: DUP_). File name: Duplication.Exact.txt"
print "\t\t" + str(len(almost_dup_table)) + " anchors for  " + str(cov90_dup_count) + "  duplications with overlaped region >90%. File name: Duplication.cov90.txt"
print "\t\t\tIncludes:"
print "\t\t\t\t"  + str(sam_chr_cov90_anchors)+ " anchors for  " + str(sam_chr_cov90_dup) + " >90% duplications on same chromosome (prefix: DUP_)"
print "\t\t\t\t"  + str(len(almost_dup_table) - sam_chr_cov90_anchors)+ " anchors for  " + str(cov90_dup_count - sam_chr_cov90_dup) + " >90% duplications on different chromosome (prefix: dup_)"

print "\t3. Output SVs (no indels; Indels will be called later) based on my clean unique anchors:"
print "\t\t" + str(bnd_count) + " BND (Transcloation)"
print "\t\t" + str(inv_count) + " INV (Inversion)"
print "\t\t" + str(exact_dup_count) + " exact full-length DUP (Inversion) on same chromosome"
print "\t\t" + str(sam_chr_cov90_dup) + " >90% length DUP (Inversion) on same chromosome"
print "\t\t" + str(cov90_dup_count - sam_chr_cov90_dup) + " >90% length DUP (Inversion) on different chromosome"


"""Step 2.5 Output final sam_chr_SV_By_chr"""

with open(Prefix + ".Unique_anchor.sam_chr.Synteny_check.megaSV_flanking", "w") as output_file:
    output_list = ['Ref_chr','Ref_start','Ref_end','Ref_size','Qry_chr','Qry_start','Qry_end','Qry_size',"Ref_Order","Qry_Order", \
                    "Synteny","Note","Synteny_Block","Note_block","Synteny_Block_count","Note_block_count","Distance_to_neigbour_Qry","SV_calling","MegaSV",\
                    "left_R_flanking","left_Q_flanking","right_R_flanking","right_Q_flanking"]
    output = "\t".join(map(str, output_list)) + "\n"
    output_file.write(output)
    for chr in sorted_chr_list:
        for cells in sam_chr_SV_By_chr[chr]:
            if len(cells) == 17:
                output = "\t".join(map(str, cells)) + "\n"
                output_file.write(output)
            else:
                flanking = revised_SV_flanking[cells[17]]
                output = "\t".join(map(str, cells + flanking[2:])) + "\n"
                output_file.write(output)










#srf
