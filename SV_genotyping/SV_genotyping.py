#!/usr/bin/env python

import sys
from operator import itemgetter
import re
import subprocess

INPUT = sys.argv[1]
sample = sys.argv[2]
Reference = sys.argv[3]
Query = sys.argv[4]
Ref_bam = sys.argv[5]
Qry_bam = sys.argv[6]
Mismath = sys.argv[7]
Mismath_per = float(Mismath)/100


''' Function '''
def Mismatch_constraint(reads):
    kept_reads = []
    if len(reads) > 1:
        for line in reads:
            fields = line.split("\t")
            read_len = len(fields[9])
            if fields[5] != "*":
                NM = int(re.findall(r"NM:i:\d+",line)[0].split(":")[-1])
                if NM <= Mismath_per * read_len:
                    kept_reads.append(line)
                elif "D" in fields[5] or "I" in fields[5]:
                    cigar_strs = re.findall(r"\d+\w",fields[5])
                    for cigar_str in cigar_strs:
                        if "I" in cigar_str or "D" in cigar_str:
                            NM -= int(cigar_str[:-1])
                    if NM <= Mismath_per * read_len:
                        kept_reads.append(line)
    return kept_reads

def check_INS(chr,sta,end,size,SV_size,SV_type,Reference,Ref_bam):
    SV_GT = "NA"
    if size < 100:
        SV_GT = ins_CIGAR_check(chr,sta,end,size,SV_size,SV_type,Reference,Ref_bam)
    if SV_GT == "NA":
        SV_GT = ins_break_check(chr,sta,end,size,SV_size,SV_type,Reference,Ref_bam)
    if SV_GT == "NA":
        SV_GT = ins_gap_check(chr,sta,end,size,SV_size,SV_type,Reference,Ref_bam)

    return SV_GT

def ins_gap_check(chr,sta,end,size,SV_size,SV_type,Reference,Ref_bam):
    SV_GT = "NA"
    if max(size,SV_size) < 40:
        check_sta = min(sta,end) - 10
        check_end = max(sta,end) + 10
    elif max(size,SV_size) < 400:
        check_sta = min(sta,end) - max(size,SV_size) / 4
        check_end = max(sta,end) + max(size,SV_size) / 4
    else:
        check_sta = min(sta,end) - 100
        check_end = max(sta,end) + 100
    check_size = check_end - check_sta + 1
    region = chr + ":" + str(check_sta) + "-" + str(check_end)
    #covered = float((subprocess.check_output("samtools  view -h  " + sample + ".bam  " +  region + " | python Mismatch_constraint.py  " + Mismath +  " | samtools depth - | awk '$3>0' | wc -l", shell=True)).strip())
    covered = float((subprocess.check_output(get_depth_comand(Ref_bam,region,Mismath,0), shell=True)).strip())
    if 0 < covered < check_size:
        left_flank = chr + ":" + str(check_sta - max(size,SV_size)) + "-" + str(check_sta - 1)
        #left_cov = float((subprocess.check_output("samtools  view -h  " + sample + ".bam  " +  left_flank + " | python Mismatch_constraint.py  " + Mismath +  " | samtools depth - | awk '$3>0' | wc -l", shell=True)).strip())
        left_cov = float((subprocess.check_output(get_depth_comand(Ref_bam,left_flank,Mismath,0), shell=True)).strip())

        right_flank = chr + ":" + str(check_end + 1) + "-" + str(check_end + max(size,SV_size))
        #right_cov = float((subprocess.check_output("samtools  view -h  " + sample + ".bam  " +  right_flank + " | python Mismatch_constraint.py  " + Mismath +  " | samtools depth - | awk '$3>0' | wc -l", shell=True)).strip())
        right_cov = float((subprocess.check_output(get_depth_comand(Ref_bam,right_flank,Mismath,0), shell=True)).strip())
        if left_cov / max(size,SV_size) > 0.5 and right_cov / max(size,SV_size) > 0.5 and min(left_cov / max(size,SV_size), right_cov / max(size,SV_size)) > covered / check_size:
            SV_GT = "INS_Gap"
    return SV_GT

def get_depth_comand(Ref_bam,region,Mismath,MinDep):
    check_sta = region.replace(":","-").split("-")[1]
    check_end = region.replace(":","-").split("-")[2]
    comm_list = ["samtools  view -h   ", Ref_bam, "  ", region, " | python Mismatch_constraint.py  ", Mismath, " | samtools depth - | awk '$3> ", MinDep , "  && $2 >= ", check_sta, " && $2 <= ", check_end, " ' | wc -l  " ]
    return "".join(map(str,comm_list))


def ins_CIGAR_check(chr,sta,end,size,SV_size,SV_type,Reference,Ref_bam):
    SV_GT = "NA"
    if max(size,SV_size) < 40:
        check_sta = min(sta,end) - 10
        check_end = max(sta,end) + 10
    elif max(size,SV_size) < 400:
        check_sta = min(sta,end) - max(size,SV_size) / 4
        check_end = max(sta,end) + max(size,SV_size) / 4
    else:
        check_sta = min(sta,end) - 100
        check_end = max(sta,end) + 100
    region = chr + ":" + str(check_sta) + "-" + str(check_end)
    # Extract reads on this region
    #reads = Mismatch_constraint(subprocess.check_output("samtools view   " + sample + ".bam " +  region , shell=True).strip().split("\n"))
    reads = Mismatch_constraint(subprocess.check_output("samtools view   " + Ref_bam + " " +  region , shell=True).strip().split("\n"))
    exp_I_cigar = set()
    count_I_cigar = {}
    if len(reads) > 1:
        for i in range(size - 3, size + 4) + range(SV_size - 3, SV_size + 4):
            exp_I_cigar.add(str(i) + "I")

        for Read in reads:
            fields = Read.split("\t")
            cigar_strs = re.findall(r"\d+\w",fields[5])
            if len(set(cigar_strs).intersection(exp_I_cigar)) > 0:
                # have exp cigar
                hit_cigars = set(cigar_strs).intersection(exp_I_cigar)
                breakpoint = int(fields[3])
                for cigar_str in cigar_strs:
                    if "M" in cigar_str or "D" in cigar_str:
                        breakpoint += int(cigar_str[:-1])
                    elif cigar_str in hit_cigars:
                        # break
                        if cigar_str not in count_I_cigar:
                            count_I_cigar[cigar_str] = {}
                            count_I_cigar[cigar_str][breakpoint] = 1
                        elif breakpoint not in count_I_cigar[cigar_str]:
                            count_I_cigar[cigar_str][breakpoint] = 1
                        else:
                            count_I_cigar[cigar_str][breakpoint] += 1
    largest_count = 1
    if len(count_I_cigar) > 0:
        for cigar_str in count_I_cigar.keys():
            for breakpoint in count_I_cigar[cigar_str]:
                if count_I_cigar[cigar_str][breakpoint] > largest_count:
                    largest_count = count_I_cigar[cigar_str][breakpoint]
                    best_cigar = cigar_str
                    best_point = breakpoint
                    best_count = largest_count
    if largest_count > 2:
        SV_GT = "_".join(map(str,["Cigar", best_cigar,best_point,best_count]))
    return SV_GT


def ins_break_check(chr,sta,end,size,SV_size,SV_type,Reference,Ref_bam):
    SV_GT = "NA"
    if max(size,SV_size) < 40:
        check_sta = min(sta,end) - 10
        check_end = max(sta,end) + 10
    elif max(size,SV_size) < 400:
        check_sta = min(sta,end) - max(size,SV_size) / 4
        check_end = max(sta,end) + max(size,SV_size) / 4
    else:
        check_sta = min(sta,end) - 100
        check_end = max(sta,end) + 100
    region = chr + ":" + str(check_sta) + "-" + str(check_end)
    reads = Mismatch_constraint(subprocess.check_output("samtools view   " + Ref_bam + " " +  region , shell=True).strip().split("\n"))

    if len(reads) < 3:
        SV_GT == "Too_few_Reads"
    else:
        left_break_dict = {}
        right_break_dict = {}
        for Read in reads:
            fields = Read.split("\t")
            cigar_strs = re.findall(r"\d+\w",fields[5])
            if "H" in cigar_strs[-1] or "S" in cigar_strs[-1]:
                # right soft or hard clipped
                breakpoint = int(fields[3])
                for cigar_str in cigar_strs:
                    if "M" in cigar_str or "D" in cigar_str:
                        breakpoint += int(cigar_str[:-1])
                count_sth(breakpoint,left_break_dict)
            if "H" in cigar_strs[0] or "S" in cigar_strs[0]:
                # right soft or hard clipped
                breakpoint = int(fields[3]) - 1
                count_sth(breakpoint,right_break_dict)
        best_pair = []
        smallest_dis = 100000000000
        largest_evi_reads = 0
        if len(left_break_dict) > 0 and len(right_break_dict) > 0:
            # check the quality of breakpoints
            Good_left = []
            Good_right = []
            for left_break in left_break_dict.keys():
                if left_break_dict[left_break] > 1 and check_sta <= left_break <= check_end and left_break_dict[left_break] * 2 > Get_normal_reads_cover_point(chr,left_break,size,SV_size,SV_type,Reference,Ref_bam):
                    Good_left.append(left_break)
            for right_break in right_break_dict.keys():
                if right_break_dict[right_break] > 1  and check_sta <= right_break <= check_end and right_break_dict[right_break] * 2 > Get_normal_reads_cover_point(chr,right_break,size,SV_size,SV_type,Reference,Ref_bam):
                    Good_right.append(right_break)
            if Good_left != [] and Good_right != []:
                # have good breakpoints in both
                # get good pair
                for left_break in Good_left:
                    for right_break in Good_right:
                        evi_reads = left_break_dict[left_break] + right_break_dict[right_break]
                        if sta <= left_break <= end  or  sta <= right_break <= end:
                            obs_dis = 0
                        else:
                            obs_dis = min(abs(left_break - sta),abs(left_break - end),abs(right_break - sta),abs(right_break - end))
                        if obs_dis < smallest_dis:
                            smallest_dis = obs_dis
                            best_pair = [obs_dis,left_break,right_break]
                            largest_evi_reads = evi_reads
                            SV_GT = "BreakDistance_" + "_".join(map(str,best_pair))
                        elif obs_dis == smallest_dis and largest_evi_reads < evi_reads:
                            best_pair = [obs_dis,left_break,right_break]
                            largest_evi_reads = evi_reads
                            SV_GT = "BreakDistance_" + "_".join(map(str,best_pair))
            elif Good_left != []:
                # only have good left
                for left_break in Good_left:
                    evi_reads = left_break_dict[left_break]
                    if evi_reads > 2:
                        if sta <= left_break <= end:
                            obs_dis = 0
                        else:
                            obs_dis = min(abs(left_break - sta),abs(left_break - end))
                        if obs_dis < smallest_dis:
                            smallest_dis = obs_dis
                            largest_evi_reads = evi_reads
                            SV_GT = "No_GoodRight_point_" + str(left_break)
                        elif obs_dis == smallest_dis and evi_reads > largest_evi_reads:
                            largest_evi_reads = evi_reads
                            SV_GT = "No_GoodRight_point_" + str(left_break)

            elif Good_right != []:
                # only have good right
                for right_break in Good_right:
                    evi_reads = right_break_dict[right_break]
                    if evi_reads > 2:
                        if sta <= right_break <= end:
                            obs_dis = 0
                        else:
                            obs_dis = min(abs(right_break - sta),abs(right_break - end))
                        if obs_dis < smallest_dis:
                            smallest_dis = obs_dis
                            largest_evi_reads = evi_reads
                            SV_GT = "No_Left_point_" + str(right_break)
                        elif obs_dis == smallest_dis and evi_reads > largest_evi_reads:
                            largest_evi_reads = evi_reads
                            SV_GT = "No_Left_point_" + str(right_break)
        elif len(left_break_dict) > 0:
            for left_break in left_break_dict.keys():
                if left_break_dict[left_break] > 2 and check_sta <= left_break <= check_end and left_break_dict[left_break] * 2 > Get_normal_reads_cover_point(chr,left_break,size,SV_size,SV_type,Reference,Ref_bam):
                    evi_reads = left_break_dict[left_break]
                    if sta <= left_break <= end:
                        obs_dis = 0
                    else:
                        obs_dis = min(abs(left_break - sta),abs(left_break - end))
                    if obs_dis < smallest_dis:
                        smallest_dis = obs_dis
                        largest_evi_reads = evi_reads
                        SV_GT = "No_Right_point_" + str(left_break)
                    elif obs_dis == smallest_dis and evi_reads > largest_evi_reads:
                        largest_evi_reads = evi_reads
                        SV_GT = "No_Right_point_" + str(left_break)
        elif len(right_break_dict) > 0:
            for right_break in right_break_dict.keys():
                if right_break_dict[right_break] > 2  and check_sta <= right_break <= check_end and right_break_dict[right_break] * 2 > Get_normal_reads_cover_point(chr,right_break,size,SV_size,SV_type,Reference,Ref_bam):
                    evi_reads = right_break_dict[right_break]
                    if sta <= right_break <= end:
                        obs_dis = 0
                    else:
                        obs_dis = min(abs(right_break - sta),abs(right_break - end))
                    if obs_dis < smallest_dis:
                        smallest_dis = obs_dis
                        largest_evi_reads = evi_reads
                        SV_GT = "No_Left_point_" + str(right_break)
                    elif obs_dis == smallest_dis and evi_reads > largest_evi_reads:
                        largest_evi_reads = evi_reads
                        SV_GT = "No_Left_point_" + str(right_break)
    return SV_GT


def check_DEL(chr,sta,end,size,SV_size,SV_type,Reference,Ref_bam):
    SV_GT = "NA"
    if size < 100:
        SV_GT = del_CIGAR_check(chr,sta,end,size,SV_size,SV_type,Reference,Ref_bam)
    if SV_GT == "NA":
        # depth check fail
        depth_GT = del_depth_check(chr,sta,end,size,SV_size,SV_type,Reference,Ref_bam)
        break_GT = del_break_check(chr,sta,end,size,SV_size,SV_type,Reference,Ref_bam)
        if depth_GT == "NA" and break_GT == "NA":
            SV_GT = "NA"
        elif depth_GT == "NA":
            SV_GT = break_GT
        elif break_GT == "NA":
            SV_GT = depth_GT
        else:
            SV_GT = break_GT + "|" + depth_GT
    return SV_GT


def del_depth_check(chr,sta,end,size,SV_size,SV_type,Reference,Ref_bam):
    SV_GT = "NA"
    region = chr + ":" + str(sta) + "-" + str(end)
    if sta > size:
        left_flank = chr + ":" + str(sta - size) + "-" + str(sta - 1)
        left_size = size
    else:
        left_flank = chr + ":" + str(1) + "-" + str(sta - 1)
        left_size = sta - 1
    if end + 2 < chr_size_dict[Reference][chr]:
        right_flank = chr + ":" + str(end + 1) + "-" + str(end + size)
        right_size = size
    else:
        right_flank = chr + ":" + str(end + 1) + "-" + str(chr_size_dict[Reference][chr])
        right_size = chr_size_dict[Reference][chr] - end

    #covered = float((subprocess.check_output("samtools  view -h  " + sample + ".bam  " +  region + " | python Mismatch_constraint.py  " + Mismath +  " | samtools depth - | awk '$3>1' | wc -l", shell=True)).strip())
    covered = float((subprocess.check_output(get_depth_comand(Ref_bam,region,Mismath,1), shell=True)).strip())
    if left_size > 0:
        #left_cov = float((subprocess.check_output("samtools  view -h  " + sample + ".bam  " +  left_flank + " | python Mismatch_constraint.py  " + Mismath +  " | samtools depth - | awk '$3>1' | wc -l", shell=True)).strip())
        left_cov = float((subprocess.check_output(get_depth_comand(Ref_bam,left_flank,Mismath,1), shell=True)).strip())
    if right_size >0:
        #right_cov = float((subprocess.check_output("samtools  view -h  " + sample + ".bam  " +  right_flank + " | python Mismatch_constraint.py  " + Mismath +  " | samtools depth - | awk '$3>1' | wc -l", shell=True)).strip())
        right_cov = float((subprocess.check_output(get_depth_comand(Ref_bam,right_flank,Mismath,1), shell=True)).strip())
    if covered / size < 0.5 or covered / SV_size < 0.5:
        if left_size > 0 and right_size > 0:
            if left_cov/left_size > 0.5 and right_cov/right_size > 0.5:
                # SV region covered < 0.5
                # at least one fland regioin covered > 0.5
                SV_GT = "Cov2_" + str(int(100 * covered / size))
            elif left_cov/left_size > 0.5 or right_cov/right_size > 0.5:
                SV_GT = "One_Cov2_" + str(int(100 * covered / size))
            elif left_cov/left_size < 0.5 and right_cov/right_size < 0.5:
                SV_GT = "Too_Large"
        elif left_size > 0:
            if left_cov/left_size > 0.5:
                SV_GT = "Cov2_" + str(int(100 * covered / size))
            else:
                SV_GT = "Too_Large"
        elif right_size > 0:
            if right_cov/right_size > 0.5:
                SV_GT = "Cov2_" + str(int(100 * covered / size))
            else:
                SV_GT = "Too_Large"
    return SV_GT


def del_CIGAR_check(chr,sta,end,size,SV_size,SV_type,Reference,Ref_bam):
    SV_GT = "NA"
    region = chr + ":" + str(sta) + "-" + str(end)
    # Extract reads on this region
    reads = subprocess.check_output("samtools view   " + Ref_bam + " " +  region , shell=True).strip().split("\n")
    reads = Mismatch_constraint(reads)
    exp_D_cigar = set()
    count_D_cigar = {}
    if len(reads) > 1:
        for i in range(size - 3, size + 4) + range(SV_size - 3, SV_size + 4):
            exp_D_cigar.add(str(i) + "D")

        for Read in reads:
            fields = Read.split("\t")
            cigar_strs = re.findall(r"\d+\w",fields[5])
            if len(set(cigar_strs).intersection(exp_D_cigar)) > 0:
                # have exp cigar
                hit_cigars = set(cigar_strs).intersection(exp_D_cigar)
                breakpoint = int(fields[3])
                for cigar_str in cigar_strs:
                    if "M" in cigar_str:
                        breakpoint += int(cigar_str[:-1])
                    elif cigar_str in hit_cigars:
                        # break
                        if cigar_str not in count_D_cigar:
                            count_D_cigar[cigar_str] = {}
                            count_D_cigar[cigar_str][breakpoint] = 1
                        elif breakpoint not in count_D_cigar[cigar_str]:
                            count_D_cigar[cigar_str][breakpoint] = 1
                        else:
                            count_D_cigar[cigar_str][breakpoint] += 1
                        breakpoint += int(cigar_str[:-1])
    largest_count = 1
    if len(count_D_cigar) > 0:
        for cigar_str in count_D_cigar.keys():
            for breakpoint in count_D_cigar[cigar_str]:
                if count_D_cigar[cigar_str][breakpoint] > largest_count:
                    largest_count = count_D_cigar[cigar_str][breakpoint]
                    best_cigar = cigar_str
                    best_point = breakpoint
                    best_count = largest_count
    if largest_count > 2:
        SV_GT = "_".join(map(str,["Cigar", best_cigar,best_point,best_count]))
    return SV_GT

def del_break_check(chr,sta,end,size,SV_size,SV_type,Reference,Ref_bam):
    SV_GT = "NA"
    # left break
    # right break
    # if have both, pass
    len_dev = size / 2
    if len_dev > 100:
        len_dev = 100

    if sta < 100:
        # only one flank
        # check right flank
        right_break_dict = {}
        r_sta = end - len_dev
        r_end = end + len_dev
        right_flank = chr + ":" + str(end - len_dev) + "-" + str(end + len_dev)

        reads = Mismatch_constraint(subprocess.check_output("samtools view   " + Ref_bam + " " +  right_flank , shell=True).strip().split("\n"))

        if len(reads) > 1:
            for Read in reads:
                fields = Read.split("\t")
                cigar_strs = re.findall(r"\d+\w",fields[5])
                if "H" in cigar_strs[0] or "S" in cigar_strs[0]:
                    # right soft or hard clipped
                    breakpoint = int(fields[3]) - 1
                    count_sth(breakpoint,right_break_dict)
        Most_break = 0
        if len(right_break_dict) > 0:
            for breakpoint in right_break_dict.keys():
                if r_sta < breakpoint < r_end and right_break_dict[breakpoint] > 3 and right_break_dict[breakpoint] > Most_break:
                    Most_break = right_break_dict[breakpoint]
                    SV_GT = "Right_Break_check_" + str(breakpoint) + "_" + str(right_break_dict[breakpoint])
    elif end + 100 > chr_size_dict[Reference][chr]:
        left_break_dict = {}
        l_sta = sta - len_dev
        l_end = sta + len_dev
        left_flank = chr + ":" + str(sta - len_dev) + "-" + str(sta + len_dev)

        reads = Mismatch_constraint(subprocess.check_output("samtools view   " + Ref_bam + " " +  left_flank , shell=True).strip().split("\n"))

        if len(reads) > 1:
            for Read in reads:
                fields = Read.split("\t")
                cigar_strs = re.findall(r"\d+\w",fields[5])
                if "H" in cigar_strs[-1] or "S" in cigar_strs[-1]:
                    # right soft or hard clipped
                    breakpoint = int(fields[3])
                    for cigar_str in cigar_strs:
                        if "M" in cigar_str or "D" in cigar_str:
                            breakpoint += int(cigar_str[:-1])
                    count_sth(breakpoint,left_break_dict)
        Most_break = 0
        if len(left_break_dict) > 0:
            for breakpoint in left_break_dict.keys():
                if l_sta < breakpoint < l_end and left_break_dict[breakpoint] > 3 and left_break_dict[breakpoint] > Most_break:
                    Most_break = left_break_dict[breakpoint]
                    SV_GT = "Left_Break_check_" + str(breakpoint) + "_" + str(left_break_dict[breakpoint])
    else:
        # check left flank
        # right soft
        left_break_dict = {}
        if sta > len_dev and sta + len_dev < chr_size_dict[Reference][chr]:
            left_flank = chr + ":" + str(sta - len_dev) + "-" + str(sta + len_dev)
        elif sta + len_dev < chr_size_dict[Reference][chr]:
            left_flank = chr + ":" + str(1) + "-" + str(sta + len_dev)
        elif sta > len_dev:
            left_flank = chr + ":" + str(sta - len_dev) + "-" + str(chr_size_dict[Reference][chr])

        reads = Mismatch_constraint(subprocess.check_output("samtools view   " + Ref_bam + " " +  left_flank , shell=True).strip().split("\n"))

        if len(reads) > 1:
            for Read in reads:
                fields = Read.split("\t")
                cigar_strs = re.findall(r"\d+\w",fields[5])
                if "H" in cigar_strs[-1] or "S" in cigar_strs[-1]:
                    # right soft or hard clipped
                    breakpoint = int(fields[3])
                    for cigar_str in cigar_strs:
                        if "M" in cigar_str or "D" in cigar_str:
                            breakpoint += int(cigar_str[:-1])
                    count_sth(breakpoint,left_break_dict)
        # check right flank
        right_break_dict = {}
        if end > len_dev and end + len_dev < chr_size_dict[Reference][chr]:
            right_flank = chr + ":" + str(end - len_dev) + "-" + str(end + len_dev)
        elif end + len_dev < chr_size_dict[Reference][chr]:
            right_flank = chr + ":" + str(1) + "-" + str(end + len_dev)
        elif end > len_dev:
            right_flank = chr + ":" + str(end - len_dev) + "-" + str(chr_size_dict[Reference][chr])

        reads = Mismatch_constraint(subprocess.check_output("samtools view   " + Ref_bam + " " +  right_flank , shell=True).strip().split("\n"))

        if len(reads) > 1:
            for Read in reads:
                fields = Read.split("\t")
                cigar_strs = re.findall(r"\d+\w",fields[5])
                if "H" in cigar_strs[0] or "S" in cigar_strs[0]:
                    # right soft or hard clipped
                    breakpoint = int(fields[3]) - 1
                    count_sth(breakpoint,right_break_dict)
        # check right flank
        # Get best pair
        best_pair = []
        smallest_dev = 100000000000
        largest_evi_reads = 0

        Good_left = []
        Good_right = []
        for left_break in left_break_dict.keys():
            if left_break_dict[left_break] > 1 and left_break_dict[left_break] * 2 > Get_normal_reads_cover_point(chr,left_break,size,SV_size,SV_type,Reference,Ref_bam):
                Good_left.append(left_break)
        for right_break in right_break_dict.keys():
            if right_break_dict[right_break] > 1 and right_break_dict[right_break] * 2 > Get_normal_reads_cover_point(chr,right_break,size,SV_size,SV_type,Reference,Ref_bam):
                Good_right.append(right_break)

        if Good_left != [] and Good_right != []:
            # have good breakpoints in both
            # check get good pair
            for left_break in Good_left:
                for right_break in Good_right:
                    evi_reads = left_break_dict[left_break] + right_break_dict[right_break]
                    if right_break - left_break > size / 2 or right_break - left_break > SV_size / 2:
                        obs_size = right_break - left_break + 1
                        obs_dev = min(abs(obs_size - size), abs(obs_size - SV_size))
                        if obs_dev < smallest_dev:
                            smallest_dev = obs_dev
                            best_pair = [obs_size, left_break,right_break]
                            largest_evi_reads = evi_reads
                            SV_GT = "BreakSize_" + "_".join(map(str,best_pair))
                        elif obs_dev == smallest_dev and largest_evi_reads < evi_reads:
                            best_pair = [obs_size, left_break,right_break]
                            SV_GT = "BreakSize_" + "_".join(map(str,best_pair))
        if SV_GT == "NA":
            # No good pair found
            # check good break on one side
            # left
            smallest_left_dev = 100000000000
            largest_left_evi = 0
            if Good_left != []:
                for left_break in Good_left:
                    left_dev = abs(left_break - sta)
                    left_evi = left_break_dict[left_break]
                    if left_dev < smallest_left_dev and (left_dev < size / 2.0 or left_dev < SV_size / 2.0):
                        best_left = left_break
                        smallest_left_dev = left_dev
                        largest_left_evi = left_evi
                    elif left_dev == smallest_left_dev and left_evi > largest_left_evi and (left_dev < size / 2.0 or left_dev < SV_size / 2.0):
                        best_left = left_break
                        smallest_left_dev = left_dev
                        largest_left_evi = left_evi
            # right
            smallest_right_dev = 100000000000
            largest_right_evi = 0
            if Good_right != []:
                for right_break in Good_right:
                    right_dev = abs(right_break - end)
                    right_evi = right_break_dict[right_break]
                    if right_dev < smallest_right_dev and (right_dev < size / 2.0 or right_dev < SV_size / 2.0):
                        best_right = right_break
                        smallest_right_dev = right_dev
                        largest_right_evi = right_evi
                    elif right_dev == smallest_right_dev and right_evi > largest_right_evi and (right_dev < size / 2.0 or right_dev < SV_size / 2.0):
                        best_right = right_break
                        smallest_right_dev = right_dev
                        largest_right_evi = right_evi
            # get best as output
            if smallest_left_dev < smallest_right_dev or (smallest_left_dev == smallest_right_dev and largest_left_evi > largest_right_evi):
                # left is better
                # if right is gap, this should be ok
                if max(size,SV_size) < 50:
                    right_region = chr + ":" + str(end - max(size,SV_size)/2) + "-" + str(end + max(size,SV_size)/2)
                else:
                    right_region = chr + ":" + str(end - 25) + "-" + str(end + 25)
                right_len = len(re.sub('.*N', '', "".join(subprocess.check_output("samtools faidx " + Reference + " " +  right_region, shell=True).split("\n")[1:-1])))
                #right_cov = int((subprocess.check_output("samtools  view -h  " + sample + ".bam  " +  right_region + " | python Mismatch_constraint.py  " + Mismath +  " | samtools depth - | awk '$3>0' | wc -l", shell=True)).strip())
                right_cov = float((subprocess.check_output(get_depth_comand(Ref_bam,right_region,Mismath,0), shell=True)).strip())
                if right_cov < right_len:
                    SV_GT = "Left_Break_" + str(best_left) + "_" + str(left_break_dict[best_left])
            elif smallest_left_dev > smallest_right_dev or (smallest_left_dev == smallest_right_dev and largest_left_evi < largest_right_evi):
                # right is better
                if max(size,SV_size) < 50:
                    left_region = chr + ":" + str(sta - max(size,SV_size)) + "-" + str(sta + max(size,SV_size))
                else:
                    left_region = chr + ":" + str(sta - 25) + "-" + str(sta + 25)
                left_len = len(re.sub('.*N', '', "".join(subprocess.check_output("samtools faidx " + Reference + " " +  left_region, shell=True).split("\n")[1:-1])))
                #left_cov = int((subprocess.check_output("samtools  view -h  " + sample + ".bam  " +  left_region + " | python Mismatch_constraint.py  " + Mismath +  " | samtools depth - | awk '$3>0' | wc -l", shell=True)).strip())
                left_cov = float((subprocess.check_output(get_depth_comand(Ref_bam,left_region,Mismath,0), shell=True)).strip())
                if left_cov < left_len:
                    SV_GT = "Right_Break_" + str(best_right) + "_" + str(right_break_dict[best_right])


    return SV_GT


def Get_normal_reads_cover_point(chr,check_break,size,SV_size,SV_type,Reference,Ref_bam):
    normal_reads = 0
    # not contain cigar
    exp_cigar = set()
    for i in range(size - 3, size + 4) + range(SV_size - 3, SV_size + 4):
        if SV_type == "INS":
            exp_cigar.add(str(i) + "I")
        else:
            exp_cigar.add(str(i) + "D")
    # check region
    region = chr + ":" + str(check_break) + "-" + str(check_break+1)
    reads = Mismatch_constraint(subprocess.check_output("samtools view   " + Ref_bam + " " +  region , shell=True).strip().split("\n"))

    if len(reads) > 1:
        for Read in reads:
            fields = Read.split("\t")
            cigar_strs = re.findall(r"\d+\w",fields[5])
            breakpoint = int(fields[3])
            if len(set(cigar_strs).intersection(exp_cigar)) == 0 and not ("H" in cigar_strs[-1] or "S" in cigar_strs[-1] or "H" in cigar_strs[0] or "S" in cigar_strs[0]):
                # No exp cigar
                # Not clipped reads
                # Get "M" region
                for cigar_str in cigar_strs:
                    if "M" in cigar_str:
                        M_sta = breakpoint
                        M_end = breakpoint + int(cigar_str[:-1])
                        breakpoint = M_end
                        if M_sta < check_break - 3 and M_end > check_break + 3:
                            normal_reads += 1
                    elif "D" in cigar_str:
                        breakpoint += int(cigar_str[:-1])
    return normal_reads

def count_sth(breakpoint,break_dict):
    if breakpoint not in break_dict:
        break_dict[breakpoint] = 1
    else:
        break_dict[breakpoint] += 1

def convert_int(fields):
    for i in range(0,len(fields)):
        if fields[i].isdigit() or (fields[i].startswith("-") and fields[i][1:].isdigit()):
            fields[i] = int(fields[i])
    return fields

def check_Tandem_expansion(chr,sta,end,size,SV_type,Reference,Ref_bam):
    SV_GT = "NA"
    region = chr + ":" + str(sta) + "-" + str(end)
    depth_table = (subprocess.check_output(get_depth_data(Ref_bam,region,Mismath,1), shell=True)).strip().split("\n")
    if len(depth_table) >= 3 and len(depth_table) / float(size) > 0.8:
        mean_dep = get_mean_dep(depth_table)
        sv_size = end - sta + 1
        # check flanking size
        check_size = sv_size
        '''
        if sv_size <= 200:
            # check ecqual lenth
            check_size = sv_size
        else:
            check_size = 200'''
        # get flanking region
        if sta > check_size:
            left_region = chr + ":" + str(sta-check_size) + "-" + str(sta-1)
            left_size = check_size
        else:
            left_region = chr + ":" + str(1) + "-" + str(sta-1)
            left_size = sta-1
        if end + check_size < chr_size_dict[Reference][chr]:
            right_region = chr + ":" + str(end + 1) + "-" + str(end + check_size)
            right_size = check_size
        else:
            right_region = chr + ":" + str(end + 1) + "-" + str(chr_size_dict[Reference][chr])
            right_size = chr_size_dict[Reference][chr] - end
        # get depth table
        #print region + "\t" + left_region + "\t" + right_region
        left_table = (subprocess.check_output(get_depth_data(Ref_bam,left_region,Mismath,1), shell=True)).strip().split("\n")
        right_table = (subprocess.check_output(get_depth_data(Ref_bam,right_region,Mismath,1), shell=True)).strip().split("\n")
        if len(left_table) >= 3 and len(left_table) / float(left_size) > 0.5 and len(right_table) >= 3 and len(right_table) / float(right_size) > 0.5:
            left_dep = get_mean_dep(left_table)
            right_dep = get_mean_dep(right_table)
            if (mean_dep > left_dep * 1.5 or mean_dep > right_dep * 1.5) and mean_dep > left_dep * 1.2 and mean_dep > right_dep * 1.2:
                SV_GT = "Fold_1.5"
        elif len(left_table) >= 3 and len(left_table) / float(left_size) > 0.5:
            left_dep = get_mean_dep(left_table)
            if mean_dep > left_dep * 1.5:
                SV_GT = "Left_fold_1.5"
        elif len(right_table) >= 3 and len(right_table) / float(right_size) > 0.5:
            right_dep = get_mean_dep(right_table)
            if mean_dep > right_dep * 1.5:
                SV_GT = "Right_fold_1.5"
    return SV_GT

def get_depth_data(Ref_bam,region,Mismath,MinDep):
    check_sta = region.replace(":","-").split("-")[1]
    check_end = region.replace(":","-").split("-")[2]
    comm_list = ["samtools  view -h   ", Ref_bam, " ", region, " | python Mismatch_constraint.py  ", Mismath, " | samtools depth -d 200 - | awk '$3> ", MinDep , "  && $2 >= ", check_sta, " && $2 <= ", check_end, " ' " ]
    return "".join(map(str,comm_list))

def get_mean_dep(depth_table):
    S = 0.0
    for line in depth_table:
        cells = convert_int(line.split("\t"))
        S += cells[2]
    return S / len(depth_table)

def gap_count(Chr,Sta,End,Chr_size,Genome):
    SV_sta = min(Sta,End)
    SV_end = max(Sta,End)
    SV_size = SV_end - SV_sta
    if SV_size < 10:
        if SV_sta > 10 and SV_end + 10 < Chr_size[Chr]:
            Left_gap = gap_in_region(Chr,SV_sta-10,SV_end,Genome)
            Right_gap = gap_in_region(Chr,SV_sta,SV_end+10,Genome)
        elif SV_sta > 10:
            # So, SV_end + 10 > Chr_size[Chr]
            Left_gap = gap_in_region(Chr,SV_sta-10,SV_end,Genome)
            Right_gap = 0
        elif SV_end + 10 < Chr_size[Chr]:
            # So, SV_sta < 10
            Left_gap = 0
            Right_gap = gap_in_region(Chr,SV_sta,SV_end+10,Genome)
    else:
        if SV_size < 100:
            Diff = 10
        elif SV_size < 500:
            Diff = int(SV_size / 10.0)
        else:
            Diff = 50

        if SV_sta > Diff and SV_end + Diff < Chr_size[Chr]:
            Left_gap = gap_in_region(Chr,SV_sta-Diff,SV_sta+Diff,Genome)
            Right_gap = gap_in_region(Chr,SV_end-Diff,SV_end+Diff,Genome)
        elif SV_sta > Diff:
            # So, SV_end + 10 > Chr_size[Chr]
            Left_gap = gap_in_region(Chr,SV_sta-Diff,SV_sta+Diff,Genome)
            Right_gap = 0
        elif SV_end + Diff < Chr_size[Chr]:
            # So, SV_sta < 10
            Left_gap = 0
            Right_gap = gap_in_region(Chr,SV_end-Diff,SV_end+Diff,Genome)
    return Left_gap,Right_gap

def gap_in_region(Chr,Sta,End,Genome):
    region = Chr + ":" + str(Sta) + "-" + str(End)
    Gap = "".join(subprocess.check_output("samtools faidx " + Genome + " " +  region, shell=True).split("\n")[1:-1]).upper().count("N")
    return Gap

def Remove_bad_Evi(Ref_Evi,Ref_Left,Ref_Right):
    Evi = Ref_Evi
    if Ref_Evi.startswith("One_Cov2"):
        Evi = "NA"

    if (Ref_Left > 0 and Ref_Evi.startswith("Left_Break_")) or (Ref_Right > 0 and Ref_Evi.startswith("Right_Break_")):
        if "|" not in Ref_Evi:
            Evi = "NA"
        else:
            Evi = Ref_Evi.split("|")[1]
    return Evi

''' Chr size '''
chr_size_dict = {}
with open(Reference + ".fai") as infile:
    chr_size_dict[Reference] = {}
    for line in infile:
        cells = line.strip().split("\t")
        chr_size_dict[Reference][cells[0]] = int(cells[1])

with open(Query + ".fai") as infile:
    chr_size_dict[Query] = {}
    for line in infile:
        cells = line.strip().split("\t")
        chr_size_dict[Query][cells[0]] = int(cells[1])


''' Main program '''
num = -1
with open(INPUT) as infile, open(sample + ".GT.txt", "w") as GT_out:
    for line in infile:
        num += 1
        line = line.strip()
        cells = convert_int(line.strip().split("\t"))
        if num == 0:
            #GT_out.write("#Genotype: L for Heinz(SL); P for LA2093(SP); H for Heterozygous; U for Undetermined.\n")
            GT_out.write("#Genotype: R for homozygous Reference genotype; Q for homozygous Query genotype; H for Heterozygous genotype; U for Undetermined.\n")
            GT_out.write("SV_ID\t" + sample + "\n")
            #Evi_out.write("SV_ID\tReference\tQuery\n")
            #gap_out.write(line + "\tRef_Right\tQry_Left\tQry_Right\n")
            j = 0
            for j in range(0,len(cells)):
                if cells[j] == "SV_type":
                    SV_type_col = j
                elif cells[j] == "SV_ID":
                    SV_ID_col = j
                elif cells[j] == "Ref_chr":
                    Ref_chr_col = j
                elif cells[j] == "Ref_start":
                    Ref_start_col = j
                elif cells[j] == "Ref_end":
                    Ref_end_col = j
                elif cells[j] == "Ref_Size":
                    Ref_Size_col = j
                elif cells[j] == "Qry_chr":
                    Qry_chr_col = j
                elif cells[j] == "Qry_start":
                    Qry_start_col = j
                elif cells[j] == "Qry_end":
                    Qry_end_col = j
                elif cells[j] == "Qry_Size":
                    Qry_Size_col = j
        else:
            SV_ID = cells[SV_ID_col]
            SV_type = cells[SV_type_col]
            '''check evidence from alignments on reference genome'''
            chr = cells[Ref_chr_col]
            sta = min(cells[Ref_start_col],cells[Ref_end_col])
            end = max(cells[Ref_start_col],cells[Ref_end_col])
            SV_size = abs(cells[Qry_Size_col] - cells[Ref_Size_col])

            if "INS" in SV_type.upper() and cells[Ref_Size_col] >= -10 and SV_size >= 10:
                New_type = "1_INS"    # for ref
                size = cells[Qry_Size_col]
                Ref_Evi = check_INS(chr,sta,end,size,SV_size,New_type,Reference,Ref_bam)

            elif ("INS" in SV_type.upper() and cells[Ref_Size_col] < -10) or SV_type == "Tandem_expansion":
                """if INS is better, use INS """
                # should have reads cover on ref
                # should have no reads cover on qry
                # breakpoints are not important
                New_type = "2_Tandem_expansion"
                size = abs(cells[Ref_Size_col])
                Ref_Evi = check_Tandem_expansion(chr,sta,end,size,New_type,Reference,Ref_bam)
                if Ref_Evi == "NA" and size < 100 and size * 2 < cells[Qry_Size_col]:
                    size = cells[Qry_Size_col]
                    Ref_Evi = check_INS(chr,sta,end,size,SV_size,New_type,Reference,Ref_bam)
                    New_type = "1_INS\t2_Tandem_expansion"
                    if Ref_Evi.startswith("BreakDistance"):
                        BreakDistance = int(Ref_Evi.split("_")[1])
                        if BreakDistance > (abs(cells[Ref_Size_col]) / 3.0) and BreakDistance > 10:
                            Ref_Evi = "NA"
            else:
                """if (1) failed; (2) ref < qry size; (3) ref size < 100bp
                      Try check_INS
                """
                # For all other SVs
                # should have no reads or ref, just like DEL
                size = abs(cells[Ref_Size_col])
                #SV_type = "DEL"
                if "DEL" in SV_type.upper() and SV_size >= 10:
                    # perfect del
                    New_type = "3_DEL"
                    Ref_Evi = check_DEL(chr,sta,end,size,SV_size,New_type,Reference,Ref_bam)
                else:
                    # like del, only consider it ref size
                    New_type = "4_SUB"
                    Ref_Evi = check_DEL(chr,sta,end,size,size,New_type,Reference,Ref_bam)

                    if Ref_Evi == "NA" and size <100 and cells[Ref_Size_col] * 2 < cells[Qry_Size_col]:
                        size = cells[Qry_Size_col]
                        Ref_Evi = check_INS(chr,sta,end,size,SV_size,New_type,Reference,Ref_bam)
                        New_type = "1_INS\t4_SUB"
                        if Ref_Evi.startswith("BreakDistance"):
                            BreakDistance = int(Ref_Evi.split("_")[1])
                            if BreakDistance > (abs(cells[Ref_Size_col]) / 3.0) and BreakDistance > 10:
                                Ref_Evi = "NA"


            '''check evidence from alignments on Query genome'''
            chr = cells[Qry_chr_col]
            sta = min(cells[Qry_start_col],cells[Qry_end_col])
            end = max(cells[Qry_start_col],cells[Qry_end_col])
            SV_size = abs(cells[Qry_Size_col] - cells[Ref_Size_col])

            if "DEL" in SV_type.upper() and cells[Qry_Size_col] >= -10 and SV_size >= 10:
                New_type = "1_INS"    # for Qry
                size = cells[Ref_Size_col]
                Qry_Evi = check_INS(chr,sta,end,size,SV_size,New_type,Query,Qry_bam)

            elif ("DEL" in SV_type.upper() and cells[Qry_Size_col] < -10) or SV_type == "Tandem_contraction":
                """if INS is better, use INS """
                # should have reads cover on qry
                # should have no reads cover on ref
                # breakpoints are not important
                New_type = "2_Tandem_expansion"
                size = abs(cells[Qry_Size_col])
                Qry_Evi = check_Tandem_expansion(chr,sta,end,size,New_type,Query,Qry_bam)
                if Qry_Evi == "NA" and size < 100 and size * 2 < cells[Ref_Size_col]:
                    size = cells[Ref_Size_col]
                    Qry_Evi = check_INS(chr,sta,end,size,SV_size,New_type,Query,Qry_bam)
                    New_type = "1_INS\t2_Tandem_expansion"
                    if Qry_Evi.startswith("BreakDistance"):
                        BreakDistance = int(Qry_Evi.split("_")[1])
                        if BreakDistance > (abs(cells[Qry_Size_col]) / 3.0) and BreakDistance > 10:
                            Qry_Evi = "NA"

            else:
                """if (1) failed; (2) ref < qry size; (3) ref size < 100bp
                      Try check_INS
                """
                # For all other SVs
                # should have no reads or qry, just like DEL
                size = abs(cells[Qry_Size_col])
                #SV_type = "DEL"    # for Qry
                if "INS" in SV_type.upper() and SV_size >= 10:
                    # perfect del
                    New_type = "3_DEL"
                    Qry_Evi = check_DEL(chr,sta,end,size,SV_size,New_type,Query,Qry_bam)
                else:
                    # like del, only consider it ref size
                    New_type = "4_SUB"
                    Qry_Evi = check_DEL(chr,sta,end,size,size,New_type,Query,Qry_bam)

                    if Qry_Evi == "NA" and size <100 and cells[Qry_Size_col] * 2 < cells[Ref_Size_col]:
                        size = cells[Ref_Size_col]
                        Qry_Evi = check_INS(chr,sta,end,size,SV_size,New_type,Query,Qry_bam)
                        New_type = "1_INS\t4_SUB"
                        if Qry_Evi.startswith("BreakDistance"):
                            BreakDistance = int(Qry_Evi.split("_")[1])
                            if BreakDistance > (abs(cells[Qry_Size_col]) / 3.0) and BreakDistance > 10:
                                Qry_Evi = "NA"

            #Evi_out.write(SV_ID + "\t" + Ref_Evi + "\t" + Qry_Evi + "\n")

            '''Genotype'''
            fail_set = set(["Too_Large","NA","INS_Gap"])
            Ref_Left,Ref_Right = gap_count(cells[Ref_chr_col],cells[Ref_start_col],cells[Ref_end_col],chr_size_dict[Reference],Reference)
            Qry_Left,Qry_Right = gap_count(cells[Qry_chr_col],cells[Qry_start_col],cells[Qry_end_col],chr_size_dict[Query],Query)
            #gap_out.write("\t".join(map(str,[line,Ref_Left,Ref_Right,Qry_Left,Qry_Right])) + "\n")

            Ref_Evi = Remove_bad_Evi(Ref_Evi,Ref_Left,Ref_Right)
            Qry_Evi = Remove_bad_Evi(Qry_Evi,Qry_Left,Qry_Right)


            if Ref_Evi  not in fail_set and Qry_Evi  not in fail_set:
                GT = "H"  # H for Heterozygous

            elif Ref_Evi  not in fail_set and Qry_Evi  in fail_set:
                #GT = "P"  # P for LA2093(SP)
                GT = "Q"  # Q for Query

            elif Ref_Evi  in fail_set and Qry_Evi  not in fail_set:
                #GT = "L"  # L for Heinz(SL)
                GT = "R"  # R for Reference

            else:
                GT = "U"  # U for Undetermined

            GT_out.write(SV_ID + "\t" + GT + "\n")







#sf
