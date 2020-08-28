#!/usr/bin/env python

import sys
import re

"""
samtools  view -h NXGL.bam NXGLv5ch02:32279286-32279287 | python Mismatch_constraint.py Mismatch_percentage
"""
Mismath = float(sys.argv[1])/100

for line in sys.stdin:
    line = line.strip()
    if line.startswith("@"):
        print line
    else:
        fields = line.split("\t")
        read_len = len(fields[9])
        if fields[5] != "*":
            NM = int(re.findall(r"NM:i:\d+",line)[0].split(":")[-1])
            if NM <= Mismath * read_len:
                print line
            elif "D" in fields[5] or "I" in fields[5]:
                cigar_strs = re.findall(r"\d+\w",fields[5])
                for cigar_str in cigar_strs:
                    if "I" in cigar_str or "D" in cigar_str:
                        NM -= int(cigar_str[:-1])
                if NM <= Mismath * read_len:
                    print line
















# s
