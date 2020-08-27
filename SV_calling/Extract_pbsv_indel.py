"""
By Lei Gao
usage: Extract_pbsv_indel.py [-h] --vcf_file VCF_FILE --ref_genome REF_GENOME
                             --query_genome QUERY_GENOME --Output OUTPUT

optional arguments:
  -h, --help            show this help message and exit
  --vcf_file VCF_FILE   pbsv output vcf file
  --ref_genome REF_GENOME
                        Reference genome file, fasta format
  --query_genome QUERY_GENOME
                        Query species genome file, fasta format.
  --Output OUTPUT       Output prefix


"""
import argparse
from pathlib import Path
from operator import itemgetter
import subprocess
import re
import os

parser = argparse.ArgumentParser()
parser.add_argument("--vcf_file", type=str, help="pbsv output vcf file", required=True, default="")
parser.add_argument("--ref_genome", type=str, help="Reference genome file, fasta format", required=True, default="")
parser.add_argument("--query_genome", type=str, help="Query species genome file, fasta format.", required=True, default="")
parser.add_argument("--Output", type=str, help="Output prefix", required=True, default="")

args = parser.parse_args()

vcf_file = args.vcf_file
ref_genome = args.ref_genome
query_genome = args.query_genome
Output_prefix = args.Output
# log
print "Input vcf file: " + vcf_file
print "Reference genome file: " + ref_genome
print "Query species' genome: " + query_genome
print "Output prefix: " + Output_prefix
print "#" * 20

# chr size
chr_size = {}
with open(ref_genome + ".fai") as fai_file:
    for line in fai_file:
        fields = line.strip().split("\t")
        chr_size[fields[0]] = int(fields[1])

''' Main Program '''

sv_count = {}
gap_indel_count = {}
sv_types = ["DEL","INS"]

for svType in sv_types:
    sv_count[svType] = 0
    gap_indel_count[svType] = 0


with open(vcf_file) as input_file, \
     open(Output_prefix + ".Indel", "w") as Indel_file, \
     open(Output_prefix + ".gap_Indel", "w") as gap_file:
    Indel_out = ["SV_id", "Type", "GT", "Chromosome", "Start", "End", "Ref_seq", "Query_seq", "Ref_size", "Query_size", "Ref_Gap%", "Query_Gap%", "Left_flank", "Right_flank"]
    Indel_file.write("\t".join(Indel_out) + "\n")
    gap_file.write("\t".join(Indel_out) + "\n")


    for line in input_file:
        if not line.startswith("#"):
            fields = line.strip().split("\t")
            filter = fields[6]
            svType = re.findall(r'SVTYPE=(\w+)',fields[7])[0]
            if (svType == "DEL" or svType == "INS") and filter == "PASS":
                svOrder = int(fields[2].split(".")[-1])
                sv_id =  svType + '{0:06}'.format(svOrder)
                GT = fields[9]
                Chromosome = fields[0]
                Start = int(fields[1])
                End = int(re.findall(r'END=(\d+)',fields[7])[0])
                Ref_seq = fields[3]
                Query_seq = fields[4]
                Ref_size = len(Ref_seq)
                Query_size = len(Query_seq)
                if "N" in fields[3] + fields[4]:
                    gap_indel_count[svType] += 1
                    ref_gap = Ref_seq.count("N")  * 100.0 / Ref_size
                    qry_gap = Query_seq.count("N")  * 100.0 / Query_size

                else:
                    # non-gap in both
                    # only keep this
                    sv_count[svType] += 1
                    ref_gap = 0
                    qry_gap = 0

                Indel_out = [sv_id, svType, GT, Chromosome, Start, End, Ref_seq, Query_seq, Ref_size, Query_size, ref_gap, qry_gap]
                if "N" in fields[3] + fields[4]:
                    gap_file.write("\t".join(map(str,Indel_out)) + "\n")
                else:
                    Indel_file.write("\t".join(map(str,Indel_out)) + "\n")

# log
print "Kept indels: " + str(sv_count["DEL"] + sv_count["INS"])
print "DEL: " + str(sv_count["DEL"])
print "INS: " + str(sv_count["INS"])

print "Indels spanning gaps: " + str(gap_indel_count["DEL"] + gap_indel_count["INS"])
print "DEL: " + str(gap_indel_count["DEL"])
print "INS: " + str(gap_indel_count["INS"])



























#srf
