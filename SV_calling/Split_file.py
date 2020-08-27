"""
By Lei Gao
usage: Split_file.py [-h] --Infile INFILE --Split_Number SPLIT_NUMBER
                     --Suffix_length SUFFIX_LENGTH --Out_prefix OUT_PREFIX

optional arguments:
  -h, --help            show this help message and exit
  --Infile INFILE       Input file
  --Split_Number SPLIT_NUMBER
                        Number of pieces you want to generate
  --Suffix_length SUFFIX_LENGTH
                        Generate suffixes of length N
  --Out_prefix OUT_PREFIX
                        prefix for output files

"""
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--Infile", type=str, help="Input file", required=True, default="")
parser.add_argument("--Split_Number", type=str, help="Number of pieces you want to generate", required=True, default="")
parser.add_argument("--Suffix_length", type=str, help="Generate suffixes of length N", required=True, default="")
parser.add_argument("--Out_prefix", type=str, help="prefix for output files", required=True, default="")

args = parser.parse_args()

Infile = args.Infile
Split_Number = int(args.Split_Number)
Suffix_length = int(args.Suffix_length)
Out_prefix = args.Out_prefix

f ="{0:0" + str(Suffix_length) + "}"

''' Creat empty output files '''
for piece in range(0,Split_Number):
    outfile_name = Out_prefix + "_" + f.format(piece)
    outfile = open(outfile_name,"w")
    outfile.close()

''' split file  '''
i = -1
with open(Infile) as infile:
    for line in infile:
        i += 1
        piece = i % Split_Number
        outfile_name = Out_prefix + "_" + f.format(piece)
        with open(outfile_name, "a+") as outfile:
            outfile.write(line)
















#srf
