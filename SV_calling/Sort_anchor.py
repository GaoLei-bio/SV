"""
By Lei Gao
python Sort_anchor.py prefix.coords.csv > prefix.Sorted_anchor.csv

"""
import sys
from operator import itemgetter

def convert_int(fields):
    for i in range(0,len(fields)):
        if fields[i].isdigit() or (fields[i].startswith("-") and fields[i][1:].isdigit()):
            fields[i] = int(fields[i])
    return fields

i = -1
my_table = []
with open(sys.argv[1]) as infile:
    for line in infile:
        cells = convert_int(line.strip().replace("ch00","chXX").split(","))
        i += 1
        if i == 0:
            print "\t".join(map(str,cells))
        else:
            my_table.append(cells)


sorted_table = sorted(my_table, key = itemgetter(6,7,0,1,2,3))

for cells in sorted_table:
    print "\t".join(map(str,cells)).replace("chXX","ch00")





#srf
