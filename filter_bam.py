import pysam
from sys import argv

input_bam=pysam.AlignmentFile(argv[1], "rb")
output=pysam.AlignmentFile(argv[2], "w", template=input_bam)

ids=open(argv[3],'r')
id_list={}
print("reading in ids")
for line in ids:
    name,is_r1=line.strip().split()
    id_list[name] = is_r1
print("done")
num_reads=len(id_list)
num=0
#id_set = set(id_list)
for line in input_bam:
    #print(line.qname)
    #print(type(line.qname))
    #print(id_list[0])
    #input(type(id_list[0]))
    #if line.qname in id_set:
    if line.qname in id_list.keys() and id_list[line.qname] == str(line.is_read1) :
        print("found read "+line.qname)
        num+=1
        output.write(line)
    if num == num_reads:
        break
print("DONE")
