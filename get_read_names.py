import pysam
from sys import argv

input_bam = pysam.AlignmentFile(argv[1], "rb")

for line in input_bam:
    print(line.qname+"\t"+str(line.is_read1))
