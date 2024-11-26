import pysam
from sys import argv

input_file=argv[1]
golden_bam=argv[2]

inbam=pysam.AlignmentFile(input_file, "rb")                                     
gbam=pysam.AlignmentFile(golden_bam, "rb")

output_file=open(input_file.replace(".bam","")+".mapstats.txt",'w')

num_mapped_only_to_right_chrom=0
num_mapped_perfectly=0
num_mapped_to_wrong_chrom=0
inbam_total = 0
reverse = 0
gbam_reads = {}
inbam_reads = {}

for gread in gbam:
    gbam_reads[gread.query_name] = [gread.reference_name,gread.reference_start,gread.is_read1,gread.query_sequence]
    
golden_depth = len(gbam_reads.keys())
number_perf_mapped=0
for read in inbam:
    if read.query_name in gbam_reads.keys():
        num_perf_mapped_per_gread=0
        gread = gbam_reads[read.query_name]

        if read.query_sequence == gread[3] and read.is_read1 == gread[2]:

            if read.reference_name == gread[0]:
                    
                if read.reference_start == gread[1]:
                    num_mapped_perfectly+=1
                else:
                    num_mapped_only_to_right_chrom+=1
            else:
                print(read.reference_name+" "+gread[0])
            inbam_total+=1

print("Number of golden reads:"+str(golden_depth))
print("Number of real reads:"+str(inbam_total))
print("Perfectly Mapped:"+str(num_mapped_perfectly))
print("mapped to wrong chromosome:"+str(num_mapped_to_wrong_chrom))
print("mapped only to right chromosome:"+str(num_mapped_only_to_right_chrom))

if inbam_total == 0:
    percentage = "N/A"
else:
    percentage = num_mapped_perfectly/golden_depth
print(str(golden_depth)+','+str(num_mapped_perfectly)+','+str(percentage))
