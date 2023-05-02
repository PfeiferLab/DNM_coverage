import random
from sys import argv
from Bio import Seq
from Bio import SeqIO
import math
import scipy.stats as st
import argparse
seq = ['A','G','C','T']
num_list = [0,1,2,3]

###############################################################
####Mean and SD are from Tatsumoto nonMV, filtered variants####
####mean=0.4945519#############################################
####sd=0.05927348##############################################
###############################################################


parser = argparse.ArgumentParser(description="Creates a VCF file with a mutation rate equivalent to the mu input in positions bounded by the input FASTA.")
parser.add_argument("-i","--input-fasta",dest="input_file",help="The input FASTA file.")
parser.add_argument("-u","--mutation-rate",dest="mutation_rate",help="The input FASTA file.")
parser.add_argument("-o","--output-VCF",dest="output_file",help="Output VCF file.")
args = parser.parse_args()

input_file = SeqIO.parse(open(args.input_file),'fasta')
mutation_rate = float(args.mutation_rate)
mutation_file = open(args.output_file,'w')

ab_mean=0.50
ab_sd=0
contig_lines=[]
mutation_lines=[]
for item in input_file:
    genome_1_seq = item.seq
    header = item.id
#print(str(genome_1_seq))
    mut_pos = []
    mean = mutation_rate * len(genome_1_seq)

    RoundUp = bool(round(random.random()))

    if random.random() < 0.5:
        num_mutations = int(math.ceil(random.normalvariate(mean,math.sqrt(mean))))
    else:
        num_mutations = int(math.floor(random.normalvariate(mean,math.sqrt(mean))))
    percentile = st.norm.cdf(num_mutations,loc=mean,scale=math.sqrt(mean))
    mutation = ""
    contig_lines.append("##contig=<ID="+header+",length="+str(len(genome_1_seq))+'>\n')
    mutations = random.sample(range(0,len(genome_1_seq)-1),k=num_mutations)
    mutations.sort()
    print(mutations)
    for pos in mutations:
        mutation = ""
        while mutation == "" or mutation.upper() == genome_1_seq[pos].upper():
            #print("ENTERED")
            mutation = random.choice(seq)
            #print(mutation)
        ancestral_allele = genome_1_seq[pos].upper()
        #genome_1_seq = genome_1_seq[:pos]+mutation+genome_1_seq[pos+1:]
        mutation_lines.append(str(header)+"\t"+str(pos+1)+"\t.\t"+ancestral_allele+"\t"+mutation+"\t1000\tPASS\tAF="+str(random.normalvariate(ab_mean,ab_sd))+"\tGT\t0/0\t0/0\t0/1\n")
#print(str(genome_1_seq))

mutation_file.write("##fileformat=VCFv4.2\n")                               
mutation_file.write("##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency among genotypes, for each ALT allele, in the same order as listed\">\n")
mutation_file.write("##INFO=<ID=ST,Number=A,Type=Float,Description=\"mean:"+str(mutation_rate)+"\tmu_used:"+str(st.norm.ppf(percentile,loc=mean/len(genome_1_seq),scale=math.sqrt(mean/len(genome_1_seq))))+"\tmu_calculation:"+str(num_mutations/len(genome_1_seq))+"\tPercentile:"+str(percentile)+"\">\n")
mutation_file.write("##FORMAT=<ID=GT,Number=A,Type=Float,Description=\"Genotype\""+">\n")
for line in contig_lines:
    mutation_file.write(line)
mutation_file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tparent_1\tparent_2\tchild\n")
for line in mutation_lines:
    mutation_file.write(line)
