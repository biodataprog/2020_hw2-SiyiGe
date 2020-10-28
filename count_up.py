gff = "Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.37.gff3.gz"
fasta = "Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna.chromosome.Chromosome.fa.gz"

import os, gzip, itertools, csv, re


# this is code which will parse FASTA files
# define what a header looks like in FASTA format
def isheader(line):
    return line[0] == '>'


def aspairs(f):
    seq_id = ''
    sequence = ''
    for header, group in itertools.groupby(f, isheader):
        if header:
            line = next(group)
            seq_id = line[1:].split()[0]
        else:
            sequence = ''.join(line.strip() for line in group)
            yield seq_id, sequence


if not os.path.exists(gff):
    os.system(
        "curl -O ftp://ftp.ensemblgenomes.org/pub/bacteria/release-45/gff3/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.37.gff3.gz")

if not os.path.exists(fasta):
    os.system(
        "curl -O ftp://ftp.ensemblgenomes.org/pub/bacteria/release-45/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/dna/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna.chromosome.Chromosome.fa.gz")

with gzip.open(gff, "rt") as fh:
    # now add code to process this
    gff = csv.reader(fh, delimiter="\t")
    for row in gff:
        if row[0].startswith("#"):
            continue
        print(row[3], row[6])

# count up and print the number of genes

n=0
total_length=0
gff2 = "Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.37.gff3.gz"
with gzip.open(gff2, "rt") as fh:
    gff = csv.reader(fh, delimiter="\t")
    for row in gff:
        if row[0].startswith("Chromosome") and row[2]=="gene":
           total_length += int(row[4])-int(row[3])

        if row[0].startswith("###"):
            n += 1
print("The number of genes is:{}".format(n))
print("The total length of the genes is:{}".format(total_length))


# Compute the total length of the genes (length is the END - START)

import itertools
import gzip
import sys
import re


# based on post here
# https://drj11.wordpress.com/2010/02/22/python-getting-fasta-with-itertools-groupby/

# define what a header looks like in FASTA format
def isheader(line):
    return line[0] == '>'


# this function reads in fasta file and returns pairs of data
# where the first item is the ID and the second is the sequence
# it isn't that efficient as it reads it all into memory
# but this is good enough for our project
def aspairs(f):
    seq_id = ''
    sequence = ''
    for header, group in itertools.groupby(f, isheader):
        if header:
            line = next(group)
            seq_id = line[1:].split()[0]
        else:
            sequence = ''.join(line.strip() for line in group)
            yield seq_id, sequence



# here is my program
# get the filename from the cmdline
with gzip.open(fasta, "rt") as fh:
# with open(filename, "r") as f:
    pairs = aspairs(fh)
    seqs = dict(pairs)
#    seqs = dict(aspairs(f))

# iterate through the sequences
n = 0
# for seqid in seqs:
#    print("id is ",seqid, "seq is ",seqs[seqid])

# for k, v in seqs.items():
#     print("id is ", k, "seq is", v)
#     n += 1

# print(n, "sequences")

print("The total length of the genome is:{}".format(len(seqs["Chromosome"])))

# Print out the percentage of the genome which is coding
print("The percentage of the genome which is coding:{:.2f}%".format(total_length/len(seqs["Chromosome"])*100))