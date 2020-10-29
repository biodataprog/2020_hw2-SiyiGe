#!/usr/bin/env python3

import os, gzip, itertools, csv

# this is code which will parse FASTA files
# define what a header looks like in FASTA format
def isheader(line):
    return line[0] == '>'

def aspairs(f):
    seq_id = ''
    sequence = ''
    for header,group in itertools.groupby(f, isheader):
        if header:
            line = next(group)
            seq_id = line[1:].split()[0]
        else:
            sequence = ''.join(line.strip() for line in group)
            yield seq_id, sequence

url1="ftp://ftp.ensemblgenomes.org/pub/bacteria/release-45/fasta/bacteria_0_collection/salmonella_enterica_subsp_enterica_serovar_typhimurium_str_lt2/cds/Salmonella_enterica_subsp_enterica_serovar_typhimurium_str_lt2.ASM694v2.cds.all.fa.gz"
url2="ftp://ftp.ensemblgenomes.org/pub/bacteria/release-45/fasta/bacteria_0_collection/mycobacterium_tuberculosis_h37rv/cds/Mycobacterium_tuberculosis_h37rv.ASM19595v2.cds.all.fa.gz"
file1="Salmonella_enterica_subsp_enterica_serovar_typhimurium_str_lt2.ASM694v2.cds.all.fa.gz"
file2="Mycobacterium_tuberculosis_h37rv.ASM19595v2.cds.all.fa.gz"

if not os.path.exists(file1):
    os.system("curl -O %s"%(url1))

if not os.path.exists(file2):
    os.system("curl -O %s"%(url2))

with gzip.open(file1,"rt") as fh:
    seqs = aspairs(fh)

    for seq in seqs:
        seqname  = seq[0]
        seqstring= seq[1]
       # print(seqname, " first 10 bases are ", seqstring[0:10])

# The total number of genes in each species.

n=0
with gzip.open(file1, "rt") as fh:

    file1_2 = csv.reader(fh, delimiter="\t")
    for row in file1_2:
        if row[0].startswith(">"):
            n+=1
print("The total number of genes in Salmonella enterica is:{}".format(n))

i=0
with gzip.open(file2, "rt") as fh2:

    file2_2 = csv.reader(fh2, delimiter="\t")
    for row in file2_2:
        if row[0].startswith(">"):
            i+=1
print("The total number of genes in Mycobacterium tuberculosis is:{}".format(i))

# Total length of these gene sequences for each file
a=0
with gzip.open(file1,"rt") as fh:
    seqs = aspairs(fh)
    for seq in seqs:
        seqstring= seq[1]
        a+=len(seqstring)
print("The total length of genes in Salmonella enterica file is:{}".format(a))

b=0
with gzip.open(file2,"rt") as fh:
    seqs = aspairs(fh)
    for seq in seqs:
        seqstring= seq[1]
        b+=len(seqstring)
print("The total length of genes in Mycobacterium tuberculosis file is:{}".format(b))

# The G+C percentage for the whole dataset (eg the frequency of G + the frequency of C)
c=0
with gzip.open(file1,"rt") as fh:
    seqs = aspairs(fh)
    for seq in seqs:
        seqstring = seq[1]
        for base in seqstring:
            if base == "C" or base =="G":
                c+=1
print("The G+C percentage for the whole dataset of Salmonella enterica is :{:.2f}%".format(c/b*100))

d=0
with gzip.open(file2,"rt") as fh:
    seqs = aspairs(fh)
    for seq in seqs:
        seqstring = seq[1]
        for base in seqstring:
            if base == "C" or base =="G":
                d+=1
print("The G+C percentage for the whole dataset of Mycobacterium tuberculosis is :{:.2f}%".format(d/b*100))

# Total number codons in each genome.
print("The total number of codons in Salmonella enterica is :{:.0f}".format(a/3))
print("The total number of codons in Mycobacterium tuberculosis is :{:.0f}".format(b/3))

# Print out table with three columns: Codon, Frequency in Sp1, Frequency in Sp2
codons_dict_1 = {}
with gzip.open(file1, "rt") as fh:
    seqs = aspairs(fh)
    for seq in seqs:
        seqstring = seq[1]
        for index in range(0, len(seqstring), 3):
            codon1 = seqstring[index:index+3]
            if codon1 in codons_dict_1.keys():
                codons_dict_1[codon1]+=1
            else:
                codons_dict_1[codon1]=1

codons_dict_2 = {}
with gzip.open(file2, "rt") as fh:
    seqs = aspairs(fh)
    for seq in seqs:
        seqstring = seq[1]
        for index in range(0, len(seqstring), 3):
            codon2 = seqstring[index:index+3]
            if codon2 in codons_dict_2.keys():
                codons_dict_2[codon2]+=1
            else:
                codons_dict_2[codon2]=1

print("Codon", "Frequency in Sp1", "Frequency in Sp2")


for codon in codons_dict_1:
    print("{}".format(codon),"\t","{:.2f}%".format(codons_dict_1[codon]/(a/3)*100), "\t", "{:.2f}%".format(codons_dict_2[codon]/(b/3)*100))


