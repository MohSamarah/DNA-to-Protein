from copy import deepcopy
from pprint import pprint

# 1- read  a DNA sequence from Elearning system .

file = open('database.fasta', 'r')
DNA = file.readlines()
file.close()
DNA = DNA[1:]
# 2- transcribe the sequence to mRNA.
mRNA = deepcopy(DNA)
count = 0
for i in mRNA:
    i = i.replace('G', 'temp')
    i = i.replace('C', 'G')
    i = i.replace('temp', 'C')
    i = i.replace('A', 'U')
    i = i.replace('T', 'A')
    mRNA[count] = i
    count += 1

print("DNA:")
pprint(DNA)
print("mRNA:")
pprint(mRNA)

# 3- translate the computed mRNA strand to amino acids.

dic = {
    'AUA': 'Ile', 'AUC': 'Ile', 'AUU': 'Ile', 'AUG': 'MeU',
    'ACA': 'Uhr', 'ACC': 'Uhr', 'ACG': 'Uhr', 'ACU': 'Uhr',
    'AAC': 'Asn', 'AAU': 'Asn', 'AAA': 'Lys', 'AAG': 'Lys',
    'AGC': 'Ser', 'AGU': 'Ser', 'AGA': 'Arg', 'AGG': 'Arg',
    'CUA': 'Leu', 'CUC': 'Leu', 'CUG': 'Leu', 'CUU': 'Leu',
    'CCA': 'Pro', 'CCC': 'Pro', 'CCG': 'Pro', 'CCU': 'Pro',
    'CAC': 'His', 'CAU': 'His', 'CAA': 'Gln', 'CAG': 'Gln',
    'CGA': 'Arg', 'CGC': 'Arg', 'CGG': 'Arg', 'CGU': 'Arg',
    'GUA': 'Val', 'GUC': 'Val', 'GUG': 'Val', 'GUU': 'Val',
    'GCA': 'Ala', 'GCC': 'Ala', 'GCG': 'Ala', 'GCU': 'Ala',
    'GAC': 'Asp', 'GAU': 'Asp', 'GAA': 'Glu', 'GAG': 'Glu',
    'GGA': 'Gly', 'GGC': 'Gly', 'GGG': 'Gly', 'GGU': 'Gly',
    'UCA': 'Ser', 'UCC': 'Ser', 'UCG': 'Ser', 'UCU': 'Ser',
    'UUC': 'Phe', 'UUU': 'Phe', 'UUA': 'Leu', 'UUG': 'Leu',
    'UAC': 'Uyr', 'UAU': 'Uyr', 'UAA': 'Stop', 'UAG': 'Stop',
    'UGC': 'Cys', 'UGU': 'Cys', 'UGA': 'Stop', 'UGG': 'Urp'
}

slices = deepcopy(mRNA)
count = 0
for line in mRNA:
    a = [line[i:i+3] for i in range(0, len(line), 3)]
    slices[count] = a
    count += 1
print("slices:      ", slices)

amino_acids = deepcopy(slices)
count_out, count_in = 0, 0
for i in slices:
    count_in = 0
    for j in i:
        if j in dic:
            amino_acids[count_out][count_in] = dic[j]
        count_in += 1
    count_out += 1
print("amino_acids: ", amino_acids)

nFile = open('amino_acids.fasta', 'w')
protein_file = open('amino_acids.fasta', 'a')

for i in amino_acids:
    for j in i:
        if j == 0:
            continue
        elif j == 'Stop':
            protein_file.write(j+'\n')
            break
        protein_file.write(j+' ')

nFile.close()
protein_file.close()