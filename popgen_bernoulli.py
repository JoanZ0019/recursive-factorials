import re
import sys
from Bio import SeqIO

class Gene_sample:
    sampleID = ""
    date_of_sequencing = ""
    phenotype = ""
    nucleotide_sequence = ""

    def __init__(self, sampleID, date_of_sequencing, phenotype, nucleotide_sequence):
        self.sampleID = sampleID
        self.date_of_sequencing = date_of_sequencing
        self.phenotype = phenotype
        self.nucleotide_sequence = nucleotide_sequence

    def __str__(self):
        return '{0} {1} {2} {3}'.format(self.sampleID, self.date_of_sequencing, self.phenotype, self.nucleotide_sequence)

def get_phenotype(current_sequence, count_orange, count_blue):
        codon = current_sequence[9:12]
        phenotype = ''
        blue_list = list()
        if codon == 'CGA' or codon == 'CGT' or codon == 'CGG' or codon == 'CGC' or codon =='AGA' or codon == 'AGG' :
            phenotype = 'orange'
            count_orange += 1    
        elif codon == 'AGC' or codon == 'AGT'  or codon == 'TCA' or codon == 'TCG' or codon == 'AGT' or codon == 'AGC':
            phenotype = 'blue'
            count_blue += 1
            blue_list.append(Gene_sample(sampleID, date_of_sequencing, phenotype, nucleotide_sequence))
        else :
            return None
        return phenotype, count_orange, count_blue, blue_list
         

with open(sys.argv[1]) as file:

    count_orange = 0
    count_blue = 0
    genetic_sequence_list = list()
    blue_list = list()
    seq_record = list(SeqIO.parse(file, "fasta"))
    total_number_of_individuals = len(seq_record)
    for i in range(0, total_number_of_individuals):
        header = seq_record[i].id
        sampleID = header.split("_")[-2]
        date_of_sequencing = header.split("_")[-1]
        nucleotide_sequence = seq_record[i].seq
        phenotype, count_orange, count_blue = get_phenotype(nucleotide_sequence, count_orange, count_blue)
        genetic_sequence_list.append(Gene_sample(sampleID, date_of_sequencing, phenotype, nucleotide_sequence))
    print(count_orange)


def main():
    for i in genetic_sequence_list:
        print(i)
    for a in blue_list:
        print(a)
if __name__ == '__main__':
    main()    
