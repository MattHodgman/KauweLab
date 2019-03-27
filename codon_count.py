import sys
import csv
import pandas as pd
from Bio import SeqIO
            
# modifies the nested dict with percent of total for each codon
def calc_perc(sample):
    for seq_record in SeqIO.parse(sample, "fasta"):
        size = len(seq_record.seq) / 3
        for gene in genesDict:
            if gene == seq_record.id[3:11]:
                for codon in genesDict[gene].keys():
                    if codon != 'Gene':
                        genesDict[gene][codon] = genesDict[gene][codon] / size

# check to see if the codon has invalid characters                   
def invalid_codon(codon):
    for c in codon:
        if c != 'A' and c != 'T' and c != 'G' and c != 'C':
            return True
        elif c == '\n':
            return True
        
       
def get_isoform_number(gene):
    n = 0
    
    for key in genesDict:
        if key.startswith(gene):
            n += 1

    return n

# ignores individual genes and instead gets total counts for the entire sample
def sample_count(sample):
    CodonsDict['Sample'] = sample[17:]; # get sample name
    
    # identify sample population
    population = ""
    df_pop = pd.read_csv('sample_population_key.csv')
    for i in range(len(df_pop)):
        if df_pop['Sample'][i] == sample[17:]:
            population = df_pop['Population'][i]
    
    # parse the sample
    for seq_record in SeqIO.parse(sample, "fasta"):
        total_size = 0;
        
        # get the gene name        
        gene = seq_record.id[3:seq_record.id.index(";")]
        
        # check if isoform
        if gene in genesDict:
            gene = gene + '(' + str(get_isoform_number(gene)) + ')'
            print(gene)
            
        # add gene to dict
        CodonsDict['Gene'] = gene;
        genesDict[gene] = CodonsDict.copy()
        
        # set gene size (number of codons)
        seq_size = len(seq_record.seq) - (len(seq_record.seq) % 3)
        total_size += seq_size / 3
        index = 0
        
        # loop through codons
        while index < seq_size:
            codon = seq_record.seq[index:index+3] # get the next codon
            init_index = index
            
            # check if codon is valid
            while invalid_codon(codon):
                index += 1
                codon = seq_record.seq[index:index+3] # try the next codon
            
            # adjust total_size for invalid codons
            total_size -= ((index - init_index) - ((index - init_index) % 3)) / 3
           
            # increment codon count
            genesDict[gene][codon] += 1
                
            # increment index to adjust reading frame
            index += 3
    
        # add population and size data to dict
        genesDict[gene]['Population'] = population
        genesDict[gene]['Size'] = total_size

genesDict = {}
columns = ['Sample', 'Population', 'Gene', 'Isoform', 'Size', 
           'TTT', 'TTC',
           'TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG',
           'TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC',
           'TAT', 'TAC',
           'TAA', 'TAG', 'TGA',
           'TGT', 'TGC',
           'TGG',
           'CCT', 'CCC', 'CCA', 'CCG',
           'CAT', 'CAC',
           'CAA', 'CAG',
           'CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG',
           'ATT', 'ATC', 'ATA',
           'ATG',
           'ACT', 'ACC', 'ACA', 'ACG',
           'AAT', 'AAC',
           'AAA', 'AAG',
           'GTT', 'GTC' ,'GTA', 'GTG',
           'GCT', 'GCC', 'GCA', 'GCG',
           'GAT', 'GAC',
           'GAA', 'GAG',
           'GGT', 'GGC', 'GGA', 'GGG']

CodonsDict = { 
    'Sample': '', 'Population': '', 'Gene': '', 'Isoform': 0, 'Size': 0, 
    'TTT': 0, 'TTC': 0, 'TTA': 0, 'TTG': 0, 'CTT': 0, 
    'CTC': 0, 'CTA': 0, 'CTG': 0, 'ATT': 0, 'ATC': 0, 
    'ATA': 0, 'ATG': 0, 'GTT': 0, 'GTC': 0, 'GTA': 0, 
    'GTG': 0, 'TAT': 0, 'TAC': 0, 'TAA': 0, 'TAG': 0, 
    'CAT': 0, 'CAC': 0, 'CAA': 0, 'CAG': 0, 'AAT': 0, 
    'AAC': 0, 'AAA': 0, 'AAG': 0, 'GAT': 0, 'GAC': 0, 
    'GAA': 0, 'GAG': 0, 'TCT': 0, 'TCC': 0, 'TCA': 0, 
    'TCG': 0, 'CCT': 0, 'CCC': 0, 'CCA': 0, 'CCG': 0, 
    'ACT': 0, 'ACC': 0, 'ACA': 0, 'ACG': 0, 'GCT': 0, 
    'GCC': 0, 'GCA': 0, 'GCG': 0, 'TGT': 0, 'TGC': 0, 
    'TGA': 0, 'TGG': 0, 'CGT': 0, 'CGC': 0, 'CGA': 0, 
    'CGG': 0, 'AGT': 0, 'AGC': 0, 'AGA': 0, 'AGG': 0, 
    'GGT': 0, 'GGC': 0, 'GGA': 0, 'GGG': 0}


i = 0
for sample in sys.argv:
    if i != 0:
        sample_count(sample)
    i += 1
    
# create csv file
with open(sample + ".csv", "w") as f:
    w = csv.DictWriter(f, columns)
    w.writeheader()
    for sample in genesDict:
        w.writerow(genesDict[sample])
