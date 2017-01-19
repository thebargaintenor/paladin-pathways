#! /usr/bin/env python3
"""
Perform taxonomy callback 
"""
#import numpy as np
import subprocess
import csv
def get_uniprot_id(ec, filename):
    """
    Pull uniprot ID from csv output of paladin using Enzyme Comission number
    ec: enzyme comission number (str)
    filename: path to csv output of paladin (str)
    """
    acc = []
    with open(filename) as csvfile:
        f = csv.reader(csvfile)
        for line in f:
            if line[1].strip() == ec.strip():
                acc.append(line[0])
        return acc
    
def uniprot_coords(uniprot_id, filename):
    """
    Pull the header for the contig file from the .sam output file of paladin
    uniprot_id: a uniprot_id (str)
    filename: path to sam file (str)
    """
    lines = []
    index = 0
    with open(filename) as f:
        for line in f:
            if not line[0] == '@':  ## if not a header
                line = line.split()
                if uniprot_id in line[2]: ## and matching the uniprot_id
                    lines.append(line[0].split(':')[-1])  ## append the contig header
    return lines


def find_seq(loh, filename):
    """                                                                                                            
    Return a list of seqeunces from fastq file (filename) corresponding
    to the headers in loh
    loh: list of headers (list(string))
    filename: contig.fasta file to search for the list of headers of                              
    Returns:                                                                                                  
    los: dictionary with the contig headers as keys and sequences as values     
    """
    flag = False
    header = ''
    los = {}
    with open(filename) as f:
        for line in f:
            if flag:
                los[header] = line.rstrip()
                flag = False
                header = ''
            else:
                if line[0] == '@':
                    if line[1:-3] in loh:
                        header = line.rstrip()
                        flag = True
    return los

def mkquery(los):
    with open('blast_query', 'w') as blast_query:
        for item in los.items():
            header = item[0]
            seq = item[1]
            print('>' + header[1:], file= blast_query)
            print(seq, file = blast_query)
def blaster():
    command = ['blastx', '-query',  'blast_query', '-db', '/usr/local/share/paladin/uniref90.fasta',  '-outfmt', '6', '-out', 'blastout', '-num_threads', '32', '-num_alignments', '1', '-max_hsps', '1']
    blasting = subprocess.run(command, stdout=subprocess.PIPE)

if __name__== "__main__":
    import time
    uid = get_uniprot_id('1.4.1.13', '/home/work/Documents/pathways/paladin-pathways/Pathways_Output/B.japonicum_100_test.csv')
    loh = uniprot_coords(uid[0], '/home/work/Documents/pathways/paladin-pathways/Paladin_Output/B.japonicum_100.sam')
    st1 = time.time()
    los = find_seq(loh, "../B.japonicum_100.1.fq")
    #print(time.time()-st1)
    st1 = time.time()
    #los2 = find_seq2(loh, "../B.japonicum_100.1.fq")
    mkquery(los)
    blaster()
