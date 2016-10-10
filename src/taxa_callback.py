#! /usr/bin/env python3
"""
Perform taxonomy callback 
"""
import numpy as np
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
    lines = []
    index = 0
    with open(filename) as f:
        for line in f:
            if not line[0] == '@':
                line = line.split()
                if uniprot_id in line[2]:
                    lines.append(line[0].split(':')[-1])
    return lines

def find_seq(loh, filename):
    """ 
    Return a list of seqeunces from fastq file (filename) corresponding 
    to the headers in loh
    """
    flag = False
    header = ''
    los = {}
    with open(filename) as f:
        lines = f.readlines(100000)
        while lines:
            for line in lines:
                if flag:
                    los[header] = line.rstrip()
                    flag = False
                    header = ''
                else:
                    if line[0] == '@':
                        if line[1:-3] in loh:
                            header = line.rstrip()
                            flag = True
            lines = f.readlines(100000)
    return los
        
if __name__== "__main__":
   uid = get_uniprot_id('1.4.1.13', '/home/work/Documents/pathways/paladin-pathways/Pathways_Output/B.japonicum_100_test.csv')
   loh = uniprot_coords(uid[0], '/home/work/Documents/pathways/paladin-pathways/Paladin_Output/B.japonicum_100.sam')
   los = find_seq(loh, "../B.japonicum_100.1.fq")
