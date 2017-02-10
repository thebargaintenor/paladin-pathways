#! /usr/bin/env python3
"""
Perform taxonomy callback
"""
import subprocess
import csv
import sys


def get_uniprot_id(ec, filename):
    """
    Pull uniprot ID from csv output of paladin using Enzyme Comission number
    ec: enzyme comission number (str)
    filename: path to csv output of paladin (str)
    """
    acc = []
    with ope.san(filename) as csvfile:
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
    with open(filename) as f:
        for line in f:
            if not line[0] == '@':  # if not a header
                line = line.split()
                if uniprot_id in line[2]:  # and matching the uniprot_id
                    lines.append(line[0].split(':')[-1])
                    # append the contig header
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
    """
    
    """
    with open('blast_query', 'w') as blast_query:
        for item in los.items():
            header = item[0]
            seq = item[1]
            print('>' + header[1:], file=blast_query)
            print(seq, file=blast_query)


def blaster():
    command = ['blastx', '-query',  'blast_query',
               '-db', '/usr/local/share/paladin/uniref90.fasta',
               '-outfmt', '6 qacc sacc pident length mismatch gapopen \
               qstart qend sstart send evalue bitscore ssciname',
               '-out', 'blastout',
               '-num_threads', '32',
               '-num_alignments', '1',
               '-max_hsps', '1']
    subprocess.run(command, stdout=subprocess.PIPE)

def get_taxa(blast_results="blastout", db="/usr/local/share/paladin/uniref90.fasta"):
    uids = []
    with open(blast_results) as br:
        for line in br:
            uids.append(line.split()[1])
    uids = set(uids)
    taxa_dict = {}
    with open(db) as database:
        for line in database:
            if line[0] == '>':
                words = line.split()
                uniref_id = words[0][1:]
                if uniref_id in uids:
                    start = line.find("Tax")
                    stop = line[start:].find(' ')
                    taxa_dict[uniref_id] = line[start:stop + 1]
    return taxa_dict

def output_taxa(taxa_dict, f):
    for item in taxa_dict.items():
        print(" ".join(item), file=f)
if __name__ == "__main__":
    enzyme_code = sys.argv[1]
#   with open(enzyme_code) as ec_file:
#       enzyme_code = ec_file.readline().rstrip()
    pathways_out = sys.argv[2]
    paladin_out = sys.argv[3]
    reads = sys.argv[4]
    uid = get_uniprot_id(enzyme_code, pathways_out)
    loh = uniprot_coords(uid[0], paladin_out)
    los = find_seq(loh, reads)
    mkquery(los)
    blaster()
    taxa_dict = get_taxa()
    output_taxa(taxa_dict, "taxonomy")
