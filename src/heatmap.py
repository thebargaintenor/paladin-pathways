#! /usr/bin/env python3
"""
Create heatmap plots from files.

Usage heatmap.py files
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import sys


def parse_files(list_of_files, binary=False):
    genomes = []
    genome_names = []
    en = []
    for file_path in list_of_files:  # open each of the files
        with open(file_path) as f:
            f = f.readlines()
        enzymes = {}
        genome_names.append(f[0].rstrip().split(',')[2])
        # pull out the naming info for labeling
        for line in f[1:]:
            words = line.rstrip().split(',')
            enzyme_name = "".join(words[1:-1])
            enzymes.update({enzyme_name: int(float(words[-1]))})
            en.append(enzyme_name)
        genomes.append(enzymes)
    en = list(set(en))
    en_disp = []

    for name in en:
        en_disp.append(name.replace('"', ""))
    copymat = []
    for enz in en:
        copyn = []
        for genome in genomes:
            try:
                copyn.append(genome[enz])
            except:
                copyn.append(0)
        copymat.append(copyn)
    copymat = np.asarray(copymat, dtype=int)
    if binary:
        copymat = copymat / copymat
    return copymat, genome_names, en_disp


def render(copymat, genome_names, en_disp, outname='hist.png', colormap=cm.Reds, figsize = (19.2, 10.8)):
    plt.figure(figsize=figsize)
    plt.imshow(1*(copymat), interpolation='none', cmap=colormap, aspect='auto')
    plt.xticks(np.arange(len(genome_names)), genome_names, rotation=20)
    plt.yticks(range(len(en_disp)), en_disp)
    cb = plt.colorbar()
    cb.set_label('copy number')
    plt.tight_layout()
    plt.savefig(outname)


if __name__ == '__main__':
    copymat, genome_names, en_disp = parse_files(sys.argv[1:])
    render(copymat, genome_names, en_disp)
