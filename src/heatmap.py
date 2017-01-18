#! /usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import sys
files = sys.argv[1:]
genomes = []
genome_names = []
en = []
for file_path in files:
    with open(file_path) as f:
        f = f.readlines()
    enzymes = {}
    genome_names.append(f[0].rstrip().split(',')[2])
    for line in f[1:]:
        words = line.rstrip().split(',')
        enzyme_name = "".join(words[1:-1])
        print(enzyme_name)
        enzymes.update({enzyme_name:int(float(words[-1]))})
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
plt.figure(figsize=(19.2,10.8))
plt.imshow(1*(copymat        ), interpolation='none', cmap=cm.Reds, aspect='auto')
plt.xticks(np.arange(len(genome_names)), genome_names, rotation=20)
plt.yticks(range(len(en)), en_disp)

cb = plt.colorbar()
cb.set_label('log copy number')
plt.tight_layout()
plt.show()
