# PALADIN _Pathways_

The goal of this project was to create an extension of [PALADIN](https://github.com/twestbrookunh/paladin) to elaborate on the pathway output produced by PALADIN Align. In order to create an extension for the PALADIN program to output functional pathway data, we compared several pathway databases or size, accessibility, level of completeness, and inclusion of both common and uncommon bacterial pathways. The analytic scripts for the Pathways pipeline were written in both Python and bash to produce a working bioinformatics tool. Pathways was applied to a sample dataset from PCE contaminated drinking water and analyzed for the chloroalkane and chloroalkene degradation pathway. The results show the enrichment of the enzymes needed for chloroalkene degradation through two visualization tools and demonstrated the potential for the Pathways tool to search for community function within WMS data. Our PALADIN Pathways tool allows the user to investigate whole metagenomes for well-characterized protein pathways and determine the relative completeness. It is our hope that PALADIN Pathways will provide novel insights into microbial communities and streamline functional metagenomic studies.

## INSTALLATION
### Dependencies
PALADIN _Pathways_ has only a handful of dependencies enumerated in `install.sh` needed to support the python scripts.  These are all handled with the pip package manager, and may require root/administrator privileges for installation on your workstation.

Overall, the process is as follows (for Linux or macOS):

```
git clone https://www.github.com/thebargaintenor/paladin-pathways.git
cd paladin-pathways
sudo bash install.sh # sudo may not be necessary
```

## INSTRUCTIONS
Coming soon to a repo near you!

In the meantime, please have the usage instructions from `pathways.py`:

```
Usage: pathways [-k pathway] [-p paladin_tsv] [-o output] [-c counts] [-v]
    Script runs entire pathways pipeline, filtering a read report from PALADIN
    to calculate the mathmatical completeness of a given metabolic pathway.

    -k pathway      - KEGG ID for pathway
    --kegg
    -p paladin_tsv  - path to PALADIN output report
    --paladin
    -o output       - output file name
    --output
    -c counts       - supplement data file name
    -v              - verbose mode (print more info to stdout)
    --verbose

    EXAMPLE: python3 pathways.py -k 00625 -p ./data/report.tsv -o results.csv -c counts.dat 
```

Output will also provide the names of any files produced during the process.  Included is the completeness results, JSON-formatted pathway markup, and a log file that contains the output from main script.  Output redirection is fine too.

```
python3 pathways.py -k 00625 -p ./data/report.tsv -o results.csv -c counts.dat -v > foo.log 
```

## ACKNOWLEDGEMENTS

In addition to PALADIN (without which we'd essentially have no project), the _Pathways_ dev team would also like to give a shout-out to a few other projects that made our lives easier.  (And possibly yours as well!)

* [D3.js](https://d3js.org/)
* [palette.js](https://github.com/google/palette.js)
* [d3heatmap](https://github.com/google/palette.js)
* [dataset](https://dataset.readthedocs.io/): databases for lazy people
* [Requests](https://requests.readthedocs.io/): HTTP for Humans
* And of course, [KEGG](http://www.kegg.jp/) for providing all the pathway annotations and KGML

Thanks!