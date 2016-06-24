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