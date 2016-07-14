#!/usr/bin/env python3

# pathways.py
#
# 2016-06-20 tdm - script created
# 2016-07-08 tdm - added support for reactions with multiple enzyme dependencies 
#
# ----
#
# PATHWAYS
#   Script runs entire pathways pipeline, filtering a read report from PALADIN
#   to calculate the mathmatical completeness of a given metabolic pathway.
#
# ----
#
# EXAMPLE: python3 pathways.py -k 00625 -p ./data/report.tsv -o results.csv -c counts.dat

# library imports
import json
import sys
import time
import csv

# external script imports
import paladin_postprocess
import kegg
from ec_lookup import lookup

# -- GLOBAL VARIABLES --

usageMsg = '''

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

    EXAMPLE: python3 pathways.py -k 00625 -p ./data/report.tsv -o results.csv -c counts.dat --verbose

'''

# -- FUNCTIONS --

def usage():
    '''If we have missing arguments or argument is -h print some "usage" information and quit.'''
    if len(sys.argv) < 1 or sys.argv[1] == "-h":  # command
        print(usageMsg)
        exit(0)

def log(f, line, verbose = False):
    if verbose:
        print(line)

    with open(f, 'at') as handle:
        handle.write(line + '\n')

def main():
    '''Main'''
    usage()

    kegg_id = None
    paladin_report = None
    output_csv = None
    count_csv = None
    verbose = False

    # check args
    for i in range(1, len(sys.argv)):
        if sys.argv[i] == '--kegg' or sys.argv[i] == '-k':
            kegg_id = sys.argv[i + 1]
        elif sys.argv[i] == '--paladin' or sys.argv[i] == '-p':
            paladin_report = sys.argv[i + 1]
            if not paladin_report.endswith('.tsv'):
                paladin_report += '.tsv'
        elif sys.argv[i] == '--output' or sys.argv[i] == '-o':
            output_csv = sys.argv[i + 1]
            if not output_csv.endswith('.csv'):
                output_csv += '.csv'
        elif sys.argv[i] == '--counts' or sys.argv[i] == '-c':
            count_csv = sys.argv[i + 1]
            if not count_csv.endswith('.csv'):
                count_csv += '.csv'
        if sys.argv[i] == '--verbose' or sys.argv[i] == '-v':
            verbose = True

    suffix = time.strftime('_%m%d%y_%H%M')
    log_file = 'pathways_{0}{1}.log'.format(kegg_id, suffix)

    log(log_file, 'Fetching pathway information for KEGG ID ' + kegg_id, verbose)

    # retrieve dictionary of pathway information
    pathway = kegg.get(kegg_id)

    if pathway:
        pathway = pathway['pathway']
        files_created = []

        # make name for filtered output
        filtered_paladin_csv = paladin_report[:-4] + suffix + '.csv'

        log(log_file, 'Running TSV post-process...', verbose)
        paladin_postprocess.run(paladin_report, filtered_paladin_csv, verbose)

        files_created.append(filtered_paladin_csv)

        enzymes = []
        enzyme_counts = {}
        enzyme_names = {}

        # initialize data structures
        for entry in pathway['entry']:
            if entry['@type'] == 'enzyme':
                # assume ALL listed refs are needed when listed in @name
                ecs = [ec[3:] for ec in entry['@name'].split(' ')]
                for ec in ecs:
                    if not ec in enzyme_names:
                        enzymes.append(ec)
                        enzyme_names[ec] = ''
                        enzyme_counts[ec] = 0

        # initialize list of enzyme IDs (for simpler use with generators)
        for enzyme in enzymes:
            enzyme_counts[enzyme] = 0
            enzyme_names[enzyme] = ''

        log(log_file, 'Calculating pathway completeness...', verbose)

        # get file handle for filtered PALADIN report
        with open(filtered_paladin_csv, 'r') as paladin_file:
            paladin_records = csv.reader(paladin_file, delimiter=',')

            # open output file for immediate writing in same pass as reading paladin report
            with open(output_csv, 'w') as csvfile:
                writer = csv.writer(csvfile, quoting=csv.QUOTE_MINIMAL, lineterminator='\n')
                writer.writerow(['uniprot', 'brenda', 'gene_name', 'organism', 'count', 'abundance'])

                # iterate over incoming report
                for row in paladin_records:
                    matches = tuple(e for e in enzymes if is_match(e, row[1]))

                    # I don't want to calculate matches twice, so incur memory impact of tuple
                    if len(matches) > 0:
                        writer.writerow(row)  # gib whole row plox.
                        for match in matches:  # there could be more than one given wildcards
                            enzyme_counts[match] += int(row[4])  # use count from file
                            if row[2] and not enzyme_names[match]:
                                enzyme_names[match] = row[2]

        files_created.append(output_csv)

        # attempt to look up representative names of enzymes absent in pathway
        # (as we don't have names when no annotation exists)
        missing = [e for e in enzymes if enzyme_counts[e] == 0]
        log(log_file, 'Looking up representative names of missing enzymes...', verbose)
        extra_names = lookup(missing, 'kegg')
        for e in extra_names.keys():
            if len(extra_names[e]) > 0:
                enzyme_names[e] = extra_names[e][0]

        # generators are quick - use more generators
        log_entries = []
        log_entries.append('Present in pathway: [{0}]'.format(', '.join(e for e in enzymes if enzyme_counts[e] > 0)))
        log_entries.append('Missing from pathway: [{0}]'.format(', '.join(missing)))

        completeness = sum(1 for e in enzymes if enzyme_counts[e] > 0) / len(enzymes) * 100
        log_entries.append('Path completion: {0:.1f}%'.format(completeness))

        for item in log_entries:
            if item:
                log(log_file, item, True)

        # example for counts
        if count_csv:
            with open(count_csv, 'wt') as count_file:
                writer = csv.writer(count_file, lineterminator='\n')
                writer.writerow(['ec', 'gene_name', 'count'])
                for e in enzymes:
                    writer.writerow([e, enzyme_names[e], enzyme_counts[e]])

            with open(count_csv[:-4] + '.log', 'wt') as count_log_file:
                count_log_file.writelines([item + '\n' for item in log_entries])

            files_created.append(count_csv)
            files_created.append(count_csv[:-4] + '.log')

        pathway_json = 'pathway_{0}{1}.json'.format(kegg_id, suffix)
        with open(pathway_json, 'wt') as pjfile:
            json.dump(pathway, pjfile, indent=4)
            files_created.append(pathway_json)

        log(log_file, '\nFiles created:', True)
        for f in files_created:
            log(log_file, f, True)

        print('\nDone.')
    else:
        log(log_file, "ERROR: No pathway found with ID: " + kegg_id, True)
        exit(0)

def is_match(known, potential):
    """Detects if potential EC reference matches the known one, accounting for '-' as wildcard"""
    dash = known.find('-')
    if dash > -1:
        return potential.startswith(known[:dash])
    else:
        return known == potential

# actual script entry point
if __name__ == "__main__":
    main()
