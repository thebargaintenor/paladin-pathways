#!/usr/bin/env python3

# visbuilder.py
#
# 2016-04-29 tdm - script created
#
# ----
#
# VISUALIZATION BUILDER
#   Script loads a user-created build configuration file and assembles the
#   CSV dataset needed for pathway visualizer based on data aggregated from
#   files specified by the configuration.
#
# ----

# -- GLOBAL VARIABLES --

import sys
import os
import json
import csv

usageMsg = '''

Usage: visbuilder.py CONFIG
    Script loads a user-created build configuration file and assembles the
    CSV dataset needed for pathway visualizer based on data aggregated from
    files specified by the configuration.

    CONFIG      - path to JSON-formatted build configuration file

    EXAMPLE: python visbuilder.py configuration.json

'''


# -- FUNCTIONS --

def usage():
    """If we have missing arguments or the first argument is -h, print some "usage" information and quit."""
    if len(sys.argv) < 2 or sys.argv[1] == '-h' or sys.argv[1] == '--help':
        print(usageMsg)
        exit(0)


def is_match(known, potential):
    """Detects if potential EC reference matches the known one, accounting for - as wildcard"""
    dash = known.find('-')
    if dash > -1:
        return potential.startswith(known[:dash])
    else:
        return known == potential


def main():
    usage()

    # load build configuration manifest
    config_path = sys.argv[1]
    with open(config_path, 'rt') as config_file:
        config = json.load(config_file)  # this will remain in scope... WORRY NOT.
        config_file.close()

        # load pathway data into dictionary using path in config
        with open(config['pathway'], 'rt') as path_file:
            dictionary = json.load(path_file)
            pathway = dictionary['pathway']
            if 'name' in config:
                pathway['@description'] = config['name']
            path_file.close()

            print('Detected output directory:', config['output']['path'])

            # this is the master list of EC references for iteration purposes
            enzyme_ids = []

            enzyme_names = {}
            enzyme_counts = {}

            sets = []
            enzymes_by_set = {}

            print('Building manifest...')

            # initialize data structures
            for entry in pathway['entry']:
                if entry['@type'] == 'enzyme':
                    ec = entry['graphics']['@name']
                    if not ec in enzyme_names:
                        enzyme_ids.append(ec)
                        enzyme_names[ec] = []
                        enzyme_counts[ec] = 0

            # *** BY ARBITRARY CONVENTION, all output files must go in the same folder.
            #     This is for simplicity of packaging, so deal with it.

            # As of python 3, makedirs can just handle itself if the target path
            # already exists.  If this is a problem, uncomment the check here:
            # eh... I'll do it anyway just for being verbose
            print('Checking output path...')

            path = config['output']['path']
            if not os.path.exists(path):
                print('Creating output directory...')
                os.makedirs(path, exist_ok=True)

            # commence actual aggregation
            print('Starting data aggregation...')
            for entry in config['data']:
                with open(entry['path'], 'rt') as datasource_file:
                    print('Aggregating ' + entry['path'] + '...')
                    reader = csv.reader(datasource_file)
                    lines_processed = 0
                    for record in reader:
                        lines_processed += 1
                        if lines_processed % 1000 == 0:
                            print(lines_processed, 'lines processed')

                        if record[0] != 'name':
                            ec = record[1]
                            # check all enzymes for match (line may match multiple)
                            for enzyme in enzyme_ids:
                                if is_match(enzyme, ec):
                                    if record[2]:
                                        if len(enzyme_names[enzyme]) < 5:
                                            # prune EC from description (since it's already stored)
                                            desc_end = record[2].find('(EC')
                                            enzyme_names[enzyme].append(record[2][:desc_end].rstrip())
                                            # enzyme_names[enzyme].append(record[2])
                                    #print(record)
                                    enzyme_counts[enzyme] += int(record[4])

                    print(lines_processed, 'lines processed')

                    # compute mathematical (euclidean) completeness of pathway
                    completeness = sum(1 for e in enzyme_ids if enzyme_counts[e] > 0) / len(enzyme_ids)

                    # add names to master dictionary
                    enzymes_by_set[entry['name']] = enzyme_names.copy()

                    # done with set, add output row
                    csv_set = [entry['name'], completeness]
                    for enzyme in enzyme_ids:
                        csv_set.append(enzyme_counts[enzyme])
                        enzyme_counts[enzyme] = 0
                        enzyme_names[enzyme] = [] # clear list for next pass
                    sets.append(csv_set)

            # dump aggregated sets to file
            print('Writing aggregated data ')
            with open(path + config['output']['dataset'], 'wt') as dataset_file:
                writer = csv.writer(dataset_file, quoting=csv.QUOTE_MINIMAL, lineterminator='\n')
                # header row
                headers = []
                headers.append('name')
                headers.append('completeness')
                headers.extend(enzyme_ids)
                writer.writerow(headers)
                for row in sets:
                    writer.writerow(row)

            # build manifest
            print('Building manifest...')
            manifest = {}
            manifest['pathway'] = pathway['@number']
            manifest['src'] = config['output']['pathway']
            manifest['dataset'] = config['output']['dataset']

            # dump manifest to JSON file
            with open(path + config['output']['manifest'], 'wt') as manifest_file:
                manifest_file.write(json.dumps(manifest, indent=4))

            # finally, update pathway data and dump to file in output dir
            print('Writing pathway definition...')
            pathway['enzymes'] = enzymes_by_set # break down by series
            with open(path + manifest['src'], 'wt') as src_file:
                json.dump(dictionary, src_file, indent=4)
    print('Done.')

if __name__ == "__main__":
    main()