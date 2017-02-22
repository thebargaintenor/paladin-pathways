#! /usr/bin/env python3
import argparse
import shlex
import plugins.core
import json
import time
import csv
import paladin_postprocess
import kegg
import requests
import dataset
import xmltodict
import os
import sys
'''
pathways.py

DESCRIPTION:
Script runs entire pathways pipeline, filtering a read report from PALADIN
to calculate the mathmatical completeness of a given metabolic pathway.

USAGE:
python3 pathways.py -k 00625
                    -p ./data/report.tsv
                    -o results.csv
                    -c counts.dat
'''


# Plugin connection definition
def pluginConnect(passDefinition):
    passDefinition.name = "pathways"  # Plugin name shown in plugin list (should match filename so user knows what to type eg @@pluginName) 
    passDefinition.description = "elaborate on the pathway output produced by PALADIN Align"  # Plugin description shown in  plugin list
    passDefinition.versionMajor = 1  # Plugin version shown in plugin list
    passDefinition.versionMinor = 0
    passDefinition.versionRevision = 0
    passDefinition.dependencies = []  # If plugin depends on other plugins, specify their plugin name here

    #passDefinition.callbackInit = templateInit # Reference plugin initialization method here (run once at startup).  Not required.
    passDefinition.callbackMain = pathwaysMain # Reference plugin main method here.  Will receive plugin arguments.  Required.



# library imports

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



def log(f, line, verbose=False):
    if verbose:
        plugins.core.sendOutput(line, 'stdout')

    with open(f, 'at') as handle:
        handle.write(line + '\n')


"""
EC_LOOKUP
"""


def lookup(ec_list, database):
    enzyme_names = {}
    if database == 'kegg':
        groups = group(ec_list, 10)
        for g in groups:
            plugins.core.sendOutput(" ".join(['Fetching',
                                              str(len(g)),
                                              'EC references...']), 'stdout')
            url = 'http://rest.kegg.jp/list/' +\
                  '+'.join('ec:' + ec for ec in g if '-' not in ec)
            plugins.core.sendOutput(url, "stdout")

            r = requests.get(url)
            if r.status_code == 200:  # 200 is no error, so proceed happily
                lines = r.text.split('\n')
                for line in lines:
                    if line:
                        record = line.split('\t')
                        enzyme_names[record[0][3:]] = record[1].split('; ')
            else:
                plugins.core.sendOutput(" ".join(["HTTP error",
                                                  r.status_code]), "stderr")
                plugins.core.sendOutput("Download process halting", "stderr")
                return  # don't bother requesting anything else
    return enzyme_names


def group(lst, n):
    for i in range(0, len(lst), n):
        val = lst[i:i + n]
        if len(val) == n or len(val) == len(lst):
            yield val

"""
KEGG
"""

def kegg_get(pathway_id, overwrite=False):
    '''Load KGML from either the local database or KEGG
    (caching the new KGML)'''
    path = os.path.dirname(__file__)
    # I stubbornly insist on continued windows compatibility
    if sys.platform == 'win32':
        path += '\\'
    else:
        path += '/'
    # database is stored in folder with scripts so that changing the working
    # directory doesn't require a new cache
    # (plus it can be shared with multiple users)
    db = dataset.connect('sqlite:///' + path + 'kgmlcache.db')
    kgml_table = db.get_table('kgml')

    if not kgml_table:
        record = download_kgml(pathway_id)
        kgml_table.insert(record, ['id'])
        kgml_table.create_index(['id'], 'pathway_id')
        kgml = record['kgml']
    elif overwrite:
        record = download_kgml(pathway_id)
        kgml_table.upsert(record, ['id'])
        kgml = record['kgml']
    else:
        record = kgml_table.find_one(id=pathway_id)
        if record:
            kgml = record['kgml']
        else:
            record = download_kgml(pathway_id)
            kgml_table.insert(record)
            kgml = record['kgml']

    if kgml:
        return json.loads(kgml)  # return none if no pathway json found


def download_kgml(pathway_id):
    '''
    Download KGML for the pathway and return a record ready for DB storage
    '''
    # but really, get the KGML, build dictionary from it, add data from
    # KEGG dat file (which has different stuff in it, because reasons)
    # I need that data stuffed into my JSON for the visualization, as
    # the regular KGML has no names for the compounds in the chart

    url = 'http://rest.kegg.jp/get/ec' + pathway_id
    plugins.core.sendOutput(": ".join(['Retrieving data from', url]), 'stdout')
    r = requests.kegg_get(url)

    # retrieve DAT version of pathway first (Because.)
    datfile = None
    if r.status_code == 200:  # 200 is no error, so proceed happily
        datfile = r.text
        if not datfile:
            plugins.core.sendOutput(" ".join(['No data for pathway',
                                              pathway_id,
                                              'found on KEGG.']), 'stdout')
    else:
        plugins.core.sendOutput(" ".join(['HTTP Error',
                                          r.status_code]) +
                                "\nDownload process halting", "stderr")
        return  # don't bother requesting anything else

    url += '/kgml'
    plugins.core.sendOutput(" ".join(['Retrieving data from ', url]), "stdout")
    r = requests.kegg_get(url)

    # retrieve KGML of pathway
    kgml = None
    if r.status_code == 200:  # 200 is no error, so proceed happily
        kgml = r.text
        if not kgml:
            plugins.core.sendOutput(" ".join(['No data for pathway',
                                              pathway_id,
                                              'found on KEGG.']), "stdout")
    else:
        plugins.core.sendOutput(" " .join(['HTTP Error',
                                           r.status_code]) +
                                'Download process halting', "stdout")
        return  # don't bother requesting anything else

    # convert KGML to dictionary for JSON output,
    # but add a dictionary mapping KEGG IDs to compound names
    kgml_dict = xmltodict.parse(kgml)

    # extract compounds and add them to dictionary
    if kgml_dict['pathway']:
        kgml_dict['pathway']['compounds'] = extract_compounds(datfile)

    # create database record from amended dictionary
    kgml = json.dumps(kgml_dict, indent=4)
    record = dict(id=pathway_id,
                  name=kgml_dict['pathway']['@title'],
                  kgml=kgml)
    return record


def extract_compounds(data):
    '''
    build a dictionary of compound IDs and names from KEGG pathway DAT file
    '''
    compounds = {}
    in_compound_section = False
    for line in data.split('\n'):
        section = line[:12].rstrip()
        if section == "COMPOUND":
            in_compound_section = True
        elif section and section != "COMPOUND":
            in_compound_section = False

        # extract compound information from record
        if in_compound_section and line:
            cid = line[12:20].rstrip()
            cname = line[20:].rstrip()
            compounds[cid] = cname

    return compounds


"""
MAIN
"""


def main_pathways(arguments):
    # Parse arguments
    argParser = argparse.ArgumentParser(
                        description='PALADIN Pipeline Plugins: pathways',
                        prog='pathways')
    argParser.add_argument(['-k', '--kegg'],
                           help='KEGG ID for pathway',
                           required=True,
                           dest="kegg")
    argParser.add_argument(['-p', '--paladin'],
                           help='path to PALADIN output report',
                           required=True,
                           dest="paladin")
    argParser.add_argument(['-o', '--output'],
                           help='output filename',
                           required=True,
                           dest="output")
    argParser.add_argument("-c",
                           help='suppliment data filename',
                           required=True)
    argParser.add_argument(['-v', '--verbose'],
                           help='verbose mode',
                           action='store_true',
                           dest="verbose")
    arguments = vars(argParser.parse_known_args(arguments)[0])
    kegg_id = arguments["kegg"]
    paladin_report = arguments["paladin"]
    output_csv = arguments["output"]
    count_csv = arguments["c"]
    verbose = arguments["verbose"]
    suffix = time.strftime('_%m%d%y_%H%M')
    log_file = 'pathways_{0}{1}.log'.format(kegg_id, suffix)
    log(log_file, 'Fetching pathway information for KEGG ID ' +
        kegg_id, verbose)

    # retrieve dictionary of pathway information
    pathway = kegg_get(kegg_id)

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
                    if ec not in enzyme_names:
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

            # open output file for immediate writing
            # in same pass as reading paladin report
            with open(output_csv, 'w') as csvfile:
                writer = csv.writer(csvfile,
                                    quoting=csv.QUOTE_MINIMAL,
                                    lineterminator='\n')
                writer.writerow(['uniprot',
                                 'brenda',
                                 'gene_name',
                                 'organism',
                                 'count',
                                 'abundance'])

                # iterate over incoming report
                for row in paladin_records:
                    matches = tuple(e for e in enzymes if is_match(e, row[1]))
                    # I don't want to calculate matches twice,
                    # so incur memory impact of tuple
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
        plugins.core.sendOutput("\nDone.", 'stdout')
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





def pathwaysMain(passArguments):
    argParser = argparse.ArgumentParser(
                        description='PALADIN Pipeline Plugins: pathways',
                        prog='pathways')
    argParser.add_argument('--version',
                           action='version',
                           version='Paladin Pathways 1.0.0')
    argParser.add_argument("modules",
                           nargs="*",
                           help="Paladin pathways step to run")
    args = argParser.parse_known_args(shlex.split(passArguments))
    modules = set(args[0].modules)
    passArguments = args[1]
    if "main" in modules:
        main_pathways(passArguments)
