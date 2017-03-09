#! /usr/bin/env python3
import argparse
import shlex
import plugins.core
import plugins.taxonomy
import json
import time
import csv
import requests
import dataset
import xmltodict
import re
import tempfile
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import subprocess
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
    passDefinition.name = "pathways_plugin"  # Plugin name shown in plugin list (should match filename so user knows what to type eg @@pluginName) 
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


def kegg_get(pathway_id, overwrite=False, kegg_db=None):
    '''Load KGML from either the local database or KEGG
    (caching the new KGML)'''
    #path = os.path.dirname(__file__)
    if kegg_db:
        path = kegg_db
    else:
        path = tempfile.mktemp() + 'kgmlcache.db'
    # I stubbornly insist on continued windows compatibility
    #if sys.platform == 'win32':
    #    path += '\\'
    #else:
    #    path += '/'
    # database is stored in folder with scripts so that changing the working
    # directory doesn't require a new cache
    # (plus it can be shared with multiple users)
    db = dataset.connect('sqlite:///' + path)
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


def download_kgml(pathway_id, return_kgml=False):
    '''
    Download KGML for the pathway and return a record ready for DB storage
    '''
    # but really, get the KGML, build dictionary from it, add data from
    # KEGG dat file (which has different stuff in it, because reasons)
    # I need that data stuffed into my JSON for the visualization, as
    # the regular KGML has no names for the compounds in the chart

    url = 'http://rest.kegg.jp/get/ec' + pathway_id
    plugins.core.sendOutput(": ".join(['Retrieving data from', url]), 'stdout')
    r = requests.get(url)

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
    r = requests.get(url)

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
    if return_kgml:
        return kgml
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
HEATMAP
"""


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


def render(copymat, genome_names, en_disp,
           outname='hist.png', colormap=cm.Reds,
           figsize=(19.2, 10.8)):
    plt.figure(figsize=figsize)
    plt.imshow(1*(copymat), interpolation='none', cmap=colormap, aspect='auto')
    plt.xticks(np.arange(len(genome_names)), genome_names, rotation=20)
    plt.yticks(range(len(en_disp)), en_disp)
    cb = plt.colorbar()
    cb.set_label('copy number')
    plt.tight_layout()
    plt.savefig(outname)


"""
POSTPROCESS
"""


def dump_records(records, f):
    '''Dumps records as lines into the given file'''
    f.writelines(records)


def format_value(input):
    if ',' in input:
        return '"' + input.replace('"', '""') + '"'
    else:
        return input


def run_postprocess(input_path, output_path, verbose):
    '''Run postprocess job (now accessible from other scripts)'''
    ec_pattern = re.compile("(?<=\\(EC )([0-9]+\\.[0-9\\-]" +
                            "+\\.[0-9\\-]+\\.[0-9\\-]+)(?=\\))")

    lines_processed = 0
    # output buffer (to decrease disk write frequency)
    records = []

    # get start of long op
    start = time.clock()

    with open(output_path, 'wt') as outfile:
        with open(input_path, 'rt') as infile:
            # iterate over tsv file
            outfile.write('uniprot,brenda,gene_name,organism,count,abundance\n')
            lines = infile.readlines(100000)
            while lines:
                for line in lines:
                    lines_processed += 1
                    if line:
                        # split fields in record
                        fields = line.split('\t')
                        if len(fields) > 5:
                            # look for EC reference
                            m = ec_pattern.search(fields[5])
                            if m:
                                # trim EC annotation from description (already in its own column)
                                ecidx = fields[5].find('(EC')
                                if ecidx > -1:
                                    desc = fields[5][:ecidx].rstrip()
                                else:
                                    desc = fields[5]

                                # add new CSV record for this item
                                records.append('{0},{1},{2},{3},{4},{5}\n'.format(
                                    fields[3], m.group(1), format_value(desc), format_value(fields[4]), fields[0], fields[1]))
                                if len(records) > 1000:
                                    # dump buffer to file
                                    dump_records(records, outfile)
                                    del records[:]

                    if verbose and lines_processed % 10000 == 0:
                        plugins.core.sendOutput('{0} lines processed.'.format(lines_processed), "stdout")

                lines = infile.readlines(100000)

            # flush output buffer
            if len(records) > 0:
                dump_records(records, outfile)
                del records[:]
                if verbose:
                    plugins.core.sendOutput('{0} lines processed.'.format(lines_processed), "stdout")

    # stop time
    end = time.clock()
    plugins.core.sendOutput('Post-process operation finished in {0:.2f} seconds.'.format(end - start), "stdout")

"""
MAIN
"""


def postprocess(arguments):
    argParser = argparse.ArgumentParser(
                        description='PALADIN Pipeline Plugins: pathways',
                        prog='pathways')
    argParser.add_argument('--verbose', "-v",
                           help='verbose mode',
                           action='store_true')
    argParser.add_argument("--i-paladin-tsv",
                           help="path to tsv output of paladin",
                           required=True)
    argParser.add_argument('--output', "-o",
                           help='output path',
                           required=True)
    arguments = vars(argParser.parse_known_args(arguments)[0])
    input_path = arguments["i_paladin_tsv"]
    input_filename = input_path.split("/")[-1].split(".")[0]
    output_path = "".join([arguments["output"],
                           input_filename,
                           "csv"])
    verbose = arguments["verbose"]
    run_postprocess(input_path, output_path, verbose)


def main_pathways(arguments):
    # Parse arguments
    argParser = argparse.ArgumentParser(
                        description='PALADIN Pipeline Plugins: pathways',
                        prog='pathways')
    argParser.add_argument('--kegg', "-k",
                           help='KEGG ID for pathway',
                           required=True)
    argParser.add_argument('--paladin', "-p",
                           help='path to PALADIN output report',
                           required=True)
    argParser.add_argument('--output', "-o",
                           help='output filename',
                           required=True)
    argParser.add_argument("-c",
                           help='suppliment data filename',
                           required=True)
    argParser.add_argument('--verbose', "-v",
                           help='verbose mode',
                           action='store_true')
    argParser.add_argument('--kegg_db',
                           help='kegg_database',
                           required=False)
    
    arguments = vars(argParser.parse_known_args(arguments)[0])
    kegg_id = arguments["kegg"]
    paladin_report = arguments["paladin"]
    output_csv = arguments["output"] + '/pathways.csv'
    count_csv = arguments["c"]
    verbose = arguments["verbose"]
    suffix = time.strftime('_%m%d%y_%H%M')
    if "kegg_db" in arguments:
        kegg_db = arguments["kegg_db"]
    else:
        kegg_db = None
    log_file = 'pathways_{0}{1}.log'.format(kegg_id, suffix)
    log(log_file, 'Fetching pathway information for KEGG ID ' +
        kegg_id, verbose)

    # retrieve dictionary of pathway information
    pathway = kegg_get(kegg_id, kegg_db=kegg_db)

    if pathway:
        pathway = pathway['pathway']
        files_created = []

        # make name for filtered output
        filtered_paladin_csv = paladin_report[:-4] + suffix + '.csv'

        log(log_file, 'Running TSV post-process...', verbose)
        run_postprocess(paladin_report, filtered_paladin_csv, verbose)

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


def get_taxa(blast_results="blastout",
             db="/usr/local/share/paladin/uniref90.fasta"):
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


def taxa_callback(passArguments):
    argParser = argparse.ArgumentParser(
                        description='PALADIN Pipeline Plugins: pathways',
                        prog='pathways')
    argParser.add_argument('--paladin', "-p",
                           help='path to PALADIN output report',
                           required=True)
    argParser.add_argument('--output', '-o',
                           help='output directory',
                           required=True)
    argParser.add_argument("--i-enzyme_code",
                           help='enzyme code',
                           required=True)
    argParser.add_argument("--i-reads",
                           help='raw reads',
                           required=True)
    arguments = vars(argParser.parse_known_args(passArguments)[0])
    enzyme_code = arguments["i_enzyme_code"]
    pathways_out = arguments["output"]
    paladin_out = arguments["paladin"]
    reads = arguments["i_reads"]
    uids = set(get_uniprot_id(enzyme_code, pathways_out))
    paladin_entries = plugins.core.PaladinEntry.getEntries(paladin_out, 0)
    filtered_entries = {}
    for entry in paladin_entries.items():
        key = entry[1].id
        if key in uids:
            filtered_entries[entry[0]] = entry[1]
    taxa = plugins.taxonomy.getSpeciesLookup(filtered_entries)
    print(taxa)   
    """
    loh = uniprot_coords(uid[0], paladin_out)
    los = find_seq(loh, reads)
    mkquery(los)
    blaster()
    taxa_dict = get_taxa()
    output_taxa(taxa_dict, "taxonomy")
    """

def heatmap(passArguments):
    argParser = argparse.ArgumentParser(
                        description='PALADIN Pipeline Plugins: pathways',
                        prog='pathways')
    argParser.add_argument("--i-heatmap-folder",
                           help="path to folder containing files to heatmap",
                           required=True)
    argParser.add_argument('--output', "-o",
                           help='output path',
                           required=True)
    arguments = vars(argParser.parse_known_args(passArguments)[0])
    infiles = glob.glob(arguments["i_heatmap_folder"] +
                        "/*.csv")
    copymat, genome_names, en_disp = parse_files(infiles)
    render(copymat, genome_names, en_disp,
           outname=arguments["output"] + "/heatmap.png")


def visbuilder(passArguments):
# load build configuration manifest
    argParser = argparse.ArgumentParser(
                        description='PALADIN Pipeline Plugins: pathways',
                        prog='pathways')
    argParser.add_argument("--kegg_db",
                           help="kegg_db file",
                           required=True)
    arguments = vars(argParser.parse_known_args(passArguments)[0])
    config_path = arguments["kegg_db"]
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
                    #ec = entry['graphics']['@name']
                    # because a reaction could now have several enzymes
                    ecs = [ec[3:] for ec in entry['@name'].split(' ')]
                    for ec in ecs:
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
    argParser.add_argument("-l",
                           action="store_true")
    args = argParser.parse_known_args(shlex.split(passArguments))
    modules = set(args[0].modules)
    passArguments = args[1]
    if args[0].l:
        mods = ["main",
                "postprocess",
                "heatmap",
                "taxa_callback",
                "visbuilder"]
        plugins.core.sendOutput("\n".join(mods), "stdout")
    if "main" in modules:
        main_pathways(passArguments)
    if "postprocess" in modules:
        postprocess(passArguments)
    if "heatmap" in modules:
        heatmap(passArguments)
    if "taxa_callback" in modules:
        taxa_callback(passArguments)
    if "visbuilder" in modules:
        visbuilder(passArguments)
