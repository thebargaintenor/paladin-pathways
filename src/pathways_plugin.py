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
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import subprocess
from Bio.KEGG.KGML import KGML_parser
from Bio.Graphics.KGML_vis import KGMLCanvas
from PyPDF2 import PdfFileWriter, PdfFileReader, PdfFileMerger
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
    if kegg_db:
        path = kegg_db
    else:
        path = tempfile.mktemp() + 'kgmlcache.db'
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
                           help='Pathways output folder',
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
                        description='PALADIN Pipeline Plugins: pathways main',
                        prog='pathways')
    argParser.add_argument('--kegg', "-k",
                           help='KEGG ID for pathway',
                           required=True)
    argParser.add_argument('--paladin', "-p",
                           help='path to PALADIN output report',
                           required=True)
    argParser.add_argument('--output', "-o",
                           help='Pathways output folder',
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
    output_csv = arguments["output"] + '/pathways.tsv'
    count_csv = arguments["output"] + '/pathways_counts.tsv'
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
                                    lineterminator='\n',
                                    delimiter='\t')
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
                writer = csv.writer(count_file, lineterminator='\n', delimiter="\t")
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


def taxa_callback(passArguments):
    argParser = argparse.ArgumentParser(
                        description='PALADIN Pipeline Plugins: pathways taxa_callback',
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
    arguments = vars(argParser.parse_known_args(passArguments)[0])
    enzyme_code = arguments["i_enzyme_code"]
    pathways_out = arguments["output"]
    paladin_out = arguments["paladin"]
    uids = set(get_uniprot_id(enzyme_code, pathways_out))
    paladin_entries = plugins.core.PaladinEntry.getEntries(paladin_out, 0)
    filtered_entries = {}
    for entry in paladin_entries.items():
        key = entry[1].id
        if key in uids:
            filtered_entries[entry[0]] = entry[1]
    taxa = plugins.taxonomy.getSpeciesLookup(filtered_entries)
    
    
def heatmap(passArguments):
    argParser = argparse.ArgumentParser(
                        description='PALADIN Pipeline Plugins: pathways heatmap',
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


def visualize_counts(passArguments):
    argParser = argparse.ArgumentParser(
                        description='PALADIN Pipeline Plugins: pathways visualize counts',
                        prog='pathways')
    argParser.add_argument('--output', "-o",
                           help='output path',
                           required=True)
    argParser.add_argument('--kegg', "-k",
                           help='kegg id',
                           required=True)
    arguments = vars(argParser.parse_known_args(passArguments)[0])
    count_file = glob.glob(arguments["output"] + "/*count*")[0]
    pathways_counts = {}
    counts = []
    with open(count_file) as cf:
        cf = cf.readlines()
    for line in cf[1:]:
        words = line.split(",")
        pathways_counts[words[0]] = float(words[2])
        counts.append(float(words[2]))
    max_count = np.max(counts)
    kgml = download_kgml(arguments["kegg"], True)
    pathway = KGML_parser.read(kgml)
    cmap = cm.inferno
    for element in pathway.entries.items():
        key = element[0]
        e_object = element[1]
        if e_object.name[3:] in pathways_counts:
            element_count = pathways_counts[e_object.name[3:]]
            if element_count == 0:
                normalized_count = 1
            else:
                normalized_count = np.log10(element_count)/np.log10(max_count)
            acc = 0
            for graphic in e_object.graphics:
                new_rgba_color = cmap(normalized_count)
                rgb_color = tuple([int(255*val) for val in new_rgba_color[0:3]])
                rgb_string = '#' +\
                             ''.join([hex(val)[2:] for val in rgb_color]).upper()
                e_object.graphics[acc].bgcolor = rgb_string
                acc += 1
        else:
            acc = 0
            for graphic in e_object.graphics:
                e_object.graphics[acc].bgcolor = "#b2aba7".upper()
                acc += 1
        pathway.entries[key] = e_object
    canvas = KGMLCanvas(pathway, fontsize=12, import_imagemap=True,
                        margins=(0.1, 0.02), label_orthologs=False,
                        label_maps=False, label_compounds=False,
                        label_reaction_entries=True)
    canvas.draw("fab_map_new_colours.pdf")
    fig = plt.figure(figsize=(1, 5))
    ax1 = fig.gca()
    cb1 = mpl.colorbar.ColorbarBase(ax1, cmap=cmap,
                                    orientation='vertical')
    cb1.set_label("log(count)/max(log(count))")
    plt.tight_layout()
    plt.savefig(arguments["output"] + "/colorbar.pdf")
    watermark = PdfFileReader(open(arguments["output"] + "/colorbar.pdf", "rb"))
    output_file = PdfFileWriter()
    input_file = PdfFileReader(open(arguments["output"] + "/fab_map_new_colours.pdf", "rb"))
    page_count = input_file.getNumPages()
    for page_number in range(page_count):
        input_page = input_file.getPage(page_number)
        input_page.mergePage(watermark.getPage(0))
        output_file.addPage(input_page)
    with open(arguments["output"] + "/document-output.pdf", "wb") as outputStream:
        output_file.write(outputStream)


def visualize_taxa(passArguments):
    argParser = argparse.ArgumentParser(
                        description='PALADIN Pipeline Plugins: pathways',
                        prog='pathways')
    argParser.add_argument('--output', "-o",
                           help='output path',
                           required=True)
    argParser.add_argument('--kegg', "-k",
                           help='kegg id',
                           required=True)
    arguments = vars(argParser.parse_known_args(passArguments)[0])
    count_file = glob.glob(arguments["output"] + "/pathways.tsv")[0]
    pathways_counts = {}
    with open(count_file) as cf:
        cf = cf.readlines()
    for line in cf:
        words = line.split(",")
        pathways_counts[words[1]] = words[-3]
    kgml = download_kgml(arguments["kegg"], True)
    pathway = KGML_parser.read(kgml)
    cmap = cm.inferno
    locs = []
    tax_list = []
    for element in pathway.entries.items():
        key = element[0]
        e_object = element[1]
        if e_object.name[3:] in pathways_counts:
            tax_list.append(pathways_counts[e_object.name[3:]])
    tax_list = list(set(tax_list))
    taxa = {}.fromkeys(pathways_counts.values())
    inds = np.linspace(0, 1, len(tax_list))
    acc = 0
    for key in tax_list:
        taxa[key] = inds[acc]
        acc += 1
    for element in pathway.entries.items():
        key = element[0]
        e_object = element[1]
        if e_object.name[3:] in pathways_counts:
            locs.append(taxa[pathways_counts[e_object.name[3:]]])
            normalized_count = taxa[pathways_counts[e_object.name[3:]]]
            acc = 0
            for graphic in e_object.graphics:
                new_rgba_color = cmap(normalized_count)
                rgb_color = tuple([int(255*val) for val in new_rgba_color[0:3]])
                rgb_string = '#' +\
                             ''.join([hex(val)[2:] for val in rgb_color]).upper()
                e_object.graphics[acc].bgcolor = rgb_string
                acc += 1
        else:
            acc = 0
            for graphic in e_object.graphics:
                e_object.graphics[acc].bgcolor = "#b2aba7".upper()
                acc += 1
        pathway.entries[key] = e_object
    canvas = KGMLCanvas(pathway, fontsize=12, import_imagemap=True,
                        margins=(0.1, 0.02), label_orthologs=False,
                        label_maps=False, label_compounds=False,
                        label_reaction_entries=True)
    canvas.draw(arguments["output"] + "/taxa_map.pdf")
    fig = plt.figure()
    ax = fig.gca()
    for value in tax_list:
        item = (value, taxa[value])
        plt.plot(0, 0, "s", label=item[0], color=cmap(item[1]))
    plt.xlim((-100, -99))
    ax.set_axis_off()
    plt.legend()
    plt.tight_layout()
    plt.savefig(arguments["output"] + "legend.pdf", bbox_inches='tight')
    merger = PdfFileMerger()
    for pdf in [arguments["output"] + "/taxa_map.pdf", arguments["output"] + "/legend.pdf"]:
        merger.append(open(pdf, 'rb'))
    with open(arguments["output"] + '/final_taxa_map.pdf', 'wb') as fout:
        merger.write(fout)


def barplot_vis(passArguments):
    argParser = argparse.ArgumentParser(
                        description='PALADIN Pipeline Plugins: pathways',
                        prog='pathways')
    argParser.add_argument('--output', "-o",
                           help='output path',
                           required=True)
    arguments = vars(argParser.parse_known_args(passArguments)[0])
    pathways_outfile = arguments["output"] + "/pathways.tsv"
    with open(pathways_outfile) as pout:
        pout = pout.readlines()
    brendas = []
    gene_names = []
    organisms = []
    counts = []
    for line in pout[1:]:
        try:
            uniprot, brenda, gene_name, organism, count, abundance = line.rstrip().split(",")
            brendas.append(brenda)
            gene_names.append(gene_name)
            organisms.append(organism)
            counts.append(count)
        except:
            pass
    brenda_bins = {}
    organism_bins = {}
    for i in range(len(brendas)):
        brenda = brendas[i]
        organism = organisms[i]
        if brenda not in brenda_bins:
            brenda_bins[brenda] = [i]
        else:
            brenda_bins[brenda].append(i)
        if organism not in organism_bins:
            organism_bins[organism] = [i]
        else:
            organism_bins[organism].append(i)
    # Bin by brenda
    cmap = cm.autumn
    xloc = np.linspace(0, 1, len(organism_bins))
    org_colors = {}
    acc = 0
    for key in organism_bins.keys():
        org_colors[key] = cmap(xloc[acc])
        acc += 1
    xloc = np.linspace(0, 1, len(brenda_bins))
    brenda_colors = {}
    acc = 0
    for key in brenda_bins.keys():
        brenda_colors[key] = cmap(xloc[acc])
        acc += 1
    acc = 0
    plt.figure(figsize=(30, len(brenda_bins)/2), dpi=300)
    labels = []
    for item in brenda_bins.items():
        brenda = item[0]
        indicies = item[1]
        bcount = []
        borg = []
        borg_org = []
        for index in indicies:
            borg.append(org_colors[organisms[index]])
            borg_org.append(organisms[index])
            bcount.append(int(counts[index]))
        total_count = np.sum(bcount)
        bcount = np.asarray(bcount)
        ratio = bcount / total_count
        acci = 1
        for i in range(len(ratio)):
            val = ratio[i]
            ratio[i] = acci
            acci = acci - val
        labels.append(brenda)
        pos = ratio * 0 + acc
        plt.barh(pos, ratio, align="center", height=.94, color=borg)
        if len(ratio) > 1:
            for i in range(len(ratio) - 1):
                if ratio[i] - ratio[i + 1] > len(borg_org[i])/300:
                    plt.text(ratio[i + 1], acc, " " + borg_org[i],
                             fontsize=10, color="k", fontweight='bold',
                             verticalalignment="center")
        if ratio[-1] > len(borg_org[-1])/300:
            plt.text(0, acc, " " + borg_org[-1], fontsize=10,
                     color="k", fontweight='bold',
                     verticalalignment="center")
        acc += 1
    plt.xlabel("ratio of counts: max counts per EC code")
    plt.yticks(np.arange(len(labels)), labels)
    for item in org_colors.items():
        plt.plot(-1, -1, "s", color=item[1], label=item[0])
    plt.xlim((0, 1))
    plt.ylim((-.5, len(brenda_bins) - .5))
    plt.savefig(arguments["output"] + "/brenda_bar.png", bbox_inches="tight", transparent=True)
    plt.figure(figsize=(15, len(organism_bins)/5), dpi=300)
    acc = 0
    labels = []
    for item in organism_bins.items():
        organism = item[0]
        indicies = item[1]
        ocount = []
        oorg = []
        oorg_org = []
        for index in indicies:
            oorg.append(brenda_colors[brendas[index]])
            oorg_org.append(brendas[index])
            ocount.append(int(counts[index]))
        total_count = np.sum(ocount)
        ocount = np.asarray(ocount)
        ratio = ocount / total_count
        acci = 1
        for i in range(len(ratio)):
            val = ratio[i]
            ratio[i] = acci
            acci = acci - val
        labels.append(organism)
        pos = 0 * ratio + acc
        plt.barh(pos, ratio, align="center", color=oorg)
        if len(ratio) > 1:
            for i in range(len(ratio) - 1):
                if ratio[i] - ratio[i + 1] > len(oorg_org[i])/300*2:
                    plt.text(ratio[i + 1], acc, " " + oorg_org[i],
                             fontsize=10, color="k", fontweight='bold',
                             verticalalignment="center")
        if ratio[-1] > len(oorg_org[-1])/300*2:
            plt.text(0, acc, " " + oorg_org[-1], fontsize=10,
                     color="k", fontweight='bold',
                     verticalalignment="center")
        acc += 1
    plt.yticks(np.arange(len(labels)), labels)
    for item in brenda_colors.items():
        plt.plot(-1, -1, "s", color=item[1], label=item[0])
    plt.xlim((0, 1))
    plt.ylim((-.5, len(organism_bins) - .5))
    plt.savefig(arguments["output"] + "/organism_bar.png", bbox_inches='tight', transparent=True)
    

def piechart_vis(passArguments):
    argParser = argparse.ArgumentParser(
                        description='PALADIN Pipeline Plugins: pathways',
                        prog='pathways')
    argParser.add_argument('--output', "-o",
                           help='output path',
                           required=True)
    argParser.add_argument("--piechart_brenda",
                           help="brenda id's to plot piecharts for list seperated by commas (no spaces!)",
                           required=False,
                           default="")
    argParser.add_argument("--piechart_organism",
                           help="organisms to plot piecharts for, list seperated by commas, (no spaces!)",
                           required=False,
                           default="")
    arguments = vars(argParser.parse_known_args(passArguments)[0])
    pathways_outfile = arguments["output"] + "/pathways.tsv"
    with open(pathways_outfile) as pout:
        pout = pout.readlines()
    brendas = []
    gene_names = []
    organisms = []
    counts = []
    for line in pout[1:]:
        try:
            uniprot, brenda, gene_name, organism, count, abundance = line.rstrip().split(",")
            brendas.append(brenda)
            gene_names.append(gene_name)
            organisms.append(organism)
            counts.append(count)
        except:
            pass
    brenda_bins = {}
    organism_bins = {}
    for i in range(len(brendas)):
        brenda = brendas[i]
        organism = organisms[i]
        if brenda not in brenda_bins:
            brenda_bins[brenda] = [i]
        else:
            brenda_bins[brenda].append(i)
        if organism not in organism_bins:
            organism_bins[organism] = [i]
        else:
            organism_bins[organism].append(i)
    # Bin by brenda
    cmap = cm.autumn
    xloc = np.linspace(0, 1, len(organism_bins))
    org_colors = {}
    acc = 0
    for key in organism_bins.keys():
        org_colors[key] = cmap(xloc[acc])
        acc += 1
    xloc = np.linspace(0, 1, len(brenda_bins))
    brenda_colors = {}
    acc = 0
    for key in brenda_bins.keys():
        brenda_colors[key] = cmap(xloc[acc])
        acc += 1
    for brenda in arguments["piechart_brenda"].split(","):
        if brenda != "":
            plt.figure()
            indicies = brenda_bins[brenda]
            bcount = []
            borg = []
            borg_org = []
            for index in indicies:
                borg.append(org_colors[organisms[index]])
                borg_org.append(organisms[index])
                bcount.append(int(counts[index]))
            bcount = np.asarray(bcount)
            ratio = bcount / np.sum(bcount)
            explode = np.fmin(np.ones_like(ratio) * .2, -np.log10(ratio)/10)
            explode = np.logspace(-3, -.7, len(ratio), base=2)
            plt.pie(bcount,  labels=borg_org, autopct='%1.1f%%', explode=explode,
                    shadow=True, startangle=90, colors=borg)
            plt.gca().axis('equal')
            plt.savefig(arguments["output"] + "/pie_" + brenda + ".png",
                        bbox_inches="tight", transparent=True)
    for organism in arguments["piechart_organism"].split(","):
        if organism != "":
            plt.figure()
            indicies = organism_bins[organism]
            ocount = []
            oorg = []
            oorg_org = []
            for index in indicies:
                oorg.append(brenda_colors[brendas[index]])
                oorg_org.append(brendas[index])
                ocount.append(int(counts[index]))
            ocount = np.asarray(ocount)
            ratio = ocount / np.sum(ocount)
            explode = np.fmin(np.ones_like(ratio) * .2, -np.log10(ratio)/10)
            explode = np.logspace(-3, -.7, len(ratio), base=2)
            plt.pie(ocount, labels=oorg_org, autopct='%1.1f%%', explode=explode,
                    shadow=True, startangle=90, colors=oorg)
            plt.gca().axis('equal')
            plt.savefig(arguments["output"] + "/pie_" +
                        ''.join(e for e in organism if e.isalnum())+".png",
                        bbox_inches="tight", transparent=True)


def pathwaysMain(passArguments):
    argParser = argparse.ArgumentParser(
                        description='PALADIN Pipeline Plugins: pathways',
                        prog='pathways')
    argParser.add_argument('--version',
                           action='version',
                           version='Paladin Pathways 1.1.0')
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
                "visualize_counts",
                "visualize_taxa",
                "barplot_vis",
                "piechart_vis"]
        plugins.core.sendOutput("\n".join(mods), "stdout")
    if "main" in modules:
        main_pathways(passArguments)
    if "postprocess" in modules:
        postprocess(passArguments)
    if "heatmap" in modules:
        heatmap(passArguments)
    if "taxa_callback" in modules:
        taxa_callback(passArguments)
    if "visualize_counts" in modules:
        visualize_counts(passArguments)
    if "visualize_taxa" in modules:
        visualize_taxa(passArguments)
    if "barplot_vis" in modules:
        barplot_vis(passArguments)
    if "piechart_vis" in modules:
        piechart_vis(passArguments)
