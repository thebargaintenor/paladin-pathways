#!/usr/bin/env python3

# paladin_postprocess.py
#
# 2016-04-11 tdm - script created
# 2016-06-05 tdm - broke out postprocess job into new method for use by other scripts
#
# ----
#
# PALADIN POSTPROCESS
#   Script parses a specified report file from PALADIN as input and
#   creates output CSV file with each line holding a UniProt ID, an
#   enzyme commission (EC) reference, the name of the enzyme from
#   TrEMBL entry, the organism providing the representative sequence,
#   a read count, and the abundance of the read in the sample
#
# ----

import sys
import re
import time
import csv

# -- GLOBAL VARIABLES --

usageMsg = '''

Usage: paladin_postprocess -i report -o output
    Script parses a specified report file from PALADIN as input and
    creates output CSV file with each line holding a UniProt ID, an
    enzyme commission (EC) reference, the name of the enzyme from
    TrEMBL entry, the organism providing the representative sequence,
    a read count, and the abundance of the read in the sample

    --input
    -i report   - input report file path
    --output
    -o output   - output file path
    --verbose
    - v         - print status messages

    EXAMPLE: python paladin_postprocess.py -i paladin_report_sample.tsv -o output.csv

'''

# -- FUNCTIONS --

# usage()
# If we have missing arguments or the first argument is -h
# print some "usage" information and quit.
def usage():
    if len( sys.argv ) < 4 or sys.argv[ 1 ] == "-h" or sys.argv[ 1 ] == "--help":
        print( usageMsg )
        exit( 0 )

def dump_records(records, f):
    '''Dumps records as lines into the given file'''
    f.writelines(records)

def format_value(input):
    if ',' in input:
        return '"' + input.replace('"', '""') + '"'
    else:
        return input


def main():
    usage()

    # some starting defaults
    input_path = ''
    output_path = ''
    verbose = False

    # check args
    for i in range(0, len(sys.argv)):
        if sys.argv[i] == '--input' or sys.argv[i] == '-i':
            input_path = sys.argv[i + 1]
            if not input_path.endswith('.tsv'):
                input_path += '.tsv'
        elif sys.argv[i] == '--output' or sys.argv[i] == '-o':
            output_path = sys.argv[i + 1]
            if not output_path.endswith('.csv'):
                output_path += '.csv'
        elif sys.argv[i] == '--verbose' or sys.argv[i] == '-v':
            verbose = True

    # run job
    run(input_path, output_path, verbose)

def run(input_path, output_path, verbose):
    '''Run postprocess job (now accessible from other scripts)'''
    ec_pattern = re.compile("(?<=\\(EC )([0-9]+\\.[0-9\\-]+\\.[0-9\\-]+\\.[0-9\\-]+)(?=\\))")

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
                        print('{0} lines processed.'.format(lines_processed))

                lines = infile.readlines(100000)

            # flush output buffer
            if len(records) > 0:
                dump_records(records, outfile)
                del records[:]
                if verbose:
                    print('{0} lines processed.'.format(lines_processed))

    # stop time
    end = time.clock()
    print('Post-process operation finished in {0:.2f} seconds.'.format(end - start))

if __name__ == "__main__":
    main()