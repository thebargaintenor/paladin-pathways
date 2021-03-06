#!/usr/bin/env python3

import requests

# http://www.kegg.jp/kegg/xml/docs/
# http://www.kegg.jp/kegg/rest/keggapi.htm
def lookup(ec_list, database):
    enzyme_names = {}
    if database == 'kegg':
        groups = group(ec_list, 10)
        for g in groups:
            print('Fetching', len(g), 'EC references...')
            url = 'http://rest.kegg.jp/list/' + '+'.join('ec:' + ec for ec in g if not '-' in ec)
            print(url)

            r = requests.get(url)
            if r.status_code == 200:  # 200 is no error, so proceed happily
                lines = r.text.split('\n')
                for line in lines:
                    if line:
                        record = line.split('\t')
                        enzyme_names[record[0][3:]] = record[1].split('; ')
            else:
                print('HTTP Error', r.status_code)
                print('Download process halting')
                return  # don't bother requesting anything else
    return enzyme_names

def group(lst, n):
    for i in range(0, len(lst), n):
        val = lst[i:i + n]
        if len(val) == n or len(val) == len(lst):
            yield val

def main():
    print(lookup(['1.97.1.8', '1.14.13.69', '4.5.1.3', '3.8.1.3', '1.14.12.11', '1.2.98.1'], 'kegg'))

# actual script entry point
if __name__ == "__main__":
    main()