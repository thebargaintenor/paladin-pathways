#!/usr/bin/env python3

# kegg.py
#
# 2016-06-05 tdm - script created
# 2016-08-24 tdm - SQLite DB now in script directory only
#
# ----
#
# KEGG
#   This script is ONLY a library for interacting with the KGML cache and downloading
#   path information.  If you want a default run behavior, write it.  Ideally, keep
#   that in your own script and just import this one.
#
# ----
#
# EXAMPLE: pathway = json.loads(kegg.get(00625))

# -- IMPORTS --
# for working with SQLite
import dataset
# this is much nicer with the Requests library
# https://requests.readthedocs.org/en/master/
# pip install requests
# (needed for nicer HTTP communication)
import requests

# for loading between python dictionaries and XML
import xmltodict

import os
import sys
# for loading between python dictionaries and JSON strings
import json

def get(pathway_id, overwrite=False):
    '''Load KGML from either the local database or KEGG (caching the new KGML)'''
    
    path = os.path.dirname(__file__)
    # I stubbornly insist on continued windows compatibility
    if sys.platform == 'win32':
        path += '\\'
    else:
        path += '/'
    
    # database is stored in folder with scripts so that changing the working
    # directory doesn't require a new cache (plus it can be shared with multiple users)
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
        return json.loads(kgml) # return none if no pathway json found

def download_kgml(pathway_id):
    '''Download KGML for the pathway and return a record ready for DB storage'''
    # but really, get the KGML, build dictionary from it, add data from
    # KEGG dat file (which has different stuff in it, because reasons)
    # I need that data stuffed into my JSON for the visualization, as
    # the regular KGML has no names for the compounds in the chart

    url = 'http://rest.kegg.jp/get/ec' + pathway_id
    print('Retrieving data from ', url)
    r = requests.get(url)

    # retrieve DAT version of pathway first (Because.)
    datfile = None
    if r.status_code == 200:  # 200 is no error, so proceed happily
        datfile = r.text
        if not datfile:
            print('No data for pathway', pathway_id, 'found on KEGG.')
    else:
        print('HTTP Error', r.status_code)
        print('Download process halting')
        return  # don't bother requesting anything else

    url += '/kgml'
    print('Retrieving data from ', url)
    r = requests.get(url)

    # retrieve KGML of pathway
    kgml = None
    if r.status_code == 200:  # 200 is no error, so proceed happily
        kgml = r.text
        if not kgml:
            print('No data for pathway', pathway_id, 'found on KEGG.')
    else:
        print('HTTP Error', r.status_code)
        print('Download process halting')
        return  # don't bother requesting anything else

    # convert KGML to dictionary for JSON output,
    # but add a dictionary mapping KEGG IDs to compound names
    kgml_dict = xmltodict.parse(kgml)

    # extract compounds and add them to dictionary
    if kgml_dict['pathway']:
        kgml_dict['pathway']['compounds'] = extract_compounds(datfile)

    # create database record from amended dictionary
    kgml = json.dumps(kgml_dict, indent=4)
    record = dict(id=pathway_id, name=kgml_dict['pathway']['@title'], kgml=kgml)
    return record

def extract_compounds(data):
    '''build a dictionary of compound IDs and names from KEGG pathway DAT file'''
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