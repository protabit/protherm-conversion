# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from optparse import OptionParser
import sys, json, os, re, io
import pickle
from utils.query import pdb_query, uniprot_query, pubmed_bulk_query, search_for_pubmedid
from utils.sequence import valid_protein_sequence, valid_mut_desc, remove_extra_chains, identify_all_potential_shifts
from utils.parsing import parse_protherm_input_to_entry, group_entries_into_studies 
#######################################################################################
#Author: Connie Wang
#Date: 4/24/2018
#Description: Parser for ProTherm Entriesinto ProtaBank Studies
#Usage: ./convert_protherm_to_protabank protherm_data.dat
#######################################################################################

def collect_study_pmids(studies):
    """
    collect all the pmids from the studies to query
     pubmed in bulk TODO add pdb?
    """
    all_pmids = []
    for id,study in studies.items():
        # collect the pmids
        if study.collect_publications() != 1:
            print str(study.id), "does not have 1 pubmed id"
            print study.publications
        all_pmids+=study.pubmed_ids
    literature_dict = pubmed_bulk_query(list(set(all_pmids)))
    if options.dump:
        json.dump(literature_dict, open('database_info/literature_dict.json','w'), indent=4)
    return literature_dict

def collect_study_pdb_ids(studies):
    all_pdbids=[]
    for id,study in studies.items():
        if not study.proteins_collected:
            study.collect_proteins()
        for protein in study.proteins:
            if protein[1] != "" and protein[1] != None:
                all_pdbids.append(protein[1])
            if protein[2] != "" and protein[2] != None:
                all_pdbids.append(protein[2])
    pdb_dict={}
    for id in list(set(all_pdbids)):
        pdb_dict[id]=pdb_query(id)
    if options.dump:
        json.dump(pdb_dict, open('database_info/pdb_dict.json','w'), indent=4)
    return pdb_dict


def collect_study_uniprot_ids(studies):
    all_uniprotids=[]
    for id,study in studies.items():
        if not study.proteins_collected:
            study.collect_proteins()
        for protein in study.proteins:
            if protein[0] != "" and re.search(r'\((.*?)\)',protein[0]) != None:
                uniprot_id = re.search(r'\((.*?)\)',protein[0]).group(1)
                all_uniprotids.append(uniprot_id)
    uniprot_dict={}
    for id in list(set(all_uniprotids)):
        uniprot_dict[id]=uniprot_query(id)
    if options.dump:
        json.dump(uniprot_dict, open('database_info/uniprot_dict.json','w'), indent=4)
    return uniprot_dict



#######################################################################################
#Main
#######################################################################################
if __name__=='__main__':
    use = "Usage: %prog [options] protherm_input_file"
    parser = OptionParser(usage = use)
    parser.add_option("-v", "--verbose", dest="verbose", action="store_true", default=False, help="Set mode to verbose.")
    parser.add_option("-d", "--dump", dest="dump", action="store_true", default=False, help="Dump intermediate steps to restart if there is an error.")
    parser.add_option("-r", "--refetch", dest="refetch", action="store_true", default=False, help="Refetch outside database information")
    parser.add_option("-s", "--load_studies", dest="load_studies", action="store_true", default=False, help="Start from dumped studies steps.")
    parser.add_option("-e", "--load_entries", dest="load_entries", action="store_true", default=False, help="Start from dumped protherm entries step.")
    parser.add_option("-f", "--filename", dest="write", metavar="FILE", help="write output to FILE")
    parser.add_option("-t", "--testing", dest="testing", action="store_true", default=False, help="Test creation do not write or submit")
    (options, args) = parser.parse_args()

    protherm_input_file=args[0]

    protherm_entries={}
    studies={}
    if options.load_studies:
        studies = load_files_from_dump('study','studies')
    elif options.load_entries:
        protherm_entries = load_files_from_dump('entry', 'entries')
        studies = group_entries_into_studies(protherm_entries, options.dump)
    else:
        protherm_entries = parse_protherm_input_to_entry(protherm_input_file, options.dump)
        studies = group_entries_into_studies(protherm_entries, options.dump)


    if os.path.isfile('database_info/literature_dict.json') and not options.refetch:
        literature_dict = json.load(open('database_info/literature_dict.json','r'))
    else:
        literature_dict = collect_study_pmids(studies)

    if os.path.isfile('database_info/pdb_dict.json') and not options.refetch:
        pdb_dict = json.load(open('database_info/pdb_dict.json','r'))
    else:
        pdb_dict = collect_study_pdb_ids(studies)

    if os.path.isfile('database_info/uniprot_dict.json') and not options.refetch:
        uniprot_dict = json.load(open('database_info/uniprot_dict.json','r'))
    else:
        uniprot_dict = collect_study_uniprot_ids(studies)

    submitter=2
    success=open('submitted_studies.txt', 'w')
    fail=open('not_submitted_studies.txt','w') 
    for id,study in studies.items():
        if not study.has_all_entries():
            fail.write("Study "+str(id)+" not submitted due to missing entries\n")
        else:
            if options.verbose: print "Creating ProtaBank study from ProTherm Study", id
            study_dict = study.create_protabank_study(literature_dict, pdb_dict, uniprot_dict, submitter)
            with open('protabank_studies/study'+str(id)+'.json', 'w') as fp:
                json.dump(study_dict, fp, indent=4)
                if len(study.flags) == 0:
                    success.write(str(id)+"\n")
                else:
                    fail.write("Study "+str(id)+" not submitted due to flags "+";".join(study.flags)+"\n")
    success.close()
    fail.close()
