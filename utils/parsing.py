# -*- coding: utf-8 -*-
from __future__ import unicode_literals
import re, json, io, pickle
from protherm import SequenceStructure, ExperimentalCondition, ThermoData, Literature, ProThermEntry, ProThermStudy

seq_expected_attributes = ['PROTEIN', 'SOURCE', 'LENGTH', 'MOL-WEIGHT', 'PIR_ID', 'SWISSPROT_ID', 'E.C.NUMBER', 'PMD.NO', 'PDB_wild', 'PDB_mutant', 'MUTATION', 'MUTATED_CHAIN', 'NO_MOLECULE', 'SEC.STR.', 'ASA']
exp_expected_attributes = ['T', 'pH', 'BUFFER_NAME', 'BUFFER_CONC', 'ION_NAME', 'ION_CONC', 'ADDITIVES', 'PROTEIN_CONC', 'MEASURE', 'METHOD']
data_expected_attributes = ['dG_H2O', 'ddG_H2O', 'dG', 'ddG', 'Tm', 'dTm', 'dHvH', 'dHcal', 'm', 'Cm', 'dCp','STATE', 'REVERSIBILITY', 'ACTIVITY', 'ACTIVITY_Km', 'ACTIVITY_Kcat', 'ACTIVITY_Kd']
lit_expected_attributes = ['KEY_WORDS', 'REFERENCE', 'AUTHOR', 'REMARKS', 'RELATED_ENTRIES']

class ComplexEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, (list, dict, str, unicode, int, float, bool, type(None))):
            return json.JSONEncoder.default(self, obj)
        elif hasattr(obj,'toJSON'):
            return obj.toJSON()
        return {'_python_object': pickle.dumps(obj)}


def collect_value(current_value, value):
    """helper function for parse_lines"""
    if value == None or value == "":
        return current_value
    elif current_value == "":
        return value
    elif isinstance(current_value,list):
        return current_value.append(value)
    else:
        return [current_value, value]

def parse_lines(lines, expected_attributes, object):
    # the first 16 characters assign the category or are blank and continuation of prev
    values = dict([(category,"") for category in expected_attributes])
    category=""
    value=""
    for line in lines:
        if len(line[:15].strip()) > 0: # is not blank
            category=line[:15].strip() # set the new category
            category = re.sub('_[0-9]','',category)
            if category not in expected_attributes:
                print "Unexpected category in line"
                print line, category, expected_attributes
                return None
            values[category]=collect_value(values[category],line[15:].strip())
        else: # if category is blank put in previous category
            values[category]=collect_value(values[category],line[15:].strip())
    return object(*[values[category] for category in expected_attributes])

def parse_entry(entry_lines,id_num):
    """
    Parse the entry lines into a ProThermEntry object
    """
    # split the lines into a dictionary with seq_struc:"", experiment:"", data:"", literature:""
    header_dict = {'Sequence and structural information':(seq_expected_attributes, SequenceStructure), 'Experimental condition':(exp_expected_attributes, ExperimentalCondition), 'Thermodynamic data':(data_expected_attributes, ThermoData), 'Literature':(lit_expected_attributes, Literature)}
    header = ""
    section=[]
    objects={}
    for line in entry_lines:
        if '*****' in line.decode('utf-8'):
            if header in header_dict.keys() and section != []:
                obj = parse_lines(section, header_dict[header][0], header_dict[header][1])
                objects[header]=obj
            header = line.strip('\n* ')
            section=[]
        else:
            section.append(line)
    obj = parse_lines(section[:-1], header_dict[header][0], header_dict[header][1]) # last line is //
    objects[header]=obj

    return ProThermEntry(id_num, objects['Sequence and structural information'],objects['Experimental condition'], objects['Thermodynamic data'],objects['Literature'])

def parse_protherm_input_to_entry(protherm_input_file, dump=False):
    """
    Parse the provided protherm input file into ProTherm entry objects
    """
    protherm_entries={}
    entry=[]
    id=0
    with io.open(protherm_input_file, mode="r", encoding="utf-8") as fp:
        for line in fp:
            if line[:3] == "NO.": # started a new entry
                if entry != [] and id != 0:
                    protherm_entries[id]=parse_entry(entry, id)
                id = int(line[4:])
                entry=[line]
            else:
                entry.append(line)
    # grab the last one too!
    if len(entry) > 1:
        protherm_entries[id]=parse_entry(entry, id)


    if dump:
        for id, entry in protherm_entries.items():
            json.dump(entry, open('entries/entry_'+str(id)+'.json','w'), cls=ComplexEncoder, indent=4)

    return protherm_entries

def group_entries_into_studies(protherm_entries, dump=False):
    """
    group protherm entries into studies through related_entries information
    """
    study_entry_map={} # {entry_id:study pairs}
    study_id=0
    studies={}

    for id,entry in protherm_entries.items():
        if id not in study_entry_map.keys():
            study = ProThermStudy(study_id)
            study.add_entry(entry)
            studies[study_id]= study
            study_entry_map[id]=study_id
            for rel_en in entry.get_related_entries():
                if rel_en not in study_entry_map.keys():
                    study.add_entry(protherm_entries[rel_en])
                    study_entry_map[rel_en] = study_id
                else:
                    print "trying to add existing key"
            study_id+=1
    if dump:
        for id,study in studies.items():
            json.dump(study, open('studies/study_'+str(id)+'.json','w'),cls=ComplexEncoder, indent=4)

    return studies

