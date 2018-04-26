# -*- coding: utf-8 -*-
from __future__ import unicode_literals
import re, json
from collections import defaultdict
from query import search_for_pubmedid
from sequence import remove_extra_chains, valid_protein_sequence, valid_mut_desc, ACCEPTEDAA, identify_all_potential_shifts

unit_converter={'kJ/mol':'kJ/mol',
        'kJ/mole':'kJ/mol',
        'kJ/mol/K':'kJ/mol·K',
        'kJ/K mol':'kJ/mol·K',
        'kJ/molK':'kJ/mol·K',
        'kJ/mol K':'kJ/mol·K',
        'kJ/mol/M':'kJ/mol·M',
        'kJ/mol/k':'kJ/mol·K',
        'kJ/K/mol':'kJ/mol·K',
        'kJ/Kmol':'kJ/mol·K',
        'kJ/molK':'kJ/mol·K',
        'kJ/m/K':'kJ/mol·K',
        'kJ/M/mol':'kJ/mol·M',
        'kJmol/M':'kJ/mol·M',
        'kJ/K':'kJ/K',
        'J/K/g':'J/K/g',
        'J/g':'J/g',
        'kcal':'kcal',
        'kcal/mol':'kcal/mol',
        'kcal/mole':'kcal/mol',
        'kal/mol':'kcal/mol',
        'kcla/molK':'kcal/mol·K',
        'kcal/mol/K':'kcal/mol·K',
        'kcal/mol deg':'kcal/mol·K',
        'kcal/mol/M':'kcal/mol·M',
        'kcal/mole/K':'kcal/mol·K',
        'kcal/mol/deg':'kcal/mol·K',
        'Kcal/mol/deg':'kcal/mol·K',
        'kcal/K mol':'kcal/mol·K',
        'kcal/mol M':'kcal/mol·M',
        'kcal/mo/Ml':'kcal/mol·M',
        'Kcal/K mol':'kcal/mol·K',
        'cal/K/mol':'cal/mol·K',
        '1/s':'1/s',
        '1/min':'1/min',
        'mmol/mg/min':'mmol/mg/min',
        '%':'%',
        '(relative)':'(relative)',
        'micro molar':'uM',
        'micro m':'uM',
        'microM':'uM',
        'micro M':'uM',
        '(x 10**6 /M)':'uM', #incorrectly described in paper
        '(micro Molar)':'uM',
        'mM':'mM',
        'nM':'nM',
        'M':'M',
        '1/M':'1/M',
        'mg/mL':'mg/mL',
        'units/mg':'units/mg',
        '(units/mg)':'units/mg',
        'U/mg':'U/mg',
        '1/mg':'1/mg',
        'K':'K',
        'C':'°C',
        'k':'K',

}
#######################################################################################
# Classes
#######################################################################################
class SequenceStructure(object):
    """
    Represents the Sequence and Structural information part of the ProTherm entry

    Has fields PROTEIN, SOURCE, LENGTH, MOL-WEIGHT, PIR_ID, SWISSPROT_ID, E.C.NUMBER, PMD.NO PDB_wild, PDB_mutant, MUTATION, MUTATED_CHAIN, NO_MOLECULE, SEC.STR., ASA

    Information will be parsed into PROTEIN, and SEQUENCE objects in ProtaBank
    """
    def __init__(self, protein, source, length, mol_weight, pir_id, swissprot_id, ec_no, pmd_no, pdb_wild, pdb_mutant, mutation, mutated_chain, no_mol, sec_str, asa):
        self.protein = protein
        self.source = source
        self.length= length
        self.mol_weight = mol_weight
        self.pir_id = pir_id
        self.swissprot_id = swissprot_id
        self.ec_no = ec_no
        self.pmd_no = pmd_no
        self.pdb_wild = pdb_wild.split(',')[0]
        self.pdb_mutant = pdb_mutant.split(',')[0]
        self.mutation = mutation
        self.mutated_chain = mutated_chain
        self.no_mol = no_mol
        self.sec_str = sec_str
        self.asa = asa
        self.pdb_sequence=""

    def __str__(self):
        return self.__dict__()

    def toJSON(self):
        return self.__dict__

    def create_protein_dict(self):
        if self.pdb_mutant != "":
            pdb_dict = pdb_query(self.pdb_mutant)
        else:
            pdb_dict = pdb_query(self.pdb_wild)
        if self.pdb_sequence=="":
            self.pdb_sequence = pdb_dict['id_sequence']
        return pdb_dict

    def identify_PIR_ID(pir_id):
        """Use the PIR ID to identify the protein.
        PIR no longer updated, integrated with uniprot
        https://pir.georgetown.edu/cgi-bin/pir_psd_get.pl?id=PSBOA"""
        return ""

    def fetch_by_swissprot_id(swissprot_id):
        """ expected format is NAME (UNIPROT ID)

        first try using uniprot ID
        """
        uniprot_id = re.search(r'\((.*?)\)',swissprot_id).group(1)
        return uniprot_query(uniprot_id)

    def fetch_sequence_by_id(self, pdb_dict, uniprot_dict, proteins=[]):
        if self.pdb_sequence=="":
            prot_dict={}
            if self.pdb_wild!="":
                prot_dict = pdb_dict[self.pdb_wild]
            elif self.pdb_mutant != "":
                prot_dict = pdb_dict[self.pdb_mutant]
            elif self.swissprot_id !="" and re.search(r'\((.*?)\)',self.swissprot_id) is not None:
                uniprot_id = re.search(r'\((.*?)\)',self.swissprot_id).group(1)
                prot_dict = uniprot_dict[uniprot_id]
            else: # check other proteins for options
                for protein in proteins:
                    swissprot_id, pdb = protein
                    prot_dict={}
                    if pdb!="":
                        prot_dict = pdb_dict[pdb]
                        break
                    elif swissprot_id !="" and re.search(r'\((.*?)\)',swissprot_id) is not None:
                        uniprot_id = re.search(r'\((.*?)\)',swissprot_id).group(1)
                        prot_dict= uniprot_dict[uniprot_id]
                        break
                    else:
                        continue
            if 'id_sequence' not in prot_dict.keys() or prot_dict['id_sequence'] == None:
                # print 'id_sequence missing from protein dict', prot_dict
                self.pdb_sequence=""
            else:
                self.pdb_sequence = remove_extra_chains(prot_dict['id_sequence'])
        return self.pdb_sequence

    def mutation_in_protabank_format(self):
        """
        Convert the protherm mutation format to the protabank formation

        protherm format is described as
        'Details about the mutation: residue in wild type, residue number and residue in mutant protein (e.g. A 123 G). In the case of insertion or deletion mutations, all inserted/deleted residues appear along with the preceding residue that corresponds to the residue number at which the mutation occurs (for example, A 378 AVL refers that the two amino acid residues "V" and "L" are inserted at 378th position where the residue "A" is present in the wild type, and it is represented conversely for the deletion, viz, AVL 378 A).

        The mutation residue number is given as in the reference. In some cases, the residue number may not match with PIR and PDB sequences. In such cases, the residue numbers corresponding to PIR and PDB sequences are also given in the parenthesis.'
        """

        simple_mutation = r'(?P<wt>[A-Z]) *(?P<pos>[0-9]+) *(?P<mut>[A-Z]) *$'
        insertion_mutation = r'(?P<wt>[A-Z]) *(?P<pos>[0-9]+) *(?P<mut>[A-Z]+) *$'
        deletion_mutation = r'(?P<wt>[A-Z]+) *(?P<pos>[0-9]+) *(?P<mut>[A-Z]) *$'
        delete_all_mutation = r'(?P<wt>[A-Z]+) *(?P<pos>[0-9]+) *del *$'
        mismatch_with_pdb = r'(?P<wt>[A-Z]) *(?P<pos>[0-9]+) *(?P<mut>[A-Z]) *\(PDB: (?P<pdb_wt>[A-Z]) *(?P<pdb_pos>[0-9]+) *(?P<pdb_mut>[A-Z]).*\) *$' # NOTE i'm ignoring the PIR sequence number because will pull sequence from the PDB
        cdr_mutation = r'(?P<wt>[A-Z]) *(?P<pos>[0-9]+)(?P<letter>[a-zA-Z]) *(?P<mut>[A-Z]) *$'
        translated_mutants=[]
        if self.mutation.strip() == 'wild' or self.mutation.strip() =='wild*' or self.mutation.strip() == 'wild**':
            return "WT"

        mutants = self.mutation.split(',')
        for mutant in mutants: # loop through the mutants
            if re.search(simple_mutation, mutant): # check if its a simple mutation
                simple = re.search(simple_mutation, mutant)
                if simple.group('wt') in ACCEPTEDAA and simple.group('mut') in ACCEPTEDAA:
                     translated_mutants.append(simple.group('wt')+simple.group('pos')+simple.group('mut'))
                else:
                    print "Matched simple mutation, but unexpected residue", mutant
                    return None
            elif re.search(insertion_mutation, mutant):
                insertion = re.search(insertion_mutation, mutant)
                if insertion.group('wt') in ACCEPTEDAA and all(res in ACCEPTEDAA for res in insertion.group('mut')):
                    translated_mutants.append(insertion.group('wt')+insertion.group('pos') + '(' + insertion.group('mut')+')')
                else:
                    print "Matched insertion mutation, but unexpected residues", mutant
                    return None
            elif re.search(deletion_mutation, mutant):
                deletion = re.search(deletion_mutation,mutant)
                if all(res in ACCEPTEDAA for res in deletion.group('wt')) and  deletion.group('mut') in ACCEPTEDAA:
                    for res in deletion.group('wt')[:-1]:
                        translated_mutants.append(res+deletion.group('pos')+'()')
                    translated_mutants.append(deletion.group('wt')[-1]+deletion.group('pos')+deletion.group('mut'))
                else:
                    print "Matched deletion  mutation, but unexpected residues", mutant
                    return None
            elif re.search(delete_all_mutation, mutant):
                deletion = re.search(delete_all_mutation,mutant)
                if all(res in ACCEPTEDAA for res in deletion.group('wt')):
                    for res in deletion.group('wt'):
                        translated_mutants.append(res+deletion.group('pos')+'()')
                print "Matched deletion all mutation, but unexpected residues", mutant
                return None
            elif re.search(cdr_mutation, mutant):
                cdr_num =  re.search(cdr_mutation, mutant)
                print 'cdr!'
                letter_index=dict([(letter,i) for i,letter in enumerate("abcdef")])
                if cdr_num.group('wt') in ACCEPTEDAA and cdr_num.group('mut') in ACCEPTEDAA and cdr_num.group('letter').lower() in letter_index.keys():
                    pos = int(cdr_num.group('pos')) + letter_index[cdr_num.group('letter').lower()]
                    translated_mutants.append(cdr_num.group('wt')+str(pos)+cdr_num.group('mut'))
                else:
                    print "Matched cdr mutation, but unexpected residue", mutant
                    return None
            elif re.search(mismatch_with_pdb, mutant): # automatically uses the PDB/PIR value
                pdb_mismatch = re.search(mismatch_with_pdb, mutant)
                if pdb_mismatch.group('pdb_wt') in ACCEPTEDAA and pdb_mismatch.group('pdb_mut') in ACCEPTEDAA:
                     translated_mutants.append(pdb_mismatch.group('pdb_wt')+pdb_mismatch.group('pdb_pos')+pdb_mismatch.group('pdb_mut'))
                else:
                    print "Matched pdb mismatch mutation, but unexpected residue", mutant
                    return None
            else:
                print "Did not match any RE", mutant
                return None
        return "+".join(translated_mutants)

class ExperimentalCondition(object):
    """
    Represents the Experimental Condition  part of the ProTherm entry

    Has fields temp, pH, BUFFER_NAME, BUFFER_CONC, ION_NAME[], ION_CONC[], ADDITIVES, PROTEIN_CONC, MEASURE, METHOD

    Information will be parsed into ASSAY objects in ProtaBank

    """

    def __init__(self, temp, pH, buffer_name, buffer_conc, ion_name, ion_conc, additives, prot_conc, measure, method):
        self.temp = temp
        self.pH = pH
        self.buffer_name=buffer_name
        self.buffer_conc = buffer_conc
        self.ion_name=ion_name
        self.ion_conc=ion_conc
        self.additives=additives
        self.prot_conc = prot_conc
        self.measure=measure
        self.method=method

    def __str__(self):
        return str(self.__dict__)

    def toJSON(self):
        return self.__dict__

    def create_assay_dict(self, orig_property):
        """
        Create experimental assays for the all the listed data"""
        property_dict = {
            'dG_H2O':{'category':'Stability', 'property':'ΔG', 'units':'kcal/mol'},
            'ddG_H2O':{'category':'Stability', 'property':'ΔG', 'units':'kcal/mol'},
            'dG_H20':{'category':'Stability', 'property':'ΔG', 'units':'kcal/mol'},
            'ddG_H20':{'category':'Stability', 'property':'ΔG', 'units':'kcal/mol'},
            'dG':{'category':'Stability', 'property':'ΔG', 'units':'kcal/mol'},
            'ddG':{'category':'Stability', 'property':'ΔG', 'units':'kcal/mol'},
            'Tm':{'category':'Stability', 'property':'Tm', 'units':'°C'},
            'dTm':{'category':'Stability', 'property':'Tm', 'units':'°C'},
            'dHvH':{'category':'Stability', 'property':'ΔH', 'units':'kcal/mol'},
            'dHcal':{'category':'Stability', 'property':'ΔH', 'units':'kcal/mol'},
            'm':{'category':'Stability', 'property':'m', 'units':'kcal/mol·M'},
            'Cm':{'category':'Stability', 'property':'Cm', 'units':'M'},
            'dCp':{'category':'Stability', 'property':'ΔCp', 'units':'kcal/mol·K'},
            'activity':{'category':'Activity', 'property':'Relative Activity', 'units':"%"},
            'activity_km':{'category':'Activity', 'property':'Km', 'units':"mM"},
            'activity_kcat':{'category':'Activity', 'property':'kcat', 'units':"1/s"},
            'activity_kd':{'category':'Activity', 'property':'Kd', 'units':"mM"},
            }# returns category, technique, property,unit for each prothem category

        technique_dict={
            'Fl':'Fluorescence Spectroscopy',
            'Fluorescence':'Fluorescence Spectroscopy',
            'CD':'Circular Dichroism (CD)',
            'CD(far-UV)':'Circular Dichroism (CD)',
            'CD(near-UV)':'Circular Dichroism (CD)',
            'CD (far-UV)':'Circular Dichroism (CD)',
            'CD (near-UV)':'Circular Dichroism (CD)',
            'far-UV CD':'Circular Dichroism (CD)',
            'near-UV CD':'Circular Dichroism (CD)',
            'DSC':'Differential Scanning Calorimetry (DSC)',
            'DSMC':'Differential Scanning Calorimetry (DSC)', # guessing typo
            'SEC':'Size Exclusion Chromatography',
            'Abs':'Absorbance',
            'Absorption':'Absorbance',
            'Absorbance':'Absorbance',
            'Refraction':'Refraction',
            'EPR':'Electron Paramagnetic Resonance (EPR)',
            'ESR':'Electron Spin Resonance',
            'NMR':'NMR',
            'HPLC':'Chromatography (HPLC, TLC)',
            'FTIR':'Fourier-transform Infrared Spectroscopy (FTIR)',
            'UV':'UV Spectroscopy',
            'UV spectroscopy':'UV Spectroscopy',
            'SAXS':'Small-angle X-ray Scattering (SAXS)',
            'Anisotropy':'Fluorescence Polarization/Anisotropy',
            'ANS':'Fluorescence Polarization/Anisotropy',
            'Fluorescence (ANS)':'Fluorescence Polarization/Anisotropy',
            'ANS binding':'Fluorescence Polarization/Anisotropy',
            'optical':'Optical Density',
            'Optical':'Optical Density',
            'Activity':'Activity',
            'activity':'Activity',
            'Isothermal denaturation':'Isothermal Titration Calorimetry (ITC)',
            'Hydrogen exchange':'Hydrogen Exchange',
            'NMR Hydrogen exchange':'Hydrogen Exchange',
            'NMR amide hydrogen exchange':'Hydrogen Exchange',
            'CD + Fluorescence':'Circular Dichroism (CD); Fluorescence Spectroscopy',
            'Fluorescence, CD':'Circular Dichroism (CD); Fluorescence Spectroscopy',
            'CD, Fluorescence':'Circular Dichroism (CD); Fluorescence Spectroscopy',
            'Light scattering':'Dynamic Light Scattering',
            'Light-scattering':'Dynamic Light Scattering',
            'Fluorescence (Trp)':'Tryptophan Fluorescence',
            'Gel electrophoresis':'Gel Electrophoresis',
        }
        denaturation_dict={
            'Thermal':'Thermal Denaturation',
            'Urea':'Urea Denaturation',
            'GdnHCl':'Guanidinium Denaturation',
        }

        property = orig_property.replace("H20", "H2O")
        new_assay_dict = {
            'name':property,
            'source':'Exp',
            'category':'',
            'units':'',
            'technique':'',
            'property':'',
            'buffers':property_str(self.buffer_name) + ": "+ property_str(self.buffer_conc),
            'details': "Additives "+property_str(self.additives),
            'ionic': property_str(self.ion_name) + ": " + property_str(self.ion_conc)
        }

        if property in property_dict.keys():
            new_assay_dict['category']=property_dict[property]['category']
            new_assay_dict['units']=property_dict[property]['units']
            new_assay_dict['property']=property_dict[property]['property']
        if self.measure in technique_dict.keys():
            new_assay_dict['technique']=technique_dict[self.measure]
        else:
            new_assay_dict['technique']=self.measure
        if self.method in denaturation_dict.keys() and new_assay_dict['category'] == 'Stability' :
            new_assay_dict['technique']+=";"+denaturation_dict[self.method]
        if self.temp != "":
            try:
                new_assay_dict['temp'] = str(float(self.temp)) + " C"
            except:
                new_assay_dict['temp'] = self.temp
        if self.pH !="":
            try:
                new_assay_dict['pH'] = str(float(self.pH))
            except:
                new_assay_dict['pH'] = self.pH
        if self.prot_conc !="":
            new_assay_dict['prot_conc'] = self.prot_conc

        # specify initial and final for dG
        if property  in ['dG', 'dG_H2O', 'dG_H20']:
            new_assay_dict['initial_state'] = "folded"
            new_assay_dict['final_state'] = "unfolded"
        return new_assay_dict

class ThermoData(object):
    """Represents the Thermodynamic data part of the ProTherm entry

    Has fields dG_H20, ddG_H2O, dG, ddG, Tm, dTm, dHvH, dHcal, m, Cm, dCp, STATE, REVERSIBILITY, ACTIVITY, ACTIVITY_Km, ACTIVITY_Kcat, AC  ACTIVITY_Kd

    Information will be parsed into DATUM and ASSAY objects in ProtaBank
    """

    def __init__(self, dG_H20, ddG_H2O, dG, ddG, Tm, dTm, dHvH, dHcal, m, Cm, dCp, state, reversibility, activity, activity_km, activity_kcat, activity_kd):
        self.dG_H20 = dG_H20
        self.ddG_H2O = ddG_H2O
        self.dG=dG
        self.ddG=ddG
        self.Tm = Tm
        self.dTm = dTm
        self.dHvH = dHvH
        self.dHcal = dHcal
        self.m =m
        self.Cm = Cm
        self.dCp = dCp
        self.state = state
        self.reversibility = reversibility
        self.activity = activity
        self.activity_km = activity_km
        self.activity_kcat = activity_kcat
        self.activity_kd = activity_kd

    def __str__(self):
        return str(self.__dict__)

    def toJSON(self):
        return self.__dict__

    def non_null_results(self):
        """returns a dictionary of property name and values that are not empty"""
        non_null={}
        for property, value in self.__dict__.items():
            if isinstance(value,list):
                print 'Property has multiple values', property, value
            elif property in ['reversibility', 'state']:
                # skip these for now
                pass
            elif value.strip() != "":
                non_null[property]=value
        return non_null
    # def create_datum_dict(self,property, sequence, assay_id):
    #     """expects a property amoung the ThermoData fields, i.e. dG
    #     passes in a sequence_dict which has the full sequence in sequence and potentially a mut_desc
    #
    #     optionally passes in an assay_dict if it already exists
    #     """
    #
    #     datum_dict = {'mut_desc':sequence_dict['mut_desc'], 'sequence':sequence_dict['sequence'], 'result':getattr(self,property)}
    #     datum_dict['units'] = assay_dict['units']
    #     datum_dict['assay'] = assay_id
    #     return datum_dict

class Literature(object):
    """Represents the Literature data part of the ProTherm entry

    Has fields REFERENCE, PMID, AUTHOR, RELATED_ENTRIES
    (ignoring keywords and remarks...)

    Information will be parsed into PUBLICATION objects in ProtaBank
    """
    def __init__(self, key_words, reference, author, remarks, related_entries):
        self.key_words=""
        if key_words is not None and isinstance(key_words,list):
            self.key_words = ",".join([word for word in key_words if word is not None])
        elif key_words is not None:
            self.key_words = key_words
        self.reference = reference
        self.author = author
        if isinstance(author,list):
            self.author=" ".join(author)
        self.remarks=""
        if remarks is not None and isinstance(remarks,list):
            self.remarks = ",".join([word for word in remarks if word is not None])
        elif remarks is not None:
            self.remarks = remarks
        self.related_entries = []
        if isinstance(related_entries,list):
            self.related_entries = related_entries
        else:
            for entry in related_entries.strip(', ').split(','):
                if entry.strip() != "":
                    try:
                        self.related_entries.append(int(entry))
                    except:
                        print "Could not identify related entry", entry
                        print self.reference, related_entries

    def toJSON(self):
        return self.__dict__

    def __str__(self):
        return str(self.__dict__)

    def __eq__(self,other):
        return (self.reference == other.reference and self.author == other.author)

    def __hash__(self):
        return hash("Literature(%s,%s)"%(self.reference,self.author))

    def get_pmid(self):
        if 'PMID:' in self.reference and len(self.reference) > (self.reference.find('PMID:') + 5):
            return self.reference.split('PMID:')[1].strip()
        else:
            search = re.search(r'^([A-Z ]+)([0-9]+)[., ]*([0-9]+)-([0-9]+).*\(([0-9]+)\).*', self.reference)
            if search:
                lib_dict = {'journal':search.group(1), 'volume':search.group(2), 'pages':search.group(3)+"-"+search.group(4), 'year':search.group(5)}
                lib_dict['authors'] = self.author
                return search_for_pubmedid(lib_dict)
            else:
                print "NO RE MATCH FOUND!", self.reference
            return None

    def generate_dict(self):
        search = re.search(r'^([A-Z ]+)([0-9]+)[., ]*([0-9]+)-([0-9]+).*\(([0-9]+)\).*', self.reference)
        if search:
            pub_dict = {'journal':search.group(1), 'volume':search.group(2), 'pages':search.group(3)+"-"+search.group(4), 'year':search.group(5)}
            pub_dict['authors'] = self.author
            return pub_dict
        else:
            return None

class ProThermEntry(object):
    """ Represents a ProTherm Entry comprisied of
    SequenceStructure Info, ExperimentalCondition, ThermoData, Literature

    Information will be parsed into STUDY objects in ProtaBank

    """
    def __init__(self, id, sequence_structure, experimental_condition, data, literature):
        self.id=id
        self.sequence_structure=sequence_structure
        self.experimental_condition=experimental_condition
        self.data=data
        self.literature=literature

    def __str__(self):
        return "Entry No." + str(self.number) + str(self.sequence_structure.__dict__) + str(self.experimental_condition.__dict__)+str(self.data.__dict__)+str(self.literature.__dict__)

    def toJSON(self):
        return self.__dict__

    def get_related_entries(self):
        if self.literature is not None:
            return self.literature.related_entries
        else:
            print str(self)
            return None

class ProThermStudy(object):
    """
    Organizes ProTherm Entries that are explicitly related entries into a STUDY
    which corresponds to protabank studies
    """
    def __init__(self, id):
        self.id = id
        self.entries = []

        # keep track of processing
        self.publications_collected = False
        self.proteins_collected = False
        self.assay_data_collected = False
        # protabank objects
        self.assays = {}
        self.sequences = {}
        self.data = {}
        self.libraries = {} # libraries organized by the starting sequence
        self.proteins = set()
        self.publications=set()
        self.pubmed_ids=[]
        self.extra_details="" # append to description/abstract
        self.related_entry_ids=set()

        self.flags=[] # put flags in here to help identify issues
    def __str__(self):
        return str(self.__dict__)

    def toJSON(self):
        return self.__dict__

    def add_entry(self, entry):
        self.entries.append(entry)

    def has_all_entries(self):
        """Checks that all entries listed in any entries related entries are assigned to this study"""
        for entry in self.entries:
            self.related_entry_ids.update(entry.get_related_entries())
            self.related_entry_ids.add(entry.id)
        if set(self.related_entry_ids) == set([entry.id for entry in self.entries]):
            return True
        else:
            print 'Study', id, 'does not have all entries', set(self.related_entry_ids), set([entry.id for entry in self.entries])
            return False

    def collect_publications(self):
        """
        Collect all the literature results from the entries and return the set of publications"""
        for entry in self.entries:
            self.publications.add(entry.literature)

        for publication in self.publications:
            pmid= publication.get_pmid()
            if pmid not in self.pubmed_ids and pmid is not None and pmid != "":
                self.pubmed_ids.append(pmid)
                self.extra_details+=publication.remarks + "\n" + publication.key_words + "\n"
        self.publications_collected = True
        return len(self.pubmed_ids)

    def collect_proteins(self):
        for entry in self.entries:
            if entry.sequence_structure.pdb_wild != "" or entry.sequence_structure.swissprot_id !="":
                self.proteins.add((entry.sequence_structure.swissprot_id,entry.sequence_structure.pdb_wild))
            if entry.sequence_structure.pdb_wild != "" or entry.sequence_structure.pdb_mutant !="":
                self.proteins.add((entry.sequence_structure.swissprot_id,entry.sequence_structure.pdb_mutant))
        self.proteins_collected = True

    def collect_assay_datum_sequence(self,pdb_dict, uniprot_dict):
        """
        Assays are created from entry experimental condition and thermo data
        Collect assays from entry.experimental and save a unique set"""
        next_assay_id=1
        datum_id=1
        library_id=1
        for entry in self.entries:
            sequence = entry.sequence_structure.fetch_sequence_by_id(pdb_dict, uniprot_dict, self.proteins)
            mut_desc = entry.sequence_structure.mutation_in_protabank_format()
            if not valid_protein_sequence(sequence):
                self.flags.append('Not Valid Protein Sequence: ' + str(sequence))
            valid, new_mut_desc = valid_mut_desc(mut_desc, sequence)
            if not valid:
                self.flags.append('Not Valid mut_desc: ' + str(new_mut_desc) + " " + str(sequence))
            library = 0
            if sequence in self.libraries.keys():
                library =self.libraries[sequence]
            else:
                self.libraries[sequence] = library_id
                library = library_id
                library_id +=1
            for property, value in entry.data.non_null_results().items():
                assay = entry.experimental_condition.create_assay_dict(property)
                parsed_value = parse_result(value, assay)
                assay_id = next_assay_id
                if assay not in self.assays.values():
                    self.assays[assay_id] = assay
                    next_assay_id+=1
                else:
                    for id, a_dict in self.assays.items():
                        if a_dict == assay:
                            assay_id = id
                            break
                datum = {'mut_desc': new_mut_desc, 'assay':assay_id, 'result':parsed_value, 'library':library, 'id':datum_id}
                self.data[datum_id] = datum
                datum_id+=1
        self.assay_data_collected = True

    # def find_index_shift(self, starting_sequence, mut_desc_set, allowed_shifts):
    #     for chain in starting_sequence.split("/"):
    #         found_shift=False
    #         proposed_shift = 0
    #         for mut_desc in mut_desc_set:
    #             chain_id = chain[:2] if ':' in chain else "" #chainid and :
    #             chain_sequence=chain[(chain.find(':')+1):]
    #             new_shift = identify_potential_shift(mut_desc, chain_sequence, allowed_shifts)
    #             if not found_shift and new_shift is None: # nothing found reject chain
    #                 break
    #             elif not found_shift and new_shift: # set the proposed shift
    #                 proposed_shift=new_shift
    #                 found_shift=True
    #             elif found_shift and (new_shift is None or proposed_shift != new_shift): #not matching reject
    #                 found_shift = False
    #                 break
    #         if found_shift:
    #             return proposed_shift
    #     return None
    def find_index_shift(self, starting_sequence, mut_desc_set, allowed_shifts):
        for chain in starting_sequence.split("/"):
            chain_id = chain[:2] if ':' in chain else "" #chainid and :
            chain_sequence=chain[(chain.find(':')+1):]
            potential_shifts=range(-1*len(chain_sequence)+1, len(chain_sequence))
            for mut_desc in mut_desc_set:
                potential_shifts = identify_all_potential_shifts(mut_desc, chain_sequence, potential_shifts)
            if potential_shifts == None or len(potential_shifts)==0:
                break
            if len(potential_shifts)>0:
                print potential_shifts
            if len(potential_shifts)==1:
                return (potential_shifts[0], chain_id) # need to return chain_id!
        return (None,None)

    def give_assay_descriptive_names(self):
        assay_by_name=defaultdict(list)
        for assay in self.assays.values():
            assay_by_name[assay['name']].append(assay)

        for name, assays in assay_by_name.items():
            if len(assays) > 1:
                conditions = set(assays[0].items())
                for assay in assays:
                    common_conditions = conditions.intersection(assay.items())
                for assay in assays:
                    differences = list(set(assay.items())-common_conditions)
                    assay['name']+=" "+", ".join([diff[0]+":"+diff[1] for diff in differences if diff[0] !="technique"])
                    #ensure that assay names do not exceed 100char

    def shorten_assay_names(self):
        for assay in self.assays.values():
            assay['name'] = assay['name'][:100]

    def assays_have_same_mutants(self, aid1, aid2):
        mutdesc1=[]
        mutdesc2=[]
        for datum in self.data.values():
            if datum['assay'] == aid1 and datum['mut_desc'] !='WT' :
                mutdesc1.append(datum['mut_desc'])
            if datum['assay'] == aid2 and datum['mut_desc'] !='WT':
                mutdesc2.append(datum['mut_desc'])
        return mutdesc1 == mutdesc2

    def organize_derived_assays(self):
        derivedAssays={}
        derived_names={'dTm':'Tm', 'ddG':'dG', 'ddG_H20':'dG_H20', 'ddG_H2O':'dG_H2O'}
        for did, assay in self.assays.items():
            for quantity in derived_names.keys():
                if assay['name'].find(quantity) == 0:
                    derivedAssays[did] = (quantity,assay['name'])

        for did, (quantity, name) in derivedAssays.items():
            for eid, eassay in self.assays.items():
                if eid != did and eassay['name'].find(derived_names[quantity]) == 0 and eassay['name'].replace(derived_names[quantity], quantity) == name and self.assays_have_same_mutants(did, eid):
                    derivedAssays[did] = eid
                    break
        for key,value in derivedAssays.items():
            if isinstance(value,int):
                self.assays[key]['source'] = 'Der'
                self.assays[key]['related_assays'] = [value]

    def create_protabank_study(self, literature_dict, pdb_dict, uniprot_dict, submitter_id=1): # default submitter is admin
        """
        Create a json that can be submitted to protabank api

        literature_dict holds a dictionary of PMID:{pub_dict} pairs from a bulk pubmed query because individual queries take too long.
        """
        response = {'title':"No Title", 'description':'No Abstract', 'depositor':submitter_id, 'publication_set':[], 'protein_set':[], 'assay_set':[], 'library_set':[], 'datum_set':[], 'geneexpression_set':[], 'database_source':'PRO'}

        if not self.publications_collected:
            self.collect_publications()
        if not self.proteins_collected:
            self.collect_proteins()
        if not self.assay_data_collected:
            self.collect_assay_datum_sequence(pdb_dict, uniprot_dict)
            self.give_assay_descriptive_names()
            self.shorten_assay_names()
            self.organize_derived_assays()

        if len(self.pubmed_ids) == 1 and self.pubmed_ids[0] in literature_dict.keys():
            pub_dict = literature_dict[self.pubmed_ids[0]]
            response['title']=pub_dict['title']
            response['description']=pub_dict['abstract']
            response['description']+=" Study holds ProTherm entries: " + ", ".join([str(entry.id) for entry in self.entries])
            if self.extra_details != "":
                response['description']+=" Extra Details: "+self.extra_details
            pub_dict['id'] = 1
            pub_dict['year'] = pub_dict['pubdate']# reassign this
            response['publication_set'].append(pub_dict)
        elif len(self.pubmed_ids) > 1 or len(self.publications) > 1:
            self.flags.append('Multiple publications')
            response['description']+=" Study holds ProTherm entries: " + ",".join([str(entry.id) for entry in self.entries])
        else: # no pubmed id so extract from publications
            self.flags.append('No PubMed ID found, need title')
            publication = self.publications.pop()
            pub_dict = publication.generate_dict()
            if pub_dict != None:
                pub_dict['id'] = 1
                response['publication_set'].append(pub_dict)
                response['description']+=" Study holds ProTherm entries: " + ",".join([str(entry.id) for entry in self.entries])
        for id, protein in enumerate(self.proteins):
            swissprot_id, pdb = protein
            prot_dict={}
            if pdb!="":
                prot_dict = pdb_dict[pdb]
            elif swissprot_id !="" and re.search(r'\((.*?)\)',swissprot_id) is not None:
                uniprot_id = re.search(r'\((.*?)\)',swissprot_id).group(1)
                prot_dict= uniprot_dict[uniprot_id]
                if 'id_pdb_id' in prot_dict:
                    prot_dict['id_pdb_id']=prot_dict['id_pdb_id'][:249]
            prot_dict['id_id'] = id
            if 'id_sequence' in prot_dict:
                prot_dict['id_sequence'] = remove_extra_chains(prot_dict['id_sequence'])
            no_id_dict = dict([(key[3:], item) for key,item in prot_dict.items()]) # query has extra id_ in keys
            if len(no_id_dict.keys()) > 1:
                response['protein_set'].append(no_id_dict)
        for id, assay_dict in self.assays.items():
            assay_dict['id'] = id
            response['assay_set'].append(assay_dict)

        for starting_sequence, lid in self.libraries.items():
            lib_dict={}
            lib_dict['id'] = lid
            lib_dict['libinput_set']=[]
            lib_dict['starting_sequence'] = starting_sequence
            lib_dict['description'] = "Mutations for sequence " + starting_sequence
            lib_dict['syntax'] = 'WTnumMUT'
            lib_dict['seq_type'] = 0
            lib_dict['start_ind'] = 1
            for prot_dict in response['protein_set']:
                if 'sequence' in prot_dict.keys() and prot_dict['sequence'] == starting_sequence:
                    lib_dict['protein'] = prot_dict['id']
                    break
            response['library_set'].append(lib_dict)

            # lib_dict['protein']=# find the right protein?

        for datum_dict in self.data.values():
            response['datum_set'].append(datum_dict)

        if len(self.flags) >0:
            error_set=defaultdict(list)
            for flag in self.flags:
                if "Not Valid mut_desc:" in flag:
                    if len(flag.replace("Not Valid mut_desc:","").split()) == 2:
                        mutdesc, sequence=flag.replace("Not Valid mut_desc:","").split()
                        error_set[sequence].append(mutdesc)
            for starting_sequence,mut_desc_set in error_set.items():
                # check that all the data are in error
                lid = self.libraries[starting_sequence]
                if len(mut_desc_set) == len(set([d['mut_desc'] for d in self.data.values() if d['library'] ==lid and d['mut_desc'] != 'WT'])):
                    shift, chain_id = self.find_index_shift(starting_sequence,mut_desc_set,[-1,1])
                    if shift:
                        print "implementing shift", self.id, shift, starting_sequence, mut_desc_set
                        response['library_set'][0]['start_ind']-=shift
                        for d in self.data.values():
                            if d['library'] ==lid:
                                d['mut_desc']=chain_id+d['mut_desc']
                        fixed_flags=[]
                        for i,flag in enumerate(self.flags):
                            if len(flag.replace("Not Valid mut_desc:","").split())==2:
                                mutdesc, sequence=flag.replace("Not Valid mut_desc:","").split()
                            if sequence == starting_sequence and mutdesc in mut_desc_set:
                                fixed_flags.append(i)
                        for flag in sorted(fixed_flags, reverse=True):
                            del self.flags[flag]
        return response

def create_entry_object(obj):
    """
    Create the python object from the dictionary
    """
    entry_objects = {
        'sequence_structure':(SequenceStructure,['protein', 'source', 'length', 'mol_weight', 'pir_id', 'swissprot_id', 'ec_no', 'pmd_no', 'pdb_wild', 'pdb_mutant', 'mutation', 'mutated_chain', 'no_mol', 'sec_str', 'asa']),
        'experimental_condition':(ExperimentalCondition,['temp', 'pH', 'buffer_name', 'buffer_conc', 'ion_name', 'ion_conc', 'additives', 'prot_conc', 'measure', 'method']),
        'data':(ThermoData,['dG_H20', 'ddG_H2O', 'dG', 'ddG', 'Tm', 'dTm', 'dHvH', 'dHcal', 'm', 'Cm', 'dCp', 'state', 'reversibility', 'activity', 'activity_km', 'activity_kcat', 'activity_kd']),
        'literature':(Literature,['key_words', 'reference', 'author', 'remarks', 'related_entries']),
    }
    new_entry_list=[obj['id']]
    for key in ['sequence_structure', 'experimental_condition', 'data', 'literature']:
        value = entry_objects[key]
        params = [obj[key][param] for param in value[1]]
        new_obj = value[0](*params)
        new_entry_list.append(new_obj)
    return ProThermEntry(*new_entry_list)

def create_study_object(study_dict):
    """
    Create the study object from the dictionary
    get the id and create all the entries
    """
    new_study = ProThermStudy(study_dict['id'])
    for entry in study_dict['entries']:
        new_entry = create_entry_object(entry)
        new_study.add_entry(new_entry)
    return new_study

def load_files_from_dump(prefix, dirname):
    """
    Load the data from dumped files
    """
    objs = {}
    for filename in os.listdir(dirname):
        f = re.search(r'^'+prefix+'_(?P<id>[0-9]+).json', filename)
        id=int(f.group('id'))
        obj = json.load(open(dirname+'/'+filename,'r'))
        if prefix == 'entry':
            objs[id] = create_entry_object(obj)
        elif prefix == 'study':
            objs[id] = create_study_object(obj)
    return objs

def property_str(property):
    if isinstance(property, list):
        return ",".join([p.strip() for p in property])
    elif property is not None:
        return str(property.strip())
    else:
        return ""

def parse_result(result, assay_dict):
    """Parse the protherm result to strip out units and assay conditions and return the result only. Change the units and add any conditions to the assay details of the assay_dict"""
    unit=r'(?P<result>[0-9\.-]+) *(?P<unit>[\w \/\%\(\)\*]+)'
    temp=r'\((?P<val>[0-9]+) degree C\)'
    at_temp=r'(?P<unit>[\w \/\%\(\)\*]+)*at (?P<val>[0-9]+) degree C'
    note=r'(?P<unit>[\w \/\%\(\)\*]+)*\((?P<note>[\w -]+)\)'
    if len(result)==0:
        return result
    else:
        try:
            a=float(result.replace(',',""))
            return result.replace(',','')
        except:
            unit_match = re.match(unit, result)
            if unit_match:
                temp_match = re.match(temp, unit_match.group('unit'))
                at_temp_match = re.match(at_temp, unit_match.group('unit'))
                note_match = re.match(note, unit_match.group('unit'))
                if unit_match.group('unit').strip() in unit_converter.keys():
                    assay_dict['units'] = unit_converter[unit_match.group('unit').strip()]
                elif temp_match:
                    assay_dict['temp'] = temp_match.group('val')+ " degree C"
                elif at_temp_match:
                    if at_temp_match.group('unit') and at_temp_match.group('unit').strip() in unit_converter.keys():
                        assay_dict['units'] = unit_converter[at_temp_match.group('unit').strip()]
                    assay_dict['temp'] = at_temp_match.group('val')+ " degree C"
                elif note_match:
                    if note_match.group('unit') and note_match.group('unit').strip() in unit_converter.keys():
                        assay_dict['units'] = unit_converter[note_match.group('unit').strip()]
                    assay_dict['details']+= " ProTherm noted: "+note_match.group('note')
                return unit_match.group('result')
            return result
