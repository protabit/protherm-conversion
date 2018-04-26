# -*- coding: utf-8 -*-
from __future__ import unicode_literals
import re
ACCEPTEDAA=['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']

def valid_protein_sequence(seq):
    """
    seq -- is a string of characters, checks if the sequence is allowed
    """
    if seq is None or seq=="":
        return False
    for chain in seq.split("/"):
        chain_sequence=chain[(chain.find(':')+1):]
        for letter in chain_sequence:
            if letter not in ACCEPTEDAA:
                return False
    return True

def valid_mut_desc(mut_desc, seq):
    """
    checks if the mut_desc is a valid mutant description for the sequence seq
    returns (is_valid?, mut_desc) where the new mut_desc may have a chain_id assigned
    """
    md_re=r'(?P<wt>[A-Z])(?P<pos>[0-9]+)(?P<mut>[A-Z\(\)]+)$'
    if mut_desc == 'WT' or mut_desc=='Wild':
        return (True, mut_desc)
    if mut_desc == None or mut_desc=='' or seq == None:
        return (False, mut_desc)


    for chain in seq.split("/"):
        chain_id = chain[:2] if ':' in chain else "" #chainid and :
        chain_sequence=chain[(chain.find(':')+1):]

        mut_desc_for_chain=True
        # TODO need to split mut_desc by chain??
        for mut_desc_single in mut_desc.split('+'):
            search = re.search(md_re, mut_desc_single)
            if search is None:
                return (False, mut_desc)
            WT=search.group('wt')
            num=int(search.group('pos'))-1 # 0 index
            mut=search.group('mut')

            if num >= len(chain_sequence):
                mut_desc_for_chain=False
                break # not this chain
            elif chain_sequence[num]!=WT:
                mut_desc_for_chain=False
                break
        if mut_desc_for_chain: # found chain
            return (True, chain_id+mut_desc)
    return (False, mut_desc)

def remove_extra_chains(sequence):
    """
    If it has multiple chains that are identical, only return one chain
    """
    if not valid_protein_sequence(sequence):
        return sequence
    elif '/' in sequence:
        chains=set()
        for chain in sequence.split('/'):
            chains.add(chain[(chain.find(':')+1):])
        if len(chains)==1:
            return chains.pop()
    return sequence

def identify_all_potential_shifts(mut_desc, chain_sequence,allowed_shifts):
        """
        identify all potential shifts
        """
        md_re=r'(?P<wt>[A-Z])(?P<pos>[0-9]+)(?P<mut>[A-Z\(\)]+)$'
        saved_shifts=[]
        for shift in allowed_shifts:
            proposed_shift=shift
            for mut_desc_single in mut_desc.split('+'):
                search = re.search(md_re, mut_desc_single)
                if search is None:
                    return None
                WT=search.group('wt')
                num=int(search.group('pos'))-1 # 0 index
                mut=search.group('mut')
                if (num+proposed_shift >= len(chain_sequence) or num+proposed_shift < 0 or chain_sequence[num+proposed_shift] != WT):
                    proposed_shift =0
                    break
            if proposed_shift !=0: # found a shift that works
                saved_shifts.append(proposed_shift)
        return saved_shifts
