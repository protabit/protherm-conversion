# ProTherm conversion to ProtaBank studies

We have been given permission to deposit all of the data collected in ProTherm (https://doi.org/10.1093/nar/gkh082) into ProtaBank (https://www.protabank.org/). These scripts automate the process of converting ProTherm entry data into the ProtaBank study format using the following steps: 

* All ProTherm entries are grouped into ProtaBank studies by literature citation where each publication corresponds to a ProtaBank study. 
* All Proteins are extracted from the provided PDB and SWISSPROT\_IDs in the ProTherm sequence and structure section of the entry. Proteins without either of these ids have not been converted, see *not\_submitted\_studies.txt* for all studies with missing information entries.
* All assay information is extracted from the experimental condition section of the ProTherm entry. All conditions are assumed to apply for all measurements. Units are assumed to be the ProTherm standard unless otherwise specified in the Thermodynamic data section.
*  Mutational data is associated with a sequence by extracting the protein sequence from the provided PDB or uniport ID and ensuring that the specified mutation matches the WT#NUM description in MUTATION. Entries without either of these ids ,or where no match to the wild type sequence after attempting any possible starting index shifts have not been converted, see *not\_sumbitted\_studies.txt* for all studies with missing information.

### Prerequisites

Install needed python packages

```
pip install -r requirements.txt
```

Unzip the ProTherm data

```
unzip ProTherm.dat.gz
```

### Comments
Modifications to the original ProTherm.dat file

*Deleted 2PGK from PDB ids because PDB entry has all Xs and gives no sequence                                        

*Deleted 1XBE because structure is withdrawn replaced with Q7YUF0\_TRYCR (Q7YUF0)

*Deleted 94.3 from ProTherm.dat because it is not a pdb code

*Placed 1OHV info in 1GTX entry because 1GTX is obsolete and superseeded by 1OHV

*Placed 1R51 info in 1UOX  entry because 1UOX is obsolete and superseeded by 1R51

*Edited ProTherm to repalce 166H (does not exist) with 166L which matches the mutant description

*Edited ProTherm to repalce RF5V (does not exist) with 3F5V which matches the mutant description
## Authors

* **Connie Wang** 

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details

## Acknowledgments

* Thanks to ProTherm for giving us permission to enter their data into ProtaBank
