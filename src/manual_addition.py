
# manual addition
from src.create_metabolite_dict import get_charged_formula
import re

#metabolite_dict = manually_add_to_metabolite_dict(metabolite_dict,'ethanol')
#reaction_dict = manually_add_to_reaction_dict(reaction_dict)

def manually_add_to_metabolite_dict(metabolite_dict,compound):  
    metabolite_dict[compound] = {}
    print(' //////////// \n')
    print('Manually adding ',compound,' to metabolite dictionary')
    print(' \\    If any info is unknown, enter "n"')
    print(' \\    If any info is not available, enter "na"')
    metabolite_dict[compound]['kegg_id'] = input('Input KEGG ID: ')
    neutral_formula = [input('Input neutral formula: ')]
    metabolite_dict[compound]['neutral_formula'] = neutral_formula
    metabolite_dict[compound]['kegg_reactions'] = []
    print('Please incude compartment (i.e. _c, _e)')
    metabolite_dict[compound]['bigg_id'] = input('Input BIGG ID: ')
    metabolite_dict[compound]['bigg_charges'] = []
    charge = int(input('Input charge: '))
    metabolite_dict[compound]['formula'] = get_charged_formula(neutral_formula,charge)
    metabolite_dict[compound]['bigg_formula'] = metabolite_dict[compound]['formula']
    metabolite_dict[compound]['charge'] = charge
    metabolite_dict['metabolites'].append(compound)
    return metabolite_dict

def manually_add_to_reaction_dict(bigg_rxn,reaction_dict):
    print(' /////////// \n')
    print('Manually adding reaction to reaction dictionary')
    print(' \\    If any info is unknown, enter "n"')
    print(' \\    If any info is not available, enter "na"')
    print(' \n WARNING: All metabolites in reaction must be manually added if not already done so')
    print(' \n ')
    name = input('Input reaction name: ')
    print(' \n Please use BIGG format for naming reaction ID')
    bigg_rxn = input('Input reaction ID: ')
    reaction_dict['reactions'].append(bigg_rxn)
    reaction_dict[bigg_rxn] = {}
    reaction_dict[bigg_rxn]['bigg_rxn'] = bigg_rxn
    reaction_dict[bigg_rxn]['kegg_rxn'] = input('Input KEGG ID: ')
    print(' \n Please use BIGG format for reaction string')
    reaction_dict[bigg_rxn]['bigg_string'] = input('Input reaction string: ')
    reaction_dict[bigg_rxn]['bigg_name'] = name
    reaction_dict[bigg_rxn]['bigg_mets'] = []
    reaction_dict[bigg_rxn]['kegg_name'] = name
    reaction_dict[bigg_rxn]['kegg_mets'] = []
    reaction_dict[bigg_rxn]['kegg_string'] = []
    reaction_dict[bigg_rxn]['kegg_definition'] = ''
    brenda = input('Input Brenda EC Number (separate w comma)')
    reaction_dict[bigg_rxn]['BRENDA_enzymes'] = [i for i in re.split(',',brenda)]
    print('Use full names for nodes')
    node1 = input('Input node1 for reaction: ')
    node2 = input('Input node2 for reaction: ')
    reaction_dict[bigg_rxn]['nodes'] = [node1,node2]
    return reaction_dict
