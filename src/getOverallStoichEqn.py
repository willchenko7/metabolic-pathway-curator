# get the overall stoichiometric equation of a path
from src.parseRxnString import parseRxnString
import numpy as np

#path = 'ethanol_pyruvate_0'
#precursor = 'pyruvate'
#ov_string = write_ov_path_eqn(path,precursor,pathway_dict,reaction_dict,metabolite_dict)

def write_ov_path_eqn(path,precursor,pathway_dict,reaction_dict,metabolite_dict):
    print('Making overall path for ... ',path)
    [ov_mets, ov_stoich] = get_pathway_stoich(path,precursor,pathway_dict,reaction_dict,metabolite_dict)
    ov_string = write_rxn_string_from_mets_and_stoich(ov_mets,ov_stoich)
    return ov_string

def get_mets_in_path(path,pathway_dict,reaction_dict):
    ov_mets = []
    rxns = pathway_dict[path]['rxns']
    for rxn in rxns:
        rxn_string = reaction_dict[rxn]['bigg_string']
        [mets,stoich,rev] = parseRxnString(rxn_string)
        ov_mets = list(set(ov_mets+mets))
    return ov_mets

def get_pathway_stoich(path,precursor,pathway_dict,reaction_dict,metabolite_dict):
    ov_mets = get_mets_in_path(path,pathway_dict,reaction_dict)
    ov_stoich = list(np.zeros((len(ov_mets),), dtype=int))
    rxns = pathway_dict[path]['rxns']
    current_node = precursor
    for rxn in rxns:
        rxn_string = reaction_dict[rxn]['bigg_string']
        [mets,stoich,rev] = parseRxnString(rxn_string)
        current_node_index = mets.index(metabolite_dict[current_node]['bigg_id'])
        current_node_stoich = stoich[current_node_index]
        nodes = reaction_dict[rxn]['nodes']
        next_node = get_next_node(current_node,rxn,reaction_dict)
        if current_node_stoich > 0 and rev == 1:
            stoich = [i*-1 for i in stoich]
        for met in mets:
            met_stoich = stoich[mets.index(met)]
            ov_met_index = ov_mets.index(met)
            ov_stoich[ov_met_index] = ov_stoich[ov_met_index] + met_stoich
        current_node = next_node
    return ov_mets, ov_stoich

def write_rxn_string_from_mets_and_stoich(mets,stoich):
    [reactant_indices, product_indices] = get_half_rxn_indices(stoich)
    reactant_string = write_half_rxn(mets,stoich,reactant_indices)
    product_string = write_half_rxn(mets,stoich,product_indices)
    rxn_string = reactant_string + ' => ' + product_string
    return rxn_string

def write_half_rxn(mets,stoich,indices):
    half_rxn = ''
    last_index = indices[len(indices)-1]
    for index in indices:
        s = stoich[index]
        if s < 0:
            s = s*-1
        met = mets[index]
        if s == 1:
            met_str = met
        else:
            met_str = str(s) + ' ' + met
        if index == last_index:
            half_rxn = half_rxn + met_str
        else:
            half_rxn = half_rxn + met_str + ' + '
    return half_rxn

def get_half_rxn_indices(stoich):
    i = 0
    reactant_indices = []
    product_indices = []
    for s in stoich:
        if s < 0:
            reactant_indices.append(i)
        elif s > 0:
            product_indices.append(i)
        i = i + 1
    return reactant_indices, product_indices
        
def get_next_node(current_node,rxn,reaction_dict):
    nodes = reaction_dict[rxn]['nodes']
    current_node_index = nodes.index(current_node)
    if current_node_index == 0:
        next_node_index = 1
    elif current_node_index == 1:
        next_node_index = 0
    next_node = nodes[next_node_index]
    return next_node
