# add pathways
from src.parseRxnString import parseRxnString
from src.makeMetabolicGraph import makeMetabolicGraph
from src.findAllPathsInMetabolicGraph import findAllPathsBetweenTwoMetabolites
from src.retrieve_data_from_local_db import add_db_pathways_between_mets_to_pathway_dict
import re
import os
from src.makeMetabolicGraph import create_node_reaction_dict


def manually_add_to_pathway_dict(pathway_dict,pathway_name,product,pathway_rxns,precursor,metabolite_dict):
    product_id = metabolite_dict[product]['bigg_id'][:-2]
    precursor_id = metabolite_dict[precursor]['bigg_id'][:-2]
    pathway_dict['pathways'].append(product_id)
    pathway_dict[product_id] = {}
    pathway_dict[product_id]['id'] = product_id
    pathway_dict[product_id]['name'] = pathway_name
    pathway_dict[product_id]['rxns'] = pathway_rxns
    pathway_dict[product_id]['product_id'] = product_id
    pathway_dict[product_id]['precursor_ids'] = precursor_id
    return pathway_dict

def add_to_pathway_dict(pathway_dict,precursor,target,reaction_dict,metabolite_dict,mmdb_path):
    node_reaction_dict = create_node_reaction_dict(reaction_dict,metabolite_dict)
    metabolic_graph = makeMetabolicGraph(node_reaction_dict,reaction_dict,metabolite_dict)
    pathways = findAllPathsBetweenTwoMetabolites(precursor,target,metabolic_graph,reaction_dict)
    [pathway_dict,reaction_dict,metabolite_dict] = add_db_pathways_between_mets_to_pathway_dict(target,precursor,mmdb_path,pathway_dict,metabolite_dict,reaction_dict)
    i = get_pathway_naming_index(pathway_dict,target,precursor,metabolite_dict)
    paths_in_dict = get_all_dict_path_rxns_between_mets(target,precursor,pathway_dict,metabolite_dict)
    for pathway in pathways:
        for path in pathway:
            if path in paths_in_dict:
                do_nothing = 0
            else:
                pathway_name = target + '_' + precursor + '_' + str(i)
                pathway_dict['pathways'].append(pathway_name)
                pathway_dict[pathway_name] = {}
                pathway_dict[pathway_name]['id'] = pathway_name
                pathway_dict[pathway_name]['name'] = pathway_name
                pathway_dict[pathway_name]['rxns'] = path
                pathway_dict[pathway_name]['product_id'] = metabolite_dict[target]['bigg_id']
                pathway_dict[pathway_name]['precursor_ids'] = metabolite_dict[precursor]['bigg_id']
                i = i + 1
    return pathway_dict,reaction_dict,metabolite_dict

def get_pathway_naming_index(pathway_dict,target,precursor,metabolite_dict):
    paths = get_all_dict_paths_between_mets(target,precursor,pathway_dict,metabolite_dict)
    if paths == []:
        current_indices = [-1]
    else:    
        current_indices = []
        for path in paths:
            current_indices.append(int(re.split('_',path)[2]))
    return max(current_indices) + 1

def get_all_dict_paths_between_mets(target,precursor,pathway_dict,metabolite_dict):
    precursor_id = metabolite_dict[precursor]['bigg_id']
    target_id = metabolite_dict[target]['bigg_id']
    paths = []
    for path in pathway_dict['pathways']:
        if pathway_dict[path]['product_id'] == target_id and pathway_dict[path]['precursor_ids'] == precursor_id:
            paths.append(path)
    return paths

def get_all_dict_path_rxns_between_mets(target,precursor,pathway_dict,metabolite_dict):
    paths = get_all_dict_paths_between_mets(target,precursor,pathway_dict,metabolite_dict)
    path_rxns = []
    for path in paths:
        path_rxns.append(pathway_dict[path]['rxns'])
    return path_rxns


        


