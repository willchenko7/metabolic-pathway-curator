# add metadata for pathways
from src.getOverallStoichEqn import write_ov_path_eqn
from src.balance_check import findMetNameFromID

def add_to_metadata(path,metadata,pathway_dict,reaction_dict,metabolite_dict):
    precursor_id = pathway_dict[path]['precursor_ids']
    precursor = findMetNameFromID(precursor_id,metabolite_dict)
    ov_string = write_ov_path_eqn(path,precursor,pathway_dict,reaction_dict,metabolite_dict)
    num_rxns = len(pathway_dict[path]['rxns'])
    target_id = pathway_dict[path]['product_id']
    target = findMetNameFromID(target_id,metabolite_dict)
    metadata[path] = {}
    metadata[path]['precursor'] = precursor
    metadata[path]['target'] = target
    metadata[path]['number_of_rxns'] = num_rxns
    metadata[path]['ov_string'] = ov_string
    metadata['pathways'].append(path)
    return metadata

def add_all_paths_to_metadata(metadata,pathway_dict,reaction_dict,metabolite_dict):
    for path in pathway_dict['pathways']:
        metadata = add_to_metadata(path,metadata,pathway_dict,reaction_dict,metabolite_dict)
    return metadata
    
    
