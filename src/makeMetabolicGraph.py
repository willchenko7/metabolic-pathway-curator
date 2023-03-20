from src.parseRxnString import parseRxnString

#metabolic_graph = makeMetabolicGraph(reaction_dict,metabolite_dict)

def findRxnsThatConsumeMet(met, reaction_dict, metabolite_dict):
    consuming_rxns = []
    met_id = metabolite_dict[met]['bigg_id']
    for rxn in reaction_dict['reactions']:
        rxn_string = reaction_dict[rxn]['bigg_string']
        [rxn_mets, rxn_stoich, rxn_rev] = parseRxnString(rxn_string)
        if rxn_rev == 1:
            if met_id in rxn_mets:
                consuming_rxns.append(rxn)
        elif rxn_rev == 0:
            if met_id in rxn_mets:
                if rxn_stoich[rxn_mets.index(met_id)] < 0:
                    consuming_rxns.append(rxn)
    return consuming_rxns

def findAdjacentNodes(met,node_reaction_dict,reaction_dict,metabolite_dict):
    adjacent_nodes = []
    consuming_rxns = findRxnsThatConsumeMet(met,node_reaction_dict, metabolite_dict)
    for rxn in consuming_rxns:
        rxn_nodes = reaction_dict[rxn]['nodes']
        adjacent_nodes.append([i for i in rxn_nodes if i != met][0])
    return list(set(adjacent_nodes))

def getListofNodes(reaction_dict):
    nodes = []
    for rxn in reaction_dict['reactions']:
        rxn_nodes = reaction_dict[rxn]['nodes']
        nodes.extend(rxn_nodes)
    return list(set(nodes))

def makeMetabolicGraph(node_reaction_dict,reaction_dict,metabolite_dict):
    nodes = getListofNodes(reaction_dict)
    metabolic_graph = {}
    for node in nodes:
        adjacent_nodes = findAdjacentNodes(node,node_reaction_dict,reaction_dict,metabolite_dict)
        metabolic_graph[node] = adjacent_nodes
    return metabolic_graph

def create_node_reaction_dict(reaction_dict,metabolite_dict):
    node_reaction_dict = {}
    node_reaction_dict['reactions'] = []
    for rxn in reaction_dict['reactions']:
        node_rxn_string = write_node_rxn_string(rxn,reaction_dict,metabolite_dict)
        node_reaction_dict['reactions'].append(rxn)
        node_reaction_dict[rxn] = {}
        node_reaction_dict[rxn]['bigg_string'] = node_rxn_string
    return node_reaction_dict

def write_node_rxn_string(rxn,reaction_dict,metabolite_dict):        
    rxn_string = reaction_dict[rxn]['bigg_string']
    [mets,stoich,rev] = parseRxnString(rxn_string)
    node1 = reaction_dict[rxn]['nodes'][0]
    node2 = reaction_dict[rxn]['nodes'][1]
    node1_id = metabolite_dict[node1]['bigg_id']
    node2_id = metabolite_dict[node2]['bigg_id']
    node1_ind = mets.index(node1_id)
    node2_ind = mets.index(node2_id)
    node_mets = [node1_id,node2_id]
    node_stoich = [stoich[node1_ind],stoich[node2_ind]]
    pos = [idx for idx, val in enumerate(node_stoich) if val > 0][0]
    neg = [idx for idx, val in enumerate(node_stoich) if val < 0][0]
    reactant = node_mets[neg]
    product = node_mets[pos]
    arrow = get_rxn_arrow_from_rev(rev)
    node_rxn_string = reactant + ' ' + arrow + ' ' + product
    return node_rxn_string

def get_rxn_arrow_from_rev(rev):
    if rev == 1:
        arrow = '<=>'
    elif rev == 0:
        arrow = '=>'
    return arrow




