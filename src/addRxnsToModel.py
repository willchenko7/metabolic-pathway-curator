
# computing elementary modes
# Step 1: add all metabolites and reactions in local dict to a cobra model
# step 2: define precursor, product
# step 3: compute elementary modes
# step 4: find all pathways, add them to pathway dict

from cobra import Model, Reaction, Metabolite
from src.parseRxnString import parseRxnString

#model = addAllRxnsFromRxnDict(reaction_dict,metabolite_dict)
#model = addExchangeRxns(model,'pyruvate','ethanol',metabolite_dict)

def addAllRxnsFromRxnDict(reaction_dict,metabolite_dict):
    # add metabolites and reactions to a cobra model
    model = Model('local_network')
    for rxn in reaction_dict['reactions']:
        nodes = reaction_dict[rxn]['nodes']
        node1 = nodes[0]
        node2 = nodes[1]
        id1 = metabolite_dict[node1]['bigg_id']
        id2 = metabolite_dict[node2]['bigg_id']
        rxn_string = reaction_dict[rxn]['bigg_string']
        [mets,stoich,rev] = parseRxnString(rxn_string)
        node1_i = mets.index(id1)
        node2_i = mets.index(id2)
        stoich1 = stoich[node1_i]
        stoich2 = stoich[node2_i]
        [ub,lb] = getBoundsFromRev(rev)
        reaction = Reaction(rxn)
        reaction.lower_bound = lb  
        reaction.upper_bound = ub
        n1 = Metabolite(id1)
        n2 = Metabolite(id2)
        reaction.add_metabolites({
            n1:stoich1,
            n2:stoich2
        })
        model.add_reactions([reaction])
    return model

def getBoundsFromRev(rev):
    if rev == 1:
        ub = 1000
        lb = -1000
    elif rev == 0:
        ub = 1000
        lb = -1000
    return ub, lb
        
def addExchangeRxns(model,precursor,target,metabolite_dict):
    #uptake rxn
    p_id = metabolite_dict[precursor]['bigg_id']
    p = Metabolite(p_id)
    uptake_rxn = Reaction('EX_'+p_id)
    uptake_rxn.lower_bound = -10
    uptake_rxn.upper_bound = 0
    uptake_rxn.add_metabolites({p:-1})
    model.add_reactions([uptake_rxn])
    #excretion rxn
    t_id = metabolite_dict[target]['bigg_id']
    t = Metabolite(t_id)
    exc_rxn = Reaction('EX_'+t_id)
    exc_rxn.lower_bound = 0
    exc_rxn.upper_bound = 1000
    exc_rxn.add_metabolites({t:-1})
    model.add_reactions([exc_rxn])
    return model

