
from src.parseRxnString import parseRxnString
import json

def export2json(filename,reaction_dict,metabolite_dict):
    escher = {
        "metabolites": create_met_json_array(metabolite_dict),
        "reactions":create_rxn_json_array(reaction_dict),
        "genes":[],
        "id":"test",
        "name":"j",
        "compartments":{
            "c":"cytosol",
            "e":"extracellular"
            },
        "version":"1"
        }
    with open(filename, "w") as outfile: 
        json.dump(escher,outfile,indent=0)
    return escher

def create_rxn_json_array(reaction_dict):
    json_rxn_array = []
    for rxn in reaction_dict['reactions']:
        rxn_string = reaction_dict[rxn]['bigg_string']
        ind_rxn_dict = {
            "id":rxn,
            "name":reaction_dict[rxn]['bigg_name'],
            "metabolites":create_rxn_met_dict(rxn_string),
            "lower_bound":from_rev_get_bounds(rxn_string)[0],
            "upper_bound":from_rev_get_bounds(rxn_string)[1]
            }
        json_rxn_array.append(ind_rxn_dict)
    return json_rxn_array

def create_rxn_met_dict(rxn_string):
    [mets,stoich,rev] = parseRxnString(rxn_string)
    metabolites = {}
    for i in range(0,len(mets)):
        met = mets[i]
        s = stoich[i]
        metabolites[met] = s
    return metabolites
            
def from_rev_get_bounds(rxn_string):
    [mets,stoich,rev] = parseRxnString(rxn_string)
    if rev == 0:
        lb = 0
        ub = 1000
    elif rev == 1:
        lb = -1000
        ub = 1000
    return lb,ub

def create_met_json_array(metabolite_dict):
    json_met_array = []
    for met in metabolite_dict['metabolites']:
        ind_met_dict = {
        "id":metabolite_dict[met]['bigg_id'],
        "name":met,
        "compartement":metabolite_dict[met]['bigg_id'][-1],
        "charge":metabolite_dict[met]['charge'],
        "formula":metabolite_dict[met]['formula'],
        "annotation":{}
        }
        json_met_array.append(ind_met_dict)
    return json_met_array

    
