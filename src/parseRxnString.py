# parse a reaction into mets, stoich, and reversibility
import re

#rxn_string = reaction_dict['PDH']['bigg_string']
#[mets,stoich,rev] = parseRxnString(rxn_string)

def parseRxnString(rxn_string):
    [reactants, products, rev] = splitRxn(rxn_string)
    [reactant_stoich, reactant_mets] = parseHalfRxn(reactants)
    [product_stoich, product_mets] = parseHalfRxn(products)
    reactant_stoich = [i*-1 for i in reactant_stoich]
    mets = reactant_mets + product_mets
    stoich = reactant_stoich + product_stoich
    return mets, stoich, rev

def splitRxn(rxn_string):
    if '<=>' in rxn_string:
        rev = 1
        reactants = re.split('<=>',rxn_string)[0]
        products = re.split('<=>',rxn_string)[1]
    else:
        rev = 0
        reactants = re.split('=>',rxn_string)[0]
        products = re.split('=>',rxn_string)[1]
    return reactants, products, rev

def parseHalfRxn(half_rxn):
    half_rxn_stoich = []
    half_rxn_mets = []
    half_rxn = re.split('\+',half_rxn)
    for r in half_rxn:
        r = re.split(' ',r)
        r = [i for i in r if i]
        if len(r) == 1:
            r_stoich = 1
            r_met = r[0]
            half_rxn_stoich.append(r_stoich)
            half_rxn_mets.append(r_met)
        elif len(r) == 2:
            r_stoich = float(r[0])
            r_met = r[1]
            half_rxn_stoich.append(r_stoich)
            half_rxn_mets.append(r_met)
    return half_rxn_stoich, half_rxn_mets
