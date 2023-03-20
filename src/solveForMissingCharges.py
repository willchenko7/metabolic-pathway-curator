
# calculate missing charge
from src.parseRxnString import parseRxnString
from src.balance_check import find_missing_metabolites
from src.balance_check import findMetNameFromID
from src.create_metabolite_dict import get_charged_formula


#metabolite_dict = solveForMissingCharges(metabolite_dict,reaction_dict)

def solveForMissingChargesinRecursiveLoop(metabolite_dict,reaction_dict):
    initial_missing_mets = find_missing_metabolites(metabolite_dict)
    metabolite_dict = solveForMissingCharges(metabolite_dict,reaction_dict)
    later_missing_mets = find_missing_metabolites(metabolite_dict)
    if initial_missing_mets != later_missing_mets:
        metabolite_dict = solveForMissingChargesinRecursiveLoop(metabolite_dict,reaction_dict)
    else:
        do_nothing = 0
    return metabolite_dict

def findRxnsFromMet(met, reaction_dict, metabolite_dict):
    involved_rxns = []
    met_id = metabolite_dict[met]['bigg_id']
    for rxn in reaction_dict['reactions']:
        rxn_string = reaction_dict[rxn]['bigg_string']
        [rxn_mets, rxn_stoich, rxn_rev] = parseRxnString(rxn_string)
        if met_id in rxn_mets:
            involved_rxns.append(rxn)
    return involved_rxns

def solveForMissingCharges(metabolite_dict,reaction_dict):
    missing_mets = find_missing_metabolites(metabolite_dict)
    for met in missing_mets:
        met_name = findMetNameFromID(met,metabolite_dict)
        involved_rxns = findRxnsFromMet(met_name,reaction_dict,metabolite_dict)
        solvable_rxns = findSolvableRxns(involved_rxns,reaction_dict,metabolite_dict)
        if len(solvable_rxns) > 0:
            rxn = solvable_rxns[0]
            missing_charge = calculateMissingCharge(rxn,met,reaction_dict,metabolite_dict)
            metabolite_dict[met_name]['charge'] = missing_charge
            neutral_formula = metabolite_dict[met_name]['neutral_formula']
            metabolite_dict[met_name]['formula'] = get_charged_formula(neutral_formula,missing_charge)
    return metabolite_dict
        
def calculateMissingCharge(rxn,missing_met,reaction_dict,metabolite_dict):
    [mets,stoich,rev] = parseRxnString(reaction_dict[rxn]['bigg_string'])
    cum = 0
    met_i = mets.index(missing_met)
    for met in mets:
        if met != missing_met:
            m_i = mets.index(met)
            met_name = findMetNameFromID(met,metabolite_dict)
            cum = cum + stoich[m_i]*metabolite_dict[met_name]['charge']
    missing_charge = -1*(cum/stoich[met_i])
    return missing_charge
        
def isRxnSolvable(rxn,reaction_dict,metabolite_dict):
    [mets,stoich,rev] = parseRxnString(reaction_dict[rxn]['bigg_string'])
    solve_check = 0
    for met in mets:
        met_name = findMetNameFromID(met,metabolite_dict)
        if metabolite_dict[met_name]['charge'] == 'NA':
            solve_check = solve_check + 1
    if solve_check > 1:
        solvable = False
    else:
         solvable = True
    return solvable

def findSolvableRxns(involved_rxns,reaction_dict,metabolite_dict):
    solvable_rxns = []
    for rxn in involved_rxns:
        solvable = isRxnSolvable(rxn,reaction_dict,metabolite_dict)
        if solvable is True:
            solvable_rxns.append(rxn)
    return solvable_rxns

        
    
