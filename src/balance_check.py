# check
from src.parseRxnString import parseRxnString
from molmass import Formula

#[unbalanced_rxns, charge_check] = mass_charge_balance_check(reaction_dict,metabolite_dict)


def find_missing_metabolites(metabolite_dict):
    missing_mets = []
    for met in metabolite_dict['metabolites']:
        if metabolite_dict[met]['charge'] == 'NA':
            missing_mets.append(metabolite_dict[met]['bigg_id'])
    return missing_mets

def charge_balance_check(reaction_dict,metabolite_dict):
    # possibly unused function
    unbalanced_rxns = []
    charge_check = []
    missing_mets = find_missing_metabolites(metabolite_dict)
    for rxn in reaction_dict['reactions']:
        rxn_charge_check = check_for_rxn(rxn,reaction_dict,metabolite_dict,missing_mets)
        charge_check.append(rxn_charge_check)
        if rxn_charge_check is False:
            unbalanced_rxns.append(rxn)
    return unbalanced_rxns, charge_check

def charge_check_for_rxn(rxn,reaction_dict,metabolite_dict):
    missing_mets = find_missing_metabolites(metabolite_dict)
    [mets, stoich, rev] = parseRxnString(reaction_dict[rxn]['bigg_string'])
    met_check =  any(item in missing_mets for item in mets)
    if met_check is True:
        charge_check = False
    else:
        charge_sum = rxnChargeBalance(rxn,reaction_dict,metabolite_dict)
        if charge_sum == 0:
            charge_check = True
        else:
            charge_check = False
    return charge_check

def rxnChargeBalance(rxn,reaction_dict,metabolite_dict):
    charge_sum = 0
    [mets, stoich, rev] = parseRxnString(reaction_dict[rxn]['bigg_string'])
    for met in mets:
        met_i = mets.index(met)
        met_id = met[:-2]
        met_name = findMetNameFromID(met,metabolite_dict)
        charge = metabolite_dict[met_name]['charge']
        charge_sum = charge_sum + charge*stoich[met_i]
    return charge_sum

def findMetNameFromID(met_id,metabolite_dict):
    met_name = []
    for name in metabolite_dict['metabolites']:
        if metabolite_dict[name]['bigg_id'] == met_id:
            met_name.append(name)
    return met_name[0]

def mass_check_for_rxn(rxn,reaction_dict,metabolite_dict):
    missing_mets = find_missing_metabolites(metabolite_dict)
    [mets, stoich, rev] = parseRxnString(reaction_dict[rxn]['bigg_string'])
    met_check =  any(item in missing_mets for item in mets)
    if met_check is True:
        mass_check = False
    else:
        mass_sum = check_rxn_mass_sum(rxn,reaction_dict,metabolite_dict)
        if mass_sum <= 0.01 and mass_sum >= -0.01:
            mass_check = True
        else:
            mass_check = False
    return mass_check
    
def check_rxn_mass_sum(rxn,reaction_dict,metabolite_dict):
    [mets, stoich, rev] = parseRxnString(reaction_dict[rxn]['bigg_string'])
    mass_sum = 0
    for met in mets:
        met_i = mets.index(met)
        #print(met)
        met_name = findMetNameFromID(met,metabolite_dict)
        met_formula = metabolite_dict[met_name]['formula']
        #print(met_formula)
        mf = Formula(met_formula)
        met_mass = mf.mass
        mass_sum = mass_sum + met_mass*stoich[met_i]
    return mass_sum

def mass_check_reaction_dict(reaction_dict,metabolite_dict):
    for rxn in reaction_dict['reactions']:
        print('Checking mass balance for ',rxn)
        mass_check = mass_check_for_rxn(rxn,reaction_dict,metabolite_dict)
        reaction_dict[rxn]['mass_balance'] = mass_check
    return reaction_dict

def charge_check_reaction_dict(reaction_dict,metabolite_dict):
    for rxn in reaction_dict['reactions']:
        print('Checking charge balance for ',rxn)
        charge_check = charge_check_for_rxn(rxn,reaction_dict,metabolite_dict)
        reaction_dict[rxn]['charge_balance'] = charge_check
    return reaction_dict

def find_unbalanced_rxns(reaction_dict):
    unbalanced_rxns = []
    for rxn in reaction_dict['reactions']:
        if rxn != 'reactions':
            if reaction_dict[rxn]['charge_balance'] == False or reaction_dict[rxn]['mass_balance'] == False:
                unbalanced_rxns.append(rxn)
    return unbalanced_rxns
    
