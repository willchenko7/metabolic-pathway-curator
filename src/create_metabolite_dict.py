import requests
from bs4 import BeautifulSoup
import re
import numpy as np
from src.retrieve_data_from_local_db import get_list_of_mets_in_db
from src.retrieve_data_from_local_db import retrieve_metabolite_info_from_db
import os

#compound = 'pyruvate'
#metabolite_dict = add_to_metabolite_dict(metabolite_dict,compounds,bigg2kegg_db)

def add_to_metabolite_dict(metabolite_dict,compounds,bigg2kegg_db,mmdb_path):
    db_mets = get_list_of_mets_in_db(mmdb_path)
    for compound in compounds:
        if compound in db_mets:
            print(', \n ',compound,' already has an entry in your local database.')
            print('Would you like to retrieve information from your database?')
            user_input = input('Enter "y" for yes or "n" for no: ')
            if user_input == "y":
                metabolite_dict = retrieve_metabolite_info_from_db(compound,mmdb_path,metabolite_dict)
            elif user_input == "n":
                [kegg_id, kegg_formula, kegg_reactions, kegg_dblinks, bigg_id, bigg_charges, bigg_formula, charge, charged_formula, bigg2kegg_db] = get_compound_information(compound,bigg2kegg_db)
                metabolite_dict[compound] = {}
                metabolite_dict[compound]['kegg_id'] = kegg_id
                metabolite_dict[compound]['neutral_formula'] = kegg_formula
                metabolite_dict[compound]['kegg_reactions'] = kegg_reactions
                metabolite_dict[compound]['bigg_id'] = bigg_id
                metabolite_dict[compound]['bigg_charges'] = bigg_charges
                metabolite_dict[compound]['bigg_formula'] = bigg_formula
                metabolite_dict[compound]['formula'] = charged_formula
                metabolite_dict[compound]['charge'] = charge
                metabolite_dict[compound]['kegg_dblinks'] = kegg_dblinks
                metabolite_dict['metabolites'].append(compound)
        else:
            [kegg_id, kegg_formula, kegg_reactions, kegg_dblinks, bigg_id, bigg_charges, bigg_formula, charge, charged_formula, bigg2kegg_db] = get_compound_information(compound,bigg2kegg_db)
            metabolite_dict[compound] = {}
            metabolite_dict[compound]['kegg_id'] = kegg_id
            metabolite_dict[compound]['neutral_formula'] = kegg_formula
            metabolite_dict[compound]['kegg_reactions'] = kegg_reactions
            metabolite_dict[compound]['bigg_id'] = bigg_id
            metabolite_dict[compound]['bigg_charges'] = bigg_charges
            metabolite_dict[compound]['bigg_formula'] = bigg_formula
            metabolite_dict[compound]['formula'] = charged_formula
            metabolite_dict[compound]['charge'] = charge
            metabolite_dict[compound]['kegg_dblinks'] = kegg_dblinks
            metabolite_dict['metabolites'].append(compound)
    return metabolite_dict,bigg2kegg_db
    
def get_compound_information(compound,bigg2kegg_db):
    kegg_compound = get_kegg_compound(compound)
    [kegg_id,kegg_formula,kegg_reactions, kegg_name,kegg_dblinks] = get_kegg_compound_information(kegg_compound)
    bigg_ids = [s for s in bigg2kegg_db if kegg_id in s[1]]
    if bigg_ids == []:
        print(' \n ///// WARNING ///////')
        print('No BiGG metabolites found for the selected KEGG compound!')
        #user_input = input("Would you like to create a BIGG metabolite for this KEGG compound? \n \t Type 'y' for Yes and 'n'for No: ")
        #if user input is 'y':
        [bigg_id,bigg_charges,bigg_formula,bigg2kegg_db] = create_BIGG_met_from_KEGG_cpd(compound,kegg_id,kegg_formula,bigg2kegg_db)
        charge = bigg_charges[0]
        charged_formula = get_charged_formula(kegg_formula,charge)
    else:
        bigg_id = decide_compartment(bigg_ids)[0]
        [bigg_charges,bigg_formula] = get_bigg_metabolite_information(bigg_id)
        charge = define_metabolite_charge(bigg_charges)
        charged_formula = get_charged_formula(kegg_formula,charge)
    return kegg_id, kegg_formula, kegg_reactions, kegg_dblinks, bigg_id, bigg_charges, bigg_formula, charge, charged_formula, bigg2kegg_db

def create_BIGG_met_from_KEGG_cpd(compound,kegg_id,kegg_formula,bigg2kegg_db):
    bigg_formula = kegg_formula[0]
    print('Create a metabolite ID for ',compound)
    bigg_id = input('Enter ID: ')
    compartment = input('Enter compartment: ')
    bigg_id = bigg_id +'_'+compartment
    print(' \n  \n If charge is not known, eneter "NA"')
    bigg_charge = input('Enter compound charge: ')
    bigg_charges = [bigg_charge]
    bigg2kegg_db.append([bigg_id,kegg_id])
    return bigg_id, bigg_charges, bigg_formula, bigg2kegg_db

def get_charged_formula(kegg_formula,charge):
    neutral_formula = kegg_formula[0]
    if charge == 'NA':
        charged_formula = neutral_formula
    else:
        elements = " ".join(re.split("[^a-zA-Z]*", neutral_formula))
        elements = elements.replace(' ','')
        [h_neutral_coeff,h_start_index,h_end_index] = get_element_coeff(neutral_formula,'H')
        if h_neutral_coeff == 0:
            charged_formula = neutral_formula
        elif h_neutral_coeff == 1:
            charged_formula = neutral_formula.replace('H','')
        else:
            h_charged_coeff = int(h_neutral_coeff + charge)
            charged_formula = neutral_formula[0:h_start_index] + str(h_charged_coeff) + neutral_formula[h_end_index:len(neutral_formula)]
    return charged_formula

def get_element_coeff(formula,element):
    elements = " ".join(re.split("[^a-zA-Z]*", formula))
    elements = elements.replace(' ','')
    if elements == element or element not in elements:
        coeff = 0
        start_index = 0
        end_index = 0
    else:
        element_index = elements.index(element)
        element_index_in_formula = formula.index(element)
        start_index = element_index_in_formula + 1
        if element_index == len(elements)-1:
            end_index = len(formula)
            coeff = formula[start_index:end_index]
        else:
            element_after_h = elements[element_index + 1]
            end_index = formula.index(element_after_h)
            coeff = formula[start_index:end_index]
        if coeff == '':
            coeff = 1
    return int(coeff), start_index, end_index
       
def define_metabolite_charge(bigg_charges):
    print(' \n DEFINE METABOLITE CHARGE ')
    print('For reference, the charges from BIGG are .... ')
    print(bigg_charges)
    print(' ......')
    print('If charge is not known, enter "NA", and it can be calculated later with ... ')
    print(' ....')
    user_input = input('    Enter metabolite charge: ')
    if user_input == 'NA':
        charge = user_input
    else:
        charge = int(user_input)
    return charge

def get_bigg_metabolite_information(bigg_id):
    met = bigg_id[:-2]
    base = 'http://bigg.ucsd.edu/api/v2/universal/metabolites/'
    url = base + met
    response = requests.get(url)
    tsv_data = response.text
    data = re.split(',',tsv_data)
    bigg_formula = [s for s in data if "formulae" in s]
    bigg_formula = re.split(': ',bigg_formula[0])[1]
    bigg_formula = bigg_formula.lstrip('["')
    bigg_formula = bigg_formula.rstrip('"]')
    if bigg_formula == '':
        bigg_formula = 'NA'
    charge_index = find_index_that_contains_substring_in_list("charges",data)[0]
    name_index = find_index_that_contains_substring_in_list("name",data)[0]
    charges = data[charge_index:name_index]
    charges[0] = charges[0].replace(' "charges": [','')
    charges[-1] = charges[-1].replace(']','')
    if charges == ['']:
        charges = 'NA'
    else:
        charges = [int(s) for s in charges]
    return charges,bigg_formula

def decide_compartment(bigg_ids):
    bigg_id = []
    for ids in bigg_ids:
        if ids[0][-1] is 'c':
            bigg_id.append(ids[0])
    return bigg_id

def get_kegg_compound_information(kegg_compound):
    cpd_id = kegg_compound[0][-6:]
    base = "http://rest.kegg.jp/get/"
    url = base + cpd_id
    response = requests.get(url)
    tsv_data = response.text
    data = re.split('\n|\t',tsv_data)
    d = np.asarray(data)
    name = [s for s in d if "NAME" in s]
    name = re.split('NAME| |;',name[0])
    name = [s for s in name if s][0]
    formula = [s for s in d if "FORMULA" in s]
    formula = re.split('FORMULA| ',formula[0])
    formula = [s for s in formula if s]
    reaction_index = [i for i, s in enumerate(d) if 'REACTION' in s][0]
    categories = [s[0] for s in d if s]
    non_empty_rows_i = [i for i, s in enumerate(categories) if ' ' is not s]
    dblinks_index = [i for i, s in enumerate(d) if 'DBLINKS' in s][0]
    category_after_dblinks_index = non_empty_rows_i[non_empty_rows_i.index(dblinks_index)+1]
    datablinks = d[dblinks_index:category_after_dblinks_index]
    dblinks = []
    for link in datablinks:
        link = re.split('DBLINKS |  ',link)
        link = [s for s in link if s][0]
        dblinks.append(link)
    category_after_reaction_index = non_empty_rows_i[non_empty_rows_i.index(reaction_index)+1]
    #pathway_index = [i for i, s in enumerate(d) if 'PATHWAY' in s][0]
    reactions = d[reaction_index:category_after_reaction_index]
    rxns = []
    for rxn in reactions:
        rxn = re.split('REACTION| ',rxn)
        rxn = [s for s in rxn if s]
        rxns = rxns + rxn
    return cpd_id, formula, rxns, name, dblinks

def get_kegg_compound(compound):
    base = "http://rest.kegg.jp/find/compound/"
    url = base + compound
    response = requests.get(url)
    soup = BeautifulSoup(response.text, 'html.parser')
    tsv_data = response.text
    data = re.split('\n|\t',tsv_data)
    d = np.asarray(data)
    clean_data = clean_requests(d)
    display_compounds = display_possible_compounds(clean_data)
    display_compounds = np.asarray(display_compounds)
    display_compounds = np.vstack(display_compounds)
    print("The following compounds were returned from your input. Please select the compound wish to select.")
    print(display_compounds)
    user_input = input("Enter number of desired compound: ")
    kegg_compound = clean_data[int(user_input)]
    return kegg_compound

def display_possible_compounds(clean_data):
    i = 0
    display_compounds = []
    for compound in clean_data:
        name = compound[1]
        display_compounds.append([i,name])
        i = i + 1
    return display_compounds

def new_index(i,ii,iii,n_cols):
    if i % n_cols == 0:
        ii = ii + 1
        iii = 0
    else:
        ii = ii
        iii = iii+1
    return ii,iii

def init_matrix(m, n):
    """Creates a m by n matrix filled with zeros."""
    return [[0]*n for i in range(m)]

def clean_requests(d):
    n_cols = 2
    n_rows = int((len(d)-1)/n_cols)
    ii = 0
    iii = 0
    clean_data = init_matrix(n_rows,n_cols)
    clean_data[ii][iii] = d[0]
    for i in range(1,len(d)-1):
        c_val = d[i]
        [ii,iii] = new_index(i,ii,iii,n_cols)
        clean_data[ii][iii] = c_val
    return clean_data


def find_index_that_contains_substring_in_list(substring,list):
    index  = []
    for i in range(0,len(list)):
        j = list[i]
        if substring in j:
            index.append(i)
    return index
