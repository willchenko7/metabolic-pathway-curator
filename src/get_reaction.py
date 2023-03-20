# complete rxn from metabolites in metabolite dict
import requests
from bs4 import BeautifulSoup
import re
import numpy as np
from src.create_metabolite_dict import get_kegg_compound_information
from src.create_metabolite_dict import get_bigg_metabolite_information
from src.create_metabolite_dict import define_metabolite_charge
from src.create_metabolite_dict import get_element_coeff
from src.create_metabolite_dict import get_charged_formula
from src.retrieve_data_from_local_db import get_list_of_met_bigg_ids_in_db
from src.retrieve_data_from_local_db import find_met_name_from_bigg_id_in_db
from src.retrieve_data_from_local_db import retrieve_metabolite_info_from_db


#compounds = ["ethanol","acetaldehyde"]
#[selected_kegg_rxn, selected_bigg_rxn, bigg_string, bigg_name, bigg_mets, kegg_name, kegg_mets, kegg_string, metabolite_dict] = get_reaction(compounds)

def display_compound_length_warning(compounds):
    if len(compounds) != 2:
        print('\n WARNING! \n')
        print('This code only works for finding reactions between TWO compounds \n')
    return 

def get_rxn_info(compounds,metabolite_dict,bigg2kegg_rxn_db, bigg2kegg_db,mmdb_path):
    r1 = metabolite_dict[compounds[0]]['kegg_reactions']
    r2 = metabolite_dict[compounds[1]]['kegg_reactions']
    common_rxns = [value for value in r1 if value in r2]
    display_rxns = display_common_rxns(common_rxns)
    display_rxns = np.vstack(display_rxns)
    print("The following KEGG reactions were returned. Please select the correct reaction.")
    print(display_rxns)
    user_input = input("Enter the number of the desired KEGG reaction: ")
    kegg_rxn = common_rxns[int(user_input)]
    [kegg_name, kegg_mets, kegg_string, kegg_ec,kegg_eqn] = get_kegg_rxn_info(kegg_rxn)
    bigg_ids = [s for s in bigg2kegg_rxn_db if kegg_rxn in s[1]]
    if bigg_ids == []:
        [bigg_rxn,bigg_name,bigg_mets,bigg_string,metabolite_dict] = create_BIGG_rxn_from_KEGG_rxn(kegg_rxn,kegg_name,kegg_mets,kegg_eqn,bigg2kegg_db,metabolite_dict,mmdb_path)  
    else:
        display_rxn_strings = display_bigg_rxn_strings(bigg_ids)
        print("The following BIGG reactions were returned. Please select the correct reaction.")
        print(display_rxn_strings)
        print('//////////////////')
        print('NOTE: IF DESIRED RXN IS NOT IN LIST, ENTER: "n"')
        print('//////////////////')
        user_input = input("Enter the number of the desired BIGG reaction: ")
        if user_input is 'n':
            [bigg_rxn,bigg_name,bigg_mets,bigg_string,metabolite_dict] = create_BIGG_rxn_from_KEGG_rxn(kegg_rxn,kegg_name,kegg_mets,kegg_eqn,bigg2kegg_db,metabolite_dict,mmdb_path)  
        else:
            bigg_rxn = bigg_ids[int(user_input)][0]
            [bigg_string, bigg_name, bigg_mets] = get_bigg_rxn_info(bigg_rxn)
            metabolite_dict = add_rxn_mets_to_metabolite_dict(metabolite_dict,bigg_mets,bigg2kegg_db,mmdb_path)
    return bigg_rxn,kegg_rxn,bigg_string,bigg_name,bigg_mets,kegg_name,kegg_mets,kegg_string,kegg_eqn, kegg_ec
   
def add_to_reaction_dict(reaction_dict,compounds,metabolite_dict,bigg2kegg_rxn_db, bigg2kegg_db,mmdb_path):
    display_compound_length_warning(compounds)
    [bigg_rxn,kegg_rxn,bigg_string,bigg_name,bigg_mets,kegg_name,kegg_mets,kegg_string,kegg_eqn,kegg_ec] = get_rxn_info(compounds,metabolite_dict,bigg2kegg_rxn_db, bigg2kegg_db,mmdb_path)
    reaction_dict['reactions'].append(bigg_rxn)
    reaction_dict[bigg_rxn] = {}
    reaction_dict[bigg_rxn]['bigg_rxn'] = bigg_rxn
    reaction_dict[bigg_rxn]['kegg_rxn'] = kegg_rxn
    reaction_dict[bigg_rxn]['bigg_string'] = bigg_string
    reaction_dict[bigg_rxn]['bigg_name'] = bigg_name
    reaction_dict[bigg_rxn]['bigg_mets'] = bigg_mets
    reaction_dict[bigg_rxn]['kegg_name'] = kegg_name
    #reaction_dict[bigg_rxn]['kegg_mets'] = kegg_mets
    reaction_dict[bigg_rxn]['kegg_definition'] = kegg_string
    reaction_dict[bigg_rxn]['kegg_eqn'] = kegg_eqn
    reaction_dict[bigg_rxn]['BRENDA_enzymes'] = kegg_ec
    reaction_dict[bigg_rxn]['nodes'] = compounds
    return reaction_dict, metabolite_dict

def create_BIGG_rxn_from_KEGG_rxn(kegg_rxn,kegg_name,kegg_mets,kegg_eqn,bigg2kegg_db,metabolite_dict,mmdb_path):
    print(' \n ///// WARNING ///////')
    print('No BiGG reactions found for the selected KEGG reaction!')
    print('Creating a BIGG-like reaction ...')
    bigg_name = kegg_name
    bigg_mets = get_bigg_mets_from_kegg_mets(kegg_mets,bigg2kegg_db)
    metabolite_dict = add_rxn_mets_to_metabolite_dict(metabolite_dict,bigg_mets,bigg2kegg_db,mmdb_path)
    bigg_string = kegg2bigg_rxn_string(kegg_eqn,kegg_mets,bigg_mets)
    print('Please name the following BIGG reaction for',kegg_name)
    bigg_rxn = input('Enter BIGG ID:')
    return bigg_rxn,bigg_name,bigg_mets,bigg_string,metabolite_dict

def kegg2bigg_rxn_string(kegg_eqn,kegg_mets,bigg_mets):
    bigg_string = kegg_eqn
    for i in range(0,len(kegg_mets)):
        km = kegg_mets[i]
        bm = bigg_mets[i]+'_c'
        bigg_string = bigg_string.replace(km,bm)
    return bigg_string

def get_bigg_mets_from_kegg_mets(kegg_mets,bigg2kegg_db):
    bigg_mets = []
    for km in kegg_mets:
        bigg_ids = [s[0] for s in bigg2kegg_db if km in s[1]]
        if bigg_ids == []:
            print('///// WARNING //////')
            print('The KEGG compound,',km,'could be identified with any BIGG metabolite.')
            bigg_mets.append('NA')
        else:
            bigg_id = decide_compartment(bigg_ids)[0][:-2]
            bigg_mets.append(bigg_id)
    return bigg_mets
    
def decide_compartment(bigg_ids):
    bigg_id = []
    for ids in bigg_ids:
        if ids[-1] is 'c':
            bigg_id.append(ids)
    return bigg_id

def create_local_metabolite_conversion_table(metabolite_dict):
    local_met_conversion_table = []
    for met in metabolite_dict['metabolites']:
        local_met_conversion_table.append([met,metabolite_dict[met]['kegg_id'],metabolite_dict[met]['bigg_id']])
    return np.vstack(local_met_conversion_table)

def add_rxn_mets_to_metabolite_dict(metabolite_dict,bigg_mets,bigg2kegg_db,mmdb_path):
    bigg_ids_in_met_dict = get_bigg_ids_in_metabolite_dict(metabolite_dict)
    mets = list(set(bigg_mets) - set(bigg_ids_in_met_dict))
    db_met_ids = get_list_of_met_bigg_ids_in_db(mmdb_path)
    db_met_ids_wo_c = remove_compartments_from_a_list_of_bigg_ids(db_met_ids)
    for met in mets:
        if met in db_met_ids_wo_c:
            print(', \n ',met,' already has an entry in your local database.')
            print('Would you like to retrieve information from your database?')
            user_input = input('Enter "y" for yes or "n" for no: ')
            if user_input == "y":
                met_name = find_met_name_from_bigg_id_in_db(met+'_c',mmdb_path)
                metabolite_dict = retrieve_metabolite_info_from_db(met_name,mmdb_path,metabolite_dict)
            elif user_input == "n":
                print('adding',met,'to the metabolite dictionary')
                [bigg_charges,bigg_formula] = get_bigg_metabolite_information(met+'_c')
                kegg_id = [s for s in bigg2kegg_db if s[0] == met+'_c'][0][1]
                [cpd_id, kegg_formula, kegg_rxns, name,kegg_dblinks] = get_kegg_compound_information([kegg_id])
                charge = define_metabolite_charge(bigg_charges)
                charged_formula = get_charged_formula(kegg_formula,charge)
                metabolite_dict['metabolites'].append(name)
                metabolite_dict[name] = {}
                metabolite_dict[name]['kegg_id'] = kegg_id
                metabolite_dict[name]['neutral_formula'] = kegg_formula
                metabolite_dict[name]['kegg_reactions'] = kegg_rxns
                metabolite_dict[name]['bigg_id'] = met+'_c'
                metabolite_dict[name]['bigg_charges'] = bigg_charges
                metabolite_dict[name]['bigg_formula'] = bigg_formula
                metabolite_dict[name]['charge'] = charge
                metabolite_dict[name]['kegg_dblinks'] = kegg_dblinks
                metabolite_dict[name]['formula'] = charged_formula
        else:
            print('adding',met,'to the metabolite dictionary')
            [bigg_charges,bigg_formula] = get_bigg_metabolite_information(met+'_c')
            kegg_id = [s for s in bigg2kegg_db if s[0] == met+'_c'][0][1]
            [cpd_id, kegg_formula, kegg_rxns, name,kegg_dblinks] = get_kegg_compound_information([kegg_id])
            charge = define_metabolite_charge(bigg_charges)
            charged_formula = get_charged_formula(kegg_formula,charge)
            metabolite_dict['metabolites'].append(name)
            metabolite_dict[name] = {}
            metabolite_dict[name]['kegg_id'] = kegg_id
            metabolite_dict[name]['neutral_formula'] = kegg_formula
            metabolite_dict[name]['kegg_reactions'] = kegg_rxns
            metabolite_dict[name]['bigg_id'] = met+'_c'
            metabolite_dict[name]['bigg_charges'] = bigg_charges
            metabolite_dict[name]['bigg_formula'] = bigg_formula
            metabolite_dict[name]['charge'] = charge
            metabolite_dict[name]['kegg_dblinks'] = kegg_dblinks
            metabolite_dict[name]['formula'] = charged_formula
    return metabolite_dict

def get_bigg_ids_in_metabolite_dict(metabolite_dict):
    metabolites = metabolite_dict['metabolites']
    bigg_ids_in_met_dict = []
    for met in metabolites:
        bigg_ids_in_met_dict.append(metabolite_dict[met]['bigg_id'][:-2])
    return bigg_ids_in_met_dict

def get_kegg_rxn_info(kegg_rxn):
    # get name, mets, string, enzyme class
    base = "http://rest.kegg.jp/get/"
    url = base + kegg_rxn
    response = requests.get(url)
    tsv_data = response.text
    data = re.split('\n',tsv_data)
    name = re.split('NAME|  ',data[1])
    name = [s for s in name if s][0]
    defn_index =  [i for i, s in enumerate(data) if 'DEFINITION' in s][0]
    defn = data[defn_index]
    defn = re.split('DEFINITION  ',defn)
    defn = [s for s in defn if s][0]
    eqn_index = [i for i, s in enumerate(data) if 'EQUATION' in s][0]
    eqn = data[eqn_index]
    eqn = re.split('EQUATION|  ',eqn)
    eqn = [s for s in eqn if s][0]
    mets = eqn.replace(' <=> ',',')
    mets = mets.replace(' + ',',')
    mets = re.split(',| ',mets)
    mets = [i for i in mets if len(i) == 6]
    ec = [s for s in data if 'ENZYME' in s][0]
    ec = re.split('ENZYME| ',ec)
    ec = [s for s in ec if s]
    return name, mets, defn, ec, eqn

def get_bigg_rxn_info(selected_bigg_rxn):
    # get name, string, metabolites, EC numbers
    base = 'http://bigg.ucsd.edu/api/v2/universal/reactions/'
    url = base + selected_bigg_rxn
    response = requests.get(url)
    tsv_data = response.text
    data = re.split('"reaction_string": |"metabolites": |"database_links": |"old_identifiers: ',tsv_data)
    rxn_string = re.split(',',data[1])[0]
    rxn_string = rxn_string.strip('"')
    rxn_string = rxn_string.replace('&#8652;','<=>')
    mets = re.split('"bigg_id": |, "name": ',data[2])
    name = re.split('"name": |, "pseudoreaction":',data[3])[1]
    name = name.strip('"')
    metabolites = []
    for num in range(1,len(mets)-1,2):
        metabolites = metabolites + [mets[num].strip('"')]
    return rxn_string, name, metabolites

def display_bigg_rxn_strings(bigg_ids):
    base = 'http://bigg.ucsd.edu/api/v2/universal/reactions/'
    i = 0
    display_rxn_strings = []
    for b in bigg_ids:
        url = base + b[0]
        response = requests.get(url)
        tsv_data = response.text
        data = re.split('"reaction_string":',tsv_data)
        rxn_string = re.split(',',data[1])
        display_rxn_strings.append([i,b[0],rxn_string[0]])
        i = i + 1
    return np.vstack(display_rxn_strings)

def display_common_rxns(common_rxns):
    base = "http://rest.kegg.jp/get/"
    i = 0
    display_rxns = []
    for rxn in common_rxns:
        url = base + rxn
        response = requests.get(url)
        tsv_data = response.text
        data = re.split('\n',tsv_data)
        defn = [s for s in data if "DEFINITION" in s]
        defn = re.split('DEFINITION ',defn[0])
        defn = [s for s in defn if s][0]
        display_rxns.append([i,defn])
        i = i + 1
    return display_rxns

def remove_compartments_from_a_list_of_bigg_ids(list_of_bigg_ids):
    ids_wo_c = []
    for bigg_id in list_of_bigg_ids:
        ids_wo_c.append(bigg_id[:-2])
    return ids_wo_c
        


    
