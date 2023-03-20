# retrieve data from database

from src.balance_check import findMetNameFromID
#from src.create_pathways import  get_all_dict_paths_between_mets
import sqlite3
from src.parseRxnString import parseRxnString
import os

def get_db_paths_between_mets(mmdb_path,target,precursor,metabolite_dict):
    db_path = os.path.join(mmdb_path,r"local_db")
    os.chdir(db_path)
    conn = sqlite3.connect('mmdb_local.db')
    selectTable = "SELECT * FROM " + "pathways"
    cursor = conn.execute(selectTable)
    rows = cursor.fetchall()
    precursor_id = metabolite_dict[precursor]['bigg_id']
    target_id = metabolite_dict[target]['bigg_id']
    db_paths = []
    for row in rows:
        if precursor_id == row[2] and target_id == row[1]:
            db_paths.append(row[0])
    return db_paths

def get_db_pathway_rxns_from_id(mmdb_path,path_id):
    db_path = os.path.join(mmdb_path,r"local_db")
    os.chdir(db_path)
    conn = sqlite3.connect('mmdb_local.db')
    selectTable = "SELECT * FROM " + "pathway_reactions"
    cursor = conn.execute(selectTable)
    rows = cursor.fetchall()
    for row in rows:
        if row[0] == path_id:
            path_rxns = []
            for rxn in row:
                if rxn != path_id and rxn != None:
                    path_rxns.append(rxn)
    return path_rxns

def add_db_pathways_between_mets_to_pathway_dict(target,precursor,mmdb_path,pathway_dict,metabolite_dict,reaction_dict):
    db_paths = get_db_paths_between_mets(mmdb_path,target,precursor,metabolite_dict)
    dict_paths = get_all_dict_paths_between_mets(target,precursor,pathway_dict,metabolite_dict)
    for path in db_paths:
        if path not in dict_paths:
            [metabolite_dict,reaction_dict,pathway_dict] = get_all_info_from_db_path(path,pathway_dict,reaction_dict,metabolite_dict,mmdb_path)
    return pathway_dict,reaction_dict,metabolite_dict

def retrieve_metabolite_info_from_db(met,mmdb_path,metabolite_dict):
    db_path = os.path.join(mmdb_path,r"local_db")
    os.chdir(db_path)
    conn = sqlite3.connect('mmdb_local.db')
    selectTable = "SELECT * FROM " + "metabolites"
    cursor = conn.execute(selectTable)
    rows = cursor.fetchall()
    for row in rows:
        if row[0] == met:
            metabolite_dict[met] = {}
            metabolite_dict[met]['kegg_id'] = row[4]
            metabolite_dict[met]['neutral_formula'] = row[5]
            metabolite_dict[met]['kegg_reactions'] = retrieve_kegg_reactions_for_met(met,mmdb_path)
            metabolite_dict[met]['bigg_id'] = row[1]
            metabolite_dict[met]['bigg_formula'] = row[6]
            metabolite_dict[met]['formula'] = row[2]
            metabolite_dict[met]['charge'] = row[3]
            metabolite_dict['metabolites'].append(met)
    return metabolite_dict

def retrieve_kegg_reactions_for_met(met,mmdb_path):
    db_path = os.path.join(mmdb_path,r"local_db")
    os.chdir(db_path)
    conn = sqlite3.connect('mmdb_local.db')
    selectTable = "SELECT * FROM " + "kegg_reactions"
    cursor = conn.execute(selectTable)
    rows = cursor.fetchall()
    indices = get_db_indices(met,"kegg_reactions",mmdb_path)
    kegg_reactions = []
    for index in indices:
        for rxn in rows[index]:
            if rxn != met and rxn != None:
                kegg_reactions.append(rxn)
    return kegg_reactions

def get_db_indices(entry_name,table,mmdb_path):
    db_path = os.path.join(mmdb_path,r"local_db")
    os.chdir(db_path)
    conn = sqlite3.connect('mmdb_local.db')
    selectTable = "SELECT * FROM " + table
    cursor = conn.execute(selectTable)
    rows = cursor.fetchall()
    i = 0
    indices = []
    for row in rows:
        if row[0] == entry_name:
            indices.append(i)
        i = i + 1
    return indices

def retrive_brenda_reactions(rxn,mmdb_path):
    db_path = os.path.join(mmdb_path,r"local_db")
    os.chdir(db_path)
    conn = sqlite3.connect('mmdb_local.db')
    selectTable = "SELECT * FROM " + "brenda_reactions"
    cursor = conn.execute(selectTable)
    rows = cursor.fetchall()
    brenda_rxns = []
    for row in rows:
        if row[0] == rxn:
            for br in row:
                if br != rxn and br != None:
                    brenda_rxns.append(br)
    return brenda_rxns

def retrieve_rxn_info_from_db(rxn,reaction_dict,mmdb_path):
    db_path = os.path.join(mmdb_path,r"local_db")
    os.chdir(db_path)
    conn = sqlite3.connect('mmdb_local.db')
    selectTable = "SELECT * FROM " + "reactions"
    cursor = conn.execute(selectTable)
    rows = cursor.fetchall()
    for row in rows:
        if row[0] == rxn:
            reaction_dict['reactions'].append(rxn)
            reaction_dict[rxn] = {}
            reaction_dict[rxn]['bigg_rxn'] = row[0]
            reaction_dict[rxn]['kegg_rxn'] = row[3]
            reaction_dict[rxn]['bigg_string'] = row[2]
            reaction_dict[rxn]['bigg_name'] = row[1]
            reaction_dict[rxn]['kegg_name'] = row[4]
            reaction_dict[rxn]['kegg_definition'] = row[5]
            reaction_dict[rxn]['BRENDA_enzymes'] = retrive_brenda_reactions(rxn,mmdb_path)
            reaction_dict[rxn]['nodes'] = [row[6], row[7]]
            reaction_dict[rxn]['mass_balance'] = row[8]
            reaction_dict[rxn]['charge_balance'] = row[9]
    return reaction_dict

def get_all_info_from_db_path(path,pathway_dict,reaction_dict,metabolite_dict,mmdb_path):
    if path not in pathway_dict['pathways']:
        pathway_dict = add_path_from_db_to_dict(path,pathway_dict,mmdb_path)
        pathway_rxns = get_db_pathway_rxns_from_id(mmdb_path,path)
        for rxn in pathway_rxns:
            if rxn not in reaction_dict['reactions']:
                reaction_dict = retrieve_rxn_info_from_db(rxn,reaction_dict,mmdb_path)
                rxn_string = reaction_dict[rxn]['bigg_string']
                [mets,stoich,rev] = parseRxnString(rxn_string)
                for met_id in mets:
                    list_of_ids = get_list_of_met_bigg_ids_in_metabolite_dict(metabolite_dict)
                    if met_id not in list_of_ids:
                        met = find_met_name_from_bigg_id_in_db(met_id,mmdb_path)
                        metabolite_dict = retrieve_metabolite_info_from_db(met,mmdb_path,metabolite_dict)
    return metabolite_dict,reaction_dict,pathway_dict

def add_all_info_of_paths_that_make_target(target,mmdb_path,metabolite_dict,reaction_dict,pathway_dict):
    target_paths = get_pathways_that_make_target(target,mmdb_path,metabolite_dict)
    for path in target_paths:
        [metabolite_dict,reaction_dict,pathway_dict] = get_all_info_from_db_path(path,pathway_dict,reaction_dict,metabolite_dict,mmdb_path)
    return metabolite_dict,reaction_dict,pathway_dict

def get_pathways_that_make_target(target,mmdb_path,metabolite_dict):
    db_path = os.path.join(mmdb_path,r"local_db")
    os.chdir(db_path)
    conn = sqlite3.connect('mmdb_local.db')
    selectTable = "SELECT * FROM " + "pathways"
    cursor = conn.execute(selectTable)
    rows = cursor.fetchall()
    target_id = metabolite_dict[target]['bigg_id']
    paths = []
    for row in rows:
        if row[1] == target_id:
            paths.append(row[0])
    return paths
        
def add_path_from_db_to_dict(path,pathway_dict,mmdb_path):
    db_path = os.path.join(mmdb_path,r"local_db")
    os.chdir(db_path)
    conn = sqlite3.connect('mmdb_local.db')
    selectTable = "SELECT * FROM " + "pathways"
    cursor = conn.execute(selectTable)
    rows = cursor.fetchall()
    for row in rows:
        if row[0] == path:
            pathway_dict['pathways'].append(path)
            pathway_dict[path] = {}
            pathway_dict[path]['id'] = path
            pathway_dict[path]['name'] = path
            pathway_dict[path]['rxns'] = get_db_pathway_rxns_from_id(mmdb_path,path)
            pathway_dict[path]['product_id'] = row[1]
            pathway_dict[path]['precursor_ids'] = row[2]
    return pathway_dict

def get_list_of_met_bigg_ids_in_metabolite_dict(metabolite_dict):
    list_of_ids = []
    for met in metabolite_dict['metabolites']:
        list_of_ids.append(metabolite_dict[met]['bigg_id'])
    return list_of_ids

def find_met_name_from_bigg_id_in_db(met_id,mmdb_path):
    db_path = os.path.join(mmdb_path,r"local_db")
    os.chdir(db_path)
    conn = sqlite3.connect('mmdb_local.db')
    selectTable = "SELECT * FROM " + "metabolites"
    cursor = conn.execute(selectTable)
    rows = cursor.fetchall()
    met_name = []
    for row in rows:
        if row[1] == met_id:
            met_name.append(row[0])
    return met_name[0]
    
def get_list_of_mets_in_db(mmdb_path):
    db_path = os.path.join(mmdb_path,r"local_db")
    os.chdir(db_path)
    conn = sqlite3.connect('mmdb_local.db')
    selectTable = "SELECT * FROM " + "metabolites"
    cursor = conn.execute(selectTable)
    rows = cursor.fetchall()
    db_mets = []
    for row in rows:
        db_mets.append(row[0])
    return db_mets

def get_list_of_met_bigg_ids_in_db(mmdb_path):
    db_path = os.path.join(mmdb_path,r"local_db")
    os.chdir(db_path)
    conn = sqlite3.connect('mmdb_local.db')
    selectTable = "SELECT * FROM " + "metabolites"
    cursor = conn.execute(selectTable)
    rows = cursor.fetchall()
    db_met_ids = []
    for row in rows:
        db_met_ids.append(row[1])
    return db_met_ids

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

def get_all_metabolites_from_db(mmdb_path,metabolite_dict):
    db_path = os.path.join(mmdb_path,r"local_db")
    os.chdir(db_path)
    conn = sqlite3.connect('mmdb_local.db')
    selectTable = "SELECT * FROM " + "metabolites"
    cursor = conn.execute(selectTable)
    rows = cursor.fetchall()
    for row in rows:
        met = row[0]
        metabolite_dict[met] = {}
        metabolite_dict[met]['kegg_id'] = row[4]
        metabolite_dict[met]['neutral_formula'] = row[5]
        metabolite_dict[met]['kegg_reactions'] = retrieve_kegg_reactions_for_met(row[0],mmdb_path)
        metabolite_dict[met]['bigg_id'] = row[1]
        metabolite_dict[met]['bigg_formula'] = row[6]
        metabolite_dict[met]['formula'] = row[2]
        metabolite_dict[met]['charge'] = row[3]
        metabolite_dict['metabolites'].append(row[0])
    return metabolite_dict

def get_all_reactions_from_db(mmdb_path,reaction_dict):
    db_path = os.path.join(mmdb_path,r"local_db")
    os.chdir(db_path)
    conn = sqlite3.connect('mmdb_local.db')
    selectTable = "SELECT * FROM " + "reactions"
    cursor = conn.execute(selectTable)
    rows = cursor.fetchall()
    for row in rows:
        rxn = row[0]
        reaction_dict['reactions'].append(row[0])
        reaction_dict[rxn] = {}
        reaction_dict[rxn]['bigg_rxn'] = row[0]
        reaction_dict[rxn]['kegg_rxn'] = row[3]
        reaction_dict[rxn]['bigg_string'] = row[2]
        reaction_dict[rxn]['bigg_name'] = row[1]
        reaction_dict[rxn]['kegg_name'] = row[4]
        reaction_dict[rxn]['kegg_definition'] = row[5]
        reaction_dict[rxn]['BRENDA_enzymes'] = retrive_brenda_reactions(row[0],mmdb_path)
        reaction_dict[rxn]['nodes'] = [row[6], row[7]]
        reaction_dict[rxn]['mass_balance'] = row[8]
        reaction_dict[rxn]['charge_balance'] = row[9]
    return reaction_dict





