# add to local db

import sqlite3
import os
import re
from datetime import date

#add_all_data_to_local_db(mmdb_path,metabolite_dict,reaction_dict,pathway_dict)

def find_names_in_table(mmdb_path,table):
    db_path = os.path.join(mmdb_path,r"local_db")
    os.chdir(db_path)
    conn = sqlite3.connect('mmdb_local.db')
    selectTable = "SELECT * FROM " + table
    cursor = conn.execute(selectTable)
    rows = cursor.fetchall()
    names = []
    for row in rows:
        names.append(row[0])
    return names

def add_all_data_to_local_db(mmdb_path,metabolite_dict,reaction_dict,pathway_dict):
    add_metabolites_to_local_db(mmdb_path,metabolite_dict)
    add_to_kegg_reactions_db(mmdb_path,metabolite_dict)
    add_reactions_to_local_db(mmdb_path,reaction_dict)
    add_brenda_reactions(mmdb_path,reaction_dict)
    add_pathways_to_local_db(mmdb_path,pathway_dict)
    add_pathway_reactions(mmdb_path,pathway_dict)
    return 

def add_pathways_to_local_db(mmdb_path,pathway_dict):
    db_path = os.path.join(mmdb_path,r"local_db")
    os.chdir(db_path)
    # create tables
    conn = sqlite3.connect('mmdb_local.db')
    c = conn.cursor()
    paths_in_db = find_names_in_table(mmdb_path,'pathways')
    for path in pathway_dict['pathways']:
        if path in paths_in_db:
            print('Not adding ',path,' to pathways table because their is already an entry with this name')
        else:
            print('Adding ',path,' to pathway table')
            target = pathway_dict[path]['product_id']
            precursor = pathway_dict[path]['precursor_ids']
            today = date.today()
            date_of_entry = today.strftime("%d/%m/%Y")
            c.execute("INSERT INTO pathways VALUES(?,?,?,?)",(path,target,precursor,date_of_entry))
            conn.commit()
    conn.close()
    return

def add_pathway_reactions(mmdb_path,pathway_dict):
    db_path = os.path.join(mmdb_path,r"local_db")
    os.chdir(db_path)
    conn = sqlite3.connect('mmdb_local.db')
    c = conn.cursor()
    qstr = create_string_of_question_marks(50)
    paths_in_db = find_names_in_table(mmdb_path,'pathway_reactions')
    for path in pathway_dict['pathways']:
        if path in paths_in_db:
            print('Not adding ',path,' to pathway_reactions table because their is already an entry with this name')
        else:
            print('Adding ',path,' to pathway_reactions table')
            db_path_rxns = create_pathway_reactions_for_db(path,50,pathway_dict)
            c.execute("INSERT INTO pathway_reactions VALUES" + qstr,db_path_rxns)
            conn.commit()
    conn.close()
    return 
        
def add_brenda_reactions(mmdb_path,reaction_dict):
    db_path = os.path.join(mmdb_path,r"local_db")
    os.chdir(db_path)
    conn = sqlite3.connect('mmdb_local.db')
    c = conn.cursor()
    rxns_in_db = find_names_in_table(mmdb_path,'brenda_reactions')
    for rxn in reaction_dict['reactions']:
        if rxn in rxns_in_db:
            print('Not adding ',rxn,' to brenda_reactions table because their is already an entry with this name')
            print('You can manually edit entries in the table with ...')
        else:
            db_brenda = create_brenda_for_db(rxn,10,reaction_dict)
            c.execute("INSERT INTO brenda_reactions VALUES(?,?,?,?,?,?,?,?,?,?)",(db_brenda))
            conn.commit()
    conn.close()
    return

def create_pathway_reactions_for_db(path,db_size,pathway_dict):
    rxns = pathway_dict[path]['rxns']
    num_none = db_size - len(rxns) - 1
    null_vector = make_null_vector(num_none)
    db_path_rxns = [path] + rxns + null_vector
    return db_path_rxns
                  
def create_brenda_for_db(rxn,db_size,reaction_dict):
    brenda = reaction_dict[rxn]['BRENDA_enzymes']
    num_none = db_size - len(brenda) - 1
    null_vector = make_null_vector(num_none)
    db_brenda = [rxn] + brenda + null_vector
    return db_brenda
    
def add_reactions_to_local_db(mmdb_path,reaction_dict):
    db_path = os.path.join(mmdb_path,r"local_db")
    os.chdir(db_path)
    # create tables
    conn = sqlite3.connect('mmdb_local.db')
    c = conn.cursor()
    rxns_in_db = find_names_in_table(mmdb_path,'reactions')
    for rxn in reaction_dict['reactions']:
        if rxn in rxns_in_db:
            print('Not adding ',rxn,' to reactions table because their is already an entry with this name')
            print('You can manually edit entries in the table with ...')
        else:
            print('Adding ',rxn,' to reaction table')
            bigg_id = reaction_dict[rxn]['bigg_rxn']
            bigg_name = reaction_dict[rxn]['bigg_name']
            bigg_string = reaction_dict[rxn]['bigg_string']
            kegg_id = reaction_dict[rxn]['kegg_rxn']
            kegg_name = reaction_dict[rxn]['kegg_name']
            kegg_defn = reaction_dict[rxn]['kegg_definition']
            node_1 = reaction_dict[rxn]['nodes'][0]
            node_2 = reaction_dict[rxn]['nodes'][1]
            mass_balance = reaction_dict[rxn]['mass_balance']
            charge_balance = reaction_dict[rxn]['charge_balance']
            today = date.today()
            date_of_entry = today.strftime("%d/%m/%Y")
            c.execute("INSERT INTO reactions VALUES(?,?,?,?,?,?,?,?,?,?,?)",(bigg_id,bigg_name,bigg_string,kegg_id,kegg_name,kegg_defn,node_1,node_2,mass_balance,charge_balance,date_of_entry))
            conn.commit()
    conn.close()
    return 
        
def add_metabolites_to_local_db(mmdb_path,metabolite_dict):
    db_path = os.path.join(mmdb_path,r"local_db")
    os.chdir(db_path)
    conn = sqlite3.connect('mmdb_local.db')
    c = conn.cursor()
    mets_in_db = find_names_in_table(mmdb_path,'metabolites')  
    for met in metabolite_dict['metabolites']:
        if met in mets_in_db:
            print('Not adding ',met,' to metabolites table because their is already an entry with this name')
            print('You can manually edit entries in the table with ...')
        else:
            print('Adding ',met,' to metabolite table')
            bigg_id = metabolite_dict[met]['bigg_id']
            kegg_id = metabolite_dict[met]['kegg_id']
            neutral_formula = metabolite_dict[met]['neutral_formula'][0]
            bigg_formula = metabolite_dict[met]['bigg_formula']
            bigg_charges = metabolite_dict[met]['bigg_charges']
            formula = metabolite_dict[met]['formula']
            charge = metabolite_dict[met]['charge']
            today = date.today()
            date_of_entry = today.strftime("%d/%m/%Y")
            c.execute("INSERT INTO metabolites VALUES(?,?,?,?,?,?,?,?)",(met,bigg_id,formula,charge,kegg_id,neutral_formula,bigg_formula,date_of_entry))
            conn.commit()
    conn.close()
    return

def add_to_kegg_reactions_db(mmdb_path,metabolite_dict):
    db_path = os.path.join(mmdb_path,r"local_db")
    os.chdir(db_path)
    # create tables
    conn = sqlite3.connect('mmdb_local.db')
    c = conn.cursor()
    db_size = 997
    qstr = create_string_of_question_marks(db_size)
    mets_in_table = find_names_in_table(mmdb_path,'kegg_reactions')
    for met in metabolite_dict['metabolites']:
        if met in mets_in_table:
            print('Not adding ',met,' to kegg_reactions table because their is already an entry with this name')
            print('You can manually edit entries in the table with ...')
        else:
            print('Adding ',met,' to kegg_reactions table')
            kegg_reactions = metabolite_dict[met]['kegg_reactions']
            if len(kegg_reactions) >= 990 and len(kegg_reactions) <= 1980:
                [db_kg_1,db_kg_2] = create_kegg_reactions_for_db_for_large_sets(met,db_size,metabolite_dict)
                c.execute("INSERT INTO kegg_reactions VALUES" + qstr,db_kg_1)
                c.execute("INSERT INTO kegg_reactions VALUES" + qstr,db_kg_2)
            elif len(kegg_reactions)>= 1980 and len(kegg_reactions) <= 2970:
                [db_kg_1,db_kg_2,db_kg_3] = create_kegg_reactions_for_db_for_very_large_sets(met,db_size,metabolite_dict)
                c.execute("INSERT INTO kegg_reactions VALUES" + qstr,db_kg_1)
                c.execute("INSERT INTO kegg_reactions VALUES" + qstr,db_kg_2)
                c.execute("INSERT INTO kegg_reactions VALUES" + qstr,db_kg_3)
            elif len(kegg_reactions) <= 990:
                db_kg = create_kegg_reactions_for_db(met,db_size,metabolite_dict)
                c.execute("INSERT INTO kegg_reactions VALUES" + qstr,db_kg)
            else:
                print('WARNING! -, ',met,' was not added bc it has too many reactions from kegg.')
            conn.commit()
    conn.close()
    return

def create_string_of_question_marks(size):
    qstr = "("
    for i in range(0,size):
        #print(i)
        if i != size-1:
            qstr = qstr + "?,"
        else:
            qstr = qstr + "?"
    return qstr + ")"

def create_kegg_reactions_for_db(met,db_size,metabolite_dict):
    kegg_reactions = metabolite_dict[met]['kegg_reactions']
    num_none = db_size - len(kegg_reactions) - 1
    null_vector = make_null_vector(num_none)
    db_kg = [met] + kegg_reactions + null_vector
    return db_kg

def create_kegg_reactions_for_db_for_very_large_sets(met,db_size,metabolite_dict):
    kegg_reactions = metabolite_dict[met]['kegg_reactions']
    div_mark_1 = int(len(kegg_reactions)*(1/3))
    div_mark_2 = int(len(kegg_reactions)*(2/3))
    kg_1 = kegg_reactions[0:div_mark_1]
    kg_2 = kegg_reactions[div_mark_1+1:div_mark_2]
    kg_3 = kegg_reactions[div_mark_2+1:len(kegg_reactions)]
    num_none_1 = db_size - len(kg_1) - 1
    num_none_2 = db_size - len(kg_2) - 1
    num_none_3 = db_size - len(kg_3) - 1
    null_vector_1 = make_null_vector(num_none_1)
    null_vector_2 = make_null_vector(num_none_2)
    null_vector_3 = make_null_vector(num_none_3)
    db_kg_1 = [met] + kg_1 + null_vector_1
    db_kg_2 = [met] + kg_2 + null_vector_2
    db_kg_3 = [met] + kg_3 + null_vector_3
    return db_kg_1,db_kg_2,db_kg_3


def create_kegg_reactions_for_db_for_large_sets(met,db_size,metabolite_dict):
    kegg_reactions = metabolite_dict[met]['kegg_reactions']
    div_mark = int(len(kegg_reactions)/2)
    kg_1 = kegg_reactions[0:div_mark]
    kg_2 = kegg_reactions[div_mark+1:len(kegg_reactions)-1]
    num_none_1 = db_size - len(kg_1) - 1
    num_none_2 = db_size - len(kg_2) - 1
    null_vector_1 = make_null_vector(num_none_1)
    null_vector_2 = make_null_vector(num_none_2)
    db_kg_1 = [met] + kg_1 + null_vector_1
    db_kg_2 = [met] + kg_2 + null_vector_2
    return db_kg_1, db_kg_2

def make_null_vector(num_none):
    null_vector = []
    for i in range(0,num_none):
        null_vector.append(None)
    return null_vector
                   
