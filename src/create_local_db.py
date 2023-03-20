# create local db

import sqlite3
import os
import re

create_local_db(mmdb_path)

def create_local_db(mmdb_path):
    # creates tables to store information from mmdb
    # WARNING: this code should only be run once, to avoid overwriting info
    # once this code has been run, you can add to the tables by using " "
    db_path = os.path.join(mmdb_path,r"local_db")
    os.chdir(db_path)
    # create tables
    conn = sqlite3.connect('mmdb_local.db')
    c = conn.cursor()
    c.execute("""CREATE TABLE metabolites (
                name text,
                bigg_id text,
                formula blob,
                charge integer,
                kegg_id text,
                neutral_formula blob,
                bigg_formula blob,
                date_of_entry text
                )""")
    c.execute("""CREATE TABLE reactions (
                bigg_id text,
                bigg_name text,
                bigg_string text,
                kegg_id text,
                kegg_name text,
                kegg_defn text,
                node_1 text,
                node_2 text,
                mass_balance text,
                charge_balance text,
                date_of_entry text
                )""")
    c.execute("""CREATE TABLE pathways (
                id text,
                target text,
                precursor text,
                date_of_entry text
                )""")
    create_kegg_reactions_db()
    create_brenda_reactions_db()
    create_pathway_reactions_db()
    conn.close()
    return

def create_kegg_reactions_db():
    conn = sqlite3.connect('mmdb_local.db')
    c = conn.cursor()
    c.execute("""CREATE TABLE kegg_reactions (
            metabolite_name text
            )""")
    for i in range(1,997):
        addColumn = "ALTER TABLE kegg_reactions ADD COLUMN r" + str(i) + " text"
        c.execute(addColumn)
    conn.close()
    return

def create_brenda_reactions_db():
    conn = sqlite3.connect('mmdb_local.db')
    c = conn.cursor()
    c.execute("""CREATE TABLE brenda_reactions (
            reaction text
            )""")
    for i in range(1,10):
        addColumn = "ALTER TABLE brenda_reactions ADD COLUMN enzyme" + str(i) + " text"
        c.execute(addColumn)
    conn.close()
    return

def create_pathway_reactions_db():
    conn = sqlite3.connect('mmdb_local.db')
    c = conn.cursor()
    c.execute("""CREATE TABLE pathway_reactions (
            pathway text
            )""")
    for i in range(1,50):
        addColumn = "ALTER TABLE pathway_reactions ADD COLUMN rxn" + str(i) + " text"
        c.execute(addColumn)
    conn.close()
    return 





