# create conv table from bigg2kegg for reaction IDs
import requests
from bs4 import BeautifulSoup
import re
import numpy as np

#bigg2kegg_rxn_db = create_bigg2kegg_rxn_db()

def create_bigg2kegg_rxn_db():
    url = 'http://bigg.ucsd.edu/static/namespace/bigg_models_reactions.txt'
    response = requests.get(url)
    tsv_data = response.text
    data = re.split('\n',tsv_data)
    bigg2kegg_rxn_db = []
    for rxn in data:
        rxn = re.split('\t',rxn)
        bigg_id = rxn[0]
        if len(rxn) > 4:
            links = rxn[4]
            links = re.split(';',links)
            kegg_link = [s for s in links if "KEGG" in s]
            kegg_id = if_kegg_is_available(kegg_link)
        else:
            kegg_id = 'NA'
        bigg2kegg_rxn_db.append([bigg_id,kegg_id])
    return bigg2kegg_rxn_db

def if_kegg_is_available(kegg_link):
    if not kegg_link:
        kegg_id = 'NA'
    else:
        kegg_id = kegg_link[0][-6:]
    return kegg_id
