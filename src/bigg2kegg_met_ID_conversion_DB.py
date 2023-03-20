# create a namespace database to interrelate bigg_id's and kegg_id's
# load BIGG metabolites
import requests
from bs4 import BeautifulSoup
import re
import numpy as np

#bigg2kegg_db = create_bigg2kegg_DB()

def create_bigg2kegg_DB():
    url = "http://bigg.ucsd.edu/static/namespace/bigg_models_metabolites.txt"
    response = requests.get(url)
    tsv_data = response.text
    data = re.split('\n',tsv_data)
    bigg2kegg_db = []
    i = 0
    for met in data:
        i = i+1
        met = re.split('\t',met)
        bigg_id = met[0]
        if len(met) > 4:
            links = met[4]
            links = re.split(';',links)
            kegg_link = [s for s in links if "KEGG" in s]
            kegg_id = if_kegg_is_available(kegg_link)
        else:
            kegg_id = 'NA'
        bigg2kegg_db.append([bigg_id,kegg_id])
    return bigg2kegg_db

def if_kegg_is_available(kegg_link):
    if not kegg_link:
        kegg_id = 'NA'
    else:
        kegg_id = kegg_link[0][-6:]
    return kegg_id
