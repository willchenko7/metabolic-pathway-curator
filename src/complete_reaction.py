import requests
from bs4 import BeautifulSoup
import re

compounds = ["ethanol","acetaldehyde"]
metabolite_dict = create_metabolite_dict(compounds,bigg2kegg_db)

def create_metabolite_dict(compounds,bigg2kegg_db):
    metabolite_dict = {}
    metabolite_dict['metabolites'] = compounds
    for compound in compounds:
        metabolite_dict[compound] = {}
        [kegg_id, kegg_formula, kegg_reactions, bigg_id, bigg_charge, bigg_formula] = get_compound_information(compound,bigg2kegg_db)
        metabolite_dict[compound]['kegg_id'] = kegg_id
        metabolite_dict[compound]['kegg_formula'] = kegg_formula
        metabolite_dict[compound]['kegg_reactions'] = kegg_reactions
        metabolite_dict[compound]['bigg_id'] = bigg_id
        metabolite_dict[compound]['bigg_charge'] = bigg_charge
        metabolite_dict[compound]['formula'] = formula
    return metabolite_dict
    
def get_compound_information(compound,bigg2kegg_db):
    kegg_compound = get_kegg_compound(compound)
    [kegg_id,kegg_formula,kegg_reactions] = get_kegg_compound_information(kegg_compound)
    bigg_ids = [s for s in bigg2kegg_db if kegg_id in s[1]]
    bigg_id = decide_compartment(bigg_ids)[0]
    [bigg_charge,bigg_formula] = get_bigg_metabolite_information(bigg_id)
    return kegg_id, kegg_formula, kegg_reactions, bigg_id, bigg_charge, bigg_formula 

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
    bigg_charge = [s for s in data if "charges" in s]
    bigg_charge = re.split(': ',bigg_charge[0])[1]
    bigg_charge = bigg_charge.lstrip('[')
    bigg_charge = bigg_charge.rstrip(']')
    return bigg_charge,bigg_formula

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
    formula = [s for s in d if "FORMULA" in s]
    formula = re.split('FORMULA| ',formula[0])
    formula = [s for s in formula if s]
    reactions = [s for s in d if "REACTION" in s]
    reactions = re.split('REACTION| ',reactions[0])
    reactions = [s for s in reactions if s]
    return cpd_id, formula, reactions

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
