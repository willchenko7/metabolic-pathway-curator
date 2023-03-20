# mmdb example script 
import os
mmdb_path = r"C:\Users\Owner\Documents\PhD_stuff\research\mmdb"
os.chdir(mmdb_path)

from src.bigg2kegg_met_ID_conversion_DB import create_bigg2kegg_DB
from src.bigg2kegg_rxn_ID_conversion import create_bigg2kegg_rxn_db
from src.create_metabolite_dict import add_to_metabolite_dict
from src.get_reaction import add_to_reaction_dict
from src.create_pathways import add_to_pathway_dict
from src.export_to_CSV import create_new_problem
from src.export_to_CSV import add_to_pathway_csv_table
from src.export_to_CSV import add_to_reaction_csv_table
from src.export_to_CSV import add_to_metabolite_csv_table

# create conversion tables
bigg2kegg_db = create_bigg2kegg_DB()
bigg2kegg_rxn_db = create_bigg2kegg_rxn_db()

# initialize metabolite and reaction dict
metabolite_dict = {}
metabolite_dict['metabolites'] = []
reaction_dict = {}
reaction_dict['reactions'] = []
pathway_dict = {}
pathway_dict['pathways'] = []

# add compounds to metabolite_dict
compounds = ["ethanol","acetaldehyde"]
metabolite_dict = add_to_metabolite_dict(metabolite_dict,compounds,bigg2kegg_db)

# add reaction
[reaction_dict, metabolite_dict] = add_to_reaction_dict(reaction_dict,compounds,metabolite_dict,bigg2kegg_rxn_db,bigg2kegg_db)

# add pathway
precursor = 'acetaldehyde'
product = 'ethanol'
pathway_rxns = ['ALCD2x']
pathway_dict = add_to_pathway_dict(pathway_dict,product,pathway_rxns,precursor,metabolite_dict)

# create a new problem
problem_name = 'ex'
create_new_problem(problem_name,mmdb_path)

# add dictionaries to csv file in problem folder
add_to_pathway_csv_table(problem_name,pathway_dict, mmdb_path)
add_to_reaction_csv_table(problem_name,reaction_dict, mmdb_path)
add_to_metabolite_csv_table(problem_name,metabolite_dict, mmdb_path)
