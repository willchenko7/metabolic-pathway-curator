# mmdb master script
import os
mmdb_path = r"C:\Users\Owner\Documents\PhD_stuff\research\mmdb"
os.chdir(mmdb_path)
from src.bigg2kegg_met_ID_conversion_DB import create_bigg2kegg_DB
from src.bigg2kegg_rxn_ID_conversion import create_bigg2kegg_rxn_db
from src.create_metabolite_dict import add_to_metabolite_dict
from src.get_reaction import add_to_reaction_dict
from src.create_pathways import manually_add_to_pathway_dict
from src.create_pathways import add_to_pathway_dict
from src.export_to_CSV import create_new_problem
from src.export_to_CSV import add_to_pathway_csv_table
from src.export_to_CSV import add_to_reaction_csv_table
from src.export_to_CSV import add_to_metabolite_csv_table
from src.solveForMissingCharges import solveForMissingChargesinRecursiveLoop
from src.balance_check import mass_check_reaction_dict
from src.balance_check import charge_check_reaction_dict
from src.add_to_local_db import add_all_data_to_local_db
from src.retrieve_data_from_local_db import add_all_info_of_paths_that_make_target
from src.addMetadataForPathways import add_all_paths_to_metadata
from src.balance_check import find_unbalanced_rxns
from src.retrieve_data_from_local_db import get_all_reactions_from_db
from src.retrieve_data_from_local_db import get_all_metabolites_from_db
from src.retrieve_data_from_local_db import get_all_info_from_db_path
from src.export2json import export2json
from src.balance_check import rxnChargeBalance
from src.balance_check import check_rxn_mass_sum

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
metadata = {}
metadata['pathways'] = []

# add compounds to metabolite_dict
compounds = ["1-propanol"]
[metabolite_dict,bigg2kegg_db] = add_to_metabolite_dict(metabolite_dict,compounds,bigg2kegg_db,mmdb_path)

compounds = ["L-Histidine","L-Histidinol"]
[reaction_dict, metabolite_dict] = add_to_reaction_dict(reaction_dict,compounds,metabolite_dict,bigg2kegg_rxn_db,bigg2kegg_db,mmdb_path)

# solve for unknown charges
metabolite_dict = solveForMissingChargesinRecursiveLoop(metabolite_dict,reaction_dict)

# check balance
reaction_dict = mass_check_reaction_dict(reaction_dict,metabolite_dict)
reaction_dict = charge_check_reaction_dict(reaction_dict,metabolite_dict)
find_unbalanced_rxns(reaction_dict)
rxnChargeBalance('HP',reaction_dict,metabolite_dict)
check_rxn_mass_sum('R5PDIP',reaction_dict,metabolite_dict)

# find and add all pathways between two metabolites
precursor = 'pyruvate'
target = 'Isobutyl acetate'
[pathway_dict,reaction_dict,metabolite_dict] = add_to_pathway_dict(pathway_dict,precursor,target,reaction_dict,metabolite_dict,mmdb_path)

# add data to local db
add_all_data_to_local_db(mmdb_path,metabolite_dict,reaction_dict,pathway_dict)

### retrieve all data from local db
# retrieve all data from a target's pathway
target = 'nopaline'
[metabolite_dict,reaction_dict,pathway_dict] = add_all_info_of_paths_that_make_target(target,mmdb_path,metabolite_dict,reaction_dict,pathway_dict)

#retrieve one path
path = 'melatonin_D-Erythrose 4-phosphate_0'
[metabolite_dict,reaction_dict,pathway_dict] = get_all_info_from_db_path(path,pathway_dict,reaction_dict,metabolite_dict,mmdb_path)


# retrieve all reactions 
reaction_dict = get_all_reactions_from_db(mmdb_path,reaction_dict)
metabolite_dict = get_all_metabolites_from_db(mmdb_path,metabolite_dict)

# add metadata
metadata = add_all_paths_to_metadata(metadata,pathway_dict,reaction_dict,metabolite_dict)


#manual metabolite addition
metabolite_dict = manually_add_to_metabolite_dict(metabolite_dict,'Isobutyl acetate')
reaction_dict = manually_add_to_reaction_dict('AATibutylace',reaction_dict)

# manually add pathways to pathway_dict
precursor = '2-oxoglutarate'
product = 'L-hydroxyproline'
pathway_rxns = ['GLUDxi','ACOTA','ORNTAC_1','ORNCD','PROAKGOX1']
pathway_name = 'L-hydroxyproline 1'
pathway_dict = add_to_pathway_dict(pathway_dict,pathway_name,product,pathway_rxns,precursor,metabolite_dict)

# add rxns to cobra model
from src.addRxnsToModel import addAllRxnsFromRxnDict
from src.addRxnsToModel import addExchangeRxns
model = addAllRxnsFromRxnDict(reaction_dict,metabolite_dict)
model = addExchangeRxns(model,'pyruvate','ethanol',metabolite_dict)


## exporting

#excel
# create a new problem
problem_name = 'basic-modules'
create_new_problem(problem_name,mmdb_path)
# add 
add_to_pathway_csv_table(problem_name,pathway_dict, mmdb_path)
add_to_reaction_csv_table(problem_name,reaction_dict, mmdb_path)
add_to_metabolite_csv_table(problem_name,metabolite_dict, mmdb_path)

# json
export2json('classic_fermentation.json',reaction_dict,metabolite_dict)
