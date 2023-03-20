# export dict's as a csv
import csv
import os

def create_new_problem(problem_name,mmdb_path):
    problem_path = os.path.join(mmdb_path,r"problems",problem_name)
    os.mkdir(problem_path)
    os.chdir(problem_path)
    create_pathway_csv_table()
    create_reaction_csv_table()
    create_metabolite_csv_table()
    return 

def create_pathway_csv_table():
    with open('pathway_table.csv','w',newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['id','name','rxns','product_id','product_tags','precursor_ids','litref','notes'])
    return

def create_reaction_csv_table():
    with open('reaction_table.csv','w',newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['id','name','rxn_str','kegg_id','bigg_id','ec_number','notes'])
    return

def create_metabolite_csv_table():
    with open('metabolite_table.csv','w',newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['id','name','formula','charge','kegg_id','bigg_id','notes'])
    return 

def add_to_pathway_csv_table(problem_name,pathway_dict, mmdb_path):
    problem_path = os.path.join(mmdb_path,r"problems",problem_name)
    os.chdir(problem_path)
    for path in pathway_dict['pathways']:
        with open('pathway_table.csv','a',newline='') as file:
            writer = csv.writer(file)
            writer.writerow([pathway_dict[path]['id'],pathway_dict[path]['name'],pathway_dict[path]['rxns'],pathway_dict[path]['product_id'][:-2],'',pathway_dict[path]['precursor_ids'][:-2],'',''])
    return

def add_to_reaction_csv_table(problem_name,reaction_dict, mmdb_path):
    problem_path = os.path.join(mmdb_path,r"problems",problem_name)
    os.chdir(problem_path)
    for rxn in reaction_dict['reactions']:
        with open('reaction_table.csv','a',newline='') as file:
            writer = csv.writer(file)
            writer.writerow([reaction_dict[rxn]['bigg_rxn'],reaction_dict[rxn]['bigg_name'],reaction_dict[rxn]['bigg_string'],reaction_dict[rxn]['kegg_rxn'],reaction_dict[rxn]['bigg_rxn'],'',''])
    return

def add_to_metabolite_csv_table(problem_name,metabolite_dict, mmdb_path):
    problem_path = os.path.join(mmdb_path,r"problems",problem_name)
    os.chdir(problem_path)
    with open('metabolite_table.csv','a',newline='') as file:
        writer = csv.writer(file)
        for met in metabolite_dict['metabolites']:
            writer.writerow([metabolite_dict[met]['bigg_id'][:-2],met,metabolite_dict[met]['formula'],metabolite_dict[met]['charge'],metabolite_dict[met]['kegg_id'],metabolite_dict[met]['bigg_id'],''])
    return
