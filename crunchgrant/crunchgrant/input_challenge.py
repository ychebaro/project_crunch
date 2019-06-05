import numpy as np 
import pandas as pd 
import matplotlib as mpl 
import matplotlib.pyplot as plt
from collections import OrderedDict
from bokeh.io import output_notebook
from bokeh.plotting import figure, show, output_file
from bokeh.models import LabelSet, Label, HoverTool, ColumnDataSource
import lxml.etree as ET
from glob import glob
from bs4 import BeautifulSoup
import requests
import requests_cache

def process_data(file_mutation, file_cna):
    # Get mutations data
    mut_df = pd.read_csv(file_mutation, sep='\t')

    # Get copy-number alterations
    cna_df = pd.read_csv(file_cna, sep='\t')

    # Make dictionnary of mutations and CNA
    gene_mutations = {col: mut_df[col].count() for col in mut_df.columns[2:]}
    gene_cna = {col: (cna_df[col].dropna() != 0).sum() for col in cna_df.columns[2:]}
    genes = list(mut_df.columns[2:])
    gene_modifs = OrderedDict()
    for k,v in gene_mutations.items():
        gene_modifs[k] = (gene_mutations[k] + gene_cna[k])/len(mut_df.index)

    gene_cna_per = [100*v/len(mut_df.index) for v in gene_cna.values()]
    gene_mutations_per = [100*v/len(mut_df.index) for v in gene_mutations.values()]
    gene_modifs_ratio = [v for v in gene_modifs.values()]

    mutations = {col: list(mut_df[col].dropna()) for col in mut_df.columns[2:]}

    return gene_cna, gene_cna_per, gene_mutations, gene_mutations_per, genes, gene_modifs_ratio, mutations

def get_gene_add_names(session, gene_names, gene_ncbi):

    genes_add_names = []
    for gene, ncbi_id in gene_ncbi.items():
        gene_add_names = [gene]
        page = session.get('https://www.ncbi.nlm.nih.gov/gene/'+str(ncbi_id))
        soup = BeautifulSoup(page.text, features="lxml")
       
        summary = soup.find("dl", attrs={'id':'summaryDl'})
        summary_dd = summary.find_all('dd')
        summary_dt = summary.find_all('dt')

        for dt, dd in zip(summary_dt, summary_dd):
            if dt.text == 'Also known as':
                names = dd.text
                names_list = names.split(';')
                for name in names_list:
                    gene_add_names.append(name.strip())
        genes_add_names.append(gene_add_names)
        
    return genes_add_names

# Parsing NSF data
def get_nsf_data(filenames, gene_add_names):

    parser = ET.XMLParser(recover=True)
    citing_grants = {}

    for file in filenames:
        with open(file) as f:
            tree = ET.parse(f, parser=parser)
            root = tree.getroot()

            check_if_conference = [elm.text for elm in root.iter('AwardTitle') if 'Conference' in elm.text]
            if len(check_if_conference) == 0 :
                awardid = [elm.text for elm in root.iter('AwardID')][0]
                firstnames_investigators = [elm.text for elm in root.iter('FirstName')]
                lastnames_investigators = [elm.text for elm in root.iter('LastName')]
                investigators = [str(firstnames_investigators[i])+' '+str(lastnames_investigators[i]) for i in range(len(firstnames_investigators))]
                emails = [elm.text for elm in root.iter('EmailAddress')]
                university = [elm.text for elm in root.iter('Name')]
                amount = [elm.text for elm in root.iter('AwardAmount')]
                for genes in gene_add_names:
                    for gene in genes:
                        for resume in root.iter('AbstractNarration'):
                            if resume is not None and resume.text is not None:
                                resume_splitted = (resume.text).split()
                                if gene in resume_splitted: 
                                    gene_index = resume_splitted.index(gene)
                                    if gene not in citing_grants:
                                        citing_grants[genes[0]] = [awardid, investigators, emails, university, amount]
                                    else:
                                        citing_grants[gene].update([awardid, investigators, emails, university, amount])
                                if 'p'+str(gene) in resume_splitted:
                                    gene_index = resume_splitted.index('p'+str(gene))
                                    if gene not in citing_grants:
                                        citing_grants[genes[0]] = [awardid, investigators, emails, university, amount]
                                    else:
                                        citing_grants[genes[0]].update([awardid, investigators, emails, university, amount])
                                if gene.lower() in resume_splitted:
                                    gene_index = resume_splitted.index(gene.lower())
                                    if gene not in citing_grants:
                                        citing_grants[genes[0]] = [awardid, investigators, emails, university, amount]
                                    else:
                                        citing_grants[genes[0]].update([awardid, investigators, emails, university, amount])

    return citing_grants


# Get data for p53 signaling
gene_cna, gene_cna_per, gene_mutations, gene_mutations_per, 
    genes, gene_modifs_ratio, mutations = process_data('mutations_p53.txt', 'cna_p53.txt')

# Map gene name to alternative names from the NCBI (mapping downloaded from HUGO)
mapping = pd.read_csv('gene_ncbi_uniprot.tsv', sep='\t')

# Scrap NCBI to get alternative names for the genes, i.e. encoded proteins
gene_ncbi = {gene: int(mapping.NCBI_Gene_ID[mapping[mapping.Approved_symbol == str(gene)].index.tolist()[0]]) for gene in genes}

cache_name = 'project_crunchbase'
MAX_RETRIES = 5
session = requests_cache.core.CachedSession(cache_name, backend='sqlite', expire_after=None)
session.mount('http://', requests.adapters.HTTPAdapter(max_retries=MAX_RETRIES))

filenames = glob("2018/*.xml")
genes_add_names = get_gene_add_names(session, genes, gene_ncbi)

citing_grants = get_nsf_data(filenames, genes_add_names)

liste_nsf_id = []
liste_nsf_pi = []
liste_nsf_emails =[]
liste_nsf_uni = []
liste_nsf_amount = []
liste_colors = []

for gene in genes:
    if gene not in citing_grants:
        liste_nsf_id.append('NA')
        liste_nsf_pi.append('NA')
        liste_nsf_emails.append('NA')
        liste_nsf_uni.append('NA')
        liste_nsf_amount.append('NA')
        liste_colors.append('blue')
    else:
        liste_nsf_id.append(','.join(citing_grants[gene][0]))
        liste_nsf_pi.append(','.join(citing_grants[gene][1]))
        liste_nsf_emails.append(','.join(citing_grants[gene][2]))
        liste_nsf_uni.append(','.join(citing_grants[gene][3]))
        liste_nsf_amount.append(','.join(citing_grants[gene][4]))
        liste_colors.append('red')

