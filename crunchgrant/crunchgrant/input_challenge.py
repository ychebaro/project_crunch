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
    mut_df = pd.read_csv("mutations_p53.txt", sep='\t')

    # Get copy-number alterations
    cna_df = pd.read_csv("cna_p53.txt", sep='\t')

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

# # Get data for cell cycle
# gene_cna, gene_cna_per, gene_mutations, gene_mutations_per, genes, gene_modifs_ratio, mutations = process_data('mutations_p53.txt', 'cna_p53.txt')

# colors = [
#      "#%02x%02x%02x" % (int(r), int(g), int(b)) for r, g, b, _ in 255*mpl.cm.viridis(mpl.colors.Normalize()(gene_modifs_ratio))
#   ]

# # Create source for plotting, useful for Labeling afterwards
# source = ColumnDataSource(data=dict(cna=gene_cna_per, 
#                                     mut=gene_mutations_per,
#                                     gene_names=genes,
#                                     gene_modif=gene_modifs_ratio,
#                                     col=colors))

# # Formatting title and labels
# p = figure(title='Gene alterations in the cell cycle pathway')
# p.xaxis[0].axis_label = '% copy-number alterations'
# p.yaxis[0].axis_label = '% mutations'

# p.scatter('cna', 'mut', source=source, alpha=0.6, radius='gene_modif', fill_color='col', line_color=None)
# labels = LabelSet(x='cna', y='mut', text='gene_names', source=source)
# p.add_layout(labels)

# # Map gene name to alternative names from the NCBI (mapping downloaded from HUGO)
mapping = pd.read_csv('gene_ncbi_uniprot.tsv', sep='\t')

# # Scrap NCBI to get alternative names for the genes, i.e. encoded proteins
# gene_ncbi = {gene: int(mapping.NCBI_Gene_ID[mapping[mapping.Approved_symbol == str(gene)].index.tolist()[0]]) for gene in genes}

cache_name = 'project_crunchbase'
MAX_RETRIES = 5
session = requests_cache.core.CachedSession(cache_name, backend='sqlite', expire_after=None)
session.mount('http://', requests.adapters.HTTPAdapter(max_retries=MAX_RETRIES))

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

# filenames = glob("2018/*.xml")
# genes_add_names = get_gene_add_names(session, genes, gene_ncbi)

# citing_grants = get_nsf_data(filenames, genes_add_names)


# liste_nsf_id = []
# liste_nsf_pi = []
# liste_nsf_emails =[]
# liste_nsf_uni = []
# liste_nsf_amount = []
# liste_colors = []

# for gene in genes:
#     if gene not in citing_grants:
#         liste_nsf_id.append('NA')
#         liste_nsf_pi.append('NA')
#         liste_nsf_emails.append('NA')
#         liste_nsf_uni.append('NA')
#         liste_nsf_amount.append('NA')
#         liste_colors.append('blue')
#     else:
#         liste_nsf_id.append(citing_grants[gene][0])
#         liste_nsf_pi.append(citing_grants[gene][1])
#         liste_nsf_emails.append(citing_grants[gene][2])
#         liste_nsf_uni.append(citing_grants[gene][3])
#         liste_nsf_amount.append(citing_grants[gene][4])
#         liste_colors.append('red')

# # Plotting NSF grant in gene representation
# p2 = figure(title='NSF grants on genes altered in cell cycle pathway altered in cancers')

# source2 = ColumnDataSource(data=dict(cna=gene_cna_per, 
#                             mut=gene_mutations_per,
#                             gene_names=genes,
#                             gene_modif=gene_modifs_ratio,
#                             liste_nsf_id=liste_nsf_id,
#                             liste_nsf_pi=liste_nsf_pi,
#                             liste_nsf_emails=liste_nsf_emails,
#                             liste_nsf_uni=liste_nsf_uni,
#                             liste_nsf_amount=liste_nsf_amount,
#                             col=liste_colors))

# tooltips = [('NSF Award ID', '@liste_nsf_id'), ('Investigators', '@liste_nsf_pi'), 
#             ('Emails', '@liste_nsf_emails'), ('University', '@liste_nsf_uni'), ('Amount', '@liste_nsf_amount')]

# p2.scatter('cna', 'mut', source=source2, alpha=0.6, fill_color='col', line_color=None)
# p2.xaxis[0].axis_label = '% copy-number alterations'
# p2.yaxis[0].axis_label = '% mutations'

# labels = LabelSet(x='cna', y='mut', text='gene_names', source=source2, text_font_size="8pt")

# p2.add_tools(HoverTool(tooltips=tooltips))
# p2.add_layout(labels)

# Get data for p53 signaling
gene_cna, gene_cna_per, gene_mutations, gene_mutations_per, genes, gene_modifs_ratio, mutations = process_data('mutations_p53.txt', 'cna_p53.txt')

# colors = [
#      "#%02x%02x%02x" % (int(r), int(g), int(b)) for r, g, b, _ in 255*mpl.cm.viridis(mpl.colors.Normalize()(gene_modifs_ratio))
#   ]

# # Create source for plotting, useful for Labeling afterwards
# source3 = ColumnDataSource(data=dict(cna=gene_cna_per, 
#                                     mut=gene_mutations_per,
#                                     gene_names=genes,
#                                     gene_modif=gene_modifs_ratio,
#                                     col=colors))

# # Formatting title and labels
# p3 = figure(title='Cancer associated gene alterations in the p53 signaling pathway')
# p3.xaxis[0].axis_label = '% copy-number alterations'
# p3.yaxis[0].axis_label = '% mutations'

# p3.scatter('cna', 'mut', source=source3, alpha=0.6, radius='gene_modif', fill_color='col', line_color=None)
# labels = LabelSet(x='cna', y='mut', text='gene_names', source=source3)
# p3.add_layout(labels)

# Scrap NCBI to get alternative names for the genes, i.e. encoded proteins
gene_ncbi = {gene: int(mapping.NCBI_Gene_ID[mapping[mapping.Approved_symbol == str(gene)].index.tolist()[0]]) for gene in genes}


filenames = glob("2018/*.xml")
gene_add_names = get_gene_add_names(session, genes, gene_ncbi)
citing_grants = get_nsf_data(filenames, gene_add_names)

print(citing_grants)

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

output_file("label2.html", title='test merde')

print(liste_nsf_amount)
print(liste_nsf_emails)
print(liste_nsf_id)
# Plotting NSF grant in gene representation
# p4 = figure(title='NSF grants on genes altered in the p53 signaling pathway altered in cancers')

# source4 = ColumnDataSource(data=dict(cna=gene_cna_per, 
#                             mut=gene_mutations_per,
#                             gene_names=genes,
#                             liste_nsf_id=liste_nsf_id,
#                             liste_nsf_pi=liste_nsf_pi,
#                             liste_nsf_emails=liste_nsf_emails,
#                             liste_nsf_uni=liste_nsf_uni,
#                             liste_nsf_amount=liste_nsf_amount,
#                             col=liste_colors))

# tooltips = [('NSF Award ID', '@liste_nsf_id'), ('Investigators', '@liste_nsf_pi'), 
#             ('Emails', '@liste_nsf_emails'), ('University', '@liste_nsf_uni'), ('Amount', '@liste_nsf_amount')]

# p4.scatter('cna', 'mut', source=source4, alpha=0.6, fill_color='col', line_color=None)
# p4.xaxis[0].axis_label = '% copy-number alterations'
# p4.yaxis[0].axis_label = '% mutations'

# labels = LabelSet(x='cna', y='mut', text='gene_names', source=source4, text_font_size="8pt")

# p4.add_tools(HoverTool(tooltips=tooltips))
# p4.add_layout(labels)

# show(p4)

# for gene in genes:
#     if gene not in citing_grants:
#         liste_nsf_id.append('NA')
#         liste_nsf_pi.append('NA')
#         liste_nsf_emails.append('NA')
#         liste_nsf_uni.append('NA')
#         liste_nsf_amount.append('NA')
#         liste_colors.append('blue')
#     else:
#         liste_nsf_id.append(citing_grants[gene][0])
#         liste_nsf_pi.append(citing_grants[gene][1])
#         liste_nsf_emails.append(citing_grants[gene][2])
#         liste_nsf_uni.append(citing_grants[gene][3])
#         liste_nsf_amount.append(citing_grants[gene][4])
#         liste_colors.append('red')

# # Plotting NSF grant in gene representation
# p2 = figure(title='NSF grants on genes altered in cell cycle pathway altered in cancers')

# source2 = ColumnDataSource(data=dict(cna=gene_cna_per, 
#                             mut=gene_mutations_per,
#                             gene_names=genes,
#                             gene_modif=gene_modifs_ratio,
#                             liste_nsf_id=liste_nsf_id,
#                             liste_nsf_pi=liste_nsf_pi,
#                             liste_nsf_emails=liste_nsf_emails,
#                             liste_nsf_uni=liste_nsf_uni,
#                             liste_nsf_amount=liste_nsf_amount,
#                             col=liste_colors))

# tooltips = [('NSF Award ID', '@liste_nsf_id'), ('Investigators', '@liste_nsf_pi'), 
#             ('Emails', '@liste_nsf_emails'), ('University', '@liste_nsf_uni'), ('Amount', '@liste_nsf_amount')]

# p2.scatter('cna', 'mut', source=source2, alpha=0.6, fill_color='col', line_color=None)
# p2.xaxis[0].axis_label = '% copy-number alterations'
# p2.yaxis[0].axis_label = '% mutations'

# labels = LabelSet(x='cna', y='mut', text='gene_names', source=source2, text_font_size="8pt")

# p2.add_tools(HoverTool(tooltips=tooltips))
# p2.add_layout(labels)
