
#### this script takes as input the annotation results of both tools
#### the annotation performance criteria (annotation coverage, etc.) are calculated
#### the values are saved in csv file
#### all the Tanoue proteomes are compared


import pandas as pd
from pathlib import Path
import pathlib
import os
import re

input = Path('C:/Users/Lucka/Documents/MUNI/5. semester/Bakalárska práca/Tanoue/Tanoue_protein_sequences')

def count_annotation_overlap(microbeannotator, eggnog):
    temp = pd.merge(microbeannotator, eggnog, how='outer', left_on='query_id', right_on='#query')
    df = temp[['query_id', '#query', 'database', 'seed_ortholog']]
    filtered_df = df.loc[(df['database'] == 'swissprot') | (df['database'] == 'kofam')]
    filtered_df2 = filtered_df.loc[filtered_df['#query'].notna() == True]
    n = len(filtered_df2.index)
    return n

def only_eggnog_annot(microbeannotator, eggnog):
    temp = pd.merge(microbeannotator, eggnog, how='outer', left_on='query_id', right_on='#query')
    df = temp[['query_id', '#query', 'database', 'seed_ortholog']]
    filtered_df = df.loc[(df['database'].notna() == False)]
    filtered_df = filtered_df[filtered_df['#query'].notna() == True]
    n = len(filtered_df.index)
    return n

def only_ma_annot(microbeannotator, eggnog):
    temp = pd.merge(microbeannotator, eggnog, how='outer', left_on='query_id', right_on='#query')
    df = temp[['query_id', '#query', 'database', 'seed_ortholog']]
    filtered_df = df.loc[(df['#query'].notna() == False)]
    filtered_df = filtered_df.loc[(filtered_df['database'] == 'kofam') | (filtered_df['database'] == 'swissprot')]
    n = len(filtered_df.index)
    return n

def count_proteins(file):
    sum = 0
    with open(file, 'r') as f:
        for line in f:
            if re.search(">", line):
                sum += 1
    return sum

def count_annotated_prot_ma(file):    # a prot is annotated if it has a KO or a swissprot annotation
    n = 0
    for i in range(len(file)): 
        if file.loc[i, 'database'] == 'kofam' or file.loc[i, 'database'] == 'swissprot':
            n += 1
    return n

def count_annotated_prot_egg(file):    # a prot is annotated if it has a seed ortholog - which it always does - count length of df col?
    n = file['seed_ortholog'].count()
    return n

def classification_ma(microbeannotator, sum):
    words = ['hypothetical', 'uncharacterized', 'domain of unknown function', 'protein of unknown function']
    hypothetical = 0
    for word in words:
         hypothetical += microbeannotator['ko_product'].str.contains(word).sum()
         hypothetical += microbeannotator['product'].str.contains(word).sum()
    unannotated = microbeannotator['product'].str.contains('No match found').sum()
    annotated = sum - hypothetical - unannotated
    return annotated, hypothetical, unannotated

def classification_eggnog(microbeannotator, eggnog, sum):
    words = ['hypothetical', 'uncharacterized', 'domain of unknown function', 'protein of unknown function',
             'Hypothetical', 'Uncharacterized', 'Domain of unknown function', 'Protein of unknown function']
    hypothetical = 0
    for word in words:
        hypothetical += eggnog['Description'].str.contains(word).sum()
    hypothetical += eggnog['Description'].value_counts()['-']
    ma = len(microbeannotator['query_id'].index)
    egg = len(eggnog['#query'].index)
    diff = ma - egg
    temp = pd.merge(microbeannotator, eggnog, how='outer', left_on='query_id', right_on='#query')
    filtered_df = temp.loc[(temp['#query'].notna() == False)]
    unannotated = len(filtered_df.index)
    annotated = sum - hypothetical - unannotated
    return annotated, hypothetical, unannotated

names = []
number_of_proteins = []  # number of proteins per proteome
matched_prot_ma = []
matched_prot_egg = []
overlap = []
eggnog_only_annot = []
ma_only_annot = []
annotated_prot_ma = []
hypothetical_prot_ma = []
unannotated_prot_ma = []
annotated_prot_egg = []
hypothetical_prot_egg = []
unannotated_prot_egg = []

for file in sorted(input.glob('*.fasta')):
    basename = pathlib.Path(file).stem
    names.append(basename)
    eggnog_file = os.path.join('C:/Users/Lucka/Documents/MUNI/5. semester/Bakalárska práca/Tanoue/EggNOG_results_Tanoue/{}.emapper.annotations'.format(basename))
    microbeannotator_file = os.path.join('C:/Users/Lucka/Documents/MUNI/5. semester/Bakalárska práca/Tanoue/MicrobeAnnotator_results_Tanoue/{}.fasta.annot'.format(basename))
    microbeannotator = pd.read_table(microbeannotator_file, sep='\t', header=0)
    eggnog = pd.read_table(eggnog_file, sep='\t', header=0, skiprows=4, skipfooter=3, engine='python')
    
    sum_of_proteins = count_proteins(file)
    number_of_proteins.append(count_proteins(file))
    matched_prot_ma.append(count_annotated_prot_ma(microbeannotator))
    matched_prot_egg.append(count_annotated_prot_egg(eggnog))
    overlap.append(count_annotation_overlap(microbeannotator, eggnog))
    eggnog_only_annot.append(only_eggnog_annot(microbeannotator, eggnog))
    ma_only_annot.append(only_ma_annot(microbeannotator, eggnog))
    annotated_ma, hypothetical_ma, unannotated_ma = classification_ma(microbeannotator, sum_of_proteins)
    annotated_egg, hypothetical_egg, unannotated_egg = classification_eggnog(microbeannotator, eggnog, sum_of_proteins)
    annotated_prot_ma.append(annotated_ma)
    hypothetical_prot_ma.append(hypothetical_ma)
    unannotated_prot_ma.append(unannotated_ma)
    annotated_prot_egg.append(annotated_egg)
    hypothetical_prot_egg.append(hypothetical_egg)
    unannotated_prot_egg.append(unannotated_egg)
    

df1 = pd.DataFrame(list(zip(number_of_proteins, matched_prot_ma, matched_prot_egg, overlap, eggnog_only_annot, ma_only_annot,
                            annotated_prot_ma, hypothetical_prot_ma, unannotated_prot_ma, annotated_prot_egg, hypothetical_prot_egg, unannotated_prot_egg)), 
                  columns =['Number of proteins', 'MA proteins with a match', 'eggNOG proteins with a match', 'Overlap of annotated proteins',
                            'eggNOG only annotated proteins', 'MA only annotated proteins', 'MA annotated proteins', 'MA hypothetical proteins', 
                            'MA unannotated proteins', 'eggNOG annotated proteins', 'eggNOG hypothetical proteins', 'eggNOG unannotated proteins'],
                            index=names)

df1.to_csv('general_statistics_Tanoue.csv')  