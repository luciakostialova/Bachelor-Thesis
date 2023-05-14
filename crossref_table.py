
### the script crossref_table.py quantifies the results from the Uniprot cross-reference analysis
### UniProtKB accessions, KEGG genes and eggNOG OGs (unique) are counted and saved as a csv

import pandas as pd
from pathlib import Path
import pathlib
import numpy as np
import os
import re

input = Path('Tanoue/UniProt crossreferencing/accessionToKeggOrEggnogFromUniprotWebsite/')

names = []
uniprotKB = []
kegg = []
eggnog = []
kegg_unique = []
eggnog_unique = []


for file in sorted(input.glob('*_eggnog.tsv')):
    basename = pathlib.Path(file).stem
    names.append(basename.split(sep='_')[0])
    df = pd.read_table(file, sep='\t', header=0)
    uniprotKB.append(len(df['Entry']))
    # explode kegg and egg columns, remove nan
    df = df.fillna('0')
    df['KEGG'] = df['KEGG'].str.split(';')
    df["eggNOG"] = df["eggNOG"].str.split(";")
    df = df.explode('eggNOG')
    df = df.explode('KEGG')
    df.replace('0', np.nan, inplace=True)
    kegg.append(df['KEGG'].count()); eggnog.append(df['eggNOG'].count())
    kegg_unique.append(df['KEGG'].nunique()); eggnog_unique.append(df['eggNOG'].nunique())


df1 = pd.DataFrame(list(zip(names, uniprotKB, kegg, kegg_unique, eggnog, eggnog_unique)), 
                  columns =['Proteome', 'UniProtKB Accessions', 'KEGG Genes', 'Unique KEGG Genes',
                            'eggNOG OGs', 'Unique eggNOG OGs'],
                            index=names)

df1.to_csv('general_statistics_uniprotCrossref.csv')  