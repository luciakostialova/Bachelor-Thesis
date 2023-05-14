#### used for retrieval of smaller proteomes from uniprot, using the uniparc database
#### all we need are proteome id numbers (UPO identifiers), which are saved in a csv file 
#### the main part of the script comes from the uniprot website!!

import requests
from pandas import *
from pathlib import Path
import os
import re

# reading CSV file
data = read_csv("proteomes.csv")
# converting column data to list
proteomes = data['x'].tolist()
# editing the data
proteomes = [i.strip() for i in proteomes if i==i]  # remove nan, whitespace
proteomes = ' '.join(proteomes).split()   # delete empty strings from list
print(len(proteomes))  # count number of ids
# counters
n = 0
p = 0
q = 0

for id in proteomes:
    url = 'https://rest.uniprot.org/uniparc/stream?format=fasta&query=({})'.format(id)  # retrieval from uniparc
    with requests.get(url, stream=True) as request:
        request.raise_for_status()
        # important to maintain the correct path to prevent errors
        path = 'C:/Users/Lucka/Documents/MUNI/5. semester/Bakalárska práca/uniparc_proteome/{}.fasta'.format(id)
        if os.path.exists(path) == False:  # check if file exists
            n += 1
            with open('{}.fasta'.format(id), 'wb') as f:
                for chunk in request.iter_content(chunk_size=2**20):
                    f.write(chunk)
        else:
            with open('{}.fasta'.format(id), 'r') as f2:
                first_line = f2.readline()
                if first_line[0] == '>':  # file already exists and contains proteome
                    p += 1
                    print('File {}.fasta already exists and is up to date.'.format(id))
                else:   # file exists but is empty - error
                    q += 1   
                    print('File {}.fasta exists but is empty.'.format(id))
print('Number of new files is:', n)
print('Number of already existing files is',p)
print('Number of empty files is', q)
