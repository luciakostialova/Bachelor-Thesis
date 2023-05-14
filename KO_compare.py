### this script will compare the 2 annotation files of 1 proteome!


import pandas as pd
import argparse
import sys
import os
import numpy as np
import pathlib


# merges the annotation files and prepares the data
def file_edit(microbeannotator_file, eggnog_file):
    microbeannotator = pd.read_table(microbeannotator_file, sep='\t', header=0)
    eggnog = pd.read_table(eggnog_file, sep='\t', header=0, skiprows=4, skipfooter=3, engine='python')
    temp = pd.merge(microbeannotator, eggnog, how='outer', left_on='query_id', right_on='#query')
    #KO_only = temp[['query_id', 'ko_number', 'KEGG_ko']]
    KO_andMore = temp[['query_id', 'ko_number', 'KEGG_ko', 'ko_product', 'Description', 'EC']]
    return KO_andMore

def isnan(x):   # function tests for NaN
    if x != x:
        return True
    else:
        return False
    
# count the match and mismatches, vector, and mismatch pairs

def match_mismatch_count(input_file):   
    match = 0; miss = 0; nf = 0; egg_NF = 0; ma_NF = 0; x = 0
    vector = []
    mismatches = []
    #mismatch_df = pd.DataFrame(columns=['query_id', 'ko_number', 'KEGG_ko', 'Description', 'ko_product'])
    for i in range(len(input_file)):  
        if isnan(input_file.loc[i, 'ko_number']):
            if input_file.loc[i, 'KEGG_ko'] == '-':
                nf += 1; vector.append(0)
                continue
            if isnan(input_file.loc[i, 'KEGG_ko']):
                nf += 1; vector.append(0)
                continue
            else:
                ma_NF += 1; vector.append(0)  # do i need this?
            continue
        # neither tools found a KO, microbeannot did not find KO, egg did
        if isnan(input_file.loc[i, 'KEGG_ko']) or input_file.loc[i, 'KEGG_ko'] == '-':
            egg_NF += 1; vector.append(0)
            continue
        # eggnog didn't find KO, ma did
        if ',' in input_file.loc[i, 'KEGG_ko']:
            x += 1; vector.append(1)
            mismatches.append((ma, egg))
            continue
        # eggnog found multiple KO's             
        whole_egg = str(input_file.loc[i, 'KEGG_ko'])
        egg = whole_egg[3:]
        ma = input_file.loc[i, 'ko_number']
        if ma != egg:
            miss += 1; vector.append(1)
            mismatches.append((ma, egg))
            continue
        # ma and eggnog found different KO's
        if ma == egg:
            match += 1; vector.append(0)
            continue
        # ma and eggnog found the same KO
    input_file['match'] = vector
    #match_mismatch_results = input_file
    results_filtered = input_file[input_file['match'] == 1]
    results_filtered = results_filtered.drop(['match'], axis=1)
    results_filtered.rename(columns = {'ko_number':'ko_microbeA', 'KEGG_ko':'ko_eggnog', 'ko_product':'ko_description_microbeA', 'Description':'description_eggnog', 'EC':'EC_eggnog'}, inplace = True)
    return [mismatches, results_filtered, vector, match, miss]

def KO_counter(input_file):  # count how many KO's did eggnog and mannot find
    egg_KO = 0; ma_KO = 0
    for i in range(len(input_file)): 
        if isnan(input_file.loc[i, 'ko_number']) == False:
            ma_KO += 1
        if isnan(input_file.loc[i, 'KEGG_ko']):
            continue
        if input_file.loc[i, 'KEGG_ko'] != '-':
            egg_KO += 1
    return [ma_KO, egg_KO]

def KO_counter_unique(input_file):  # count how many unique KO identifiers were found
    ma_KO_unique = input_file['ko_number'].nunique() - 1  # because I don't take into account the nan
    KEGG_KO = input_file.drop(['query_id', 'ko_number', 'Description', 'ko_product'], axis=1)
    KO_unique = KEGG_KO.drop_duplicates() # what about the eggnog KO tuples?
    KO_unique = KO_unique.dropna()
    egg_KOs = []
    for ind in KO_unique.index:
        if KO_unique['KEGG_ko'][ind] == '-':
            continue
        if ',' in KO_unique['KEGG_ko'][ind]:
            for j in KO_unique['KEGG_ko'][ind].split():
                egg_KOs.append(j)
            continue
        else:
            egg_KOs.append(KO_unique['KEGG_ko'][ind])
    egg_KO_unique = len(egg_KOs)
    print(egg_KO_unique, ma_KO_unique)
    return [ma_KO_unique, egg_KO_unique]

def KO_freq_counter(input_file):
    ko_number = input_file['ko_microbeA']
    KEGG_ko = input_file['ko_eggnog']
    df_ma = ko_number.value_counts() #.rename_axis('unique_values').to_frame('counts')
    df_ma_n = ko_number.value_counts(normalize=True)
    df_egg = KEGG_ko.value_counts(); df_egg_n = KEGG_ko.value_counts(normalize=True)
    df1 = pd.DataFrame(df_ma.index, columns=['MA KO'])
    df1['MA KO absolute freq'] = df_ma.values; df1['MA KO relative freq'] = df_ma_n.values
    df2 = pd.DataFrame(df_egg.index, columns=['EggNOG KO'])
    df2['EggNOG KO absolute freq'] = df_egg.values; df2['EggNOG KO relative freq'] = df_egg_n.values
    df = pd.concat([df1, df2], axis=1)
    print(df1.sum(), df2.sum())
    return df

def KO_mismatch_txt(input, filename):     # make a csv/txt file that contains the mismatch KO identifiers, without the descriptions
    with open('{}_KO_mismatches.csv'.format(filename), 'w') as file:              
        for pair in input:
            file.write("{}, {}\n".format(pair[0], pair[1]))
    print('The file {} contains the KO pairs that are a mismatch.'.format(file.name))

"""--- Main Function ---"""

def main():
    # Setup parser for arguments.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, 
        description='''This compares the KO identifiers from MicrobeAnnotator and EggNOG annotation files.\n'''
        '''It returns a text file of the mismatch KO pairs; a vector based of 0 and 1, based on whether 
        there was a match or not; and the various match and mismatch counts.\n'''
        '''Usage: ''' + sys.argv[0] + ''' -m [MicrobeAnnotator File] -e [EggNOG File]\n'''
        '''Global mandatory parameters: -m [MicrobeAnnotator File] -e [EggNOG File]\n''')
    parser.add_argument('-m', '--microbeannotator_file', dest='microbeannotator_file', 
                        action='store', required=True, help='MicrobeAnnotator annotation file')
    parser.add_argument('-e', '--eggnog_file', dest='eggnog_file', action='store', required=True,
                        help='EggNOG annotation file')
    args = parser.parse_args()
    microbeannotator_file = args.microbeannotator_file   # file extension is '.fasta.annot'
    eggnog_file = args.eggnog_file                  # file extension is '.emapper.annotations'
    basename = os.path.basename(microbeannotator_file).split('.')[0]  # pathlib.Path(eggnog_file).stem doesn't work
    # or eggnog_file, it doesn't matter
    input_file = file_edit(microbeannotator_file, eggnog_file)
    results = match_mismatch_count(input_file)
    KO_mismatches = results[0]
    KO_mismatch_txt(KO_mismatches, basename)
    return results

    
    
if __name__ == "__main__":
    main()









