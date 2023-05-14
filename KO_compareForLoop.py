### this script will use the script KO_compare.py in a for loop for all the Tanoue bacteria
### the annotation files given by MicrobeAnnotator and EggNOG are used
### it returns a text file of the KO mismatches for all of the Tanoue bacteria,
### a matrix containing the number of KO identifiers found by MicrobeAnnotator and EggNOG
### a matrix containing the match-mismatch vectors for each of the Tanoue bacteria

### the annotation files need to be in double quotes - windows has a problem with whitespace in paths
### for some reason the double quotes don't help the python script - need to save it in temp folder
### row names of vector table?

import pandas as pd
from pathlib import Path
import pathlib
import os
import KO_compare   #KO_compare2.main() microbeannotator_file eggnog_file

input = Path('C:/Users/Lucka/Documents/MUNI/5. semester/Bakalárska práca/Tanoue/Tanoue_protein_sequences')
vector_data = {}
eggnog_count = []
ma_count = []
names = []

for file2 in sorted(input.glob('*.fasta')):
    basename = pathlib.Path(file2).stem
    eggnog_file = os.path.join('C:/Users/Lucka/Documents/MUNI/5. semester/Bakalárska práca/Tanoue/Functional annotation results/EggNOG_results_Tanoue/{}.emapper.annotations'.format(basename))
    microbeannotator_file = os.path.join('C:/Users/Lucka/Documents/MUNI/5. semester/Bakalárska práca/Tanoue/Functional annotation results/MicrobeAnnotator_results_Tanoue/{}.fasta.annot'.format(basename))
    
    input_file = KO_compare.file_edit(microbeannotator_file, eggnog_file)
    results = KO_compare.match_mismatch_count(input_file)
    KO_mismatches = results[0]
    KO_compare.KO_mismatch_txt(KO_mismatches, basename)
    KO_mismatches_descrip = results[1]
    KO_mismatches_descrip.to_csv('{}_KO_mismatches_description_new.csv'.format(basename), encoding='utf-8', mode='w')
    print('The file {}_KO_mismatches_description.csv contains the KO mismatch pairs along with their description'.format(basename))
    KO_freq = KO_compare.KO_freq_counter(KO_mismatches_descrip)
    KO_freq.to_csv('{}_KO_mismatch_frequency.csv'.format(basename), encoding='utf-8', mode='w')
    ma_KO_unique, egg_KO_unique = KO_compare.KO_counter_unique(input_file)
    ma_count.append(ma_KO_unique); eggnog_count.append(egg_KO_unique); names.append(basename)
        
df3 = pd.DataFrame(list(zip(ma_count, eggnog_count)), columns =['MicrobeAnnotator KO counts', 'EggNOG KO counts'], index=names)
df3.to_csv('matrixOfUniqueKOcounts.csv', encoding='utf8', mode='w')
print('The matrixOfUniqueKOcounts.csv contains the number of KO identifiers found by MicrobeAnnotator and EggNOG.')


