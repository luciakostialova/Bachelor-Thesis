#!/bin/bash
#PBS -N eggNOG_rest2.0
#PBS -l select=8:ncpus=10:mem=100gb:scratch_local=100gb
#PBS -l walltime=100:00:00
#PBS -m ae

# this script uses a for cycle to parse the files in the input folder, which are then used as input for eggNOG-mapper
# used for a greater number of fasta files
# the directory paths need adjusting
# the script is run in Metacentrum

umask 0002

# path to the eggnog database
eggnog_db_dir='/storage/brno2/home/mariap3636/EggNOGv5'

# folder where the results will be located
output_dir='/storage/brno2/home/luciakostialova/Bachelor_thesis/EggNOG_Results'

# folder where fasta files are located
input_path='/storage/brno2/home/luciakostialova/Bachelor_thesis/input_COMPLETE'

module add conda-modules-py37

conda activate eggNOG-mapper_v2.1.9_py3.9  # deactivate after job finishes

for input in $(ls $input_path/*.fasta);do
        FILE='/storage/brno2/home/luciakostialova/Bachelor_thesis/EggNOG_Results/'`basename $input .fasta`'.emapper.hits'
        if test -f "$FILE"; then
                #echo "$FILE exists."
                continue
        fi
	name=`basename "$input" .fasta`    # `basename "$input" .fasta`
        emapper.py -i $input --data_dir $eggnog_db_dir -o $name --output_dir $output_dir --override --dmnd_ignore_warnings --evalue 0.001 --score 60 --pident 40 --query_cover 20 --subject_cover 20 --report_orthologs --tax_scope bacteria
done

conda deactivate

