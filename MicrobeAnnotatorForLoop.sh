#!/bin/bash
#PBS -N microbeAnnotator_Job
#PBS -l select=1:ncpus=2:mem=12gb:scratch_local=12gb
#PBS -l walltime=24:00:00
#PBS -m ae

# this script uses a for cycle to parse the files in the input folder, which are used as input for MicrobeAnnotator
# used for a greater number of fasta files
# the directory paths need adjusting

module add conda-modules-py37

microbeannotator_db_dir='/storage/brno2/home/mariap3636/MICROBEANNOTATOR/microbeAnnotator_dblight'

# folder where the results are located
output_dir='/storage/brno2/home/luciakostialova/Bachelor_thesis/MicrobeAnnotator_Results'

# folder where fasta files are located
input_path='/storage/brno2/home/luciakostialova/Bachelor_thesis/input_COMPLETE'

conda activate microbeannotator

for a in $(ls $input_path/*.fasta);do

	microbeannotator -i $a -d $microbeannotator_db_dir -o $output_dir -m diamond -p 1 -t 2 --light

done

conda deactivate
