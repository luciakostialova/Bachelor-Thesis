#!/bin/bash
## this script takes the eggnog annotation files, extracts the KO column and edits the data.
## in the next step the python script KO_mapper.py is ran on the new .ko files
## the directory paths may need changing

umask 0002

# folder where the results will be located
output_path='/storage/brno2/home/luciakostialova/Bachelor_thesis/Tanoue_complete/KO_mapper/EggNOG'

# folder where annotation files are located
input_path='/storage/brno2/home/luciakostialova/Bachelor_thesis/Tanoue_complete/eggNOG'

# extracts of KO column, the data is cleaned
for file in $(ls $input_path/*.emapper.annotations);do
        name=`basename "$file" .emapper.annotations`
        cut -f 12 $file | grep -v '^#\|-\|KEGG_ko' | sed 's/\,/\n/g' | cut -c 4- > $output_path'/'$name'.emapper.ko'
done


# for cycle,  use KO_mapper on every new .ko file

# the data directory has to be in the directory where KO_mapper.py is happening
# as does the script KO_mapper.py
# there is no need to export the bash variables

cd $output_path

module add python36-modules-gcc

for file2 in $(ls $output_path/*.ko); do
	name=`basename "$file2" .emapper.ko`'_eggnog'
	python KO_mapper.py -i $file2 -p $name
done

module unload python36-modules-gcc

exit

