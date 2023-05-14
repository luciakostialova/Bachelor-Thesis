#!/bin/bash
#PBS -N NCBIJob
#PBS -l select=1:ncpus=2:mem=8gb:scratch_local=8gb
#PBS -l walltime=10:00:00
#PBS -m ae

# downloads proteome seq in fasta format from ncbi
# taxonID.txt file contains the taxon ids attributed to the proteomes

module load edirect

listOfTAXONIDSfile='/storage/brno2/home/luciakostialova/Bachelor_thesis/taxonID.txt'

for taxonID in `cat $listOfTAXONIDSfile`;do
	out_file="/storage/brno2/home/luciakostialova/Bachelor_thesis/NCBI_proteomes/${taxonID/$'\r'/}.fasta"
	key="txid${taxonID/$'\r'/}[Organism:exp]"
	`esearch -db "ipg" -query $key | efetch -format fasta > $out_file`
done

