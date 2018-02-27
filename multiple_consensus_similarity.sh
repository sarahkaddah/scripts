#!/bin/bash

# Repertory: data/seq_AS_fasta/
# Arguments:
	# 1- fichier fasta
	# 2- fichier classification
	# 3- chemin vers resultat 
	# 4- prefixe fichiers

# Partie similatity	
#for i in {1..10}
#do 
#	echo "================Run $i=================="
#	../../script/get_similarity_by_family.py -s $1 -f $2 -o $3/$4_$i.similarity -n 500 -m 100 -d 0 --intra_only -v 9 &> $3/out/outsimilarity_intra_$4_$i.log &
#	echo "=================DONE==================="
#done

# Partie similatity	
for i in {1..10}
do 
	echo "================Run $i=================="
	../../script/build_consensus_by_family.py -s $1 -f $2 -o $3/$4_$i.consensus -m 500 -v 9 &> $3/out/outconsensus_$4_$i.log &
	echo "=================DONE==================="
done

 

