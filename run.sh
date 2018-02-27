#!/bin/bash

#Parametres:
	# 1-Nom de l'input (kmers52)
	# 2-Nom de l'output (.dat)
	# 3-Valeur de la LDA
	# 4-Family size: taille maxi des familles: 100 (calculer)

#../../script/classification.py $1 --classification_output_file ../../results/$2/$2_5kmers_pm095_lda$3_fs$4.dat --verbose 9 --pairmate_threshold 0.95 --timing 1 --lda_sample_size $3 --min_family_size $4

#./classification.py $1 --classification_output_file ../../results/$2/$2_5kmers_pm095.dat --verbose 9 --pairmate_threshold 0.95 --timing 1

#Version 2
#../../script/classification.py $1 --classification_output_file ../../results/$2/$2_5kmers_pm095_lda$3_fs$4.dat --verbose 9 --pairmate_threshold 0.95 --timing 1 --lda_sample_size $3 --min_family_size $4

# Version 3
../../script/classification.py $1 --classification_output_file ../../results/$2/$2_5kmers_pm095_lda$3_fs$4.dat --verbose 9 --pairmate_threshold 0.95 --timing 1 --lda_sample_size $3 --min_family_size $4
