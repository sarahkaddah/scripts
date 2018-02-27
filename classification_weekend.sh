#!/bin/bash

# Repertory: data/seq_AS_fasta/

#M. fascicularis
#../../script/classification.py Macaca_fascicularis_SRA.fasta.regions_np_L30.monomers.noN.fst_149.161_182.fst --classification_output_file ../../results/Macaca_Fascicularis/Macaca_fascicularis_SRA_kmers52_pm090_fs0.0025.dat --verbose 9 --pairmate_threshold 0.90 --timing 1 --min_family_size 0.0025 &> ../../results/Macaca_Fascicularis/out/out_classification_pm090_fs00025.log &

#C. sabaeus
#../../script/classification.py Chlorocebus_sabaeus.fa.regions.fst_np_L30.monomers.noN.fst_149.161_182.fst --classification_output_file ../../results/Chlorocebus_Sabaeus/Chlorocebus_sabaeus.kmers52_pm090_fs0.0034.dat --verbose 9 --pairmate_threshold 0.90 --timing 1 --min_family_size 0.0034 &> ../../results/Chlorocebus_Sabaeus/out/out_classification_pm090fs00034.log &

#C. solatus
#../../script/classification.py Cercopithecus_solatus_XmnI.fst --classification_output_file ../../results/Cercopithecus_Solatus/Cercopithecus_solatus_kmers52_pm090_lda75000_fs0.00095.dat --verbose 9 --pairmate_threshold 0.90 --timing 1 --lda_sample_size 75000 --min_family_size 0.00095 &> ../../results/Cercopithecus_Solatus/out/out_classification_pm090_lda75000_fs0.00095.log &

###################################22/02/2018
#Classification pour solatus avec pm=0.95 et lda=75000
#../../script/classification.py Cercopithecus_solatus_XmnI.fst.kmers52 --classification_output_file ../../results/Cercopithecus_Solatus/Cercopithecus_solatus_kmers52_pm095_lda75000_fs0.00095.dat --verbose 9 --pairmate_threshold 0.95 --timing 1 --lda_sample_size 75000 --min_family_size 0.00095 &> ../../results/Cercopithecus_Solatus/out/out_classification_pm095_lda75000_fs0.00095.log

#Classification pour Macaca Assembly avec pm=0.95 et lda=75000
#../../script/classification.py Macaca_fascicularis.fa.regions.fst_np_L30.monomers.noN.fst_149.161_182.fst.kmers52 --classification_output_file ../../results/Macaca_Fascicularis/Macaca_fascicularis_kmers52_pm095_lda75000_fs0.00051.dat --verbose 9 --pairmate_threshold 0.95 --timing 1 --lda_sample_size 75000 --min_family_size 0.00051 &> ../../results/Macaca_Fascicularis/out/out_classification_pm095_lda75000_fs0.00051.log
################################## classification_weekend.sh
#RELANCEMENT
#../../script/classification.py Cercopithecus_solatus_XmnI.fst.kmers52  --classification_output_file ../../results/Cercopithecus_Solatus/Cercopithecus_solatus_kmers52_pm095_lda75000_fs0.00095.dat --verbose 9 --pairmate_threshold 0.95 --timing 1 --lda_sample_size 75000 --min_family_size 0.00095 &> ../../results/Cercopithecus_Solatus/out/out_classification_pm095_lda75000_fs0.00095.log

#RELANCEMENT
#../../script/classification.py Macaca_fascicularis.fa.regions.fst_np_L30.monomers.noN.fst_149.161_182.fst.kmers52 --classification_output_file ../../results/Macaca_Fascicularis/Macaca_fascicularis_kmers52_pm095_lda75000_fs0.00051.dat --verbose 9 --pairmate_threshold 0.95 --timing 1 --lda_sample_size 75000 --min_family_size 0.00051 &> ../../results/Macaca_Fascicularis/out/out_classification_kmers52_pm095_lda75000_fs0.00051.log

#../../script/classification.py Macaca_fascicularis_SRA.fasta.regions_np_L30.monomers.noN.fst_149.161_182.fst.kmers52 --classification_output_file ../../results/Macaca_Fascicularis/Macaca_fascicularis_SRA_kmers52_pm090_fs0.0025.dat --verbose 9 --pairmate_threshold 0.90 --timing 1 --min_family_size 0.0025 &> ../../results/Macaca_Fascicularis/out/out_classification_pm090_fs00025_again.log &

#############################################Lancement: 23/02/2018 MOKA
#solatus
#../../script/classification.py Cercopithecus_solatus_XmnI.fst.kmers52  --classification_output_file ../../results/Cercopithecus_solatus/Cercopithecus_solatus_kmers52_pm090_lda100000_fs0.00095.dat --verbose 9 --pairmate_threshold 0.90 --timing 1 --lda_sample_size 100000 --min_family_size 0.00095 &> ../../results/Cercopithecus_solatus/out/out_classification_pm090_lda100000_fs0.00095.log
#pogonias
#../../script/classification.py Cercopithecus_pogonias_XmnI.fst.kmers52  --classification_output_file ../../results/Cercopithecus_pogonias/Cercopithecus_kmers52_pm090_lda100000_fs0.00089.dat --verbose 9 --pairmate_threshold 0.90 --timing 1 --lda_sample_size 100000 --min_family_size 0.00089 &> ../../results/Cercopithecus_pogonias/out/out_classification_pm090_lda100000_fs0.00095.log
#fascicularis_assembly	
#../../script/classification.py Macaca_fascicularis.fa.regions.fst_np_L30.monomers.noN.fst_149.161_182.fst.kmers52 --classification_output_file ../../results/Macaca_fascicularis/Macaca_fascicularis_kmers52_pm090_lda100000_fs0.00051.dat --verbose 9 --pairmate_threshold 0.90 --timing 1 --lda_sample_size 100000 --min_family_size 0.00051 &> ../../results/Macaca_fascicularis/out/out_classification_pm090_lda100000_fs0.00051.log

##############################################Lancement: 26/02/2018 MOKA
#pogonias lda=100000 et pm=0.90
../../script/classification_moka.py Cercopithecus_pogonias_XmnI.fst.kmers52  --classification_output_file ../../results/Cercopithecus_pogonias/Cercopithecus_kmers52_pm090_lda100000_fs0.00089.dat --verbose 9 --pairmate_threshold 0.90 --timing 1 --lda_sample_size 100000 --min_family_size 0.00089 &> ../../results/Cercopithecus_pogonias/out/out_classification_pm090_lda100000_fs0.00095_BIS.log
#fascicularis lda=75000 et pm=0.95
../../script/classification.py Macaca_fascicularis.fa.regions.fst_np_L30.monomers.noN.fst_149.161_182.fst.kmers52 --classification_output_file ../../results/Macaca_fascicularis/Macaca_fascicularis_kmers52_pm095_lda_75000_fs0.00051.dat --verbose 9 --pairmate_threshold 0.95 --timing 1 --lda_sample_size 75000 --min_family_size 0.00051 &> ../../results/Macaca_fascicularis/out/out_classification_kmers52_pm095_lda_75000_fs0.00051_BIS.log
