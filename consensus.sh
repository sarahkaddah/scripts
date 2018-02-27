#!/bin/bash

# Repertory: data/seq_AS_fasta/

# C. sabaeus
../../script/build_consensus_by_family_TEST.py -s Chlorocebus_sabaeus.fa.regions.fst_np_L30.monomers.noN.fst_149.161_182.fst -f ../../results/Chlorocebus_Sabaeus/Chlorocebus_sabaeus.kmers52_pm090_fs0.0034.dat -m 500 -o ../../results/Chlorocebus_Sabaeus/Chlorocebus_sabaeus.kmers52_pm090_fs0.0034.consensus -v 9 &> ../../results/Chlorocebus_Sabaeus/out/out_consensus_kmers52_pm090_fs0.0034.log

# C. solatus
../../script/build_consensus_by_family_TEST.py -s Cercopithecus_solatus_XmnI.fst -f ../../results/Cercopithecus_Solatus/Cercopithecus_solatus_kmers52_pm090_lda75000_fs0.00095.dat -m 500 -o ../../results/Cercopithecus_Solatus/Cercopithecus_solatus_kmers52_pm090_lda75000_fs0.00095.consensus -v 9 &> ../../results/Cercopithecus_Solatus/out/out_consensus_kmers52_pm090_fs0.00095.log

# C. solatus 2
../../script/build_consensus_by_family_TEST.py -s Cercopithecus_solatus_XmnI.fst -f ../../results/Cercopithecus_Solatus/Cercopithecus_solatus_kmers52_pm090_lda75000_fs0.00095.dat -m 500 -o ../../results/Cercopithecus_Solatus/Cercopithecus_solatus_kmers52_pm090_lda75000_fs0.00095_2.consensus -v 9 &> ../../results/Cercopithecus_Solatus/out/out_consensus_kmers52_pm090_fs0.00095_2.log

# C. solatus 3
../../script/build_consensus_by_family_TEST.py -s Cercopithecus_solatus_XmnI.fst -f ../../results/Cercopithecus_Solatus/Cercopithecus_solatus_kmers52_pm090_lda75000_fs0.00095.dat -m 500 -o ../../results/Cercopithecus_Solatus/Cercopithecus_solatus_kmers52_pm090_lda75000_fs0.00095_3.consensus -v 9 &> ../../results/Cercopithecus_Solatus/out/out_consensus_kmers52_pm090_fs0.00095_3.log

# C. pogonias
../../script/build_consensus_by_family_TEST.py -s Cercopithecus_pogonias_XmnI.fst -f ../../results/Cercopithecus_Pogonias/Cercopithecus_pogonias_kmers52.pm090_lda75000_fs000089.dat -m 500 -o ../../results/Cercopithecus_Pogonias/Cercopithecus_pogonias_kmers52.pm090_lda75000_fs000089.consensus -v 9 &> ../../results/Cercopithecus_Pogonias/out/out_consensus_kmers52_pm090_fs000089.log

# M. fascicularis assembly
../../script/build_consensus_by_family_TEST.py -s Macaca_fascicularis.fa.regions.fst_np_L30.monomers.noN.fst_149.161_182.fst -f ../../results/Macaca_Fascicularis/Macaca_fascicularis_pm090_lda75000_fs0.00051.dat -m 500 -o ../../results/Macaca_Fascicularis/Macaca_fascicularis_kmers52_pm090_fs0.00051.consensus -v 9 &> ../../results/Macaca_Fascicularis/out/out_consensus_kmers52_pm090_fs0.00051.log
