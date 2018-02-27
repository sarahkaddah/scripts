#!/bin/bash

# Répertoire: lancer depuis data/

#Parametres:
	# 1-Chemin/Fichier fasta à parser
	# 2-Chemin/Fichier de classification
	# 3-Chemin/output (radical)
	# 4-chemin/out/output

# Recherche motif CENP-B:
../script/get_motifs_by_family.py -s $1 -f $2 -p TTCGTTGGAARCGGGA:2:complement  -o $3.CENPB.motif -v 9 &> $5.cenpb.log &

# Recherche motif pJalpha:
../script/get_motifs_by_family.py -s $1 -f $2 -p TTCCTTTTYCACCRTAG:2:complement  -o $3.pJalpha.motif -v 9 &> $5.pjalpha.log &

#Recherche motif pKbeta:
../script/get_motifs_by_family.py -s $1 -f $2 -p CTATAGGGCCAAAGGAA:2:complement  -o $3.pKbeta.motif -v 9 &> $5.pkbeta.log &



