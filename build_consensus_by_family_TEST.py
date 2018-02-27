#! /usr/bin/python3.4

import os
import sys
import time
import argparse
import tempfile
import numpy as np
#import networkx as nx
from Bio.Align.Applications import MuscleCommandline
import subprocess
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio.Align import AlignInfo
import random

#######################################################################
### read family file 
########################################################################
def get_families(famfile, verbose=1):
    ''' 
    '''

    faminfo = np.genfromtxt(famfile, skip_header=0, names=['sequence', 'family'], dtype=['S500','i8'])
    fam_info_dict = {}
    for pairs in faminfo:
        seq,fam=pairs
        if fam not in fam_info_dict:
            fam_info_dict[fam]=[seq]
        else:
            fam_info_dict[fam].append(seq)

    if verbose > 0:
        print("family file: number of sequences: %d" % (faminfo.shape[0]))
        print("family file: number of families: %d" % (len(fam_info_dict)))
		
    return(fam_info_dict)
###################################################################################								
def build_all_consensus(families, seq_file, max_seq=None, verbose=0):
	all_consensus = []
	os.system("makeblastdb  -logfile /dev/null -in "+seq_file+" -dbtype nucl -parse_seqids")
	
	for fam  in families.keys():
		cons_name = str(fam) + "_" + str(len(families[fam]))
		if verbose > 0:
			print("Consensus for %s ..." % (cons_name))
		if max_seq != None and len(families[fam]) > max_seq: 
			cons_name = cons_name + "_" + str(max_seq)
			all_consensus.append([build_consensus(random.sample(families[fam], max_seq), seq_file, verbose), cons_name])
		else:
			all_consensus.append([build_consensus(families[fam], seq_file, verbose), cons_name])
	return(all_consensus)
###################################################################################					
def extract_sequences_to_new_file(sequences, db, output_file, mode="w"):
		if mode != "a":
			if os.path.exists(output_file):
				os.remove(output_file)
		for seq in sequences:
			blastcmd='blastdbcmd '
			blastcmd=blastcmd+' -db ' + db
			blastcmd=blastcmd+' -entry ' + '"'+seq.decode('UTF-8')+'"'
			blastcmd=blastcmd+' -dbtype nucl'
			blastcmd=blastcmd+' -outfmt %f'
			blastcmd=blastcmd+' >> '+ output_file
			os.system(blastcmd)
		return(0)
###################################################################################									
def write_all_consensus_sequences(all_consensus, output_file, mode="w"):
	seq_file_hd = open(output_file, mode)
	for seq,name in all_consensus:
		SeqIO.write(SeqRecord(Seq(seq), id=name), seq_file_hd, "fasta")
	seq_file_hd.close()
###################################################################################								
def build_consensus(cluster, seq_file, verbose=0):
	tmp_file_in = tempfile.NamedTemporaryFile(delete=False).name
	tmp_file_out = tempfile.NamedTemporaryFile(delete=False).name
	if verbose > 3:
		print("Extracting sequences (%d sequences)..." %(len(cluster)))
	extract_sequences_to_new_file(cluster, seq_file, tmp_file_in, mode="w")
	
	if verbose > 3:
		print("Aligning sequences ...")
	muscle_cline = MuscleCommandline(input=tmp_file_in, out=tmp_file_out,clwstrict=False,quiet=True)
	code = subprocess.call(str(muscle_cline), shell=(sys.platform!="win32"))
	
	if verbose > 3:
		print("Reading alignment ...")
	alignment = AlignIO.read(tmp_file_out, "fasta")
	summary_align = AlignInfo.SummaryInfo(alignment)
	cons = summary_align.dumb_consensus()
	os.remove(tmp_file_in)
	os.remove(tmp_file_out)
	return(str(cons))

###################################################################################					
parser = argparse.ArgumentParser(description='Build multiple consensus sequences from a faily file and a fasta file ...', epilog='')
parser.add_argument('-s', dest="seq_file", required=True)
parser.add_argument('-f', dest="fam_file", required=True)
parser.add_argument('-m', dest="max_seq", type=int, default=None)
parser.add_argument('-o', dest="out_file", required=True)
parser.add_argument('-v', dest="verbose",   type=int, default=0)

args = parser.parse_args()

if not os.path.exists(args.seq_file):
	print("Error: sequence file ("+args.seq_file+") not found !!")
	quit()	
if not os.path.exists(args.fam_file):
	print("Error: family file ("+args.fam_file+") not found !!")
	quit()	
	
if args.verbose > 0:
	print(" ---> Reading families ...")
families = get_families(args.fam_file, args.verbose)

if args.verbose > 0:
	print(" ---> Building consensus ...")
all_consensus = build_all_consensus(families,  args.seq_file, args.max_seq, args.verbose)

if args.verbose > 0:
	print(" ---> Writing sequences ...")
write_all_consensus_sequences(all_consensus, args.out_file, "w")
