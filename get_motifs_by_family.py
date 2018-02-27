#! /usr/bin/python3.4

# Loic Ponger <loic.ponger@mnhn.fr>

#################################################################################

#This work is licensed under the Creative Commons Attribution-ShareAlike 3.0
#Unported License. To view a copy of this license,
#visit http://creativecommons.org/licenses/by-sa/3.0/
# or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.

#################################################################################

# A computer program does what you tell it to do, not what you want it to do
#-- Greer's Law

# example : 

import os, sys, re
import numpy as np
#import matplotlib.pyplot as plt
import time
import random
from sklearn.decomposition import PCA, IncrementalPCA
import scipy
import argparse
import fastcluster
import scipy.cluster.hierarchy
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align.Applications import MuscleCommandline
from collections import Counter
import logging
import tempfile
import numpy as np
version=0.5


########################################################################
### write sequences 
########################################################################
def write_sequences(names, sequences, outfile):
	seqfile_hdl = open(outfile, "w")
	nb = 0
	for seq in names:
		nb = nb + 1
		SeqIO.write(sequences[seq], seqfile_hdl, "fasta")
	seqfile_hdl.close()
	return(nb)

########################################################################
### multiple alignment of a fasta file
########################################################################
def get_multiple_alignment(seq_file, alignment_file=None, verbose=9):
	muscle_exe = 'muscle'
	if alignment_file ==  None:
		out_file = tempfile.NamedTemporaryFile(delete=False).name
	else:
		out_file = alignment_file
		
	muscle_cline = MuscleCommandline(muscle_exe, input=seq_file, out=out_file)
	if verbose > 6:
		print("      muscle command line:")
		print(muscle_cline)
	stdout, stderr = muscle_cline()	
#	MultipleSeqAlignment = AlignIO.read(out_file, "fasta")
#	if alignment_file ==  None:
#		os.remove(out_file)
#	return(MultipleSeqAlignment)
	return(out_file)

########################################################################
### get patterns by using fuzznuc
########################################################################
def search_pattern(seq_file, pattern, mismatch=0, complement=True, verbose=0):
	regex = re.compile(r".*lcl\|")
	fuzznuc_file = tempfile.NamedTemporaryFile(delete=False).name
	if complement == False:
		os.system("fuzznuc "+ seq_file+" "+ fuzznuc_file+" -pattern "+ pattern + " -pmismatch "+ str(mismatch))
	else:
		os.system("fuzznuc "+ seq_file+" "+ fuzznuc_file+" -pattern "+ pattern + " -pmismatch "+ str(mismatch)+" -complement")
	fuzznuc = read_fuzznuc(fuzznuc_file, verbose)
	all_seq = SeqIO.to_dict(SeqIO.parse(seq_file, "fasta")).keys()
	nb0=len(fuzznuc.keys())
	nb=0
	for ss in all_seq:
		ss=re.sub(regex, "", ss)
		if ss not in fuzznuc.keys():
			nb=nb+1
			fuzznuc[ss] = []
	if verbose > 2:
		print("Number of sequences: %d" % (len(all_seq)))
		print("Number of sequences with pattern: %d" % (nb0))
		print("Number of sequences without pattern: %d" % (nb))
	if verbose > 8:
		print("Fuzznuc temporary file: %s" % (fuzznuc_file))
	else:
		os.remove(fuzznuc_file)
	return(fuzznuc)
########################################################################
### get distances form distmat
########################################################################
def read_fuzznuc(fuzznuc_file, verbose=0):
	fuzz_hdl = open(fuzznuc_file, "r")
	nb = 0
	output = {}
	for line in fuzz_hdl:
		sl=re.split("[\u200b\s]+", line, flags=re.UNICODE)		
		if len(sl) > 1 and sl[1]== "Sequence:":
			seq=sl[2]
			output[seq] = []
		if len(sl) == 8 and (sl[3] == "+" or sl[3] == "-"):
			if sl[5] == ".":
				sl[5]="0"
			output[seq].append((int(sl[1]),int(sl[2]), sl[3], int(sl[5]), sl[6]))
	return(output)
####################################################################
### write data
########################################################################
def	write_intra_data(faminfo, all_distances, out_file, verbose):
	outf = open(out_file, 'w')
	outf.write("family total_nb_seq nb_distances median mean sd\n")
	nb_seq=0
	for fam in faminfo.keys():
		tmp = []
		nb_seq=nb_seq+len(faminfo[fam])
		if fam in all_distances.keys():
			for ii in range(0, len(all_distances[fam])):
				if all_distances[fam][ii][0] != all_distances[fam][ii][1]:
					tmp.append(all_distances[fam][ii][2])
			outf.write("%s %d %d %f %f %f\n" %(str(fam), len(faminfo[fam]), len(all_distances[fam]), np.median(tmp), np.mean(tmp), np.std(tmp)))
		else:
			outf.write("%s %d NA NA NA NA\n" %(str(fam), len(faminfo[fam])))
	outf.close()
#######################################################################
### write data
########################################################################
def	write_fam_patterns(fam_patterns, out_file, verbose):
	outf = open(out_file, 'w')
	outf.write("family1 nb_seq")
	for pat in fam_patterns.keys():
		outf.write(" %s" %(pat))
	outf.write("\n")
	refpat=pat
	for fam in fam_patterns[refpat].keys():
		outf.write("%s %d" %(fam, fam_patterns[refpat][fam][1]))
		for pat in fam_patterns.keys():
			outf.write(" %f" %(fam_patterns[pat][fam][0]))
		outf.write("\n")
	outf.close()
######################################################################
### statistics about patterns into families 
########################################################################
def stat_pattern(seq_pattern, fam_info, verbose=0):
	output= {}
	for fam in fam_info.keys():
		seq_with_pattern = 0
		nb_seq = len(fam_info[fam])
		for ss in fam_info[fam]:
			if len(seq_pattern[ss]) > 0:
				seq_with_pattern = seq_with_pattern + 1
		output[fam] = (seq_with_pattern / nb_seq, nb_seq)
	return(output)
#######################################################################
### rename seq
########################################################################
def rename_sequences(sequences):
	regex = re.compile(r".*lcl\|")
	seq2={}
	for name in sequences:
		seq2[re.sub(regex, "", name)] = sequences[name]
	return(seq2)
#########write_fam_pattern##############################################################
### read family file 
########################################################################
def get_families(famfile, verbose=1):
    ''' 
    '''
    regex = re.compile(r".*lcl\|")
    faminfo = np.genfromtxt(famfile, skip_header=0, names=['sequence', 'family'], dtype=['S500','i8'])
    fam_info_dict = {}
    for pairs in faminfo:
        seq,fam=pairs
        seq=re.sub(regex, "", seq.decode('utf-8'))
        if fam not in fam_info_dict:
            fam_info_dict[fam]=[seq]
        else:
            fam_info_dict[fam].append(seq)

    if verbose > 0:
        print("family file: number of sequences: %d" % (faminfo.shape[0]))
        print("family file: number of families: %d" % (len(fam_info_dict)))
		
    return(fam_info_dict)
  
########################################################################
### Getting arguments
########################################################################

def get_args():
    ''' get the argument from the command line
        NO ARG
        RETURN : a tupple with the different args used in the programm
    '''
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description='%(prog)s can be used to ...')

    parser.add_argument("-s","--sequence_file", 
					help = "Seuqences in fasta format", 
					default = None, 
					required=True)
    parser.add_argument("-f","--family_file", 
					help = "Tabular file with family information", 
					default = None, 
					required=True)
    parser.add_argument("-o", "--output", 
						help = "Ouput file",
                        default=None, 
                        required=True)
    parser.add_argument('-p','--pattern', action='append', help='pattern(s)', required=True)
	
#                        
    parser.add_argument("--verbose","-v", 
						help ="Degree of verbosity",
                        type=int, 
                        default=0)                  
    parser.add_argument("--version", 
						help ="Print the version number", 
						action='version', 
						version='Version number: '+ str(version))
    args = parser.parse_args()
    patterns=[]
    mismatchs=[]
    complements=[]
    for ii in range(0,len(args.pattern)):
         l=args.pattern[ii].split(":")
         if len(l) != 3:
             print("Error:  pattern syntaxe should be: pattern:mismatchs:complement (%s)" % (args.pattern[ii]))
             exit(0)   
 		
         patterns.append(l[0])
         l[1]=int(l[1])
         if l[1] < 0:
            print("Error:  number of mismatchs should be [0-infinity] (%d)" % (l[1]))
            exit(0)   		 
         mismatchs.append(l[1])
         complements.append(str2bool(l[2]))
         
    return(args.sequence_file,
			args.family_file,
			args.pattern,
			patterns,
			mismatchs,
			complements,
			args.output,
			args.verbose
			)
########################################################################
def str2bool(v):
  return v.lower() in ("complement", "true", "1")
########################################################################
### Getting arguments
########################################################################

(seq_file, fam_file, full_patterns, patterns, mismatchs, complements, out_file, verbose) = get_args()






########################################################################
### Checking input/output files
########################################################################
if os.path.exists(seq_file) ==  False:
	print("Error: input file %s doesn't exist" % (seq_file))
	print("Bye bye :-(")
	quit(-1)
if os.path.exists(fam_file) ==  False:
	print("Error: input file %s doesn't exist" % (fam_file))
	print("Bye bye :-(")
	quit(-1)
if verbose > 2:
	print("Reading family information ...")
faminfo = get_families(fam_file, verbose)
seq_patterns = {}
fam_patterns = {}
for ii in range(0,len(patterns)):
	key = full_patterns[ii]
	
	if verbose > 2:
		print("searching pattern (%s-%d-%s)..." % (patterns[ii],mismatchs[ii], complements[ii]))	
	seq_patterns[key]= search_pattern(seq_file, patterns[ii], mismatchs[ii], complements[ii], verbose)
	fam_patterns[key] = stat_pattern(seq_patterns[key], faminfo, verbose)
write_fam_patterns(fam_patterns, out_file, verbose)
	
	
