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
version=1


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
### get distances form distmat
########################################################################
def get_distances_with_distmat(alignment_file, distance_file=None, distance=0, verbose=0):
	if distance_file ==  None:
		distance_file = tempfile.NamedTemporaryFile(delete=False).name
	os.system("distmat "+ alignment_file+" -filter "+" -nucmethod "+ str(distance) +" "+distance_file)
	distmat = read_distmat(distance_file, verbose)
	return(distmat)
########################################################################
### get distances form distmat
########################################################################
def get_inter_distances_with_distmat(fam1, fam2, alignment_file, distance_file=None, distance=0, verbose=0):
	if distance_file ==  None:
		distance_file = tempfile.NamedTemporaryFile(delete=False).name
	os.system("distmat "+ alignment_file+" -filter "+" -nucmethod "+ str(distance) +" "+distance_file)
	distmat = read_distmat(distance_file, verbose)
	return(distmat)


########################################################################
### get distances for one fasta file
########################################################################
def read_distmat(distmat_file, verbose):
	distmat_hdl = open(distmat_file, "r")
	for ii in range(0,7):
		distmat_hdl.readline()
	nb_seq = len(distmat_hdl.readline())
	nb = 0
	distmat = []
	distname = []
	for line in distmat_hdl:
		sl = line.split()
		for ii in range(0, len(sl)-2):		
			distmat.append([nb, ii+nb, float(sl[ii])])
		distname.append(sl[ii+1])
		nb = nb + 1
	for jj in range(0, len(distmat)):
		distmat[jj][0] = distname[distmat[jj][0]]
		distmat[jj][1] = distname[distmat[jj][1]]
	return(distmat)

########################################################################
### get distances for one fasta file
########################################################################
def get_inter_distances(fam1, fam2, fasta_file, distance, verbose):
	aligned_file = get_multiple_alignment(fasta_file, verbose=verbose)
	distances = get_inter_distances_with_distmat(fam1, fam2, aligned_file, distance_file=None, distance=0, verbose=0)
	return(distances)

########################################################################
### get distances for one fasta file
########################################################################
def get_intra_distances(fasta_file, distance, verbose):
	aligned_file = get_multiple_alignment(fasta_file, verbose=verbose)
	distances = get_distances_with_distmat(aligned_file, distance_file=None, distance=0, verbose=0)
	return(distances)
########################################################################
### get distances for all families
########################################################################
def get_intra_fam_distances(faminfo, sequences, distance, nb_seq=0, verbose=1):
	all_distances = {} 
	tmp_seqfile = tempfile.NamedTemporaryFile(delete=False).name
	for fam in faminfo.keys():
		write_sequences(faminfo[fam], sequences, tmp_seqfile)
		all_distances[fam] = get_intra_distances(tmp_seqfile, distance, verbose)
	return(all_distances)
########################################################################
### get distances for all families
########################################################################
	
def get_inter_fam_distances(faminfo, sequences, distance, nb_seq=0, verbose=1):
	all_distances = {} 
	tmp_seqfile = tempfile.NamedTemporaryFile(delete=False).name
	for fam1 in faminfo.keys():
		for fam2 in faminfo.keys():
			#if fam1 != "all" and fam2 != "all" and fam1 < fam2:
			if fam1 < fam2:
				write_sequences(faminfo[fam1]+faminfo[fam2], sequences, tmp_seqfile)
				all_distances[str(fam1)+"___"+str(fam2)] = get_inter_distances(faminfo[fam1], faminfo[fam2], tmp_seqfile, distance, verbose)
	return(all_distances)
	
	
########################################################################
### get sequences
########################################################################
def get_sequences(seq_file, verbose=1):
	all_seq = SeqIO.to_dict(SeqIO.parse(seq_file, "fasta"))
	return(all_seq)
#########################################################################
### randomly select sequences into families
########################################################################
def filter_seqfam(faminfo, nbseq=0, verbose=1):
	if nbseq == 0:
#		faminfo['all'] = list(sum(faminfo.values(), []))
		return(faminfo)
	filter_infofam = {}
	#tmp = list(sum(faminfo.values(), []))
	#if len(tmp) > nbseq:
		#filter_infofam['all'] = random.sample(tmp, nbseq)
	#else:
		#filter_infofam['all'] = tmp
		
	for fam in faminfo.keys():
		init_len = len(faminfo[fam])
		if len(faminfo[fam]) > nbseq:
			filter_infofam[fam] = random.sample(faminfo[fam], nbseq)
		else:
			filter_infofam[fam] = faminfo[fam]
		if verbose > 5:
			print("%s %d %d" %(fam, init_len, len(faminfo[fam])))
	return(filter_infofam)
#########################################################################
### randomly select sequences into families
########################################################################
def remove_small_families(faminfo, minimal=2, verbose=1):
	filter_infofam = {}
	for fam in faminfo.keys():
		if len(faminfo[fam]) >= minimal:
			filter_infofam[fam] = faminfo[fam]
	return(filter_infofam)
	
	
#######################################################################
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
def	write_inter_data(faminfo, all_distances, out_file, verbose):
	outf = open(out_file, 'w')
	outf.write("family1 total_nb_seq1 family2 total_nb_seq2 nb_dist_inter median_inter mean_inter sd_inter nb_dist_intra1 median_intra1 mean_intra1 sd_intra1 nb_dist_intra2 median_intra2 mean_intra2 sd_intra2\n")
	#nb_seq=0
	inv_faminfo = {}
	for k, v in faminfo.items():
		for ii in range(0,len(v)):
			inv_faminfo[v[ii]]=k
	
	for fam1 in faminfo.keys():
		for fam2 in faminfo.keys():
			if fam1 < fam2:
				re1=re.compile("^.*_"+str(fam1)+"$")
				re2=re.compile("^.*_"+str(fam2)+"$")
				tmp = []
				tmp1 = []
				tmp2 = []
				#nb_seq=nb_seq+len(faminfo[fam])
				kk=str(fam1) + "___" + str(fam2)
				if kk in all_distances.keys():
					for ii in range(0, len(all_distances[kk])):
						if all_distances[kk][ii][0] != all_distances[kk][ii][1]:
							if   inv_faminfo[all_distances[kk][ii][0]] == fam1 and inv_faminfo[all_distances[kk][ii][1]] == fam1:
								tmp1.append(all_distances[kk][ii][2])
							elif inv_faminfo[all_distances[kk][ii][0]] == fam2 and inv_faminfo[all_distances[kk][ii][1]] == fam2:
								tmp2.append(all_distances[kk][ii][2])
							else:
								tmp.append(all_distances[kk][ii][2])
					outf.write("%s %d %s %d %d %f %f %f %d %f %f %f %d %f %f %f\n" %(str(fam1), len(faminfo[fam1]), str(fam2), len(faminfo[fam2]), len(tmp), np.median(tmp), np.mean(tmp), np.std(tmp), len(tmp1), np.median(tmp1), np.mean(tmp1), np.std(tmp1), len(tmp2), np.median(tmp2), np.mean(tmp2), np.std(tmp2)))
				else:
					outf.write("%s %d %s %d NA NA NA NA NA NA NA NA NA NA NA NA\n" %(str(fam1), len(faminfo[fam1]), str(fam2), len(faminfo[fam2])))
	outf.close()

#######################################################################
### rename seq
########################################################################
def rename_sequences(sequences):
	regex = re.compile(r".*lcl\|")
	seq2={}
	for name in sequences:
		seq2[re.sub(regex, "", name)] = sequences[name]
	return(seq2)
#######################################################################
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
    parser.add_argument("-n", "--nb_seq", 
						help = "Maximum number of sequences used for the analysis", 
						type = int, 
						default = 0)
    parser.add_argument("-m", "--minseq", 
						help = "Minimal number of sequences in families", 
						type = int, 
						default = 2)
    parser.add_argument("-d", "--distance", help = "Distance", type = int, default = 0)
    parser.add_argument( "--intra_only", help = "Type of distances", action = 'store_true', default = False)
    parser.add_argument("--verbose","-v", 
						help ="Degree of verbosity",
                        type=int, 
                        default=0)                  
    parser.add_argument("--version", 
						help ="Print the version number", 
						action='version', 
						version='Version number: '+ str(version))
    args = parser.parse_args()
    if  args.nb_seq  < 0 :
        print("Error:  maximum number of sequences should be in  [1-infinity] (value: %d)" % (args.nb_seq))
        exit(0)   
    if  args.minseq  < 2 :
        print("Error:  Minimal number of sequences should be in  [2-infinity] (value: %d)" % (args.nb_seq))
        exit(0)   
        
             
    return(args.sequence_file,
			args.family_file,
			args.intra_only,
			args.distance,
			args.output,
			args.nb_seq,
			args.minseq,
			args.verbose
			)
########################################################################
### Getting arguments
########################################################################

(seq_file, fam_file, intra_only, distance, out_file, nbseq, minseq, verbose) = get_args()


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

if verbose > 2:
	print("Filtering families ...")
faminfo3 = filter_seqfam(faminfo, nbseq, verbose)

if verbose > 2:
	print("Removing small families ...")
faminfo2 = remove_small_families(faminfo3, minseq, verbose)

if verbose > 2:
	print("Reading sequences ...")
sequences = get_sequences(seq_file)
sequences = rename_sequences(sequences)

if intra_only == True:
	if verbose > 2:
		print("Getting intra-families distances ...")
	all_intra_distances = get_intra_fam_distances(faminfo2, sequences, distance, verbose)
	write_intra_data(faminfo, all_intra_distances, out_file, verbose)
else:
	if verbose > 2:
		print("Getting intra and inter-families distances ...")
	all_inter_distances = get_inter_fam_distances(faminfo2, sequences, distance, verbose)
	write_inter_data(faminfo, all_inter_distances, out_file, verbose)
