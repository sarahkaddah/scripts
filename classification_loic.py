#! /usr/bin/python3.4

# Florence Jornod <florence@jornod.com>
# Loic Ponger <loic.ponger@mnhn.fr>

#################################################################################

#This work is licensed under the Creative Commons Attribution-ShareAlike 3.0
#Unported License. To view a copy of this license,
#visit http://creativecommons.org/licenses/by-sa/3.0/
# or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.

#################################################################################

# A computer program does what you tell it to do, not what you want it to do
#-- Greer's Law

# example : ./classification.py ../data/nt011630.rev.monomers.fst_149.length_166_177.kmers5 -o res.txt -- --pairmate 0.90



# modify lda solver to lsqr to avoid warning message for colinearity

import os, sys
import numpy as np
import matplotlib.pyplot as plt
import time
import random
from sklearn.decomposition import PCA, IncrementalPCA
import scipy
import argparse
import fastcluster
import scipy.cluster.hierarchy
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from Bio import SeqIO
from collections import Counter
import logging

version=1

########################################################################
### Writing/saving the classification results
########################################################################
def write_classification(results, sequence_names, output_file, sort_by_sequence_name=True):
	if sort_by_sequence_name == True:
		save_results_sorted_by_sequence(output_file, results, sequence_names)
	else:
		save_results_sorted_by_cluster(output_file, results, sequence_names)

########################################################################
### The classification method
########################################################################
def classification(kmer_table, lda_sample_size, verbose=0):
	cluster_number = 0
	ii = 0
	res = {}
	queue = []
	queue.append(range(len(kmer_table)))
	#familysize = int(maxfam * len(kmer_table))
	familysize = int(maxfam)
	#filename = open(outputfile,"w")
	#filename.write("input filename :" +kmer_file+" pairmate : " + str(pairmate) + \
	 #" nb comp pca : " + str(nbcomp_pca) + " size min family : " + str(familysize)\
	 #+ " lda :"+ lda + "\n")
	
	while(queue):
	    if verbose > 1:
	    	print("___________________")
	    	print("Queue length: %d" % (len(queue)))
	    kmer_indice = queue.pop(0) # pop the first value of the queue
	    pca, expl_var_ratio = compute_pca(kmer_table[kmer_indice])
	    nb_comp_dispo = len(pca[0])
	    if verbose > 4:
	        print(" Turn %d with %d sequences" % (ii, len(pca)))
	    # to reduce computation time, we use just nbcomp_pca components if possible
	    if(len(kmer_indice) > lda_sample_size):
	        sample=range(len(kmer_indice))
	        kmer_sample = random.sample(sample, lda_sample_size)
	    else:
	        kmer_sample=range(len(kmer_indice))

	    if pca.shape[1] > nbcomp_pca:
	        pca=pca[:,range(nbcomp_pca)]
	        
	    distance, dist_time = compute_dist2(pca[kmer_sample,:])
	    if verbose > 4:
	        print("   Ratio of explicative variable: %f" % (expl_var_ratio))
	        print("   Distance time: %f s" % (dist_time))
	        print("   Number of composantes: %d s" % (nb_comp_dispo))
	        print("   Number of sequences: %d s" % (len(pca)))
	        
	    clusters = compute_clustering_fast(distance)
	    if(len(kmer_indice) > lda_sample_size):
	        diff = list( set(kmer_indice) - set(kmer_sample) )
	        clf = LinearDiscriminantAnalysis(solver='lsqr',shrinkage='auto')      
	        clf.fit(pca[kmer_sample,:], clusters)
	        clusters2 = clf.predict(pca[range(len(kmer_indice)),:])
	    else:
	        clusters2=clusters
	    mp1, mp2 = compute_pairmate2(clusters, distance)
	    group1, group2 = [],[]
	    if verbose > 4:
	        print("   Mate pair values: %f and %f" % (mp1, mp2))
	    if verbose > 0:
	        logging.info("mp "+str(mp1)+" "+str(mp2))
	    
	    if mp1 >= pairmate and mp2 >= pairmate:
	        if max(clusters2) == min(clusters2):
	            exit("LDA problem")
	        group1, group2 = create_2_groups(clusters2, kmer_indice)
	        if(len(group1)>familysize):
	            queue.append(group1)
	        else:
	            cluster_number = cluster_number+1
	            res[cluster_number] = group1
	        if(len(group2)>familysize):
	            queue.append(group2)
	        else:
	            cluster_number = cluster_number+1
	            res[cluster_number] = group2
	    else:
	        cluster_number = cluster_number+1
	        if verbose > 1:
	            logging.info("matepair test non ok")
	        res[cluster_number] = kmer_indice
	    ii += 1
	return(res)




def gen_sub(s, len_chunk):
   	for start in range(0, len(s)-len_chunk+1):
       		yield s[start:start+len_chunk]

def reading_fasta_file(input_fasta_file, verbose=0):
	handle1=open(input_fasta_file, 'r')
	seq_dict = SeqIO.to_dict(SeqIO.parse(handle1, "fasta"))
	handle1.close()
	if verbose > 0:
		print("Nb of sequences: %d" % (len(seq_dict)))
	return(seq_dict)
	
def filtering_sequences_with_N(seq_dict, verbose=0):
	select_seq = {}
	for ss1 in seq_dict:
		if seq_dict[ss1].seq.count('[N]') == 0 and seq_dict[ss1].seq.count('n') == 0:
			select_seq[ss1] = seq_dict[ss1]
	return(select_seq)

def filtering_sequences_by_length(seq_dict, minimal_length=None, maximal_length=None, verbose=0):
	if minimal_length == None and maximal_length == None:
		return(select_seq)
	select_seq = {}
	if minimal_length != None and maximal_length != None:
		for ss1 in seq_dict:
			if len(seq_dict[ss1].seq) >= minimal_length and len(seq_dict[ss1].seq) <= maximal_length:
				select_seq[ss1] = seq_dict[ss1]
	if minimal_length == None and maximal_length != None:
		for ss1 in seq_dict:
			if len(seq_dict[ss1].seq) <= maximal_length:
				select_seq[ss1] = seq_dict[ss1]
	if minimal_length != None and maximal_length == None:
		for ss1 in seq_dict:
			if len(seq_dict[ss1].seq) >= minimal_length:
				select_seq[ss1] = seq_dict[ss1]
	if verbose > 0:
		print("Nb of selected sequences: %d" % (len(seq_dict)))
	return(select_seq)
	
	
def compute_kmers_table_from_fasta_file(input_fasta_file, kmer_length, \
	output_kmer_file, kmer_relative_frequencies, \
	kmer_minimal_frequency=0.0, minimal_length=None, maximal_length=None, \
	filter_N=True, \
	verbose=0):
	seq_dict = reading_fasta_file(input_fasta_file, 0)
	if verbose > 0:
		print("Kmer frequencies, nb of sequences: %d" % (len(seq_dict)))
		
	if minimal_length != None or maximal_length != None:
		seq_dict = filtering_sequences_by_length(seq_dict, minimal_length, maximal_length, 0)
		if verbose > 0:
			print("Kmer frequencies, nb of sequences: %d (after length selection)" % (len(seq_dict)))

	if filter_N == True:
		seq_dict = filtering_sequences_with_N(seq_dict, 0)
		if verbose > 0:
			print("Kmer frequencies, nb of sequences: %d (after N filtering)" % (len(seq_dict)))

	all_words=[]
	kmers_count={}
	nb_kmers_count={}
	nb=0
	t1=time.time()
	for ss1 in seq_dict:
		nb=nb+1
		#if (nb/len(seq_dict))%10 == 0:
		#print(nb)
		kmers_count[ss1]=Counter([sub for sub in gen_sub(str(seq_dict[ss1].seq.upper()),kmer_length )])
		all_words= all_words+ list(kmers_count[ss1].keys())
		all_words=list(set(all_words))
		if kmer_relative_frequencies:
			for ii in kmers_count[ss1]:
				kmers_count[ss1][ii]=float(kmers_count[ss1][ii])/(len(seq_dict[ss1].seq)-kmer_length+1.0)
	t2=time.time()
	print(t2-t1)
	all_words=sorted(all_words)
	handle=open(output_kmer_file,'w')
	handle.write("id")
	nb_kmers_count= dict(zip(all_words , [0.0] * len(all_words)))
	if kmer_minimal_frequency > 0.0:
		if verbose > 0:
			print("Kmer frequencies, initial number of words: "+str(len(all_words))	)
	else:
		if verbose > 0:
			print("Kmer frequencies, number of words: "+str(len(all_words))	)
		
	if kmer_minimal_frequency > 0.0:
		for ss1 in seq_dict:
			for kkk in all_words:
				if kmers_count[ss1][kkk] > 0:
					nb_kmers_count[kkk] = nb_kmers_count[kkk] + 1
					
		for kkk in all_words:
			nb_kmers_count[kkk] = nb_kmers_count[kkk] / len(seq_dict)
			if nb_kmers_count[kkk] < kmer_minimal_frequency:
				if verbose > 2:
						print("Kmer frequencies, removing: "+kkk+" with frequency: "+str(nb_kmers_count[kkk]))
				all_words=list(filter(lambda a: a != kkk, all_words))	
	t3=time.time()
	print(t3-t2)
						
	if kmer_minimal_frequency > 0.0:
		if verbose > 0:
			print("Kmer frequencies, final number of words: "+str(len(all_words))	)
	for kkk in all_words:
		handle.write (" "+kkk)
	handle.write ("\n")
	
	for ss1 in seq_dict:
		handle.write(ss1)
		for kkk in all_words:
			if kkk in kmers_count[ss1]:
				handle.write (" "+str(kmers_count[ss1][kkk]))
			else:
				if kmer_relative_frequencies:
					handle.write (" 0.0")
				else:
					handle.write (" 0")
		handle.write ("\n")
	if kmer_relative_frequencies == True:
			handle.write ("frequency")
			for kkk in all_words:
				handle.write (" "+str(nb_kmers_count[kkk]))
			handle.write ("\n")
	handle.close()
	if verbose > 0:
		print("Kmer frequencies, saved in: %s "	% (output_kmer_file))
	t4=time.time()
	print(t4-t3)
	print(t4-t1)
	

	return(output_kmer_file)




























def get_kmer_table(kmer_file, verbose=1):
    ''' get the names of each sequence which correspond to the kmer_file except
        the first column. Column names are lost.
        ARG    : - kmer_file : the filename with the kmer table
                 - verbose : verbosity

        RETURN : - a 2D numpy.array with the kmer frequencies
    '''
    kfile=open(kmer_file,"r")
    line=kfile.readline()
    nb_col=len(line.split())-1
    kfile.close()
    kmer_table = np.genfromtxt(kmer_file, skip_header=1, usecols = range(1,nb_col+1))
    if verbose > 0:
        print("Kmer table, number of sequences: %d" % (kmer_table.shape[0]))
        print("Kmer table, number of words: %d" % (kmer_table.shape[1]))
    return(kmer_table)
    
def get_sequences_name(kmer_file, verbose=1):
    ''' get the names of each sequence which correspond to the first column of
        the kmer_file
        ARG    : - kmer_file : the file with the kmer table
        RETURN : - a 1D numpy.array with the sequences' name
    '''
    name_table = np.genfromtxt(kmer_file, skip_header=1, usecols=0, dtype=str)
    if verbose > 0:
         print("Kmer table, number of sequence names: %d" % (name_table.shape[0]))
    return(name_table)
    
def compute_pca(kmer_table):
    ''' compute a Principal component analysis of the kmer frequencies
        ARG    : - kmer_table : an 2D array with the kmer frequencies
        RETURN : - An 2D numpy.array with the eigenvalues
    '''
    pca = PCA()
    if verbose > 2:
        logging.info("compute pca")
    pca.fit(kmer_table)
    expl_var_ratio = sum(pca.explained_variance_ratio_[0:nbcomp_pca])
    return pca.fit_transform(kmer_table), expl_var_ratio


def compute_dist2(res_pca, method = "euclidean"):
    ''' compute the distance between each sequence using the eigenvalues of the
        PCA. By default the method used is "euclidean"
        ARG    : - res_pca : an 2D array with the eigenvalues
                 - method  : the method used to compute the distance
        RETURN : - An 2D array with pairwise distances between sequences
    '''
    t1=time.clock()
    if verbose > 2:
        logging.info("compute distances with "+ method + " method")
    d =  scipy.spatial.distance.pdist(res_pca, method)
    t2=time.clock()
    dist_time = t2 - t1
    return d, dist_time

def compute_clustering_fast(distance):
    t1=time.clock()
    c=fastcluster.ward(distance)
    t2=time.clock()
    return scipy.cluster.hierarchy.fcluster(c, 2, criterion="maxclust")

def plot_clustering( res_pca, group1, group2, output):
    ''' plot the 2 first dimension of the pca colored by cluster.
        ARG    : - res_pca        : an 2D array with the eigenvalues
                 - group1, group2 : vector of indice
                 - output         : the name of the output
        NO RETURN
    '''
    output+="clustering.png"
    if verbose > 2:
        logging.info("plot clustering")
    plt.figure()
    for ind in range(len(res_pca)):
        if ind in group1:
            plt.scatter(res_pca[ind,0],res_pca[ind,1],c="r",alpha=0.2)
        elif ind in group2:
            plt.scatter(res_pca[ind,0],res_pca[ind,1],c="b",alpha=0.2)
    plt.savefig(output)
    plt.close()

def create_2_groups(cut,kmer_indice):
    ''' Create two groups following clustering results.
        ARG    : - cut : a rpy2 IntVector with the repartition of the 2 groups
        RETURN : 2 lists one for each group
    '''
    ii = 0
    plot_output="plot_"+str(ii)
    if verbose > 1:
        logging.info("matepair test ok")
    g1,g2 = [],[]
    for x in range(len(cut)):
        (g1,g2)[cut[x] == 1].append(kmer_indice[x])
    if verbose > 0:
        logging.info("length "+str(len(g1))+" "+str(len(g2)))
    return g1, g2


def compute_pairmate(cut, distance):
    ''' compute the pairmate for the two groups.Pairmate corresponds to the
        proportion of sequences which has their nearest neighbor in a same group

        ARG    : - cut      : a rpy2IntVector with the repartition of the groups
                 - distance : An 2D array with distances between each sequence
        RETURN : - a tupple with the matepair of each group
    '''
    if verbose > 2:
        logging.info("compute pairmate")
    distance[distance<=0.00000000001] = 999
    res = []
    print(max(cut),min(cut))
    if(max(cut)==min(cut)):
        print(cut)
        return(0,0)
    for i in range(len(distance)):
        j = np.argmin(distance[i,])
        res.append(str(cut[i])+str(cut[j]))
        if i == 1:
            print(min(distance[i,]))
    #    print(j)
    mp1 = float(res.count('11'))/(res.count('11')+res.count('12'))
    mp2 = float(res.count('22'))/(res.count('21')+res.count('22'))
    return mp1,mp2,res

def compute_pairmate2(cut, distance):
    ''' compute the pairmate for the two groups.Pairmate corresponds to the
        proportion of sequences which has their nearest neighbor in a same group

        ARG    : - cut      : a rpy2IntVector with the repartition of the groups
                 - distance : An 2D array with distances between each sequence
        RETURN : - a tupple with the matepair of each group
    '''
    if verbose > 2:
        logging.info("compute pairmate")
    res = []
    n = len(cut)
    if(max(cut)==min(cut)):
        return 0,0
    distance[distance<=0.00000000001] = 999
    line = distance[0:n-1]
    line = line.tolist()
    line.insert(0,999)
    j = np.argmin(np.asarray(line)) ###
    res.append(str(cut[0])+str(cut[j]))
    line = []
    for i in range(1,n):
        line = []
        pos = i - 1
        line.append(distance[pos])
        k = 1
        while k < i:
            pos = pos + n - k -1
            line.append(distance[pos])
            k=k+1
        line.append(999)
        pos = pos +n-k
        if i != n-1 :
            line2 = distance[pos:pos+n-k-1]
            line = line + line2.tolist()
        j = np.argmin(np.asarray(line)) # douze
        res.append(str(cut[i])+str(cut[j]))
    mp1 = float(res.count('11'))/(res.count('11')+res.count('12'))
    mp2 = float(res.count('22'))/(res.count('21')+res.count('22'))
    return mp1,mp2

def save_results_sorted_by_cluster(filename, res, seq_name):
    ''' save the result in a file. Each sequence is associated to its group.
        the result is order by cluster. This function is faster than sorting by
        sequence.
        ARG    : - filename : the result filename
                 - res      : a dict with cluster as keys and the seq as value
                 - seq_name : a 1D numpy.array with the name of the sequences
        NO RETURN
    '''
    if verbose > 2:
        logging.info("save file")
    filename = open(filename,"w")
    for key in res.keys():
        for indice in res[key]:
            towrite = str(seq_name[indice])+"  "+str(key)+"\n"
            filename.write(towrite)

def save_results_sorted_by_sequence(filename, res, kmer_name):
    ''' save the result in a file. Each sequence is associated to its group.
        the result is order by sequence. Attention! This function is slower than
        sorting by cluster.
        ARG    : - filename : the result filename
                 - res      : a dict with cluster as keys and the seq as value
                 - seq_name : a 1D numpy.array with the name of the sequences
        NO RETURN
    '''
    if verbose > 2:
        logging.info("save file")
    filename = open(filename,"w")
    for i in range(len(kmer_name)):
        for key in res.keys():
            if i in res[key]:
                filename.write(kmer_name[i]+"   "+str(key)+"\n")

def get_args():
    ''' get the argument from the command line
        NO ARG
        RETURN : a tupple with the different args used in the programm
    '''
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description='%(prog)s can be used to classify multiple repeated/homologous sequences into families based on similariy comparison. \
    The similarity is measured from kmers content. \
    The classification uses a pairmate approach using coordinates of sequences into a PCA space.')
    parser.add_argument("input_file", help = "INPUT kmer/fasta file.")

    parser.add_argument("--sequence_minimal_length", help = "Filtered out sequence shorter than ...", type = int, default = None)
    parser.add_argument("--sequence_maximal_length", help = "Filtered out sequence longer than ...",  type = int, default = None)
    parser.add_argument("--kmer_length", help = "Set the length of the kmers.", type = int, default = None)
    parser.add_argument("--kmer_output_file", help = "Name of the output file storing kmer frequencies.\
                          Default filename is {input_file}.kmers{kmer_length}.", default = None)
    parser.add_argument("--kmer_relative_frequencies", 
						help = "Use relative frequencies of kmers.", 
						default = False, 
						action='store_true')
    parser.add_argument("--kmer_minimal_frequency", 
						help = "Set the minimal frequency of kmers used for the analysis.", 
						type = float, 
						default = 0.0)
                           
    parser.add_argument("--pairmate_threshold", 
						help = "Set mate pair threshold.", 
						type = float, 
						default = 0.95)
    parser.add_argument("--nbcomppca", 
						help = "Number of composante of pca used to compute distance.", 
						type=int,
                        default = 1024)
    parser.add_argument("--classification_output_file", 
						help = "Name of the output file storing the classification result. Default filename is {input_file}.classification",
                        default=None)
                        
    parser.add_argument('--sort_by_sequence_name', 
						dest='sort_by_sequence_name', 
						action='store_true', 
						default=True)
    parser.add_argument('--no_N_filtering', 
						dest='N_filtering', 
						action='store_false', 
						default=True)
    parser.add_argument('--sort_by_cluster_name', 
						dest='sort_by_sequence_name', 
						action='store_false',
						default=True)
    
    parser.add_argument("--verbose","-v", 
						help ="Degree of verbosity",
                        type=int, 
                        default=0)
    parser.add_argument("--timing","-t", help ="Print time of code",
                        type=int, default=0)
    parser.add_argument("--pca", help ="yes if you want to plot the \
                        pca else : no ", default="no")
    parser.add_argument("--lda_sample_size", help ="size of the sample used to train the lda ", type = int, default=75000)
    parser.add_argument("--min_family_size", help ="Minimal size of families, expressed as a proportion of the initial dataset.", type=float, default=0.0)
                        
    parser.add_argument("--version", help ="Print the version number", action='version', version='Version number: '+ str(version))
    args = parser.parse_args()
    if  args.min_family_size  > 1.0 or  args.min_family_size < 0.0:
        print("Error:  min_family_size  should be [0-1] (value: %f)" % (args.min_family_size))
        exit(0)        
    return(args.pairmate_threshold, args.input_file, \
    	   args.nbcomppca, \
           args.classification_output_file, \
           args.sort_by_sequence_name, args.pca, args.verbose, \
           args.lda_sample_size, \
           args.min_family_size, args.kmer_length, args.kmer_output_file, \
           args.kmer_relative_frequencies, args.kmer_minimal_frequency, \
           args.timing, args.sequence_minimal_length,  \
           args.sequence_maximal_length, args.N_filtering)


########################################################################
### Getting arguments
########################################################################

pairmate, input_file, nbcomp_pca, classification_output_file, \
sort_by_sequence_name, plot, verbose, lda_sample_size, maxfam, \
kmer_size, kmer_output_file, kmer_relative_frequencies, \
kmer_minimal_frequency, timing, \
sequence_minimal_length, sequence_maximal_length, \
N_filtering = get_args()

if timing > 0:
	t1=time.clock()

########################################################################
### Checking input/output files
########################################################################
if os.path.exists(input_file) ==  False:
	print("Error: input file %s doesn't exist" % (input_file))
	print("Bye bye :-(")
	quit(-1)

if kmer_size != None and kmer_output_file == None:
	kmer_output_file = "%s.kmers%d" % (input_file, kmer_size)
if kmer_size != None and os.path.exists(kmer_output_file) ==  True:
	os.remove(kmer_output_file)
	
if classification_output_file == None:
	classification_output_file = "%s.classifcation" % (input_file)
if os.path.exists(classification_output_file) ==  True:
	os.remove(classification_output_file)
	
########################################################################
### Computing kmer frequencies
### Needed if input file is fasta sequences
########################################################################
if timing > 0:
	t_km1=time.clock()
if kmer_size != None:
	kmer_file = compute_kmers_table_from_fasta_file(input_file, \
	kmer_size, kmer_output_file, kmer_relative_frequencies, \
	kmer_minimal_frequency, sequence_minimal_length, \
	sequence_maximal_length, N_filtering, verbose)
else:
	kmer_file = input_file
if timing > 0:
	print ("Time kmer calculation: %f s" %(time.clock()-t_km1))
	
########################################################################
### Reading kmer file
########################################################################
#size = len(kmer_indice)
## getting data
if timing > 0:
	t_kmr1=time.clock()
sequence_names = get_sequences_name(kmer_file, verbose)
kmer_table = get_kmer_table(kmer_file, verbose)
if timing > 0:
	print ("Time kmer reading: %f s" %(time.clock()-t_kmr1))

########################################################################
### Classification
########################################################################
if timing > 0:
	t_cl=time.clock()
results = classification(kmer_table, \
				lda_sample_size, \
				verbose)
if timing > 0:
	print ("Time classification: %f s" %(time.clock()-t_cl))

write_classification(results, \
				sequence_names, \
				classification_output_file, \
				sort_by_sequence_name)
if timing > 0:
	print("Time full: %f s" % (time.clock()-t1))
