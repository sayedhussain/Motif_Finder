#!/usr/bin/env python3

# This program executes the algorithm for 1 dataset 

from difflib import SequenceMatcher
import numpy as np
from Bio.Seq import Seq
import math
from Bio import motifs
import os
import time
import sys



##########
# GLOBAL #
##########

############
# FUNTIONS #
############
def genrate_icpc(motif, column)-> float:
    """ This method is generate Information content
    per column given motif(position weight matrix) and column
    :param: motif matrix
    :param: column
    :return: icpc
    """
    icpc = 0
    for row in range(len(motif)):
        a=motif[row][column]
        if a == 0:
            continue
        icpc += a * math.log(a/0.25,2)
    return icpc



def generate_motifs(possible_motifs):
    """ This function will generate the list of all possible list and returns the array
    :param: list
        possible_motifs
    :return: list
        possible_motifs_list"""
    global motifs_index
    for kx,i in enumerate(possible_motifs[:]): #looping through all the possible sequence in possible motif list
        separate_values = []
        separate_values.append(Seq(i))
        possible_motifs.append(separate_values)
        first_motif = [kx] # temporary list to store index value of each sequence
        for ix,j in enumerate(list_of_sequence): #looping to all the other rows in main input list(list_of_sequence)
            if ix == 0:
                continue
            #generate all the possible combinations for that row
            possible_seq = [list_of_sequence[ix][s:s + var] for s in range(0, len(list_of_sequence[0]) - var)]
            possible_seq = np.asarray(possible_seq) #converting it to numpyarray
            list_of_sim = [SequenceMatcher(None, i, k).ratio()for k in possible_seq] #calculating similarity between 2 strings
            list_of_sim = np.asarray(list_of_sim, dtype=np.float64)
            max_value = np.max(list_of_sim) #getting the max similiraity score
            index_value = np.argmax(list_of_sim) # getting its max similiraity score value index
            first_motif.append(index_value) #appending index value into the list to later use to create predicted sites
            separate_values.append(possible_seq[index_value])# appending the most similar sequence in list
        motifs_index.append(first_motif)# appending the motif index list for genrated for each possible sequence
        m = motifs.create(separate_values)#create motif object using biopython library
        possible_motifs_list.append(m) # appending motif to the main list

    return possible_motifs_list

########
# MAIN #
########

# Check input
if len(sys.argv) != 2:
    print ("Usage: python algorithm.py dataset_path")

path = sys.argv[1]

# creating a list to store input sequence fasta file
list_of_sequence = []
with open(path+'sequences.fa', 'r') as fobj:
    for i in fobj.readlines():
        if i.strip()[0] == '>':  # skipping the line which contains '>'
            continue
        list_of_sequence.append(i)

# reading the lenght of the motif into the var variable
with open(path+'motiflength.txt', 'r') as fobj:
    for i in fobj.readlines():
        var = int(i.strip())

# start time
start = time.time()

# generating all the possible var lenght combinations of sequence from input sequence first entry
possible_motifs = [list_of_sequence[0][i:i + var] for i in range(0, len(list_of_sequence[0]) - var)]

possible_motifs_list = []  # final list of motifs
motifs_index = []  # will be used to store the index of sequence

possible_motifs_list = generate_motifs(possible_motifs)

information_content_list = []  # this list will store total IC for the motif
# generating total IC for the motifs and then storing them in a list
for i in possible_motifs_list:
    pwm = i.counts.normalize()
    total_icpc = 0
    for col_num in range(0, var):
        total_icpc += genrate_icpc(pwm, col_num)
    information_content_list.append(total_icpc)

information_content_list = np.asarray(information_content_list, dtype=np.float64)
max_icpc = np.max(information_content_list)  # fetching maximum total IC value
max_motif_index = np.argmax(information_content_list)  # fetching index of max motif

print_list = []  # adding motif matrix value into the list to write in correct format
for i in ['A', 'C', 'G', 'T']:
    print_list.append(possible_motifs_list[max_motif_index].counts.normalize()[i])

mat = np.matrix(print_list)
mat = mat.transpose()

end = time.time()

# writing the motif matrix
with open(path+'predictedmotif.txt', 'w+') as f:
    f.write('>MOTIF0 ' + str(var) + '\n')
    for line in mat:
        np.savetxt(f, line, fmt='%.1f')
    f.write('<')
f.close()

# writing index values
with open(path+'predictedsites.txt', 'w') as f:
    for i in motifs_index[max_motif_index]:
        f.write(str(i) + '\n')
f.close()

# writing time
with open(path+'runningtime.txt', 'w+') as f:
    f.write(str(round(end-start, 3)) + '\n')
f.close()

















