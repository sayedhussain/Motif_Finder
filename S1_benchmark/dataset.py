#!/usr/bin/env python3

# This porgram generates the dataset following the handout instructions 

import errno    
import math as m
import numpy as np
import os
import random
import scipy.optimize
import sys

##########
# GLOBAL #
##########

# Alphabet
ALPHABET = ['A', 'C', 'G', 'T']
# Information Contents Per Column Values
ICPC_VALUES = [2, 1, 1.5]
# Motif Length Values
ML_VALUES = [8, 6, 7] 
# Sequence Length
SL = 500
# Secuence Count Values
SC_VALUES = [10, 5, 20] 
# Number of experiments per combination
NUM_EXP = 10

# Filenames for the results
PARAMETERS_FILE = "parameters.txt"
SEQUENCES_FILE = "sequences.fa"
SITES_FILE = "sites.txt"
MOTIF_FILE = "motif.txt"
MOTIFLENGTH_FILE = "motiflength.txt"


############
# FUNTIONS #
############

def get_random_seq(sc, sl):
    """
    Produce random sequences

    Parameters
    ----------
    arg1: int
        Number of sequences
    arg2: int 
        Sequences length
    Returns
    -------
    list: 
       list of sequences
    """
    sequences = ['0']*sc
    # Produce every sequence
    for i in range(sc):
        seq = ['0']*sl
        # Produce every sequence element 
        for j in range(sl):
            seq[j] = (ALPHABET[random.randint(0, len(ALPHABET)-1)])
        sequences[i] = seq
    return sequences

def expression(c, a, b, icpc):
    """
    Expression to evaluate using scipy module

    Parameters
    ----------
    arg1: int
        Value of c (introduced by fsolve)
    arg2: int 
        Value of b
    arg3: int
        Value of c
    arg4: int
        Value of icpc
    Returns
    -------
    int: 
       results of the computation
    """
    a_val = a*m.log2(a/0.25)
    b_val = b*m.log2(b/0.25)
    c_val = c*m.log2(c/0.25)
    tot_val = (1-a-b-c)*m.log2((1-a-b-c)/0.25)
    return a_val + b_val + c_val + tot_val - icpc

def get_entry(icpc):
    """
    Computes one entry for the motif matrix

    Parameters
    ----------
    arg1: int
        Value of icpc
    Returns
    -------
    np.array: 
       row for the motif matrix 
    """
    # Case where icpc is 2
    if icpc == 2:
        return np.array([0,0,0,1])
    # General case
    i = 0
    while True:
        # Select a and b randomly
        a = random.random()
        b = random.random()*(1-a)
        # Check if there is a solution
        value = icpc-2-a*m.log2(a)-b*m.log2(b)
        sub = (1-a-b)
        if value<sub*m.log2(sub/2.) or value>sub*m.log2(sub):
            continue
        # Solve
        try:
            c = scipy.optimize.fsolve(expression, (1-a-b)/2., args=(a, b, icpc), factor=1)
            return np.array([a,b,c,1-a-b-c])
        except Exception as e: 
            continue

def get_random_motif(ml, icpc):
    """
    Generates the motif matrix

    Parameters
    ----------
    arg1: int
        Motif lenght
    arg2: int
        Value of icpc
    Returns
    -------
    np.array: 
        motif matrix 
    """
    random_motif = np.zeros([ml, len(ALPHABET)])
    # Generate every row
    for i in range(ml):
        # Get an entry and randomly permute it
        entry = get_entry(icpc)
        random_motif[i, :] = np.random.permutation(entry)
    return random_motif

def get_motifs(random_motif, sc, ml):
    """
    Generate random motifs using the motif matrix

    Parameters
    ----------
    arg1: np.array
        Motif matrix
    arg2: int
        Secuence length
    arg3: int
        Motig length
    Returns
    -------
    list: 
       list of random generated motifs
    """
    motifs = ['0']*sc
    # For each motif
    for i in range(sc):
        motif = ['0']*ml
        # generate it randomly using the matrix
        for j in range(ml):
            motif[j] = np.random.choice(ALPHABET, p=random_motif[j])
        motifs[i] = motif
    return motifs

def plant(sequences, motifs, sl, sc, ml):
    """
    Plant randomly the motif in the random sequence

    Parameters
    ----------
    arg1: list
        List of motifs
    arg2: list
        List of sequences
    arg3: int
        Sequence length
    arg4: int
        Sequence count
    arg5: int
        Motif length
    Returns
    -------
    list: 
       sequences with motif included
    list: 
       locations of motifs
    """
    locations = [0] * sc
    for i in range(sc):
        loc = random.randint(0, sl - ml - 1)
        locations[i] = loc
        sequences[i][loc:loc+ml] = motifs[i]
    return sequences, locations

def mkdir_p(path):
    """
    Creates a directory

    Parameters
    ----------
    arg1: string
        directory path
    """
    # Create dir
    try:
        os.makedirs(path)
    # Check error or already exists
    except OSError as exc:
        # Exists
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        # Error
        else:
            raise

########
# MAIN #
########
dataset_counter = 1
for i in range(3*len(ICPC_VALUES)-2):
    # Experiments changing ICPC
    if i < len(ICPC_VALUES):
        ICPC = ICPC_VALUES[i]
        ML = ML_VALUES[0]
        SC = SC_VALUES[0]
    # Experiments changing ML
    elif i < 2*len(ICPC_VALUES)-1:
        ICPC = ICPC_VALUES[0]
        ML = ML_VALUES[i%len(ICPC_VALUES)+1]
        SC = SC_VALUES[0]
    # Experiments changing SC
    else:
        ICPC = ICPC_VALUES[0]
        ML = ML_VALUES[0]
        SC = SC_VALUES[i%(len(ICPC_VALUES)+1)]
    # Get and create the directory
    path = "./dataset"+str(dataset_counter)+"/"
    mkdir_p(path)
    # Compute 10 experiments for each combination
    for nexp in range(NUM_EXP):
        # Get random sequences
        random_seq = get_random_seq(SC, SL)
        #print("Random Sequences: ", len(random_seq))

        # Get random motif
        random_motif = get_random_motif(ML, ICPC)
        # Generate motifs
        motifs = get_motifs(random_motif, SC, ML)

        # Plant one sampled site at a random location
        sequences, locations = plant(random_seq, motifs, SL, SC, ML)

        # Get and create the directory
        path = "./dataset"+str(dataset_counter)+"/"+"dataset"+str(nexp+1)+"/"
        mkdir_p(path)
        # Write parameters
        param_file = open(path+PARAMETERS_FILE, 'w')
        param_file.write("ICPC: " + str(ICPC) + "\n")
        param_file.write("ML: " + str(ML) + "\n")
        param_file.write("SL: " + str(SL) + "\n")
        param_file.write("SC: " + str(SC) + "\n")
        param_file.close()

        # Write sequences
        seq_file = open(path+SEQUENCES_FILE, 'w')
        for i in range(SC):
            seq_file.write(">seq" + str(i) + "\n")
            seq_file.write("".join(sequences[i]) + "\n")
        seq_file.close()

        # Write sites
        sites_file = open(path+SITES_FILE, 'w')
        for i in range(len(locations)):
            sites_file.write(str(locations[i]) + "\n")
        seq_file.close()

        # Write motif
        motif_file = open(path+MOTIF_FILE, 'w')
        motif_file.write(">MOTIF"+ str(nexp) + " " + str(ML) + '\n')
        for i in range(len(random_motif)):
            motif_file.write(" ".join(str(j) for j in random_motif[i]) + '\n')
        motif_file.write("<\n")
        motif_file.close()

        # Write motif length
        motiflength_file = open(path+MOTIFLENGTH_FILE, 'w')
        motiflength_file.write(str(ML))
        motiflength_file.close()
    # Increment counter
    dataset_counter += 1