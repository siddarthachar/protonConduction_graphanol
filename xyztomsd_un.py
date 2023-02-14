#!/usr/bin/env python

import numpy as np
import argparse
import os
#from pylab import *

# User specificed value
# Length of a timestep in fs
timestep = 0.25

parser = argparse.ArgumentParser()
# Option to read file from command line:
parser.add_argument("-i", "--input", help="Name of xyz file to be read", default="centerofmass.xyz")
args = parser.parse_args()
print("Opening file {0}, in {1}".format(args.input,os.getcwd()))
# load file, using numpy libraries 
filename = args.input
f=open(filename)
lines = f.readlines()
# Get the total number of atoms from the first line of the file
total_atoms = int(lines[0])
print('Total number of atoms = ', total_atoms)
# offset to account for the number of lines to skip in the xyz file before reading atom positions
offset = 2
# Specify the atom lables of the atoms to include in the MSD calculations:
#atoms_to_include = [48,49,50,51,52]
atoms_to_include = [0]
print(f"Include {len(atoms_to_include)} atoms : ",atoms_to_include)
# counting in Python is zero-indexed, so shift everything
atoms_to_include = [x for x in atoms_to_include]
# counting off the lines, select only the atoms we care about
#  throwing away all other data
lines = [line for k, line in enumerate(lines)
         if ((k - offset) % (total_atoms + offset)) in atoms_to_include]
# splits the lines into fields, throwing away the first column
lines = [line.split()[1:] for line in lines]
# convert the data into numeric data from the string form
lines = [[float(q) for q in line] for line in lines]
# finally convert into the Numpy array object type
raw = np.array(lines)

# reshapes from a blob to indexable by
#  a[snapshot][atom][xyz]
nmol = len(atoms_to_include)
nsnapshots = int(len(raw)/nmol)
print(nsnapshots)
a = raw.reshape(nsnapshots, nmol, 3)

# calculate important bounds
#  and other odd-and-ends for counting

n_cmsds = 100    # only concurrent MSDs
#n_cmsds = 2    # only concurrent MSDs
half_width = int(nsnapshots / 2)
dtspan = half_width / n_cmsds
print(dtspan, half_width, n_cmsds)
span_starts = [x for x in range(nsnapshots)
    if x % dtspan == 0 and x + half_width <= nsnapshots]
n_msds = len(span_starts)
#print n_msds

# do calculations, using conveniences provided by Numpy
#  specifically the ability to find the difference of two
#  snapshots by simply a[k] - a[k0], with a double loop
#  to loop over time origins then snapshots in the window
# the list comprehension as written is just a blob, so
#  reshape it to something indexable, specifically
#  msd[time_origin][trajectory_time][xyz]
# then free some memory, just to play nice

msd = np.array([np.mean((a[k] - a[k0]) ** 2, 0)
    for k0 in span_starts  for k in range(k0, k0 + half_width)])
msd = msd.reshape(n_msds, half_width, 3).mean(0)

#dcf = np.array([np.mean(a[k] - a[k0], 0) ** 2
#    for k0 in span_starts
#    for k in range(k0, k0 + half_width)])
#dcf = dcf.reshape(n_msds, half_width, 3).mean(0) * nmol

del raw, a

# Compute the time in fs and put into a vector
t = np.arange(0,len(msd),timestep)
#print msd
filename2 = filename.rstrip('.xyz') + '_msd.txt'
print('MSD written to', filename2, '!')

msdtxyz = np.array([np.hstack(k) for k in zip(t,msd)])
np.savetxt(filename2, msdtxyz)

