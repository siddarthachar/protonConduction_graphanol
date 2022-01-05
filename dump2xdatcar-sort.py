## Copyright: Siddarth Achar
## University of Pittsburgh
## ska31@pitt.edu


# This script reads LAMMPS dump files and converts them to XDATCAR files. 
# It also sorts the atoms based on the atom index if LAMMPS messes up the 
# order of printing due to parallelization etc.
# NOTE: At the end of the script the user would have to manually change the heading 
# of the XDATCAR file. 

import numpy as np
import sys

lammpsdump=str(sys.argv[1])   # LAMMPS dump file
xdatcar_out=str(sys.argv[2])  # output XDATCAR file
noAtoms=int(sys.argv[3])    # Integer number of atoms in the system

def line_prepender(filename, line):
    with open(filename, 'r+') as f:
        content = f.read()
        f.seek(0, 0)
        f.write(line.rstrip('\r\n') + '\n' + content)


dump = open('prod.lammpstrj')  # Reading lammps dump file
out = open('XDATCAR','w+')     # Writing to XDATCAR file
count = 0
direct_count = 0
start = 0
tag = 0
frame = []
atm = []
for l in dump:
    if l.split()[0] == 'ITEM:':
        if l.split()[1] == 'TIMESTEP':
            count = 0
    
    if tag == 1:
        atm.append([int(l.split()[0]),float(l.split()[2]),float(l.split()[3]),float(l.split()[4])])
        count += 1
    
    if count == noAtoms:
        tag = 0
        frame.append(atm)
        atm = []
    
    
    if l.split()[0] == 'ITEM:':
        if l.split()[1] == 'ATOMS':
            direct_count += 1
            tag = 1
            count = 0
frame = np.array(frame)
print("Writing into XDATCAR...")
for i in range (len(frame)):
    fr_i = frame[i]
    fr_is = fr_i[fr_i[:,0].argsort()]
    
    out.write(f"Direct configuration=\t{i+1}\n")
    
    for j in range(len(fr_is)):
        out.write('   '+str(fr_is[j,1])+'  '+str(fr_is[j,2])+'  '+str(fr_is[j,3])+'\n')
    
    if (i+1)%100 == 0:
        print(f'{i+1} frames written.')

# Change the header of the XDATCAR file. 

prependThis =' C  O  H\n\
 1.0000000000000000\n\
    10.2929999999999993    0.0000000000000000    0.0000000000000000\n\
    -3.9039999999999981    6.7619263527488975    0.0000000000000000\n\
     0.0000000000000000    0.0000000000000000   20.0000000000000000\n\
 C   O   H\n\
  24  12  26'

line_prepender('XDATCAR',prependThis)
