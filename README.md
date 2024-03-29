# Proton Conduction Graphanol
Documentation and discussion regarding the study of proton conduction in graphanol using deep learning potentials.

## Primary Files

1. **cec.py**: Script to compute the center of excess charge (CEC) to locate a proton in graphanol.
2. **lmc.py**: Script to perform lattice Monte Carlo (LMC) proton transfer simulation for graphanol. The code uses graph models constructed with Networkx. 
3. **goh.pb**: DeePMD potential (DP) for graphanol.
4. **goh.poscar**: Sample POSCAR file of 24C graphanol.

## Miscellaneous Files

1. **z-axis-check.py**: Script that takes in XDATCAR file from a simulation and checks if all the atoms are within the accepted upper and lower z-planes. This is convenient for quickly checking the stability 2-D surfaces with dangling bonds that are prone to breaking. As for my case, I am interested in studying the proton conduction in graphanol, a lot of simulations are prone to bond breaking, which get hard to detect from mean-squared displacement calculations. 
2. **dump2xdatcar-sort.py**: Script to convert LAMMPS dump file to XDATCAR. It also sorts the atoms if LAMMPS jumbles the atom index order. 
3. **breakXDATCAR.py**: Script that takes in a XDATCAR file and then breakes it into smaller XDATCAR files of equal time spans. This is useful if you were to iteratively read individual images of a large XDATCAR file using ase.io. I will be computing something called the hopping rate using these broken XDATCAR files.
4. **hoppingRATE.py**: 

