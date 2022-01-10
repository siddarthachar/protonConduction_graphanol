# Proton Conduction Graphanol
Documentation and discussion regarding the study of proton conduction in graphanol using deep learning potentials.

## Codes
1. **z-axis-check.py**: Script that takes in XDATCAR file from a simulation and checks if all the atoms are within the accepted upper and lower z-planes. This is convenient for quickly checking the stability 2-D surfaces with dangling bonds that are prone to breaking. As for my case, I am interested in studying the proton conduction in graphanol, a lot of simulations are prone to bond breaking, which get hard to detect from mean-squared displacement calculations. 
2. **dump2xdatcar-sort.py**: Script to convert LAMMPS dump file to XDATCAR. It also sorts the atoms if LAMMPS jumbles the atom index order. 
3. **breakXDATCAR.py**: Script that takes in a XDATCAR file and then breakes it into smaller XDATCAR files of equal time spans. This is useful if you were to iteratively read individual images of a large XDATCAR file using ase.io. I will be computing something called the hopping rate using these broken XDATCAR files.
4. **hoppingRATE.py**: This code can be used to copmute the rate of proton hopping for an MD simulation. The user would need to have their files sorted in this way: {working directory}/{independent run}/XDATCAR. Each independent run should have a directory called "steps" that contains the XDATCAR separated into several parts. This process that be done by running the breakXDATCAR.py file. The code calculates a "hopping rate" by generating a connectivity matrix for each XDATCAR.x file. The number of changes in the connectivity matrix from iteration i and i-1 is accounted as the number of hops.  I have added a couple of output files: hopping.txt and hopping.png from a sample run of this code. 

