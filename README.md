# Proton Conduction Graphanol
Documentation and discussion regarding the study of proton conduction in graphanol using deep learning potentials.

## Codes
1. **z-axis-check.py**: Script that takes in XDATCAR file from a simulation and checks if all the atoms are within the accepted upper and lower z-planes. This is convenient for quickly checking the stability 2-D surfaces with dangling bonds that are prone to breaking. As for my case, I am interested in studying the proton conduction in graphanol, a lot of simulations are prone to bond breaking, which get hard to detect from mean-squared displacement calculations. 
