import numpy as np
import ase 
from ase.io import read, write
import matplotlib.pyplot as plt

def line_prepender(filename, line):
	'''
	This function is used for prepending strings to files.
	'''
    with open(filename, 'r+') as f:
        content = f.read()
        f.seek(0, 0)
        f.write(line.rstrip('\r\n') + '\n' + content)


# XDATCAR file from simulations
path = lambda x: f'/bgfs/kjohnson/ska31/1DeePMD/paper-graphanol/oneside/p+transfer/TWO_protons/multiple-24C_800K/{x}/XDATCAR'

# Trimmed XDATCAR file: containing only 20 images from the last 10000 time steps.
outpath = lambda x: f'/bgfs/kjohnson/ska31/1DeePMD/paper-graphanol/oneside/p+transfer/TWO_protons/multiple-24C_800K/{x}/XDATCAR_select'

TOTsteps = 80000  #change accordingly

extract=np.linspace(TOTsteps-10000, TOTsteps, 20)
extract=extract.astype(int)
print(extract)
for k in range(1,21):
	start_collect = False
	file_inp = open(path(k),'r')
	file_out = open(outpath(k),'w+')
	lines = file_inp.readlines()
	noAtoms = 62  #change accordingly
	count = 0
	counter_direct = 0
	for l in lines:
	    if start_collect == True:
	        file_out.write(str(l))
	        count += 1
	    if count==noAtoms:
	        start_collect = False
	        count = 0
	    sp = l.split()
	    if sp[0] == 'Direct':
	        if (int(sp[2]) == extract).any() == True :  # uncomment this once test case is done
	            counter_direct += 1
	            file_out.write(f'Direct configuration=\t{counter_direct}\n')
	            start_collect = True
	file_out.close()

	# Need to add this string to the top of the file so that it becomes an XDATCAR file.
	prependThis =' C  O  H\n\
	 1.0000000000000000\n\
	    10.2929999999999993    0.0000000000000000    0.0000000000000000\n\
	    -3.9039999999999981    6.7619263527488975    0.0000000000000000\n\
	     0.0000000000000000    0.0000000000000000   20.0000000000000000\n\
	 C   O   H\n\
	  24  12  26'
	print(k)
	line_prepender(outpath(k),prependThis)

	plane_z = [7,11.8]  # considering the z-axis plane, still need to check if the upper plane should be beyond 11.8

	# All images that have floating atoms (atoms above or below plane_z) will be saved as an xyz file called selected.xyz
	xyz_displaced = lambda x: f'/bgfs/kjohnson/ska31/1DeePMD/paper-graphanol/oneside/p+transfer/TWO_protons/multiple-24C_800K/{x}/selected.xyz'
	write_it = open(xyz_displaced(k),'w+')
	for i in range(20):
	    o=read(outpath(k),index=i)
	    pos=o.get_positions()
	#     print(i)
	    below_plane_pos = pos[pos[:,2]<plane_z[0]]
	    below_planeIND = (np.argwhere(pos[:,2]<plane_z[0])).reshape(1,-1)[0]
	    symb=o.get_chemical_symbols()
	    below_plane_sym = np.array(symb)[below_planeIND]
	    if len(below_plane_pos) > 0:
	    	print("YES")

	    above_planeIND = (np.argwhere(pos[:,2]>plane_z[1])).reshape(1,-1)[0]
	    above_plane_pos = pos[pos[:,2]>plane_z[1]]
	    above_plane_sym = np.array(symb)[above_planeIND]
	    write_it.write(f'\t {len(above_plane_sym)+len(below_plane_sym)}\n\n')
	    if len(above_plane_pos) > 0:
	    	print("YES")
	    for p,s in zip(below_plane_pos,below_plane_sym):
	        write_it.write(f'{s}\t{p[0]:1.3f}  {p[1]:1.3f}  {p[2]:1.3f}\n')
	        
	    for p,s in zip(above_plane_pos,above_plane_sym):
	        write_it.write(f'{s}\t{p[0]:1.3f}  {p[1]:1.3f}  {p[2]:1.3f}\n')
	    write_it.write('\n')
