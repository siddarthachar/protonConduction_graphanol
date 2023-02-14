import numpy as np
import ase 
from ase.io import read, write
from ase import Atoms
import matplotlib.pyplot as plt
import networkx as nx
import random
from tqdm import tqdm
from tqdm import trange
from IPython import display
import time
import networkx as nx
import os, sys
from numba import jit
class HiddenPrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')
    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout
with HiddenPrints():
    print("This wont print")

delE_rot = float(sys.argv[1])
delE_hop = float(sys.argv[2])
T = float(sys.argv[3]) # temperature
# delE_rot = 25 # Fill in kJ/mol
# delE_hop = 9 # Fill in kJ/mol

hex_loops = [
    [11,12,25,24,23,10],
    [9, 10, 23, 22, 21, 8],
    [7, 8, 21, 20, 19, 6],
    [5, 6, 19, 18, 17, 4],
    [3, 4, 17, 16, 15, 2],
    [1, 2, 15, 14, 13, 0],
    [25, 26, 40, 39, 38, 24],
    [23, 24, 38, 37, 36, 22],
    [21, 22, 36, 35, 34, 20],
    [19, 20, 34, 33, 32, 18],
    [17, 18, 32, 31, 30, 16],
    [15, 16, 30, 29, 28, 14],
    [38, 39, 53, 52, 51, 37],
    [36, 37, 51, 50, 49, 35],
    [34, 35, 49, 48, 47, 33],
    [32, 33, 47, 46, 45, 31],
    [30, 31, 45, 44, 43, 29],
    [28, 29, 43, 42, 41, 27],
    [53, 54, 68, 67, 66, 52],
    [51, 52, 66, 65, 64, 50],
    [49, 50, 64, 63, 62, 48],
    [47, 48, 62, 61, 60, 46],
    [45, 46, 60, 59, 58, 44],
    [43, 44, 58, 57, 56, 42],
    [66, 67, 81, 80, 79, 65],
    [64, 65, 79, 78, 77, 63],
    [62, 63, 77, 76, 75, 61],
    [60, 61, 75, 74, 73, 59],
    [58, 59, 73, 72, 71, 57],
    [56, 57, 71, 70, 69, 55],
    [81, 82, 95, 94, 93, 80],
    [79, 80, 93, 92, 91, 78],
    [77, 78, 91, 90, 89, 76],
    [75, 76, 89, 88, 87, 74],
    [73, 74, 87, 86, 85, 72],
    [71, 72, 85, 84, 83, 70]
]

# finding a set of O centers maximized
oxygen_hex =  {'O0': [77, 78, 91, 90, 89, 76],
 'O1': [60, 61, 75, 74, 73, 59],
 'O2': [21, 22, 36, 35, 34, 20],
 'O3': [43, 44, 58, 57, 56, 42],
 'O4': [66, 67, 81, 80, 79, 65],
 'O5': [15, 16, 30, 29, 28, 14],
 'O6': [32, 33, 47, 46, 45, 31],
 'O7': [5, 6, 19, 18, 17, 4],
 'O8': [49, 50, 64, 63, 62, 48],
 'O9': [38, 39, 53, 52, 51, 37],
 'O10': [71, 72, 85, 84, 83, 70],
 'O11': [11, 12, 25, 24, 23, 10]}


def initializedHex():
    # Initialize graph
    hexlat = nx.hexagonal_lattice_graph(m=6,n=6, periodic=False) # m=2 and n=4 is the 24C strucutre
    # relabeling these nodes. There is no reason to have the sencond dimension
    mapping = {}
    H_idx = 0
    for node in hexlat.nodes:
        mapping[node]=H_idx
        H_idx+=1
    H = nx.relabel_nodes(hexlat, mapping)

    # Mapping all the hexagons in the graph

    # removing redundant nodes, which do not form hexagons.
    nodes_as_hex = set(np.unique(list(oxygen_hex.values())))
    all_nodes = set(H.nodes())
    diff_nodes = list(all_nodes.difference(nodes_as_hex))
    H_remove = H.copy()
    H_remove.remove_nodes_from(diff_nodes)
    # pos = nx.get_node_attributes(H_remove, 'pos')
    # nx.draw(H_remove, with_labels=True, pos=pos, node_color='red')
    # plt.draw()
    # plt.savefig('graph.png', transparent=True, node_color='red')

    # adding periodic edges for the corner nodes
    H_remove.add_edges_from([(28,12),
                             (25,42),
                             (56,39),
                             (53,70),
                             (83,67),
                             (81,14),
                             (84,15),
                             (85,4),
                             (74,5),
                             (89,6), 
                             (90,21),
                             (91,10),
                             (80,11)])

    pos = nx.get_node_attributes(H_remove, 'pos')
    # nx.draw(H_remove, with_labels=True, pos=pos, node_color='red')
    # plt.draw()

    # maintaining a separate copy of H_remove
    H_mc = H_remove.copy()

    # setting all nodes to have H occupied as False
    for n in H_mc.nodes:
        H_mc.nodes[n]['H occupied'] = False

    # go through each O atom and place a H atom in them with the rule that it is not occupied.
    start = 0
    for keys, values in oxygen_hex.items(): # iterate over all oxygen atoms
        if start == 0:
            H_place_attempt = values[0] # just pick the first one
            H_mc.nodes[H_place_attempt]['H occupied'] = True
            start += 1
        else:
            for H_place_attempt in values:
                neigh = [ne for ne in H_mc.neighbors(H_place_attempt)] # finding the neighbors of the attempted placement
                if sum([H_mc.nodes[x]['H occupied'] == False for x in neigh]) == 3: # if all neighbors are not occupied
                    H_mc.nodes[H_place_attempt]['H occupied'] = True
                    break

    # Adding a proton to the system, we can just place the proton on the first O atom.

    # finding the H atom in the O atom
    hinO=[hatom for hatom in oxygen_hex['O0'] if H_mc.nodes.data()[hatom]['H occupied']==True][0]
    # finding all the shortest path lengths for each site in O0 to hinO and store the ones that are over 1 distance
    poten_proton_site = [x for x in oxygen_hex['O0'] if nx.shortest_path_length(H_mc, hinO,x) > 1]
    # next find the 
    for pps in poten_proton_site:
        neigh = [ne for ne in H_mc.neighbors(pps)] # finding the neighbors of the attempted placement
        if sum([H_mc.nodes[x]['H occupied'] == False for x in neigh]) == 3: # if all neighbors are not occupied
            H_mc.nodes[pps]['H occupied'] = True
            break
            
    return H_mc

H_mc = initializedHex()

# Locating all H atom sites 
def locate_H_sites(G):
    "Returns the nodes where there are H atoms occupied"
    return [m[0] for m in G.nodes.data() if m[1]['H occupied']==True]
def O_H_conn_dict(G):
    "Returns a dictionary of the format {'O_':[H_nodes]} for each O atom in the system"
    H_site = locate_H_sites(G)
    set_H_site = set(locate_H_sites(G))
    O_dict = {}
    for keys, values in oxygen_hex.items():
        O_dict[keys] = list(set(values).intersection(set_H_site))
    return O_dict, H_site
# locate_H_sites(H_mc) # example 
# O_H_conn_dict(H_mc) # example
def LMC(H_mc, delE_rot, delE_hop, T, CEC_track=[], move_tracking=[0,0,0,0,0,0,[],[]], displaying = False, updating = False, moves = 1000000, checking_neigh = True):
    rot_boltzmann = np.exp(-delE_rot*1000/8.314/T) # 1000 for kJ to J
    hop_boltzmann = np.exp(-delE_hop*1000/8.314/T)
    
    temp = Atoms('O',positions=[(0,0,0)], pbc=True)
    temp.set_cell([7.6755306, 8.9851905, 20])
    
    rot_accepted, rot_attempted, rot_rejected, phop_accepted, phop_attempted, phop_rejected,total_H_atoms, H_atoms_traj  = move_tracking
    for move in trange(moves):    
        O_H_connect, H_site_list = O_H_conn_dict(H_mc)

        # Checking if at any time there are two occupied H sites
        if checking_neigh:
            for node in H_mc.nodes:
                if node in H_site_list:
                    if len([x for x in H_mc.neighbors(node) if H_mc.nodes[x]['H occupied']==True]) > 0:
                        print(H_site_list)
                        print('Error')

        total_H_atoms += [len(H_site_list)]
        H_atoms_traj += [H_site_list]
        
        # Pick a random O atom
        O_pick = random.choice(list(oxygen_hex.keys())) # randomly pick an O atom

        # Rotation attempt
        if len(O_H_connect[O_pick]) == 1: 
            # If picked O atom has only 1 H neighbor, then rotate
            H_pick = O_H_connect[O_pick][0]
            rot_neigh = set([ne for ne in H_mc.neighbors(H_pick)]).intersection(set(oxygen_hex[O_pick]))
            poten_rot_sites = [x for x in rot_neigh if H_mc.nodes[x]['H occupied'] == False] # find potential sites to rotate to
            rot_attempted += 1
            if len(poten_rot_sites) == 0:
                rot_rejected += 1
            elif len(poten_rot_sites) == 2: # when both moves are possible
                rot_pick = random.choice(poten_rot_sites) # pick one site at random
                # see if the picked H site does not suffer from any other neighbors around it
                rot_neigh_neigh = set([ne for ne in H_mc.neighbors(rot_pick)]).difference(set([H_pick]))
                if sum([H_mc.nodes[x]['H occupied'] == True for x in rot_neigh_neigh]) > 0:
                    rot_rejected += 1
                else:
                    zeta = random.uniform(0,1)
                    boltzman_f = rot_boltzmann # change this to exp(-deltaE/kT)
                    if zeta < boltzman_f:
                        H_mc.nodes[H_pick]['H occupied'] = False
                        H_mc.nodes[rot_pick]['H occupied'] = True
                        rot_accepted += 1

                    else:
                        rot_rejected += 1
            else: # when only one move is possible
                rot_pick = poten_rot_sites[0]
                rot_neigh_neigh = set([ne for ne in H_mc.neighbors(rot_pick)]).difference(set([H_pick]))
                if sum([H_mc.nodes[x]['H occupied'] == True for x in rot_neigh_neigh]) > 0:
                    rot_rejected += 1
                else:
                    zeta = random.uniform(0,1)
                    boltzman_f = rot_boltzmann # change this to exp(-deltaE/kT)
                    if zeta < boltzman_f:
                        H_mc.nodes[H_pick]['H occupied'] = False
                        H_mc.nodes[rot_pick]['H occupied'] = True
                        rot_accepted += 1
                    else:
                        rot_rejected += 1

        # Proton hop attempt
        elif len(O_H_connect[O_pick]) == 2:
            H_pick = random.choice(O_H_connect[O_pick])
            proton_hop_neigh = list(set([ne for ne in H_mc.neighbors(H_pick)]).difference(set(oxygen_hex[O_pick])))[0] # where to potentially move the H neighbor
            proton_hop_neighs_neigh = list(set([ne for ne in H_mc.neighbors(proton_hop_neigh)]).difference(set([H_pick])))
            phop_attempted += 1
            if H_mc.nodes[proton_hop_neigh]['H occupied'] == True: # if it is occupied, then reject
                phop_rejected += 1
            elif sum([H_mc.nodes[x]['H occupied'] == True for x in proton_hop_neighs_neigh]) > 0: # if either of the neighbors are not empty 
                phop_rejected += 1
            else:
                zeta = random.uniform(0,1)
                boltzman_f = hop_boltzmann # ch
                if zeta < boltzman_f:
                    H_mc.nodes[H_pick]['H occupied'] = False
                    H_mc.nodes[proton_hop_neigh]['H occupied'] = True
                    phop_accepted += 1
                else:
                    phop_rejected += 1
                    
        O_H_connect, H_site_list = O_H_conn_dict(H_mc)
        CEC_track_cart = oxygen_coords[[x for x in O_H_connect.keys() if len(O_H_connect[x]) == 2][0]]
        temp.set_positions([CEC_track_cart])
        CEC_track.append(temp.get_scaled_positions())
        
        if displaying and move%1 == 0:
            options = {"edgecolors": "tab:gray", "node_size": 600, "alpha": 0.9}
            Empty_nodes = [x for x in H_mc.nodes if x not in H_site_list]
            nx.draw_networkx_nodes(H_mc, pos, nodelist=H_site_list, node_color="tab:red", **options)
            nx.draw_networkx_nodes(H_mc, pos, nodelist=Empty_nodes, node_color="tab:blue", **options)
            nx.draw_networkx_labels(H_mc, pos);
            if updating:
                display.clear_output(wait=True)
                display.display(plt.gcf())
            time.sleep(1.0)
            
            
    move_tracking = [rot_accepted, rot_attempted, rot_rejected,
                 phop_accepted, phop_attempted, phop_rejected,
                total_H_atoms, H_atoms_traj]

    CEC_track = np.array(CEC_track)
    
    # Unwrapping the CEC
    
    Natoms = 1 # 12 O atoms
    aa = np.array([7.5, 0.0, 0.0])
    bb = np.array([0.0, 8.96, 0.0])
    cc = np.array([0.0, 0.0, 20.0])
    unitcell = np.array([aa, bb, cc])
    element_names = 'O'
    scale=1

    prev_step = np.zeros((Natoms,3))
    curr_step = np.zeros((Natoms,3))
    unwrapped = np.zeros((Natoms,3))
    unwrapped_list = []
    
    for ii in range(len(CEC_track)):
        if ii == 0:
            prev_step[0] = CEC_track[0][0]

    # Read in the first configuration to start with:
    # Start the unwrapped coordinates at the origin
            unwrapped = np.copy(prev_step)
            xyzcoords = np.multiply(np.dot(prev_step,unitcell),scale)
            unwrapped_list.append(xyzcoords)

    # Write out the first step in Cartesian coordinates:
    # Multiply the vectors
    # Read in the next configuration and compute the distance

        for i in range(Natoms):
            coords = CEC_track[ii]
            for j in range(3):
                curr_step[i][j] = float(coords[0][j])
                dd = prev_step[i][j] - curr_step[i][j]
                if np.fabs(dd) > 0.5:
    # The distance traveled was larger than 1/2 the box length, so unwrap the coordinate
                    unwrapped[i][j] = unwrapped[i][j] + curr_step[i][j] + np.sign(dd)*1.0 - prev_step[i][j]
                else:
                    unwrapped[i][j] = unwrapped[i][j]+curr_step[i][j] - prev_step[i][j]

    # Set the previous step to the current step:
        prev_step = np.copy(curr_step)
        # Multiply the vectors
        xyzcoords = np.multiply(np.dot(unwrapped,unitcell),scale)
        unwrapped_list.append(xyzcoords)

    unwrapped_list = np.array(unwrapped_list)
    
    
    return H_mc, move_tracking, unwrapped_list



# Initialize Metropolis MC settings
# Now that we have our network populated, let's try to run some Monte Carlo!
rot_accepted = 0
rot_attempted = 0
rot_rejected = 0
phop_accepted = 0
phop_attempted = 0
phop_rejected = 0


# plotting settings
displaying = False
updating = False

# Checing neighbors bool
checking_neigh = True

CEC_now = 'O0' # just for now
total_H_atoms = []
H_atoms_traj = []
move_tracking = [rot_accepted, rot_attempted, rot_rejected,
                 phop_accepted, phop_attempted, phop_rejected,
                total_H_atoms, H_atoms_traj]

oxygen_coords = {'O11': [0.567, 7.576, 10.877],
'O9': [3.176, 7.544, 10.870],
'O4': [5.686, 7.577, 10.897],
'O2': [1.828, 5.362, 10.885],
'O8': [4.391, 5.292, 10.871],
'O0': [6.972, 5.374, 10.915],
'O7': [0.507, 3.079, 10.917],
'O6': [3.053, 3.105, 10.891],
'O1': [5.664, 3.078, 10.915],
'O5': [1.864, 0.808, 10.892],
'O3': [4.410, 0.871, 10.893],
'O10': [7.030, 0.808, 10.916],
}


H_mc = initializedHex()
H_mc, move_tracking, CEC_track = LMC(H_mc, delE_rot, delE_hop, T, moves=1000000)

nx.write_adjlist(H_mc,'H_mc.adjlist')

filename_title = 'cec' # 
nmol = 1
nsnapshots = int(len(CEC_track)/nmol)
print(nsnapshots)
a = CEC_track
# a = raw.reshape(nsnapshots, nmol, 3)

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

del a, CEC_track

# Compute the time in fs and put into a vector
t = np.arange(0,len(msd),0.25)
#print msd

filename2 = filename_title + '_msd.txt'
print('MSD written to', filename2, '!')

msdtxyz = np.array([np.hstack(k) for k in zip(t,msd)])
np.savetxt(filename2, msdtxyz)
