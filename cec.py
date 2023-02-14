import numpy as np
from ase.io import read, write
from ase import neighborlist
import matplotlib.pyplot as plt
import os
import sys


class CEC():
    def __init__(self, working_dir, dumptraj, dumptraj_tag, noAtoms, xdatcar_prepend, timestep=0.25, rcut=1.2):
        '''
        1.
        2.
        3.
        4.
        '''

        self.working_dir = working_dir
        self.dumptraj = dumptraj
        self.dumptraj_tag = dumptraj_tag #"LAMMPS" Or "XDATCAR"
        self.noAtoms = noAtoms
        self.ts = timestep
        self.rcut = rcut
        self.xdatcar_path = f'{self.working_dir}/XDATCAR_CEC'
        self.xdatcar_prepend = xdatcar_prepend

    def get_atom_index(self, atoms, atom_char):
        '''
        Takes in Atoms object and returns the indices of where a particular atom is
        '''
        Element_idx_bool = atoms.symbols==atom_char
        Element_idx = [i for i, x in enumerate(Element_idx_bool) if x]
        return Element_idx

    def xdatcar2array(self):
        '''
        Uses the XDATCAR and stores the coordinates into an array after sorting based on the index.
        Add this to the main CEC code and add a condition that asks for XDATCAR or lammpstrj
        '''
        dump = open(self.dumptraj)
        count = 0
        direct_count = 0
        start = 0
        tag = 0
        frame = []
        atm = []
        for l in dump:

            if tag == 1:
                coords = [float(l.split()[0]),float(l.split()[1]),float(l.split()[2])]
                atm.append(coords)
                count += 1

            if count == self.noAtoms:
                tag = 0
                frame.append(atm)
                atm = []
                count = 0
            if l.split()[0] == 'Direct':
                direct_count += 1
                tag = 1
        self.frame = np.array(frame)
        
    def dump2array(self):
        '''
        Uses the dump traj file from lammps and stores the coordinates into an array after sorting based on the index.
        '''
        dump = open(self.dumptraj)
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
            if count == self.noAtoms:
                tag = 0
                frame.append(atm)
                atm = []
            if l.split()[0] == 'ITEM:':
                if l.split()[1] == 'ATOMS':
                    direct_count += 1
                    tag = 1
                    count = 0
        frame = np.array(frame)
        frame_sort = []
        for i in range (len(frame)):
            fr_i = frame[i]
            fr_is = fr_i[fr_i[:,0].argsort()]
            frame_sort.append(fr_is)
        frame_sort = np.array(frame_sort)
        self.frame = frame_sort

    def read_first_atoms_object(self):
        self.atomsObj = read(self.dumptraj, index=0)
        return self.atomsObj
    
    def convert_to_adjacency(matrix):
        start = 0
        res = []
        lst = []
        n = len(matrix)

        for i in range(n):
            res.append(lst*n)
        while start < n:
            y = matrix[start]
            for i in range(len(y)):
                if y[i] == True:
                    res[start].append(i)
            start += 1
        return res
    
    

    def CECmain(self, write_xdatcar=True, id0_path=None):
        '''
        Uses the dump traj file from lammps and stores the coordinates into an array after sorting based on the index.
        '''
        if write_xdatcar:
            print('DUMP to CEC in XDATCAR')
            print('------------------------')
            self.xdatcar_path = f'{self.working_dir}/XDATCAR_CEC'
            print(id0_path)
            xdatcar = open(self.xdatcar_path,'w+')
            
        if id0_path is None:
            print('READING AND WRITING FIRST ATOMS OBJECT')
            print('------------------------')
            id0 = self.read_first_atoms_object()
            # only run this for the run. 
            O_idx_bool = id0.symbols=='He'
            O_idx = [i for i, x in enumerate(O_idx_bool) if x]
            id0.symbols[O_idx] = 'O'

            C_idx_bool = id0.symbols=='H'
            C_idx = [i for i, x in enumerate(C_idx_bool) if x]
            id0.symbols[C_idx] = 'C'

            H_idx_bool = id0.symbols=='Li'
            H_idx = [i for i, x in enumerate(H_idx_bool) if x]
            id0.symbols[H_idx] = 'H'

            id0_path = f'{self.working_dir}/POSCAR_0'
            print('id0_path = ', id0_path)
            write(id0_path, id0) # Assumes POSCAR
            self.id0_path = id0_path # Save to use for next iterations
        else:
            id0 = read(id0_path)
            self.id0_path = id0_path
            print('FOUND BASE ATOMS OBJECT')
            
        H_idx = self.get_atom_index(id0, 'H')
        O_idx = self.get_atom_index(id0, 'O')
        
        # Create a neighborlist for all O atoms
        conn_bool_O = id0.get_all_distances(mic=True)[O_idx][:,O_idx] < 3 # assuming O-O distance is just under 3 Angs
        conn_bool_O
        O_adj_list = {}
        conn_all=conn_bool_O[0]*O_idx
        conn_fin=conn_all[conn_all!=0]

        for i_o,o in enumerate(O_idx):
            conn_all=conn_bool_O[i_o]*O_idx
            conn_fin=conn_all[conn_all!=0]
            O_adj_list[o] = conn_fin
        
        if self.dumptraj_tag == "LAMMPS":
            self.dump2array() # generates self.frame
        elif self.dumptraj_tag == "XDATCAR":
            self.xdatcar2array()
        
        # perform a global scan to generate the first CEC. Since the proton is manually added, 
        where_CEC=0 
        CEC_O_idx = O_idx[where_CEC]
        prevCEC = CEC_O_idx
        CEC_track = []
        
        # TODO: Add in a constraint where the next CEC is less than X distance away. X is calculated as the 'min' possible O-O distance in GOH. 
        # If there are instances with more than 1 potential CECs then save those instances (just the index) 

        for i,f in enumerate(self.frame):
            if len(f.T) == 4:
                pos_f = f[:,1:]
            elif len(f.T) == 3:
                pos_f = f

            id0.set_scaled_positions(pos_f)
            if i == 0:
                # conn_bool = id0.get_all_distances(mic=True)[O_idx][:,H_idx] < self.rcut # rcut
                # conn_H = np.sum(conn_bool, axis=1)
                # where_CEC = np.argwhere(conn_H>=2) # detect CEC, ranges from 0 to total number of O atoms 
                # CEC_O_idx = O_idx[where_CEC[0,0]]

                # Finds the O that has the two closest H atoms. 
                conn_OH = id0.get_all_distances(mic=True)[O_idx][:,H_idx]
                sorted_conn_OH = np.sort(conn_OH)
                where_CEC = np.argmin(np.sum(sorted_conn_OH[:,0:2], axis=1))
                CEC_O_idx = O_idx[where_CEC]
                
            else:
                pot_CEC = []
                for o in O_adj_list[prevCEC]:
                    conn_bool = id0.get_distances(o, H_idx, mic=True) < self.rcut
                    conn_H = np.sum(conn_bool)
                    if conn_H >= 2:
                        pot_CEC.append(o)
                
                if len(pot_CEC)>=2:
                    # iterate over all potential CECs
                    # this needs to be changed. Try getting the O with the two closest H atoms instead of this. The smallest O-O distance does not make sense.
                    min_dist = np.inf
                    for pc in pot_CEC:
                        O_O_dist = id0.get_distance(prevCEC, pc)
                        if O_O_dist < min_dist:
                            min_dist = O_O_dist
                            CEC_O_idx = pc
                            
                elif len(pot_CEC) == 1:
                    CEC_O_idx = pot_CEC[0]
                    
            scaled_pos = id0.get_scaled_positions()[CEC_O_idx]
            prevCEC = CEC_O_idx
            if write_xdatcar:
                xdatcar.write(f'Direct configuration= {i+1}\n')
                xdatcar.write(f' {scaled_pos[0]:1.4f} {scaled_pos[1]:1.4f} {scaled_pos[2]:1.4f}\n')
        xdatcar.close()        
        with open(self.xdatcar_path, 'r+') as f:
            content = f.read()
            f.seek(0, 0)
            f.write(self.xdatcar_prepend.rstrip('\r\n') + '\n' + content)

    def xdatcar2unwrapped(self):
        print('XDATCAR to UNWRAPPED xyz')
        print('------------------------')
        xdatcar2unwrappedxyz_path = '/bgfs/kjohnson/ska31/1DeePMD/paper-graphanol/oneside/CEC_code/source_code/xdatcar2unwrappedxyz.py'
        os.chdir(self.working_dir)
        os.system(f'python {xdatcar2unwrappedxyz_path} -i {self.xdatcar_path}')

    def calc_msd(self):
        print('UNWRAPPED xyz to MSD')
        print('------------------------')
        calc_msd_path = '/bgfs/kjohnson/ska31/1DeePMD/paper-graphanol/oneside/CEC_code/source_code/xyztomsd_un.py'
        os.chdir(self.working_dir)
        os.system(f'python {calc_msd_path} -i coordunwrapped_CEC.xyz')
