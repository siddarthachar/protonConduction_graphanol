import numpy as np
import ase
from ase.io import read, write
import os
import time
import pandas as pd

def chunks(generator, chunk_size):
    """Yield successive chunks from a generator"""
    chunk = []

    for item in generator:
        if len(chunk) >= chunk_size:
            yield chunk
            chunk = [item]
        else:
            chunk.append(item)

    if chunk:
        yield chunk

class CEC():
    def __init__(self, working_dir, dumptraj, noAtoms,calc_msd_only=False, timestep=0.25, rcut=1.2, ):
        '''
        1.
        2.
        3.
        4.
        '''

        self.working_dir = working_dir
        self.dumptraj = dumptraj
        self.ts = timestep
        self.rcut = rcut
        self.calc_msd_only = calc_msd_only
        self.write_msd_path = f'{self.working_dir}/msd_2.dat'
        self.write_cec_path = f'{self.working_dir}/cec_2.dat'

    def get_atom_index(self, atoms, atom_char):
        '''
        Takes in Atoms object and returns the indices of where a particular atom is
        '''
        Element_idx_bool = atoms.symbols==atom_char
        Element_idx = [i for i, x in enumerate(Element_idx_bool) if x]
        return Element_idx

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
        cec = []
        get_TOT_CEC = 0
        TOT_CEC = -1
        # consider the cases where there are no cecs. store as a list. make use of the NUMBER OF ATOMS to get the # before hand
        for l in dump:
        #     if direct_count == 10000: # REMOVE THIS LATER ###### 
        #         break
            if l.split()[0] == 'ITEM:':
                if l.split()[1] == 'TIMESTEP':
                    count = 0
                    continue
                    
            if get_TOT_CEC == 1:
                TOT_CEC = int(l.split()[0])
    #             print("TOT_CEC: ", direct_count, TOT_CEC, "tag: ", tag)
                if TOT_CEC == 0:
                    cec.append([])
                    tag = 0
                get_TOT_CEC = 0

            if l.split()[0] == 'ITEM:':
                if l.split()[1] == 'NUMBER':
                    get_TOT_CEC = 1
                    continue
            
                    
            if tag == 1:
                if TOT_CEC > 0:
    #                 print(direct_count, l)
    #                 print('tag: ', tag)
                    cec.append([int(l.split()[0])-1,float(l.split()[3]),float(l.split()[4]),float(l.split()[5])]) # -1 for python
                    count += 1
                    
                
            if count == TOT_CEC:
                tag = 0
                TOT_CEC = -1
                frame.append(cec)
                cec = []
                
            if l.split()[0] == 'ITEM:':
                if l.split()[1] == 'ATOMS':
                    direct_count += 1
                    if TOT_CEC > 0:
                        tag = 1
                    count = 0
        self.frame = frame

    def CECmain(self, id0_path=None):
        if self.calc_msd_only == False:
            if id0_path is None:
                raise ValueError('Enter a POSCAR path for id0. Currently None.')
            else:
                id0 = read(id0_path)
                self.id0_path = id0_path
                print('FOUND BASE ATOMS OBJECT')

            H_idx = self.get_atom_index(id0, 'H')
            O_idx = self.get_atom_index(id0, 'O')
            
            self.dump2array()

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

            cec_list = []
            tot_pot_cec = 0
            for i,potcec in enumerate(self.frame):
                if i == 0:
                    prev_cec = potcec[0]
                    cec_list.append(prev_cec)
                else:
                    
                    if len(potcec) == 1:
                        if len(potcec[0]) == 0: # no CEC case
                            cec_list.append(prev_cec)
                            tot_pot_cec = 0
                        elif potcec[0][0] in O_adj_list[prev_cec[0]]: # only consider the next CEC if it is in the O 1-NN

                            cec_list.append(potcec[0])
                            prev_cec = potcec[0]
                            tot_pot_cec = 1
                            
                        else:
                            cec_list.append(prev_cec)
                            tot_pot_cec = 1
                    elif len(potcec) >= 2:
                        tot_pot_cec = len(potcec)
                        sameCEC = False
                        otherCEC = []
                        count_O_NN_CEC = 0
                        # Give priority to the previous CEC
                        for pot in potcec:
                            if pot[0] in O_adj_list[prev_cec[0]]:
                                count_O_NN_CEC += 1
                                if prev_cec[0] == pot[0]:
                                    sameCEC = True
                                    temp_pot = pot
                                    break # found the cec here
                                else:
                                    otherCEC.append(pot)
                        
                        if sameCEC:
                            cec_list.append(temp_pot)
                        else:
                            if len(otherCEC):
                                cec_list.append(otherCEC[0]) # append the first one. will mostly be 1. If not then change the code.
                            else:
                                cec_list.append(prev_cec)
                        if count_O_NN_CEC == 0:
                            cec_list.append(prev_cec)
            self.cec_list = np.array(cec_list)
            nsnapshots = int(len(self.cec_list))
            print('No. snapshots: ',nsnapshots)
            print((self.cec_list).shape)
            np.savetxt(f'{self.write_cec_path}', self.cec_list)
                                            
        # Compute MSD from this
        else:
            id0 = read(id0_path)
            self.id0_path = id0_path
            print('FOUND BASE ATOMS OBJECT')
            # self.cec_list = np.loadtxt(self.write_cec_path)
            df = pd.read_csv(self.write_cec_path,  header=None, delimiter=' ')
            self.cec_list = df.to_numpy()
            nsnapshots = int(len(self.cec_list))

        # unwrapping coordinates
        self.cec_list = self.cec_list[:,1:]
        Natoms = 1 # 12 O atoms
        unitcell = id0.get_cell()
        element_names = 'O'
        scale=1

        prev_step = np.zeros((Natoms,3))
        curr_step = np.zeros((Natoms,3))
        unwrapped = np.zeros((Natoms,3))
        unwrapped_list = []

        for ii in range(len(self.cec_list)):
        #     cec_i = np.expand_dims(cec_list[ii],0)
            scaled_cec = np.expand_dims(np.dot(self.cec_list[ii],np.linalg.inv(unitcell)),0)
            cart_cec_wrapped = np.expand_dims(self.cec_list[ii],0)
            
            if ii == 0:
                prev_step = scaled_cec
                unwrapped = np.copy(prev_step)
                unwrapped_list.append(cart_cec_wrapped[0])
            else:


                for i in range(Natoms):
                    curr_step[i] = scaled_cec[i]
                    for j in range(3):
                        dd = prev_step[i][j] - curr_step[i][j]
                        if np.fabs(dd) > 0.5:
                            unwrapped[i][j] = unwrapped[i][j] + curr_step[i][j] + np.sign(dd)*1.0 - prev_step[i][j]
                        else:
                            unwrapped[i][j] = unwrapped[i][j] + curr_step[i][j] - prev_step[i][j]                    
                prev_step = np.copy(curr_step)
                xyzcoords = np.dot(unwrapped,unitcell)
                unwrapped_list.append(xyzcoords[0])

        self.cec_List = np.array(unwrapped_list)

# calculate msd
        half_width = int(nsnapshots / 2)
        n_cmsds = 100
        dtspan = int(half_width / n_cmsds)
        span_starts = [x for x in range(nsnapshots)
            if x % dtspan == 0 and x + half_width <= nsnapshots]
        n_msds = len(span_starts)
        a = self.cec_list[:,1:3].reshape(nsnapshots, 1, 2) # changed to only inclue x and y 
        # msd = np.array([np.mean((a[k] - a[k0]) ** 2, 0)
        #     for k0 in span_starts  for k in range(k0, k0 + half_width)])
        # trying out generator
        msd_gen = (np.mean((a[k] - a[k0]) ** 2, 0)
            for k0 in span_starts  for k in range(k0, k0 + half_width))

        msd_3d = np.array([np.array(m, dtype=np.uint16) for i, m in enumerate(chunks(msd_gen, half_width))], dtype=np.uint16)
        msd_cc = msd_3d.mean(0)
        # self.msd = msd.reshape(n_msds, half_width, 3).mean(0)
        
        # open(self.write_msd_path, "w").write('\n'.join(f'{i*0.25/1000:.3E}\t{m[0]:.3E}\t{m[1]:.3E}\t{0.5*(m[0]+m[1])/2/(i*0.25/1000):.3E}' for i, m in enumerate(msd_gen)))
        open(self.write_msd_path, 'w').write('\n'.join(f'{i*0.25/1000:.5E}\t{0.5*(m[0]+m[1])/2/(i*0.25/1000):.3E}' for i, m in enumerate(msd_cc)))
        # self.diff = np.mean(self.msd[:,0:2],axis=1)/2/self.t 
        # print('len diff, msd shape: ',len(self.diff), (self.msd).shape)
        # np.savetxt(self.write_msd_path, [self.t, self.msd[:,0], self.msd[:,1], self.msd[:,2], self.diff])
