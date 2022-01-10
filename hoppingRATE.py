## Copyright: Siddarth Achar
## University of Pittsburgh

## This code can be used to copmute the rate of proton hopping for an MD simulation. The user would need 
## to have their files sorted in this way: {working directory}/{independent run}/XDATCAR
## each independent run should have directory called "steps" that contains the XDATCAR separated into several



import numpy as np
from ase.io import read, write
from ase import neighborlist
import matplotlib.pyplot as plt
import os

class hoppingrate:

    def __init__(self, working_dir, timestep, lenFile, TOTsteps, atoms2include, stepsINstep):
        '''
        Enter the following:
        1. Working directory: working_dir
        2. Timestep (in ps): ts
        3. Length of File: lenFile
        4. Total number of steps: TOTsteps
        5. Atom indices involved in computing the hopping rate: atoms2include 
        '''
        self.atoms2include = atoms2include
        self.ts = timestep
        self.lenFile = lenFile
        self.TOTsteps = TOTsteps
        self.working_dir = working_dir
        self.stepsINstep = stepsINstep

    def updateconnectionmatrix(ase_atoms, curr_i, curr_matrix):
        '''
        There are cases where the connectivity matrix generates repeats and this function is aimed 
        at comparing the repeated connections and not eliminating the shortest distance connection. 
        '''
        allH_index=[k[0] for k in curr_matrix]
        if curr_i[0] in allH_index:
            repeat = np.argwhere(np.array(allH_index)==curr_i[0])[0,0]
            curr_dist = ase_atoms.get_distance(curr_i[0], curr_i[1], mic=True)
            repeat_dist = ase_atoms.get_distance(curr_matrix[repeat][0], curr_matrix[repeat][1], mic=True)

            if repeat_dist >= curr_dist:
                curr_matrix.remove(curr_matrix[repeat])
                curr_matrix += [[curr_i[0], curr_i[1]]]
        else:
            curr_matrix += [[curr_i[0], curr_i[1]]]

        return curr_matrix

    def check_OH_connectivity(chemSymb, curr_i):
        '''
        Function to check if there are nearest neighbor connections where the connection isnt between 
        H and O. 
        '''
        if chemSymb[curr_i[0]] != 'H':
            print("Alert i[0] NOT H")
        if chemSymb[curr_i[1]] != 'O':
            print("Alert i[1] NOT O")

    def compHopping(self, dir_in):
        '''
        Function that calculates the hopping rate for several XDATCAR files and then calculates a trend over several time scales.
        '''
        struct_path=lambda x,y:f'{self.working_dir}/{x}/steps/XDATCAR.{y}'
        goh=read(struct_path(dir_in,0))  # Assuming that all independent runs will have the same atomic configuration
        chemSymb=goh.get_chemical_symbols()

        prev_matrix = [0]  # Defining a previous matrix that compares the connectivity matrix of iteration i with i-1
        hopping_avg = []   # Iteratively store the hopping after each iteration. 

        for k in range(0,self.TOTsteps - self.lenFile , self.lenFile ):
            hopping_step=[]   # inner temp variable
            for j in range(0,self.lenFile,self.stepsINstep):
                c=0
            #     print(f't = {j*xdatcar_steps*ts:1.2f} ps')
                goh=read(struct_path(dir_in,k), index=j)
                cutOff = neighborlist.natural_cutoffs(goh)
                neighborList = neighborlist.NeighborList(cutOff, self_interaction=False, bothways=True)
                neighborList.update(goh)
                matrix = neighborList.get_connectivity_matrix(sparse=True)
                matrix_column=[]

                for i in matrix.getH().keys():   # Uses the connectivity matrix to update the matrix
                    if i[0] in self.atoms2include and goh.get_distance(i[0], i[1], mic=True) < 1.5:
                        matrix_column = hoppingrate.updateconnectionmatrix(goh, i, matrix_column)
                        hoppingrate.check_OH_connectivity(chemSymb, i)

                matrix_column = np.array(matrix_column)
                matrix_column = matrix_column[np.argsort(matrix_column[:,0])] 

                if len(prev_matrix)>1:
                    diff_matrix=matrix_column[:,1]-prev_matrix[:,1]      # Find the difference between connectivity of i and i-1 connected O atoms
                    hops=len(diff_matrix[diff_matrix!=0])
                    hopping_step += [hops]
                prev_matrix = matrix_column
            mean_hop= np.mean(np.array(hopping_step))
            hopping_avg += [mean_hop]

        return hopping_avg

if __name__ == '__main__':
    h=hoppingrate(working_dir='.', timestep=0.00025, lenFile=200, TOTsteps=80000,
     atoms2include=np.arange(48,61), # These are the indices (from 0) of H atoms that I want to compute the connectivity matrix of.
     stepsINstep=25)

    rate_10=h.compHopping(10) # 10th independent run. 
    print(rate_10)
    
    np.savetxt('10/hopping.txt', rate_10)

    time=np.linspace(0,20,len(rate_10)) # time steps

    plt.plot(time, rate_10,alpha=0.3,color='b')
    mean_hop=np.mean(rate_10)
    plt.plot(time, mean_hop+0*np.array(rate_10),'b--',label=f'24C 1 proton 800 K : {mean_hop:1.2f}')
    plt.savefig('10/hopping.png',dpi=300)