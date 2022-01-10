## Copyright: Siddarth Achar
## University of Pittsburgh

# This function is designed to break an XDATCAR file into tiny equally sized XDATCAR files
# This is useful for file reading, especially if you use the default ase.io reader. 
# I will be using this to compute properties at a fixed set of time frames.


import numpy as np
from ase.io import read, write
from ase import neighborlist
import os

class breaking:

    def __init__(self, working_dir, timestep, lenFile, TOTsteps, noAtoms):
        '''
        Enter the following:
        1. Working directory: working_dir
        2. Timestep (in ps): ts
        3. Length of File: lenFile
        4. Total number of steps: TOTsteps
        5. Number of Atoms: Number of Atoms
        '''
        self.noAtoms = noAtoms
        self.ts = timestep
        self.lenFile = lenFile
        self.TOTsteps = TOTsteps
        self.working_dir = working_dir

    def line_prepender(filename, line):
        '''
        Function that takesin filename and a line as string and prepends that line to the file.
        '''
        with open(filename, 'r+') as f:
            content = f.read()
            f.seek(0, 0)
            f.write(line.rstrip('\r\n') + '\n' + content)

    def xdatcar_header_string(self, dir_in):
        '''
        Function that extracts the first 7 lines of an XDATCAR file (the header) and returns 
        it as a string
        '''
        xdatcar_path = f'{self.working_dir}/{dir_in}/XDATCAR'
        file_inp=open(xdatcar_path)
        lines = file_inp.readlines()
        return "".join(lines[0:7])


    def breakxdatcar(self, dir_in):
        '''
        This function breaks the XDATCAR files based on the total number of steps and 
        the length of each length file. 
        The code assumes that {working_dir} has a bunch of independent directories. Each
        independent directory has an XDATCAR file. The code is aimed at breaking the XDATCAR
        file into several files and saved into a directory called {working_dir}/{dir_in}/steps/.
        '''
        inpath = lambda x:f'{self.working_dir}/{x}/XDATCAR'   # input of each of the XDATCAR files.
        
        prependThis = breaking.xdatcar_header_string(self, dir_in)

        steps_path = f'{self.working_dir}/{dir_in}/'
        list_dir = os.listdir(steps_path)
        # print(list_dir)
        if 'steps' not in list_dir:
            os.mkdir(f'{steps_path}steps')

        outpath = lambda x, y :f'{self.working_dir}/{x}/steps/XDATCAR.{y}'
        steps=np.arange(0, self.TOTsteps, self.lenFile)
        start_collect = False
        count_atoms = 0
        counter_direct = 0
        count_step = 0
        print(f'Total number of XDATCAR sub files: {len(steps)} \nEach file contains: {self.lenFile} steps = {self.lenFile*self.ts} ps')
        file_inp = open(inpath(dir_in),'r')
        lines = file_inp.readlines()

        for l in lines:
            if start_collect == True:
                file_out.write(str(l))
                count_atoms += 1

            if count_atoms==self.noAtoms:
                start_collect = False
                count_atoms = 0

            sp = l.split()
            if sp[0] == 'Direct':
                if int(sp[2]) == steps[count_step]+1: 
                    # goes in when direct is step[i]
                    
                    if count_step > 0: # for the next step, we would have to preprend the header adn close the file
                        file_out.close()
                        breaking.line_prepender(outpath(dir_in,steps[count_step-1]),prependThis)  # count_step was updated initially
                        
                    
                    file_out = open(outpath(dir_in,steps[count_step]),'w+')  # creates an xdatcar file to start writing
                    count_step += 1 # the next time steps[i+1] is encountered, a new file is created.
                    counter_direct = 0
                
                counter_direct += 1
                file_out.write(f'Direct configuration=\t{counter_direct}\n')
                start_collect = True

if __name__ == "__main__":
    breakIT = breaking(working_dir='.',timestep=0.00025, lenFile=200, TOTsteps=80000, noAtoms=61)
    breakIT.breakxdatcar(10) # An example is file directory 10.
   
