import numpy as np
import ase
from ase.io import read, write
import os
import time
import sys

sys.path.append('/bgfs/kjohnson/ska31/1DeePMD/paper-graphanol/oneside/CEC_code')

temp = 300
i = int(os.getcwd().split('/')[-1])

xdatcar_prepend =' O\n\
      1.0000000000000000\n\
         10.2929999999999993    0.0000000000000000    0.0000000000000000\n\
         -3.9039999999999981    6.7619263527488975    0.0000000000000000\n\
         0.0000000000000000    0.0000000000000000   20.0000000000000000\n\
      O\n\
      1'

from cec import CEC

id0_path = '/bgfs/kjohnson/ska31/1DeePMD/paper-graphanol/oneside/CEC_code/sample_traj/POSCAR_0'
c24 = CEC(f'/bgfs/kjohnson/ska31/1DeePMD/paper-graphanol/oneside/Diffusion_CEC_trial-2/multi-24C-random-start/{temp}K/{i}/',
          f'/bgfs/kjohnson/ska31/1DeePMD/paper-graphanol/oneside/Diffusion_CEC_trial-2/multi-24C-random-start/{temp}K/{i}/prod.lammpstrj',
          "LAMMPS",
          61,
          xdatcar_prepend,
          0.25,
          1.2)
c24.CECmain(id0_path=id0_path)
c24.xdatcar2unwrapped()
c24.calc_msd()
