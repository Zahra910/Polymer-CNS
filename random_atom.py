##########################################################################################################################
# This code is written by ZAHRA MOHAMMADYARLOO, PhD student at Techinical Uiniversity Dresden, Germany, after October    #
#2019 for the project entitled "Co-Nonsolvency Effect on Phase Segregation of Polymer Solution" as a PhD project. This is#
# delivered to Professor Jens-Uwe Sommer as my PhD Supervisor in 31th of May 2022 before submitting the paper at Leibniz #
# Institute for Polymer Research Dresden IPF.                                                                            #
#                                                                                                                        #
# This code generates input file for the Lammps code:  "in.figure3red"                                                    #
#########################################################################################################################

import numpy as np
from math import *

# reading the positions (x,y,z) of the polymers from generated data file using AVOGADRO
f = np.genfromtxt("polydata.dat")

# number of monomer to create
np = 80 # number of polymers
nm = 60 # number of monomers
nco_sol = 0 # number of co_solvency
natoms = nm*np 
bonds = natoms-np
atom_type = 2
bond_type = 1
counter2 = 0


# size of the box size in LJ units
box_sizex = 22.000000
box_sizey = 22.000000
box_sizez = 438.070000


#number of particles needed for making wall
wall_particles = int(3*box_sizex*box_sizey)

# open a text file to write the data as an input for LAMMPS
with open('data.polymer','w') as fdata:
       # first line is a comment 
       fdata.write('data description\n\n')
       
       #specify number of atoms, bonds, angles, dihedrals
       fdata.write('{} atoms\n'.format(natoms+nco_sol+int(3*(box_sizex)*(box_sizey))))
       fdata.write('{} bonds\n\n\n'.format(bonds))

       #specify types of atoms, bonds, angles, dihedrals and impropers
       fdata.write('{} atom types\n'.format(3))
       fdata.write('{} bond types\n\n\n'.format(1))

       #specify box dimensions
       fdata.write('{}   {}  xlo  xhi\n'.format(0.0, box_sizex))
       fdata.write('{}   {}  ylo  yhi\n'.format(0.0, box_sizey))
       fdata.write('{}   {}  zlo  zhi\n\n\n'.format(0.0, box_sizez*2))

       #specify masses
       fdata.write('Masses\n\n')
       fdata.write('{} {}\n'.format(1,1))
       fdata.write('{} {}\n'.format(2,1))
       fdata.write('{} {}\n\n'.format(3,1))
       
       #atom section
       fdata.write('Atoms\n\n')

       for i in range(natoms):
          fdata.write(' {}   {}   {}   {}   {}   {}\n'.format(i+1, int(f[i][0]), 1, f[i][1], f[i][2], f[i][3]))
       for i in range(nco_sol):
          fdata.write(' {}   {}   {}   {}   {}   {}\n'.format(i+natoms+1, int(f[i+natoms][1]), f[i+natoms][0], f[i+natoms][2], f[i+natoms][3], f[i+natoms][4]))
       for i in range(len(f)-natoms-nco_sol):
          walltype = 2
          if i>=(box_sizex*box_sizey):
              walltype = 3	
          fdata.write(' {}   {}   {}   {}   {}   {}\n'.format(natoms+nco_sol+i+1, int(f[i+natoms+nco_sol][0]), walltype, f[i+natoms+nco_sol][1], f[i+natoms+nco_sol][2], f[i+natoms+nco_sol][3]))
       #bond section
       fdata.write('\nBonds\n\n')

       
       for i in range(np):
           for j in range(nm-1):
               counter1 = (((i+1)*nm)-nm+j+1)
               counter2+=1
               fdata.write('{}  {}  {}  {}\n'.format(counter2, bond_type, counter1, counter1+1))


