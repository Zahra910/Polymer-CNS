# This code is written by ZAHRA MOHAMMADYARLOO, PhD student at Techinical Uiniversity Dresden, Germany, after October 2019 for the project entitled "Co-Nonsolvency Effect on Phase Segregation of Polymer Solution" as a PhD project. This is delivered to Professor Jens-Uwe Sommer as my PhD Supervisor in 31th of May 2022 before submitting the paper at Leibniz Institute for Polymer Research Dresden IPF. 

# This code generates initial positions of linear homopolymers + membrane + wall on the top and bottom of simulation box. The output would be input for "random_atom.py" code.
########################################################################################################################
import random
from math import *

box_lenx = 22
box_leny = 22
box_lenz = 438.07
mn = 60    # number of monomers
pn = 80   # number of polymers
co_sol = 0 # number of co_solvents
bond_length = 1   # bond length for FENE potential

###################################################################################################
##############  Generating the random walk

array = [[0 for i in range(3)] for j in range(pn*mn)]  # zeros array: dimonsion is 3*len(mn*pn)

with open('polydata.dat','w') as fdata:

    for i in range(pn):
        a = 0.0 
        array[0][0] = (random.random())*(box_lenx-1)
        array[0][1] = (random.random())*(box_leny-1)
        array[0][2] = (random.random())*(box_lenz-1)
        r = 0.0
        fdata.write('{}   {}   {}   {}\n'.format(i+1, array[0][0], array[0][1], array[0][2]))

        for j in range(mn-1):
            while r <= 1:
                c = random.random()
                if c < 0.5:
                     z = array[j][2]+1
                else:
                     z = array[j][2]-1
                     if z>=box_lenz-1.5:
                          z = array[j][2]-1
                     elif z<1.5:
                          z = array[j][2]+1       
                a = random.random()
                if a < 0.5:
                    x = array[j][0]+1
                else:
                    x = array[j][0]-1
                b = random.random()
                if b < 0.5:
                    y = array[j][1]+1
                else:
                    y = array[j][1]-1
                r = float(sqrt(((x-array[j][0])**2)+((y-array[j][1])**2)+((z-array[j][2])**2)))
            array[j+1][0] = array[j][0]+((x-array[j][0])/r)
            array[j+1][1] = array[j][1]+((y-array[j][1])/r)
            array[j+1][2] = array[j][2]+((z-array[j][2])/r)

            fdata.write('{}   {}   {}   {}\n'.format(i+1, array[j+1][0], array[j+1][1], array[j+1][2]))
            r = 0.0
####################################################################################################################
####### Generating co_solvent particles
    counter1 = 1
    for i in range(co_sol):
        a = random.random()
        b = random.random()
        c = random.random()
        x = a*(box_lenx-1)
        y = a*(box_leny-1)
        z = a*(box_lenz-1)
        fdata.write('{}   {}   {}   {}\n'.format(counter1, x, y, z))
        counter1+=1

###################################################################################################################
##### Generating the membrane
    x = 0.5
    y = 0.5
    z = box_lenz+20
    counter2 = 1
    wall_type = 2
    
    for i in range(box_leny):
        y = i+0.5
        for j in range(box_lenx):
            fdata.write('{}   {}    {}   {}\n'.format(counter2, x, y, z))
            x+=1
            counter2+=1
            if x==box_lenx+0.5:
                x=0.5
###################################################################################################################
#####  Generating wall on top of box
    x = 0.5
    y = 0.5
    z = (2*box_lenz)-0.5
    wall_type = 3
    
    for i in range(box_leny):
        y = i+0.5
        for j in range(box_lenx):
            fdata.write('{}   {}    {}   {}\n'.format(counter2, x, y, z))
            x+=1
            counter2+=1
            if x==box_lenx+0.5:
                x=0.5
####################################################################################################################
#####   Generating wall on floor of box
    x = 0.5
    y = 0.5
    z = 0.5

    for i in range(box_leny):
        y = i+0.5
        for j in range(box_lenx):
            fdata.write('{}   {}    {}   {}\n'.format(counter2, x, y, z))
            x+=1
            counter2+=1
            if x==box_lenx+0.5:
                x=0.5
    

