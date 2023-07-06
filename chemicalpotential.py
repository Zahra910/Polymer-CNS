# This code is written by ZAHRA MOHAMMADYARLOO, PhD student at Techinical Uiniversity Dresden, Germany, after October 2019 
#for the project entitled "Co-Nonsolvency Effect on Phase Segregation of Polymer Solution" as a PhD project. This is 
#delivered to Professor Jens-Uwe Sommer as my PhD Supervisor in 31th of May 2022 before submitting the paper at Leibniz 
#Institute for Polymer Research Dresden IPF. 
#
#The input for this code is the output from the Lammps. This code calculates the average chemical potential and plots \mu versus time. Chemical potential is calculated as: \rho/(1-\rho) where \rho is density of CNS particles above the membrane. Thus, the average of membrane height should be extracted from "wall.py" code and entered as an input in line 37 of this code: aveLz=...;  The horizontal axis of Figure5 is extracted using this code.


import matplotlib.pyplot as plt
import numpy as np
from math import *

mu = []
time = []
count1 = 0
count2 = 0

ff = open("cosolvent_co1000_epsilon1.vmd", "r")
with open('cosolvent.dat','w') as ffdata:
	n = 910  # total number of cosolvent particles

	for x in ff:

		if count1==0 or count1==1 or count1==2 or count1==3 or count1==4 or count1==5 or count1==6 or count1==7 or count1==8:
			count1+=1
		else:
				ffdata.write('{}'.format(x))
				count2+=1
				if count2==n:
					count1=0
					count2=0
ffdata.close()
ff.close()


print("done done!\n\n\n")
aveLz = 128.08
box_lengthy = 22.00
box_lengthx = 22.00
box_lengthz = 192.00

b = np.genfromtxt("cosolvent.dat")
t = 460500000 #first time step in the vmd file of wall (first frame)
tt = 470500000 #last step in the vmd file of wall (last frame)
count4 = int((tt-t)/20000)
count3 = 0
count5 = 0
mysum = 0
count6 = 0
count7 = 0
summu = 0
cut_time = 0 # favorite starting time step/20000 for measuring the chemical potential
volume = box_lengthx*box_lengthy*(box_lengthz-1.5-aveLz) #volume of upper reservior
with open('zposition.dat','w') as hdata:  # positions of cosolvent particles which are bigger than average length of middle wall, are recorded in this file

	for i in range(count4):
		for j in range(n):	
			zn = float((b[count7][1]+b[count7][2])*box_lengthz)
			count7+=1
			if zn > aveLz:
				hdata.write('{:7.6f} \n'.format(float(zn)))
				count6+=1
		count3+=1
		t+=20000

		rho = float((float(count6)/volume)*0.52)
		c = float(rho/(1-rho))
		
		if count6 !=0:
			mymu = float(log(c))
		
			count6 = 0

			if count3>cut_time:
				time.append(t)
				mu.append(mymu)
				summu+= mymu
				count5+= 1
avemu = float(summu/count5)
print(avemu)

print("thanks Zahra\n\n")

font = {'fontname':'Times New Roman'}
plt.title('NVT Simulation, F = 0.02 , N_cosolvent = 3780, epsilon =1.0, <mu> = -3.50 (after 0.3e8 time_step) ', **font, fontsize=20)
plt.suptitle('Z direction non_periodic and wall in the middle, N/V =  ,  monomer per chain = 60', **font, fontsize=15)
plt.xlabel('time (time_step)', **font, fontsize=15)
plt.ylabel('mu', **font, fontsize=15)
plt.plot(time,mu)
plt.show()
