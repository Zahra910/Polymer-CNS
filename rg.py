# This code is written by ZAHRA MOHAMMADYARLOO, PhD student at Techinical Uiniversity Dresden, Germany, after October 2019 
#for the project entitled "Co-Nonsolvency Effect on Phase Segregation of Polymer Solution" as a PhD project. This is 
#delivered to Professor Jens-Uwe Sommer as my PhD Supervisor in 31th of May 2022 before submitting the paper at Leibniz 
#Institute for Polymer Research Dresden IPF. 
#
#The input for this code is the output from "readfilerg.py". This code calculates the Radius of Gyration od a single polymer by time and ensemble averaging. Also, it plots Rg versus time. The vertical axis of Figure 7-a is calculated using this code. 

import numpy as np
from math import *
#import matplotlib.pyplot as plt

f = np.loadtxt("xyz.dat")

N = 60 # monomer per chain
np = 160 # total number of polymers
nc = 95 #cosolvent
nw = int(3*484) # N_walls
nm = int(N*np)
RG = []
t = 0 # starting time step in vmd file
time = []
radius_gyration = []
frames = 500
deltaRR = 0

for l in range(frames):
	for k in range(np):
		for i in range(N):
			for j in range(i,N):
				deltax = f[i+(k*N)+((nm+nc+nw)*l)][1]-f[j+(k*N)+((nm+nc+nw)*l)][1]
				deltay = f[i+(k*N)+((nm+nc+nw)*l)][2]-f[j+(k*N)+((nm+nc+nw)*l)][2]
				deltaz = f[i+(k*N)+((nm+nc+nw)*l)][3]-f[j+(k*N)+((nm+nc+nw)*l)][3]
				deltaRR +=(deltax**2)+(deltay**2)+(deltaz**2)
		R = float(sqrt(deltaRR/N**2))
		RG.append(R)
		deltaRR = 0
	#print(len(RG))
	#print(RG)
	sumRG = sum(RG)
	count = len(RG)
	Rg = sumRG/count
	RG = []
	t+=20000
	time.append(t)
	radius_gyration.append(Rg)

count1 = len(radius_gyration)
sumrg = sum(radius_gyration)
ave_Rg = sumrg/count1
print('<Rg> =  ', ave_Rg)

font = {'fontname':'Times New Roman'}
plt.xlabel('time (time_step)', **font, fontsize=15)
plt.ylabel('Radius of Gyration', **font, fontsize=15)
plt.title('NVT Simulation', **font, fontsize=20)
plt.suptitle('Z direction non_periodic and wall in the middle, phi = 0.05, monomer perchain = 60', **font, fontsize=15)
plt.plot(time,radius_gyration)
plt.show() 
