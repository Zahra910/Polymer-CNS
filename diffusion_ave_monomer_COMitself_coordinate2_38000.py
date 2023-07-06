# This code is written by ZAHRA MOHAMMADYARLOO, PhD student at Techinical Uiniversity Dresden, Germany, after October 2019 
#for the project entitled "Co-Nonsolvency Effect on Phase Segregation of Polymer Solution" as a PhD project. This is 
#delivered to Professor Jens-Uwe Sommer as my PhD Supervisor in 31th of May 2022 before submitting the paper at Leibniz 
#Institute for Polymer Research Dresden IPF. 
#
#The input for this code is the output from "readfile.py". This code calculates the Mean Squared Displacement (MSD) of the monomers with respect to the center of mass (COM). We find  subdiffusion in short times and plateau in t>\tau_R Rouse time. Thus we can find Rouse time using this code. The MSD is calculated in xy plane because in z direction the simulation is non-periodic. Figure10-a,b are the results from this code.

import matplotlib.pyplot as plt
import numpy as np
from math import *


a = np.loadtxt("xyz38000.dat")
with open('msd_COMitself_allmonomers_N120_CNS38000_espilon1.dat','w') as fdata:
	nm = 9600  # total monomers
	nwall = int(3*484) # total number of wall particles 
	ncos = 36981  # number of cosolvent
	frames = 2500
	t = 0
	N = frames
	msdmonomer = []
	myMSD = []
	time = []
	xx = 0
	yy = 0
	zz = 0
	l = int(1*1) # For generating MSD of Center of Mass
	NP = 1	# number of polymers
	pn = 9600 # monomer per chain
	nn = int(N*NP)

	array1 = [[0 for i in range(3)] for j in range(frames)]  # zeros array: dimonsion is 3*frames; Array1 is used for collecting the x, y, z positions of one specific COM from different frames.
	array2 = [[0 for i in range(l)] for j in range(int(frames))]  # zeros array: dimonsion is l*(N-2) for saving MSD of individual monomers or cosolvent particles. In each colum of array2, MSD of one specific monomer is saved.
	array3 = [[0 for i in range(3)] for j in range(nn)] # zero array: dimonsion is 3*(frames*number of polymers); Array3 is used for collecting the center of mass of each polymer in each time step. In each frame, there are NP center of masses: xcm, ycm and zcm. This will be used for changing the coordinates of monomers to the center of mass coordinate.


# Here I am going to find the center of mass of polymers and save them in array3:
	for k in range(frames): # loop in frames
		for i in range(NP): #loop in number of polymers
			xcm = 0
			ycm = 0
			zcm = 0
			for j in range(pn): # loop in monomers of each polymer
				xcm+= a[j+(i*pn)+(k*(nwall+nm+ncos))][0]
				ycm+= a[j+(i*pn)+(k*(nwall+nm+ncos))][1]
				zcm+= a[j+(i*pn)+(k*(nwall+nm+ncos))][2]

			array3[i+(k*NP)][0] = float(xcm/pn)
			array3[i+(k*NP)][1] = float(ycm/pn)
			array3[i+(k*NP)][2] = float(zcm/pn)	
# Until here we have found the center of mass of each polymer in each frame and have saved in array3.
	
	for i in range(NP):	#(NP):
		for k in range(frames):
			array1[k][0] = array3[i+(k*NP)][0] # X = x_cm
			array1[k][1] = array3[i+(k*NP)][1] # Y = y_cm
			array1[k][2] = array3[i+(k*NP)][2] # Z = z_cm
			
		for n in range(1, N-1):
			for m in range(1, N-n):
				xx+= (array1[m+n-1][0]-array1[m-1][0])**2
				yy+= (array1[m+n-1][1]-array1[m-1][1])**2
				zz+= (array1[m+n-1][2]-array1[m-1][2])**2
			MSD = float((xx+yy)/(N-n))
			array2[n-1][i] = MSD
			xx = 0
			yy = 0
			zz = 0
############# In this part, I am averaging the MSD of monomers. First Ido summation in each row (mymsd), then dividing by number of msd's which is number of columns (count).

	mymsd = 0
	count = int(len(array2[0]))
	print(count)
	print(len(array2))
	for i in range(len(array2)):
		for j in range(count):
			mymsd+=array2[i][j]
		aveMSD = float(mymsd/count)
		#aveMSD = float(mymsd)
		myMSD.append(aveMSD)
		t+=20000
		time.append(t)
		mymsd = 0
		fdata.write('{:7.6f}   {:7.6f}\n'.format(float(t), float(aveMSD)))		


##################   Plot ln(MSD) as a function of ln(time)

font = {'fontname':'Times New Roman'}
plt.title('MSD as a function of time, log-log plot', **font, fontsize=20)
plt.suptitle('Z direction non_periodic and wall in the middle, N/V = , monomer per chain = 60', **font, fontsize=15)
plt.xlabel('log(time)', **font, fontsize=15)
plt.ylabel('log(<R^2>)', **font, fontsize=15)
plt.plot(time,myMSD)
plt.show()

		
