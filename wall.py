# This code is written by ZAHRA MOHAMMADYARLOO, PhD student at Techinical Uiniversity Dresden, Germany, after October 2019 
#for the project entitled "Co-Nonsolvency Effect on Phase Segregation of Polymer Solution" as a PhD project. This is 
#delivered to Professor Jens-Uwe Sommer as my PhD Supervisor in 31th of May 2022 before submitting the paper at Leibniz 
#Institute for Polymer Research Dresden IPF. 
#
#The input for this code is the output from the Lammps. This code calculates the average of height of movable membrane and plots the lz versus time. The vertical axis of Figure5 is extracted using this code.
#For Figure7-c, inverse of lz is used.

import matplotlib.pyplot as plt
import numpy as np

box_lengthx = 22.00
box_lengthy = 22.00
box_lengthz = 192.00
Lz = []
time = []

count1 = 0        
count2 = 0
n = int(box_lengthx*box_lengthy)   # total number of particles in the wall
f = open("wallmiddle_co1300_epsioln1.vmd", "r")
with open('xyzfile.dat','w') as fdata:

	for x in f:

		if count1==0 or count1==1 or count1==2 or count1==3 or count1==4 or count1==5 or count1==6 or count1==7 or count1==8:
			count1+=1
		else:
				fdata.write('{}'.format(x))
				count2+=1
				if count2==n:
					count1=0
					count2=0
fdata.close()
f.close()

print("done!\n\n\n")

a = np.genfromtxt("xyzfile.dat")
t = 0 #first time step in the vmd file (first frame)
tt = 10000000 #last step in the vmd file (last frame)
count4 = int((tt-t)/20000)
count3 = 0
count5 = 0
mysum = 0
cut_time = 0 
with open('xyz.dat','w') as adata:

	for i in range(count4):	
		zn = float((a[count3*n][1]+a[count3*n][2])*box_lengthz)
		t += 20000
		count3+=1
		if count3>cut_time:
			adata.write('{:7.6f} \n'.format(float(zn)))
			mysum += zn 
			Lz.append(zn)
			time.append(t)
			count5+=1
		
avez = mysum/count5
print(avez)
font = {'fontname':'Times New Roman'}
plt.title('NVT Simulation, F = 0.02, epsilon = , N_cosolvent = , <Lz> =   ', **font, fontsize=20)
plt.suptitle('Z direction non_periodic and wall in the middle, N/V = , monomer per chain = 60', **font, fontsize=15)
plt.xlabel('time (time_step)', **font, fontsize=15)
plt.ylabel('L_z', **font, fontsize=15)
plt.plot(time,Lz)
plt.show()
