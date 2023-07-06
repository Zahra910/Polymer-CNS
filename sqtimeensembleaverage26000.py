# This code is written by ZAHRA MOHAMMADYARLOO, PhD student at Techinical Uiniversity Dresden, Germany, after October 2019 
#for the project entitled "Co-Nonsolvency Effect on Phase Segregation of Polymer Solution" as a PhD project. This is 
#delivered to Professor Jens-Uwe Sommer as my PhD Supervisor in 31th of May 2022 before submitting the paper at Leibniz 
#Institute for Polymer Research Dresden IPF. 
#
#The input for this code is the output from the "readfile.py". This code calculates the scattering function S(q) versus q. By fitting to the results, the Flory exponent \nu can be extracted. Flory exponen is used for measuring the size of polymer. \nu can be 1/2 (ideal chain: in concentrated regime and dilute regime), 1/3 (condesed completely packed chain) and 0.588 (good solvent).
# The ouput of this code is shown in Figure8.  

import numpy as np
from math import *
import matplotlib.pyplot as plt

f = np.loadtxt("xyz.dat") # positions (x,y,z) of monomers
with open('Sq.dat','w') as fdata:
	q = []
	sq = []
	Lx = 100 # length of box in x direction
	Ly = 100 # length of box in z direction
	nx = list(range(-Lx,Lx+1,1))
	ny = list(range(-Ly,Ly+1,1))
	N = 120 # monomer per chain
	sum1 = 0
	sum2 = 0
	np = 80 # number of polymers
	sumSQt = 0
	sumSQ = 0
	nm = int(N*np) # total monomers
	count = 0
	nc = 25415 # number of co-solvent
	nwall = int(3*484)
	frames = 200


	for i in range(len(nx)):
		for j in range(len(ny)):
			qx = (2*pi/Lx)*nx[i]
			qy = (2*pi/Ly)*ny[j]
			Q = sqrt((qx**2)+(qy**2))
			q.append(Q)			
			
			for l in range(np):
				for t in range(frames):
					for k in range(N):
						qdotr = qx*f[k+(l*N)+(nm+nc+nwall)*t][1]+qy*f[k+(l*N)+(nm+nc+nwall)*t][2]
						sum1+=cos(qdotr)
						sum2+=sin(qdotr)
					SQt = (1/N)*((sum1**2)+(sum2**2))
					sumSQt+= SQt
					count+=1				
					sum1 = 0
					sum2 = 0

				SQQ = sumSQt/count
				sumSQ+=SQQ
				count = 0
				sumSQt = 0
			SQ = sumSQ/np
			sq.append(SQ)
			sumSQ = 0
				
######## Until here we have genrated two lists: q = [..,..,...] and sq = [..,..,....] but data are mixed. We need to sort the data, remove the repetative data from q list and do average for sq values which are for equal q values.

	qnew = []
	sqnew = []
	qindex = []
	sumsq = 0
	lengthq = len(q)

	while lengthq != 0:
		a = min(q)
		b = q.count(a)
		qnew.append(a)
		for i in range(len(q)):
			if q[i]==a:
				sumsq+=sq[i]
				qindex.append(i)
				
		avesq = sumsq/b
		sqnew.append(avesq)
		sumsq = 0
		for i in range(len(qindex)):
			a = qindex[i]
			q.pop(a-i)
			sq.pop(a-i)
		lengthq = len(q)
		qindex = []


########### writing  data in sq.dat file for saving
 	
	for i in range(len(qnew)):
		a = qnew[i]
		b = sqnew[i]*(a**2)
		fdata.write('{:4.12f}   {:4.12f}\n '.format(float(a), float(b)))


font = {'fontname':'Times New Roman'}
plt.xlabel('q', **font, fontsize=15)
plt.ylabel('s(q)', **font, fontsize=15)
plt.title('Structure Factor for N = 60 (single polymer)', **font, fontsize=20)
##plt.suptitle('Z direction non_periodic and wall in the middle, phi = 0.5, monomer perchain = 60', **font, fontsize=15)
plt.plot(qnew,sqnew)
plt.show() 
