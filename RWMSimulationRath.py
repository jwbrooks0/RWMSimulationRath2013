#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 10:48:25 2019

@author: john

References
----------
https://iopscience.iop.org/article/10.1088/0029-5515/53/7/073052
http://stacks.iop.org/0029-5515/45/i=4/a=010
http://stacks.iop.org/0741-3335/58/i=4/a=045001
"""

import numpy as np
import matplotlib.pyplot as plt
import johnspythonlibraries as jpl



### parameters
"""
specified in niko's paper 2013 
https://iopscience.iop.org/article/10.1088/0029-5515/53/7/073052
"""
c=0.261
c_f=0.8*np.exp(-1j*np.pi/10)
gamma_w=3.42/0.001 #3.42 ms^-1 # units of angular frequency I think...
s=1
alpha=2.2
L_c=91e-6
R_c=1
M_c=2.3e-6
m=3
r_w=0.16


### subfunctions
def AMatrix3D(Omega,gain,phase):
	A=np.zeros((3,3),dtype=complex)
	A[0,0]=(1-s+1j*alpha)/alpha*Omega
	A[0,1]=-Omega/(alpha*np.sqrt(c))
	A[0,2]=-c_f*Omega/alpha
	A[1,0]=gamma_w*np.sqrt(c)/(1-c)
	A[1,1]=-gamma_w/(1-c)
	A[1,2]=gamma_w*(1-c*c_f)/(1-c)
	A[2,0]= 1/(M_c*L_c)*gain*np.exp(1j*phase)*m/r_w*2*np.sqrt(c)/(1-c)
	A[2,1]=-1/(M_c*L_c)*gain*np.exp(1j*phase)*m/r_w*(1+c)/(1-c)
	A[2,2]=-R_c/L_c
	return A


### solve for eigenvalues for fixed gain
# setup parameter space
Omega=-8e3*np.pi*2
gain=np.array([0,1e-4,2e-4, 5e-4, 10e-4, 20e-4])*1e-4
phase=np.arange(0,np.pi*2,0.1)

# initialize
gamma1=np.zeros((len(phase),len(gain)),dtype=complex)
gamma2=np.zeros((len(phase),len(gain)),dtype=complex)
gamma3=np.zeros((len(phase),len(gain)),dtype=complex)

# solve
for i in range(len(gain)):
	for j in range(len(phase)):
		# print('%.1f, %.1e, %.3f'%(Omega,gain[i],phase[j]))
		A=AMatrix3D(Omega,gain[i],phase[j])
		gamma=np.linalg.eigvals(A)
		gamma1[j,i],gamma2[j,i],gamma3[j,i]=gamma
		

### solve for eigenvalues for fixed phase 
# setup parameter space
gainB=np.linspace(gain[0],gain[-1],100)
phaseB=np.array([0,0.5,1,1.5,-100*np.pi/180])*np.pi

# initialize
gamma1B=np.zeros((len(gainB),len(phaseB)),dtype=complex)
gamma2B=np.zeros((len(gainB),len(phaseB)),dtype=complex)
gamma3B=np.zeros((len(gainB),len(phaseB)),dtype=complex)

# solve
for i in range(len(phaseB)):
	for j in range(len(gainB)):
		print('%d, %d, %.1f, %.1e, %.3f'%(i,j,Omega,gainB[j],phaseB[i]))
		A=AMatrix3D(Omega,gainB[j],phaseB[i])
		gamma=np.linalg.eigvals(A)
		gamma1B[j,i],gamma2B[j,i],gamma3B[j,i]=gamma

		
### plotting
		
# fixed gain root-locus plots
fig,ax=plt.subplots()
for i in range(len(gain)):
	if gain[i]==0:
		marker='x'
	else:
		marker='.'
	c='k'
	label='' #r'$\gamma_1, g=%.1f$'%gain[i]
	ax.plot(np.imag(gamma1[:,i])*1e-3/(np.pi*2),np.real(gamma1[:,i])*1e-3,label=label,marker=marker,linestyle='',color=c)
	ax.plot(np.imag(gamma2[:,i])*1e-3/(np.pi*2),np.real(gamma2[:,i])*1e-3,label=label,marker=marker,linestyle='',color=c)
	ax.plot(np.imag(gamma3[:,i])*1e-3/(np.pi*2),np.real(gamma3[:,i])*1e-3,label=label,marker=marker,linestyle='',color=c)
	
# fixed phase root locus plots
for i in range(len(phaseB)):
	if phaseB[i]==0 or phaseB[i]==np.pi:
		c='r'
	elif phaseB[i]==0.5*np.pi or phaseB[i]==np.pi*1.5:
		c='b'
	else:
		c='g'
	marker='.'
	ax.plot(np.imag(gamma1B[:,i])*1e-3/(np.pi*2),np.real(gamma1B[:,i])*1e-3,label=label,marker=marker,linestyle='',color=c)
	ax.plot(np.imag(gamma2B[:,i])*1e-3/(np.pi*2),np.real(gamma2B[:,i])*1e-3,label=label,marker=marker,linestyle='',color=c)
	ax.plot(np.imag(gamma3B[:,i])*1e-3/(np.pi*2),np.real(gamma3B[:,i])*1e-3,label=label,marker=marker,linestyle='',color=c)
jpl.plot.finalizeSubplot(	ax,
						xlabel='Frequency (kHz)',
						ylabel='Growth rate (1/ms)',
						title='gains = %s'%''.join(['%.1e, '%i for i in gain]))
jpl.plot.finalizeFigure(fig,figSize=[6,3])
fig.savefig('result.png')

	