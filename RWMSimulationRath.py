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
nikos2013Parameters =	{
	"c": 0.261,
	"c_f": 0.8*np.exp(-1j*np.pi/10),
	"gamma_w": 3.42/0.001, #3.42 ms^-1 # units of angular frequency I think...,
	"s":1,
	"alpha":2.2,
	"L_c":91e-6,
	"R_c":1,
	"M_c":2.3e-6,
	"m":3,
	"r_w":0.16,
	"Omega":-8e3*np.pi*2
}

#parameters=nikos2013Parameters


### subfunctions
def Niko2013Matrix(gain,phase,params,verbose=False):
	"""
	specified in niko's paper 2013 
	https://iopscience.iop.org/article/10.1088/0029-5515/53/7/073052
	"""
	A=np.zeros((3,3),dtype=complex)
	A[0,0]=(1-params['s']+1j*params['alpha'])/params['alpha']*params['Omega']
	A[0,1]=-params['Omega']/(params['alpha']*np.sqrt(params['c']))
	A[0,2]=-params['c_f']*params['Omega']/params['alpha']
	A[1,0]=params['gamma_w']*np.sqrt(params['c'])/(1-params['c'])
	A[1,1]=-params['gamma_w']/(1-params['c'])
	A[1,2]=params['gamma_w']*(1-params['c']*params['c_f'])/(1-params['c'])
	A[2,0]= 1/(params['M_c']*params['L_c'])*gain*np.exp(1j*phase)*params['m']/params['r_w']*2*np.sqrt(params['c'])/(1-params['c'])
	A[2,1]=-1/(params['M_c']*params['L_c'])*gain*np.exp(1j*phase)*params['m']/params['r_w']*(1+params['c'])/(1-params['c'])
	A[2,2]=-params['R_c']/params['L_c']
	
	if verbose: print(A)
	
	return A


def solveAForEigenvalues(AFunc,gain,phase,parameters,verbose=False):
	""" Solve A matrix for eigenvalues for a range of phase and gain parameters """
		
	# initialize
	gamma1=np.zeros((len(phase),len(gain)),dtype=complex)
	gamma2=np.zeros((len(phase),len(gain)),dtype=complex)
	gamma3=np.zeros((len(phase),len(gain)),dtype=complex)
	
	# solve
	for i in range(len(gain)):
		for j in range(len(phase)):
			if verbose: print('%.1f, %.1e, %.3f'%(parameters['Omega'],gain[i],phase[j]))
			
			A=AFunc(gain[i],phase[j],parameters,verbose=verbose)
			gamma=np.linalg.eigvals(A)
			gamma1[j,i],gamma2[j,i],gamma3[j,i]=gamma
			
	return gamma1, gamma2, gamma3
			

### solve for eigenvalues for fixed gain
gain=np.array([0,2e-4, 5e-4, 10e-4, 20e-4])*1e-4
phase=np.arange(0,np.pi*2,0.1)		
gamma1,gamma2,gamma3=solveAForEigenvalues(Niko2013Matrix,gain,phase,nikos2013Parameters,verbose=False)
		

### solve for eigenvalues for fixed phase 
gainC=np.linspace(gain[0],gain[-1],50)
phaseC=np.array([0,0.5,1,1.5,-100*np.pi/180])*np.pi
gamma1C,gamma2C,gamma3C=solveAForEigenvalues(Niko2013Matrix,gainC,phaseC,nikos2013Parameters,verbose=False)

		
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
for i in range(len(phaseC)):
	if phaseC[i]==0 or phaseC[i]==np.pi:
		c='r'
	elif phaseC[i]==0.5*np.pi or phaseC[i]==np.pi*1.5:
		c='b'
	else:
		c='g'
	marker='.'
	ax.plot(np.imag(gamma1C[i,:])*1e-3/(np.pi*2),np.real(gamma1C[i,:])*1e-3,label=label,marker=marker,linestyle='',color=c)
	ax.plot(np.imag(gamma2C[i,:])*1e-3/(np.pi*2),np.real(gamma2C[i,:])*1e-3,label=label,marker=marker,linestyle='',color=c)
	ax.plot(np.imag(gamma3C[i,:])*1e-3/(np.pi*2),np.real(gamma3C[i,:])*1e-3,label=label,marker=marker,linestyle='',color=c)
jpl.plot.finalizeSubplot(	ax,
						xlabel='Frequency (kHz)',
						ylabel='Growth rate (rad/ms)',
						title='gains = %s'%''.join(['%.1e, '%i for i in gain]))
jpl.plot.finalizeFigure(fig,figSize=[6,3])
fig.savefig('result.png')

	