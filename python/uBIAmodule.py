# uBIA python module 

from itertools import product, combinations
import numpy as np
import math
import operator
from matplotlib import pyplot as plt
import functools
def nCr(n,r):
    f = math.factorial
    return f(n) / f(r) / f(n-r)
def list_powerset(lst):
    return functools.reduce(lambda result, x: result + [subset + [x] for subset in result], lst, [[]])
from collections import Counter
import time

#------------------------------------------------------------------------

def extract_patt(samples, hcutoff=0.05):
	"""This function extract all patterns that occur in samples, and patterns that do not occur but are above a threshold in prob"""
	#hcutoff=0.05

	M=len(samples) #number of samples
	T=np.shape(samples)[1] #dimension (time bins or number of neurons, etc)
	pi=np.sum(samples,axis=0)
	pi=[1.*i/M for i in pi] #empirical individual probabilities
	smax=np.amax(np.sum(samples,axis=1)) #max number of spikes in a sample

	sets=[]
	for xi in samples:
	    iia=[i for i,x in enumerate(xi) if x == 1]
	    sets+=list_powerset(iia)
	    
	dicset=dict()
	for x in sets:
	    tx=tuple(x)
	    if tx in dicset:
	        dicset[tx]+=1
	    else:
	        dicset[tx]=1

	# dicset: a tally list of marginal pattern that occur at least once ~ 10^5
	#print("patterns with n>0:     ",len(dicset))

	cc=hcutoff*2/(M+1) #cutoff for probability of a pattern that doesn't occur
	order_max = 2
	p0max=np.amax([np.prod(np.array(pi)[np.array(i)]) for i in combinations(range(T), order_max)])

	while (p0max > cc):
	    order_max += 1
	    p0max=np.amax([np.prod(np.array(pi)[np.array(i)]) for i in combinations(range(T), order_max)])

	order_max -= 1

	# now including patterns that do not occur, but may be relevant ~ 10^3
	for io in range(1,order_max+1):
	    for ss in combinations(range(T), io):
	        if np.prod(np.array(pi)[np.array(ss)])>cc and ss not in dicset:
	            dicset[ss]=0

	#print("plus patterns with n=0:",len(dicset))
	return dicset


#------------------------------------------------------------------------

def magn_field(dicset, M, pi, Nmax1 =1000, hcutoff1=0.05):
#	field_patt=dict()
	field_patt={}
	for pat in dicset:
		nn=dicset[pat]
		pp=np.prod(np.array(pi)[np.array(pat).astype(int)])
		hh=0.25*((nn - M*pp)**2 - M*pp*(1.-pp))
		if hh>hcutoff1:
			field_patt[pat]=[hh, nn, pp*M]

	#print("patt with h>hcutoff:   ",len(field_patt))
	sorted_patt = sorted(field_patt.items(), key=operator.itemgetter(1), reverse=True)[:Nmax1]

	return sorted_patt

#------------------------------------------------------------------------

def couplings(sortedpatt, M, pi, Nmax1):#------------------------------->>>>>>CHANGE HERE
	#Nmax1=len(sortedpatt)#------------------------------->>>>>>CHANGE HERE
	Jvw=np.zeros((Nmax1,Nmax1))

	for v in range(Nmax1-1):
		for w in range(v+1,Nmax1):
			tv=sortedpatt[v][0]
			nv=sortedpatt[v][1][1]
			pvM=sortedpatt[v][1][2]

			tw=sortedpatt[w][0]
			nw=sortedpatt[w][1][1]
			pwM=sortedpatt[w][1][2]

			Iv=set(tv) & set(tw) #intersection of patterns

			if Iv!= set([]):
			#        print "no vacio"
				Uv=set(tv) | set(tw) #union of patterns
				pin= np.prod(np.array(pi)[np.array(list(Iv)).astype(int)])
				pun= np.prod(np.array(pi)[np.array(list(Uv)).astype(int)])
				Jvw[v,w] += (1/16.)*pun*(1. - pin)*((M**2)*pun*(1.-pin)-2.*M*(nv - pvM)*(nw - pwM));

	# making it symmetric now
	Jvw=Jvw+Jvw.transpose()

	return Jvw

#------------------------------------------------------------------------

def marg_magn_h1ef(hh0, Jij, M, e_max =1., n_e =20):
	hef=np.sum(Jij,axis=0)
	Nmax1=len(Jij)
	delta_e= e_max/(M*n_e)

	#marginal magnetizations for the Nmax patterns, and for different levels of inverse regularization strength ee
	mi=np.zeros((Nmax1,n_e+1))
	#effective field (coming from interactions)
	h1ef=np.zeros((Nmax1,n_e+1))
	h1ef[:,0]=hh0

	#RECHECK
	for tt in range(n_e):
		ee=1.*delta_e*tt
    
		for v in range(Nmax1):
			if tt==0:
				h1ef[v,tt+1]=hh0[v]+2.*(ee+delta_e)*hef[v]
				mi[v,tt+1]=np.tanh((ee+delta_e)*h1ef[v, tt + 1])
			else:
				h1ef[v,tt+1]=hh0[v]+2.*(ee+delta_e)*hef[v]+1.*(ee+delta_e)*np.dot(Jij[v],mi[:,tt])
				mi[v,tt+1]=np.tanh((ee+delta_e)*h1ef[v, tt + 1])

	return mi, h1ef

#------------------------------------------------------------------------

def level_inverse_reg(hh0, hh1ef, n_e =20, cutoff =1.):
	"""The columns of h1ef correspond to the sum of the fields h0 plus the effective field for each level of inverse regularization 𝜖 (from 𝜖=0 to 𝜖=1/𝑀). We will increase 𝜖 up to the point where the effective field is of the order of h0 (taking into account our Taylor approximation), or up to 1/𝑀 (saddle point approximation for posterior), whatever happens first. Using L1 norm here."""
	nh0=np.linalg.norm(hh0,ord=1)
	list_nh1nh0=[np.linalg.norm(hh1ef[:,i]-hh0,ord=1)/nh0 for i in range(n_e+1)]

	neopt = 0
	#cutoff=1.

	nh1nh0=list_nh1nh0[neopt]
	while (neopt <=n_e and nh1nh0<cutoff):
		neopt += 1
		if neopt<=n_e:
			nh1nh0=list_nh1nh0[neopt]
	neopt -=1

	return neopt

#------------------------------------------------------------------------

def marg_magn_sorted(sortedpatt, mi, neopt):
	marg_patt={}
	NN=len(mi)#len(sortedpatt)

	for i in range(NN):
		pat=sortedpatt[i][0]
		hh=sortedpatt[i][1][0]
		nn=sortedpatt[i][1][1]
		ppM=sortedpatt[i][1][2]

		marg_patt[pat]=[mi[i,neopt],hh, nn, ppM]

	marg_sorted_patt = sorted(marg_patt.items(), key=operator.itemgetter(1), reverse=True)

	return marg_sorted_patt

#------------------------------------------------------------------------

def select_codewords(sortedpatt, mi, neopt):
	code_patt={}
	NN=len(mi)#len(sortedpatt)

	for i in range(NN):
		pat=sortedpatt[i][0]
		if 0 in pat:
			hh=sortedpatt[i][1][0]
			nn=sortedpatt[i][1][1]
			ppM=sortedpatt[i][1][2]

			code_patt[pat]=[mi[i,neopt],hh, nn, ppM]

	code_marg_sorted_patt = sorted(code_patt.items(), key=operator.itemgetter(1), reverse=True)

	return code_marg_sorted_patt

#------------------------------------------------------------------------

def select_codewords_overoccur(sortedpatt, mi, neopt):
	code_patt_over={}
	NN=len(mi)#len(sortedpatt)

	for i in range(NN):
		pat=sortedpatt[i][0]
		if 0 in pat:
			hh=sortedpatt[i][1][0]
			nn=sortedpatt[i][1][1]
			ppM=sortedpatt[i][1][2]
            
			if nn>ppM:
				code_patt_over[pat]=[mi[i,neopt],hh, nn, ppM]

	code_marg_sorted_patt = sorted(code_patt_over.items(), key=operator.itemgetter(1), reverse=True)

	return code_marg_sorted_patt

#------------------------------------------------------------------------
