# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 11:41:12 2019

@author: aawasthi
"""
import dimod
import sympy as sym
import numpy as np
from dwave.system.samplers import DWaveSampler
from dwave.system.composites import EmbeddingComposite, AutoEmbeddingComposite
from pyqubo import Binary


n1 = 7
n2 = 6
N = n1*n2
precision = max(len(bin(n1)),len(bin(n2)))-2


X = 0
Y = 0
for i in range(precision):
    X += 2.0**i*Binary('x'+str(i))
    Y += 2.0**i*Binary('y'+str(i))    

Pol = (N- X*Y)**2
model = Pol.compile((N**2) * 100).to_dimod_bqm()

tuples = list(model.linear.items())        
max_val = max(tuples, key=lambda x:x[1])[1]

# In[]#
def getVal(bestSol,precision):
    x_val = ''
    y_val = ''
    for i in range(precision):
        x_val += str(bestSol['x'+str(i)])
        y_val += str(bestSol['y'+str(i)])
    x_val = x_val[::-1]
    y_val = y_val[::-1]
    return x_val, y_val


def getEnergy(x1,y1,precision,N):
    X1 = int(x1,2)
    Y1 = int(y1,2)
    
    Pol1 = (N- (X1)*(Y1))**2
    return Pol1



# In[] #################

sampler = EmbeddingComposite(DWaveSampler())
print('Dwave_Energy \t Actual Energy \t Chain_Break_Frq')
for iter in range(10):    
    response = sampler.sample(model, chain_strength=max_val, 
                          num_reads=10, annealing_time=20, return_embedding=True)
    datum=next(response.data(['sample', 'energy', 'num_occurrences', 'chain_break_fraction']))
    chainBreak = str(100 * datum.chain_break_fraction)            
    for solution in response.data():
        x_val, y_val = getVal(solution.sample, precision)
        E = getEnergy(x_val,y_val,precision,N)
#        print(solution.sample.values(), int(solution.energy),int(x_val,2), int(y_val,2), E, c)
        print(int(solution.energy),'\t\t', E, '\t\t', chainBreak)
# In[]

# In[]

#x_val = ''
#y_val = ''
##print(bestSol)
#for i in range(precision):
#    x_val += str(bestSol[i])
#    print(y_val)
#    y_val += str(bestSol[i+precision])
#x_val = x_val[::-1]
#y_val = y_val[::-1]
#print(x_val)
#int(x_val,2), int(y_val,2)