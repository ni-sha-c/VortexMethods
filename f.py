""" 
    s = f(t,alpha)
    
"""
import numpy as np
from a0a1 import a0a1
def f(X):
    
    alph = X[0]
    t = X[1]
    #Number of fourier coefficients including
    #a0 and a1
    N = np.size(X) 
    #Finding a0 and a1
    X1 = np.insert(X,2,alph)
    [a0,a1] = a0a1(X1[2:N+1])
    
    X = np.insert(X, [2,2],[a0,a1])
    
    n = np.linspace(0,N-1,N)
   
    res = np.cos(t*n)
    
    res1 = res*(X[2:N+2])
    
    res2 = np.sum(res1)
    
    return res2