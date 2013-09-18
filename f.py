""" 
    s = f(t,alpha)
    
"""
import math
import numpy as np
def f(X):
    
    t = X[0]
    N = np.size(X) - 1
    
    n = np.linspace(0,N-1,N)
    
    res = np.cos(t*n)
    
    res1 = res*(X[1:N+1])
    
    res2 = np.sum(res1)
    
    return res2