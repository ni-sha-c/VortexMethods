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
    
    
    return res