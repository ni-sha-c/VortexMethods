"""
  Populates the RHS for every iteration.

"""
from P import P_alph
import math
import numpy as np
def pop_b(X, alpha, t, N):
    
    dt = math.pi/(N-1)
    t = 0
    b = np.zeros((1,N))
    for i in range(0,N):
        ti = t + i*dt
        
        X1 = np.insert(X,[2,2],[alpha,ti])
        b[i] = -P_alph(X1)
        
    return b
        
        
        
        
        
    
    
    