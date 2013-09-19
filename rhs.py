"""
  Populates the RHS for every iteration.

"""
from P import P_alph
import math
import numpy as np
def pop_b(X_all):
    
    N = np.size(X_all)
    dt = math.pi/(N-1)
    t = 0
    b = np.zeros((1,N))
    for i in range(0,N):
        ti = t + i*dt
        
        X_all[3] = ti
        b[i] = -P_alph(X_all)
        
    return b
        
        
        
        
        
    
    
    