"""
   Del P_alph / Del X_j 
   
"""
import numpy as np
import scipy.integrate as spi   
alph = 0.2e0
X  = [ 4234.e-4, 8488.e-4, 1999.e-4, -5.e-4, -67.e-4, \
     5.e-4, 1.e-4, -0.e-4 ]
zer = [0.e0]*12
X = X + zer
k = X[0]
W = X[1]
res = spi.quadrature(G,0.e0,1.e0,\
                    (1.e0, )