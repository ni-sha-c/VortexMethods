"""
    Driver module.
"""
import math
import numpy as np
from numpy.linalg import solve
from P import P_alph
from rhs import pop_b
alpha = 0.2e0
#X1 is X in the first iteration.
k = 4234.e-4
W = 8488.e-4
t = 0.e0 
delt = 0.008
X1 = np.array([ k, W, 1999.e-4, -5.e-4, -67.e-4, \
     5.e-4, 1.e-4, -0.e-4, 0.,  0.,  0.,  0.,  0.,  0., \
     0.,  0.,  0.,  0.,  0.,  0.])
#Number of fourier coefficients + 2 
#Also, equal to number of iterations 
N = 20     
A = np.zeros((N,N))
b = np.zeros((1,N))

st_time = 0
en_time = math.pi
dt = (en_time-st_time)/(N-1)

#Populating LHS and RHS once
for j in range(0,N):
    tj = t + j*dt
    X = X1
    #Tweak X for P's sake
    X = X1
    X = np.insert(X1,[2,2],[alpha,tj])
    P_2 = P_alph(X)
    for i in range(0,N):
        X[i] = X[i] + delt
        X = np.insert(X1,[2,2],[alpha,tj])
        P_1 = P_alph(X)
        Pxj = (P_1-P_2)/delt
        A[j][i] = Pxj
    
    b[j] = P_2 

for i in range(0,N):
    d = solve(A,b)
    X1 = X1 + d
    b = pop_b(X,alpha,t)
    
    