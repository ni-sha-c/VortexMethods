
"""
    Driver module.
"""
import math
import numpy as np
from numpy.linalg import solve
from P import P_alph
from rhs import pop_b
alpha = 1.0e0
#X1 is X in the first iteration.
k = 182.e-04
W = 3703.e-04
t = 0.e0 
delt = 0.008
X1 = np.array([ k, W, -0.27e0, \
     0.0907e0, 0.0318, -0.04e0, 0.0099e0,  0.0091,  -0.0092,  0.0013,  0.0032,  -0.0024e0, \
     0.0e0,  0.0011,  -0.0007,  0.0,  0.0003,  -0.0002])
#Number of fourier coefficients + 2 
#Also, equal to number of iterations 
N = np.size(X1)     
M = 2
A = np.zeros((N,N))
b = np.zeros((N,1))

st_time = 0
en_time = math.pi
dt = (en_time-st_time)/(N-1)

#Populating LHS and RHS once
print "populating N*N matrix..."
for j in range(0,N):
    print "************ t = ", j , "******************"
    tj = t + j*dt
    #Tweak X for P's sake
    X = np.copy(X1)
    X = np.insert(X1,[2,2],[alpha,tj])
    P_2 = P_alph(X)
    for i in range(0,N):
        print "j : ", i
        X = np.copy(X1)
        X[i] = X[i] + delt
        X = np.insert(X,[2,2],[alpha,tj])
        P_1 = P_alph(X)
        Pxj = (P_1-P_2)/delt
        A[j][i] = Pxj
    
    b[j] = P_2 
    
print "A is", A
print "b is", b
print "Solving to get better approximations..."
for i in range(0,M):
    print "n = ", i
    d = solve(A,b)
    X1 = X1 + np.transpose(d)[0]
    X2 = np.insert(X1, [2,2], [alpha,t])
    b = pop_b(X2)
    

    