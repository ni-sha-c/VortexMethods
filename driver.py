
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
M = 2
A = np.zeros((N,N))
b = np.zeros((N,1))

st_time = 0
en_time = math.pi
dt = (en_time-st_time)/(N-1)

#Populating LHS and RHS once
print "populating N*N matrix..."
for j in range(0,N):
    print "************ j = ", j , "******************"
    tj = t + j*dt
    X = X1
    #Tweak X for P's sake
    X = np.insert(X1,[2,2],[alpha,tj])
    P_2 = P_alph(X)
    for i in range(0,N):
        print "i : ", i
        X = X1
        X[i] = X[i] + delt
        X = np.insert(X,[2,2],[alpha,tj])
        P_1 = P_alph(X)
        Pxj = (P_1-P_2)/delt
        A[j][i] = Pxj
    
    b[j] = P_2 
    
print "A is", A
print "b is", b

for i in range(0,M):
    d = solve(A,b)
    X1 = X1 + d
    X2 = np.insert(X1, [2,2], [alpha,t])
    b = pop_b(X2,alpha,t,N)
    
    