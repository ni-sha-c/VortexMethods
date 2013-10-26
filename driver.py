"""
    Driver module.
"""
import math
import numpy as np
from numpy.linalg import solve
from P import X_alph
#from rhs import pop_b

alpha = 0.8e0
#X1 is X in the first iteration.
k = 0.047
W = 0.4428
t = 0.e0
delt = 0.008
X1 = np.array([k, W, -0.1761, 0.0504, 0.0178, \
            -0.0184,  0.0035, 0.0035, -0.0027, 0.0003,
            0.0008, -0.0005, 0.0000, 0.0002, -0.0001,
            0.0000, 0.0000, 0.0000])
#Number of fourier coefficients 

N = np.size(X1)     
M = 10
A = np.zeros((N,N))
b = np.zeros(N)

st_time = 0
en_time = math.pi
dt = (en_time-st_time)/(N-1)

#Forming an X_alph object
x = X_alph(alpha,X1)

#Populating LHS and RHS once
print "populating N*N matrix..."
for j in range(0,N):
    print "************ t = ", j , "******************"
    #Ensuring the delt is gone before every iteration
    x.X = np.copy(X1)
    tj = t + j*dt
    P_2 = x.P_alph(tj)
    for i in range(0,N):
        print "j : ", i
        x.X[i] = x.X[i] + delt
        P_1 = x.P_alph(tj)
        Pxj = (P_1-P_2)/delt
        A[j][i] = Pxj
    
    b[j] = -P_2 
    
b = np.transpose(b)
print "Solving to get better approximations..."
for i in range(0,M):
    print "n = ", i
    d = solve(A,b)
    x.X = np.copy(X1 + np.transpose(d))
    b = x.pop_b(N)
    X1 = np.copy(x.X)
   
#Getting a0 and a1 
[a0,a1] = x.get_a0a1()
X2 = np.insert(X1,[2,2],[a0,a1])   
print "X after", M, "iteration(s):", X2
    