"""
    Function to calculate a0 and a1
    given X
"""

import numpy as np

def a0a1(X):
    N = np.size(X)
    ot = X[2:N:2]
    a1 = -sum(ot)
    sq = 0.5e0*np.sum(np.square(X[1:N]))
    alph = X[0]
    a0 = np.sqrt(alph**2 - sq - 0.5e0*(a1**2))
    x = np.array([a0,a1])
    return x
    