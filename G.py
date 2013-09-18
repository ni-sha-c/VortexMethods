
"""
    Calculates G(r, r^, z-z^)
    
"""
import scipy.special as scsp
import math

def G(x, t_cap, f, t, k):
    
    r = 1. + f*math.cos(t)
    rcap = 1 + x*math.cos(t_cap)
    zmzcap = f*math.sin(t) - x*math.sin(t_cap)
    term1 = ((r+rcap)**2 + zmzcap**2)**0.5e0
    term1 = term1*rcap
    term2 = (1.e0 - k**2)*scsp.ellipk(k)
    term2 = term2 - scsp.ellipe(k)
    return term1*term2
    
