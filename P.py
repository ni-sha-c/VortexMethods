""" 
    P_alph(X)@t = -X_1 + 0.5*X_2(1+fcost)^2
    + 0.5/(pi*alph**2)*int_0^2*pi dt_cap *
    int_0^f_cap G(1+fcost,1+s_cap cos t_cap,
    fsint - s_cap sin t_cap)s_cap ds_cap 

"""

import numpy as np
import scipy.integrate as sci
from f import f
from G import G
import math

def P_alph(X):
    X_1 = X[0]
    X_2 = X[1]
    N = np.size(X)
    fatt = f(X[3:N])
    alph = X[2]
    t = X[3]   
    term1 = -X_1
    term2 = 0.5e0*X_2*(1+fatt*math.cos(t))**2
    term3 = (1./(2.*math.pi*alph**2))*get_integral(fatt, t, X_1)
    
    res = term1 + term2 + term3
    return res
     
 
def get_integral(fatt, t, k):
     
    n = 1000
    a = 0.e0
    b = 2*math.pi
     
    dt_cap = (b-a)/(n-1)
     
    t_cap = np.linspace(a,b,n)
    t_cap_mid = np.linspace(a+dt_cap/2,b-dt_cap/2,n-1)  
     
    H_j = in_integral(t_cap, t, fatt, k)
    H_mid = in_integral(t_cap_mid, t, fatt, k)
      
    int_simpson = (dt_cap/6.) * (2.*sum(H_j) - (H_j[0] + H_j[-1]) \
                                    + 4.*sum(H_mid))
    
    return int_simpson            
    
def in_integral(t_cap_all, t, fatt, k):
    
    count = 0
    res = np.zeros(shape=t_cap_all.size)
    
    for t_cap in t_cap_all:
        res[count] = sci.quadrature(G,0.,fatt,args=(t_cap,fatt,t,k))[0]
        count = count + 1
        
    return res
         
    
    