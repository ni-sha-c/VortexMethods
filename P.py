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
        ft = f(X[2:N])
        alph = X[2]
        t = X[3]   
        term1 = -X_1
        term2 = -0.5e0*X_2*(1+ft*math.cos(t))**2
        Xin = np.insert(X,0,ft)
        term3 = (1./(2.*math.pi*alph**2))*get_integral(Xin)
    
        res = term1 + term2 + term3
        return res
     
 
def get_integral(Xin):
     
        n = 50
        a = 0.e0
        b = 2*math.pi
     
        dt_cap = (b-a)/(n-1)
     
        t_cap = np.linspace(a,b,n)
        t_cap_mid = np.linspace(a+dt_cap/2,b-dt_cap/2,n-1)  
     
        #Number of fourier coefficients
        N = np.size(Xin) - 5
        
        Xin = np.insert(Xin, 0, N)
        arg_in_end = np.concatenate((Xin, t_cap)) 
        arg_in_mid = np.concatenate((Xin,t_cap_mid))
        
        H_j = in_integral(arg_in_end)
        H_mid = in_integral(arg_in_mid)
      
        int_simpson = (dt_cap/6.) * (2.*sum(H_j) - (H_j[0] + H_j[-1]) \
                                    + 4.*sum(H_mid))
    
        return int_simpson            
    
def in_integral(t_cap_all):
    
        count = 0
        N = t_cap_all[0]
        ft = t_cap_all[1]
        k = t_cap_all[2]
        t = t_cap_all[5]
        alph = t_cap_all[4]

        foco = np.copy(t_cap_all[6:N+6])
        t_cap_all = np.copy(t_cap_all[N+6:np.size(t_cap_all)])    
        res = np.zeros(shape=t_cap_all.size)
    
        for t_cap in t_cap_all:
            foco1 = np.insert(foco,[0,0],[alph,t_cap])
            ftcap = f(foco1)
            res[count] = sci.quadrature(G,0.,ftcap,args=(t_cap,ft,t))[0]
            count = count + 1
        
        return res
         
    
    
