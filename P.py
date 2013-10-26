import numpy as np
import scipy.integrate as sci
from f import f
from G import G_int
from a0a1 import a0a1
import math

class X_alph:
    """ 
    P_alph(X)@t = -X_1 - 0.5*X_2(1+fcost)^2
    + 0.5/(pi*alph**2)*int_0^2*pi dt_cap *
    int_0^f_cap G(1+fcost,1+s_cap cos t_cap,
    fsint - s_cap sin t_cap)s_cap ds_cap 

    """
    def __init__(self,alpha,X):
        #Setting [k,W,a2..a17]
        self.X = np.copy(X)
        self.alph = alpha
        [self.a0, self.a1] = self.get_a0a1()
        
    def get_a0a1(self):
        #Calculating a0 and a1
        N = np.size(self.X)
        X_A = np.insert(self.X[2:N],0,self.alph)
        [a,b] = a0a1(X_A)
        return np.array([a,b])

    def P_alph(self,t):
        X_1 = self.X[0]
        X_2 = self.X[1]
        Xff = np.copy(self.X)
        Xff[0] = self.a0
        Xff[1] = self.a1
        ft = f(np.insert(Xff,0,t))
        term1 = -X_1
        term2 = -0.5e0*X_2*(1+ft*math.cos(t))**2
        term3 = (1./(2.*math.pi*self.alph**2))*self.get_integral_simp(ft,t)
        res = term1 + term2 + term3
        return res
     
    def get_integral_gauss(self,ft,t):
        res = sci.quadrature(self.get_in_gauss,0.,2*math.pi,args=(ft,t), \
        vec_func=False)[0]
        return res
        
    def get_in_gauss(self,tcap,ft,t):
        X_f = np.insert(self.X,0,tcap)
        X_f[1] = self.a0
        X_f[2] = self.a1
        ftcap = f(X_f)
        #Splitting the interval
        #n = 5
        #res = np.zeros(n)
        #for i in range(0,n):
            #start_scap = i/n*ftcap
            #end_scap = (i+1)/n*ftcap       
            #res[i] = sci.quadrature(G_int,start_scap,end_scap, \
            #args=(tcap,ft,t),vec_func=False)[0]
        #return np.sum(res) 
      
        #Without splitting
        res = sci.quadrature(G_int,0.,ftcap, \
            args=(tcap,ft,t),vec_func=False)[0]
        return res
        
        
    def get_integral_trap(self,ft,t):
        st_tcap = 0.e0
        end_tcap = 2.e0*math.pi
        n = 500
        tcap_all = np.linspace(st_tcap, end_tcap, n)
        i = 0
        intgd = np.zeros(n)
        for tcap in tcap_all:
            intgd[i] = self.get_in_gauss(tcap,ft,t) 
            i = i + 1
        res = sci.trapz(intgd,tcap_all)
        return res   
    
    def pop_b(self,Nt):
        t = 0
        dt = math.pi/(Nt-1)      
        b = np.zeros(Nt)
        for i in range(0,Nt):
            ti = t + i*dt
            b[i] = -self.P_alph(ti)
        
        return np.transpose(b)
                
    def get_integral_simp(self,ft,t):
        st_tcap = 0.e0
        end_tcap = 2.e0*math.pi
        n = 500
        tcap_all = np.linspace(st_tcap, end_tcap, n)
        i = 0
        intgd = np.zeros(n)
        for tcap in tcap_all:
            intgd[i] = self.get_in_gauss(tcap,ft,t) 
            i = i + 1
        res = sci.simps(intgd,tcap_all)
        return res  
                 
    
    
    
    
