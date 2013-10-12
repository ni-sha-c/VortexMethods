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
        term3 = (1./(2.*math.pi*self.alph**2))*self.get_integral_gauss(t,ft)
        res = term1 + term2 + term3
        return res
     
    def get_integral_gauss(self,t,ft):
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
                
    def get_integral_simp(Xin):
     
        n = 1000
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
    
    def in_integral_simp(t_cap_all):
    
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
            res[count] = sci.quadrature(G_int,0.,ftcap,args=(t_cap,ft,t))[0]
            count = count + 1
        
        return res
         
    
    
