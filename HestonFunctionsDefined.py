import numpy as np
from numpy.random import standard_normal
import math as math
import cmath
 
u1=0.5
u2=-0.5
### part of the heston option pricing is here 
def  b1(Parameters,psi):
    return  Parameters[0]+   Parameters[4]-Parameters[2]*Parameters[3]
def  b2(Parameters,psi):
    return Parameters[0] + Parameters[4]


#  also define ksi, for ease of calculation
def ksi1(Parameters,psi):
    return b1(Parameters,psi) - Parameters[3] *Parameters[2]*psi*1j

def ksi2(Parameters,psi):
    return b2(Parameters,psi) - Parameters[3] *Parameters[2]*psi*1j
def d1(Parameters, psi):
    
    return cmath.sqrt (  ( ksi1(Parameters,psi) )**2 
                       - ((Parameters[2])**2)*(2*u1*psi*1j-psi**2   )
                      )
def d2(Parameters, psi):
    return cmath.sqrt (  (   ksi2(Parameters,psi))**2 
                       - ((Parameters[2])**2)*(2*u2*psi*1j-psi**2   )
                      )
    
def g1(Parameters, psi ):
    return ( 
            (ksi1(Parameters,psi) +   d1(Parameters, psi))/
            ( ksi1(Parameters,psi) -   d1(Parameters, psi)  )
            )

def g2(Parameters, psi ):
    return ( 
            (ksi2(Parameters,psi) +   d2(Parameters, psi))/
            ( ksi2(Parameters,psi) -   d2(Parameters, psi)  )
            )
def C1(Parameters,psi, tau, r_for,r_dom ):
   return ( 
           (r_dom- r_for)*psi*1j*tau
           + ( Parameters[0]*Parameters[1]/(Parameters[2]**2))*(  
                ksi1(Parameters,psi)+d1(Parameters, psi)*tau
                -2*cmath.log( (1-g1(Parameters, psi)*cmath.exp(d1(Parameters, psi)*tau) )/
                             (1-g1(Parameters, psi ))))
           )
           
def C2(Parameters,psi, tau, r_for,r_dom ):
   return ( 
           (r_dom- r_for)*psi*1j*tau
           + ( Parameters[0]*Parameters[1]/(Parameters[2]**2))*(  
                ksi2(Parameters,psi)+d2(Parameters, psi)*tau
                -2*cmath.log( (1-g2(Parameters, psi)*cmath.exp(d2(Parameters, psi)*tau) )/
                             (1-g2(Parameters, psi ))))
           ) 
def D1(Parameters,psi, tau):
    return ( (ksi1(Parameters,psi) +   d1(Parameters, psi))/(Parameters[2]**2) )*(
            (1-cmath.exp(d1(Parameters, psi)*tau) )/(
            (1-g1(Parameters, psi)*cmath.exp(d1(Parameters, psi)*tau) )))
def D2(Parameters,psi, tau):
    return ( (ksi2(Parameters,psi) +   d2(Parameters, psi))/(Parameters[2]**2) )*(
            (1-cmath.exp(d2(Parameters, psi)*tau) )/(
            (1-g2(Parameters, psi)*cmath.exp(d2(Parameters, psi)*tau) )))

#now we move on to characteristic function and gaussiche integration

def f1(Parameters,psi, tau, FX0,r_for,r_dom):
    return cmath.exp ( 
                         C1(Parameters,psi, tau, r_for,r_dom )
                      +  D1(Parameters,psi, tau)*Parameters[5]
                      + 1j*psi*cmath.log(FX0)
                      )        
def f2(Parameters,psi, tau, FX0,r_for,r_dom):
    return cmath.exp ( 
                         C2(Parameters,psi, tau, r_for,r_dom )
                      +  D2(Parameters,psi, tau)*Parameters[5]
                      + 1j*psi*cmath.log(FX0)
                      )    
# gaussian integration       
def integral1(Parameters, K,tau, FX0, r_for,r_dom, M,nodes):
    x, w = np.polynomial.legendre.leggauss(nodes)  
    psi= (x[0]+1)*0.5*M
    value = w[0]*(cmath.exp(-1j*psi*math.log(K))/(1j*psi))*( 
            f1(Parameters,psi, tau, FX0,r_for,r_dom) )
    for j in range(1, nodes): #   number of nodes
        psi = (x[j] + 1)*0.5*M
        value = value + w[j]*(cmath.exp(-1j*psi*math.log(K))/(1j*psi))*( 
            f1(Parameters,psi, tau, FX0,r_for,r_dom) )
    value = value*0.5*M
    return 0.5 + (1/math.pi)*value.real

def integral2(Parameters, K,tau, FX0, r_for,r_dom, M,nodes):
    x, w = np.polynomial.legendre.leggauss(nodes)  
    psi= (x[0]+1)*0.5*M
    value = w[0]*(cmath.exp(-1j*psi*math.log(K))/(1j*psi))*( 
            f2(Parameters,psi, tau, FX0,r_for,r_dom) )
    for j in range(1, nodes): #   number of nodes
        psi = (x[j] + 1)*0.5*M
        value = value + w[j]*(cmath.exp(-1j*psi*math.log(K))/(1j*psi))*( 
            f2(Parameters,psi, tau, FX0,r_for,r_dom) )
    value = value*0.5*M
    return    0.5 + (1/math.pi)*value.real

#   define heston price 

def HestonPrice(Parameters, K, tau, FX0, r_for, r_dom, M, nodes,phi):
    return   phi *(math.exp( -r_for*tau)*FX0*
                   ( (1-phi)/2 +phi* integral1(Parameters, K,tau, FX0, r_for,r_dom, M,nodes))
                   -K*math.exp(-r_dom*tau)*
                   ( (1-phi)/2 + phi* integral2(Parameters, K,tau, FX0, r_for,r_dom, M,nodes))  )     


################ heston simulation part is here : euler disc. is used 
     
def Heston_Sim(Parameters, FX0,K, tau,r_dom, r_for , NReps, NSteps,phi)  :
    dt= 1/100 
    FX= np.zeros(NSteps + 1)
    v= np.zeros(NSteps + 1)
    #PayOff= 0
    for i in range(NReps):
        v[0]=Parameters[5] # initial value of variane 
        FX[0]= FX0   # initial exchange rate
        for j in range(NSteps):
            W1 =  standard_normal(NSteps) # for later on
            W2 = Parameters[3]*W1 + math.sqrt(1 - Parameters[3]*Parameters[3])*standard_normal(NSteps)
            # Parameterss5 is rho
            FX[j + 1] = FX[j] + (r_dom-r_for)* FX[j] * dt +  math.sqrt(v[j]) *FX[j]* math.sqrt(dt)* W1[j]
            #   take the max of variance, only to be sure that it is not negative,  
            v[j + 1] =    v[j] +  Parameters[0]*(Parameters[1]- v[j])*dt 
            + Parameters[2]*math.sqrt(v[j])*math.sqrt(dt)*W2[j]
                         
        return FX
    #so that we generated heston fx process
    
 def Heston_Sim_Option(Parameters, FX0,K, tau,r_dom, r_for , NSim, NSteps,phi):
    Payoff = 0
    P1 = np.array(Heston_Sim(Parameters, FX0,K, tau,r_dom, r_for , NSim, NSteps,phi))
    P2 = P1 -K
    Payoff = sum(x for x in P2 if x>0) #   just took the positives of our array
    return  math.exp(-r_dom*tau )*(Payoff/NSim)
    
      
    
    
 

         
            
    














