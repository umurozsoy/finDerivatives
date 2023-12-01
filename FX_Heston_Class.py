 """ this class calculates heston price on currency options. By parameters, the order
is [Kappa, Theta, Sigma, Rho, Lambda ,v0] . Kappa is mean rev. rate
theta = the long run variance
Sigma = volatility of volatility
Rho correlation and lambda is set to zero but i included for some other reasons

for numerical integration i use gaussian and i set M=100 and nodes=100
phi 1 for calls -1 for puts

I use the formulae from "FX smile in the Heston model" by Uwe Wystup et al.

to be able to have a more compact class, i define my functions in an other file

"""
# =============================================================================
# t
 
# =============================================================================
 
import HestonFunctionsDefined as hfd
class FX_Heston_Option_Pricing(object):
    def __init__(self,Parameters,FX0, K,tau, r_dom, r_for,M,nodes,   phi ):
        self.Parameters=Parameters
        self.FX0=float(FX0)
        self.K=K
        self.tau=tau
        self.r_dom=r_dom 
        self.r_for=r_for
        self.M=M
        self.nodes=nodes
        self.phi=phi
    def value(self):
           return hfd.HestonPrice(self.Parameters, self.K, self.tau, self.FX0, self.r_for, self.r_dom, self.M, self.nodes,self.phi)


    
    
 