# -*- coding: utf-8 -*-
"""


 
""" # EUR/USD = EUR is foreign currency , 
#import numpy as np
from math import log, sqrt, exp
from scipy import stats
''' FX = rate at the time, tau = time to expiration, r_for =foreign int rate, dom=domestic interest rate
 sigma = annualized volatility of the rate, phi is 1 for calls ,-1 for puts. only for vanilla fx options'''
class GarmanKohlhagen(object):
    def __init__(self, FX,K,tau,r_for,r_dom,sigma,phi):
        self.FX=float(FX)
        self.K=K
        self.tau=tau
        self.r_for=r_for
        self.r_dom=r_dom
        self.sigma=sigma
        self.phi=phi
        # phi 1 for calls , -1 for puts 
    def  value(self):
        d1=( log(self.FX/self.K ) +(self.r_dom -self.r_for+0.5*self.sigma**2 )*self.tau    )/(self.sigma*sqrt(self.tau  ))
        d2=( log(self.FX/self.K ) +(self.r_dom -self.r_for-0.5*self.sigma**2 )*self.tau    )/(self.sigma*sqrt(self.tau  ))
        value= self.phi*(self.FX*exp(-self.r_for*self.tau)*stats.norm.cdf(self.phi*d1,0.0,1.0) -self.K*exp(-self.r_dom*self.tau)*stats.norm.cdf(self.phi*d2,0.0,1.0        ))
        return  value
    def delta(self):
         d1=( log(self.FX/self.K ) +(self.r_dom -self.r_for+0.5*self.sigma**2 )*self.tau    )/(self.sigma*sqrt(self.tau  ))
         delta=self.phi*exp(-self.r_for*self.tau)*stats.norm.cdf(self.phi*d1,0.0,1.0)
         return delta
    def gamma(self):
         d1=( log(self.FX/self.K ) +(self.r_dom -self.r_for+0.5*self.sigma**2 )*self.tau    )/(self.sigma*sqrt(self.tau  ))
         gamma=exp(-self.r_for*self.tau)*( (stats.norm.pdf(d1))/(self.FX*self.sigma*sqrt(self.tau)))       
         return gamma
    def vega(self):
        d1=( log(self.FX/self.K ) +(self.r_dom -self.r_for+0.5*self.sigma**2 )*self.tau  )/(self.sigma*sqrt(self.tau  ))
        vega= self.FX* exp(-self.r_for*self.tau)*sqrt(self.tau)*stats.norm.pdf(d1)
        return vega 
    
        