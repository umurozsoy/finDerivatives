import HestonFunctionsDefined as hfd

class FX_Heston_Simulation(object):
    def __init__(self, Parameters, FX0,   tau, r_dom, r_for , NReps, NSteps):
        self.Parameters=Parameters
        self.FX0=float(FX0)
        self.tau=tau # or we could define capital T , if thats the case we need
        self.r_dom=r_dom 
        self.r_for=r_for
        self.NReps=NReps
        self.NSteps=NSteps
        
    def FX_Heston_Paths(self):
        return hfd.Heston_Sim(self.Parameters, self.FX0,self.K, self.tau,self.r_dom,self.r_for , self.NReps, self.NSteps,self.phi)

    def FX_Heston_Sim_Pricing(self):
        return hfd.Heston_Sim_Option(self.Parameters,self. FX0,self. K, self.tau,self.r_dom,self. r_for ,self. NReps, self.NSteps,self.phi)