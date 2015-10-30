from physics import mass
from ROOT import TLorentzVector

def missing_mass_calc(gamma_mu, pi_mu):
	reaction = TLorentzVector(0.0,0.0,0.0,0.0)
	p_mu = TLorentzVector(0.0,0.0,0.0,mass['PROTON'])
	reaction = (gamma_mu + p_mu - pi_mu)
	return reaction.M()