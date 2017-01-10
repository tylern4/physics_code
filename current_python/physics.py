import numpy as np
from ROOT import TLorentzVector,TVector3
from constants import *

def all_mom(p,cx,cy,cz):
	calc = lambda _p, _cos: _p * _cos
	return calc(p,cx),calc(p,cy),calc(p,cz)

def fvec(_px,_py,_pz,_mass):
	_vec = TVector3(_px,_py,_pz)
	_4vec = TLorentzVector()
	_4vec.SetVectM(_vec,_mass)
	return _4vec

def fourvec(_p,_cx,_cy,_cz,_mass):
	_px,_py,_pz = all_mom(_p,_cx,_cy,_cz)
	return fvec(_px,_py,_pz,_mass)

Square = lambda x: x**2
append = lambda _arr,_val: np.append(_arr,_val)
array = lambda _:np.array([])

#Beam 4 vector
e_mu = TLorentzVector(0.0,0.0, np.sqrt(Square(E1D_E0)-Square(mass['ELECTRON'])), E1D_E0)
_p_target = TLorentzVector(0, 0, 0, mass['PROTON'])

#Calcuating Q^2 
# q^mu^2 = (e^mu - e^mu')^2 = -Q^2
def Q2_calc(_e_mu, _e_mu_prime):
	_q_mu = (_e_mu - _e_mu_prime)
	return -_q_mu.Mag2()

#	Calcualting W
#	Gotten from s channel [(gamma - P)^2 == s == w^2]
#	Sqrt[M_p^2 - Q^2 + 2 M_p gamma]
def W_calc(_e_mu, _e_mu_prime):
	_q_mu = (_e_mu - _e_mu_prime)
	return (_p_target + _q_mu).Mag()

def xb_calc(_Q2, _E_prime):
	_g = E1D_E0 - _E_prime
	_xb = (_Q2/(2 * mass['PROTON'] * _g))
	return _xb

def xb_calc(_e_mu, _e_mu_prime):
	_Q2 = _Q2_calc(_e_mu,_e_mu_prime)
	_q = _e_mu - _e_mu_prime
	return (_Q2/ (2 * (_q.Dot(target))))