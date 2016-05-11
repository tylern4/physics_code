from ROOT import TLorentzVector
from ROOT import TVector3

mass = {'ELECTRON':0.000511, 'PIP':0.13957, 'PIM':0.13957, 'PROTON':0.93827, 'NEUTRON':0.939565, 0:0, 22:0, 'K+':0.493667, 'K-':0.493667, 45:0, 47:0, 49:0}
masses = {11:0.000511, 211:0.13957, -211:0.13957, 2212:0.93827, 2112:0.939565, 0:0, 22:0, 321:0.493667, -321:0.493667, 45:0, 47:0, 49:0}
getM2 = lambda x: [masses[pid]**2 for pid in x]

branches = ['npart',
			'gpart', 
			'id',
			'evntid',
			'W',
			'Q2',
			'MM',
			'm',
			#'evtype',
			#'evntclas',
			#'evthel',
			#'evntclas2',
			#'q_l',
			#'t_l',
			#'tr_time',
			#'rf_time1',
			#'rf_time2',
			
			#'stat',   
			#'dc',   
			#'cc',   
			#'sc',   
			#'ec',   
			#'lec',   
			'p',   
			'q',   
			'b',   
			'cx',   
			'cy',   
			'cz',   
			'vx',   
			'vy',   
			'vz'#,   
			
			#'dc_part',		
			#'dc_sect',   
			#'dc_trk',   
			#'dc_stat',   
			#'dc_xsc',   
			#'dc_ysc',   
			#'dc_zsc',   
			#'dc_cxsc',   
			#'dc_cysc',   
			#'dc_czsc',   
			#'dc_xec',   
			#'dc_yec',   
			#'dc_zec',   
			#'dc_thcc',   
			#'dc_c2',   
			
			#'ec_part',
			#'ec_stat',   
			#'ec_sect',   
			#'ec_whol',   
			#'etot',   
			#'ec_ei',   
			#'ec_eo',   
			#'ec_t',   
			#'ec_r',   
			#'ech_x',   
			#'ech_y',   
			#'ech_z',   
			#'ec_c2',   
			
			#'sc_part',
			#'sc_sect',   
			#'sc_hit',   
			#'sc_pd',   
			#'sc_stat',   
			#'edep',   
			#'sc_t',   
			#'sc_r',   
			#'sc_c2',   
			
			#'cc_part',
			#'cc_sect',   
			#'cc_hit',   
			#'cc_segm',   
			#'nphe',   
			#'cc_t',   
			#'cc_r',   
			#'cc_c2'
			] #This sets the branch names, should be able to add all of them to load everything or just use the ones I need for now

#Calcuating Q^2 
# q^mu^2 = (e^mu - e^mu')^2 = -Q^2
def Q2_calc(e_mu, e_mu_prime):
	q_mu = TLorentzVector()
	q_mu = (e_mu - e_mu_prime)
	return -q_mu.Mag2()

#	Calcualting W
#	Gotten from s channel [(gamma - P)^2 == s == w^2]
#	Sqrt[M_p^2 - Q^2 + 2 M_p gamma]
def W_calc(e_mu, e_mu_prime):
	q_mu = TLorentzVector()
	q_mu = (e_mu - e_mu_prime)
	p_mu_3 = TVector3(0,0,0)
	p_mu = TLorentzVector()
	p_mu.SetVectM(p_mu_3,mass['PROTON'])
	return (p_mu + q_mu).Mag()

def xb_calc(Q2, E_prime):
	gamma = E1D_E0-E_prime
	xb = (Q2/(2 * mass['PROTON'] * gamma))
	return xb

#overload with 4 vectors instaed of otehr calculations
def xb_calc(e_mu, e_mu_prime):
	Q2 = Q2_calc(e_mu,e_mu_prime)
	q = TLorentzVector()
	q = e_mu - e_mu_prime
	target = TLorentzVector(0, 0, 0, mass['PROTON'])
	return (Q2/ (2 * (q.Dot(target))))