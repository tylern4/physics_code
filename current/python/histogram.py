import ROOT
from ROOT import TH1D, TH2D

bins = 500
p_min = 0.0
p_max = 5.0
hname = ""
htitle = ""
w_min = 0
w_max = 3.25
q2_min = 0.0
q2_max = 10.0

WvsQ2_hist = TH2D("WvsQ2_hist", "W vs Q^{2}", bins, w_min, w_max,
                  bins, q2_min, q2_max)
W_hist = TH1D("W", "W", bins, w_min, w_max)
Q2_hist = TH1D("Q2", "Q^{2}", bins, q2_min, q2_max)
E_prime_hist = TH1D("E_prime", "Scattered Electron Energy", bins, 0.0, 5.0)
Q2_vs_xb = TH2D("Q2_vs_xb", "Q^{2} vs x_{b}", bins, 0.1, 0.6, bins, 1.0, 3.5)
WvsQ2_proton = TH2D("WvsQ2_proton", "W vs Q^{2} P", bins, w_min,
                    w_max, bins, q2_min, q2_max)
W_proton = TH1D("W_proton", "W P", bins, w_min, w_max)
Q2_proton = TH1D("Q2_proton", "Q^{2} P", bins, q2_min, q2_max)
WvsQ2_pion = TH2D("WvsQ2_pion", "W vs Q^{2} #pi^{+} only", bins,
                  w_min, w_max, bins, q2_min, q2_max)
W_pion = TH1D("W_pion", "W #pi^{+} only", bins, w_min, w_max)
Q2_pion = TH1D("Q2_pion", "Q^{2} #pi^{+} only", bins, q2_min, q2_max)
WvsQ2_single_pi = TH2D("WvsQ2_single_pi", "W vs Q^{2} #pi^{+}",
                       bins, w_min, w_max, bins, q2_min, q2_max)
W_single_pi = TH1D("W_single_pi", "W #pi^{+}", bins, w_min, w_max)
Q2_single_pi = TH1D("Q2_single_pi", "Q^{2} #pi^{+}", bins, q2_min, q2_max)
WvsQ2_single_proton = TH2D("WvsQ2_single_proton", "W vs Q^{2} P", bins, w_min, w_max, bins,
                           q2_min, q2_max)
W_single_proton = TH1D("W_single_proton", "W P", bins, w_min, w_max)
Q2_single_proton = TH1D("Q2_single_proton", "Q^{2} P", bins, q2_min, q2_max)


W_Q2 = {'WvsQ2_hist': WvsQ2_hist,
        'W_hist': W_hist,
        'Q2_hist': Q2_hist,
        'E_prime_hist': E_prime_hist,
        'Q2_vs_xb': Q2_vs_xb,
        'WvsQ2_proton': WvsQ2_proton,
        'W_proton': W_proton,
        'Q2_proton': Q2_proton,
        'WvsQ2_pion': WvsQ2_pion,
        'W_pion': W_pion,
        'Q2_pion': Q2_pion,
        'WvsQ2_single_pi': WvsQ2_single_pi,
        'W_single_pi': W_single_pi,
        'Q2_single_pi': Q2_single_pi,
        'WvsQ2_single_proton': WvsQ2_single_proton,
        'W_single_proton': W_single_proton,
        'Q2_single_proton': Q2_single_proton
        }


b_min = 0.1
b_max = 1.2
MomVsBeta_hist = TH2D("MomVsBeta", "Momentum versus #beta", bins,
                      p_min, p_max, bins, b_min, b_max)
MomVsBeta_hist_pos = TH2D("MomVsBeta_pos", "Momentum versus #beta Positive", bins, p_min,
                          p_max, bins, b_min, b_max)
MomVsBeta_hist_neg = TH2D("MomVsBeta_neg", "Momentum versus #beta Negative", bins, p_min,
                          p_max, bins, b_min, b_max)
Mom = TH1D("Momentum", "Momentum", bins, 0, 2.5)
Energy_hist = TH1D("Energy_hist", "Energy_hist", bins, 0.0, 2.5)
MomVsBeta_proton_ID = TH2D("MomVsBeta_proton_ID", "Momentum versus #beta p", bins, p_min,
                           p_max, bins, b_min, b_max)
MomVsBeta_Pi_ID = TH2D("MomVsBeta_Pi_ID", "Momentum versus #beta #pi^{+}", bins, p_min,
                       p_max, bins, b_min, b_max)
MomVsBeta_proton_Pi_ID = TH2D("MomVsBeta_proton_Pi_ID", "Momentum versus #beta P #pi^{+}",
                              bins, p_min, p_max, bins, b_min, b_max)

P_B = {'MomVsBeta_hist': MomVsBeta_hist,
       'MomVsBeta_hist_pos': MomVsBeta_hist_pos,
       'MomVsBeta_hist_neg': MomVsBeta_hist_neg,
       'Mom': Mom,
       'Energy_hist': Energy_hist,
       'MomVsBeta_proton_ID': MomVsBeta_proton_ID,
       'MomVsBeta_Pi_ID': MomVsBeta_Pi_ID,
       'MomVsBeta_proton_Pi_ID': MomVsBeta_proton_Pi_ID
       }

bins_MM = 200
mm_min = 0.0
mm_max = 3.0

Mass = TH1D("Mass", "Mass", 600, 0, 6)
Missing_Mass = TH1D("Missing_Mass", "Missing Mass", bins_MM, mm_min, mm_max)
Missing_Mass_square = TH1D(
    "Missing_Mass_square", "Missing Mass square", bins_MM, mm_min, mm_max * mm_max)

MissMass = {'Mass': Mass,
            'Missing_Mass': Missing_Mass,
            'Missing_Mass_square': Missing_Mass_square
            }

N_SIGMA = 3
Dt_min = -10
Dt_max = 10
num_points = 20
# delta_t_hist[3][num_points]
bin_width = (p_max - p_min) / num_points

sc_sector_num = 6
sc_paddle_num = 48
# delta_t_sec_pad_hist[3][sc_sector_num][sc_paddle_num]
delta_t_mass_P = TH2D("delta_t_mass_P", "#Deltat assuming mass of proton", bins, p_min,
                      p_max, bins, Dt_min, Dt_max)
delta_t_mass_P_PID = TH2D(
    "delta_t_mass_P_PID", "#Deltat assuming mass of proton with PID proton",
    bins, p_min, p_max, bins, Dt_min, Dt_max)

delta_t_mass_PIP = TH2D("delta_t_mass_PIP", "#Deltat assuming mass of #pi^{+}", bins,
                        p_min, p_max, bins, Dt_min, Dt_max)
delta_t_mass_PIP_PID = TH2D("delta_t_mass_PIP_PID",
                            "#Deltat assuming mass of #pi^{+} with PID #pi^{+}", bins, p_min,
                            p_max, bins, Dt_min, Dt_max)

delta_t = {'delta_t_mass_P': delta_t_mass_P,
           'delta_t_mass_P_PID': delta_t_mass_P_PID,
           'delta_t_mass_PIP': delta_t_mass_PIP,
           'delta_t_mass_PIP_PID': delta_t_mass_PIP_PID
           }

ec_min = 0
ec_max = 1
ec_sampling_fraction = TH2D(
    "EC_sampling_fraction", "EC_sampling_fraction", bins, p_min, p_max, bins, ec_min, ec_max)

ec = {'EC_sampling_fraction': ec_sampling_fraction}

"""
bins_CC = 50
CC_min = 0
CC_max = 250
L_R_C - ""
sector = 6
segment = 18
PMT = 3
cc_hist[6][18][3]
cc_hist_allSeg[6][3]
ndims_cc_sparse = 4
bins_cc_sparse[ndims_cc_sparse] = {sector, segment, PMT, bins_CC}
xmin_cc_sparse[ndims_cc_sparse] = {0.0, 0.0, -2.0, CC_min}
xmax_cc_sparse[ndims_cc_sparse] = {sector + 1.0, segment + 1.0, 1.0,
                                   CC_max}
x_cc_sparse[ndims_cc_sparse]
cc_sparse = THnSparseD("cc_sparse", "Histogram", ndims_cc_sparse,
                       bins_cc_sparse, xmin_cc_sparse, xmax_cc_sparse)

theta_min = 0
theta_max = 90
phi_min = -360 / 2.0
phi_max = 360 / 2.0

sector_num = 6
fid_sec_hist[sector_num]
fid_hist = TH2D("fid", "fid", bins, phi_min, phi_max, bins,
                theta_min, theta_max)
"""


histo = {}
for d in (W_Q2, P_B, MissMass, delta_t, ec):
    histo.update(d)
