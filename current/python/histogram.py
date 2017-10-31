import ROOT
from ROOT import TH1D, TH2D, TF1, gStyle, TCanvas
from math import sqrt, log
import cppyy
import numpy as np
cppyy.load_reflection_info("H10.so")

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

fid_sec_hist = []
for sec in range(sector_num):
    hname = "fid_sec%d" % (sec + 1)
    htitle = "fid_sec%d" % (sec + 1)

    fid_sec_hist.append(TH2D(hname, htitle, bins, min_phi[sec],
                             max_phi[sec], bins, theta_min, theta_max))

"""


def fitGaus(hist, min_value, max_value):
    c1 = TCanvas()
    gaus = "[0]*TMath::Gaus(x,[1],[2],1)"
    fitFunc = TF1("fitFunc", gaus, min_value, max_value)
    fitFunc.SetLineColor(2)
    par_max = hist.GetMaximum() if hist.GetMaximum() else 0
    par_mean = hist.GetMean() if hist.GetMean() else 0
    fitFunc.SetParameter(0, par_max)
    fitFunc.SetParameter(1, par_mean)
    fitFunc.SetParameter(2, 1)
    fitFunc.SetParNames("height", "mean", "FWHM")

    hist.Fit("fitFunc", "qM0+", "", min_value, max_value)

    par_mean = fitFunc.GetParameter(
        "mean") if fitFunc.GetParameter("mean") else 0

    par_FWHM = fitFunc.GetParameter(
        "FWHM") if fitFunc.GetParameter("FWHM") else 0

    fitFunc.SetParameter(0, par_max)
    fitFunc.SetParameter(1, par_mean)
    fitFunc.SetParameter(2, par_FWHM)
    hist.Fit("fitFunc", "qM+", "", min_value, max_value)
    mean = fitFunc.GetParameter("mean")
    FWHM = fitFunc.GetParameter("FWHM")
    sigma = fitFunc.GetParameter("FWHM") / (2 * sqrt(2 * log(2)))
    gStyle.SetOptFit(1111)


def gaussian(x, mu, sig, const):
    return const * 1 / (sig * np.sqrt(2 * np.pi)) * np.exp(-(x - mu)**2 / 2 * sig**2)


def fitGaus_py(hist, min_value, max_value):
    try:
        from scipy.optimize import curve_fit
        from root_numpy import hist2array
        import matplotlib.pyplot as plt
    except:
        return 0

    t = hist2array(hist, return_edges=True)
    xdata = (t[1][0][:-1] + t[1][0][1:]) / 2.0
    ydata = t[0]

    gaus_guess = np.array([0.9, 0.1, 1])

    gaus_1, gaus_cov_1 = curve_fit(
        gaussian, xdata, ydata, p0=gaus_guess, maxfev=80000)

    plt.plot(xdata, gaussian(xdata, gaus_1[0], gaus_1[1], gaus_1[2]))
    plt.xlim((np.min(xdata), np.max(xdata)))
    plt.legend(loc=0)
    plt.xlabel(r'Mass (GeV)', fontsize=20)
    plt.ylabel(r'Counts (#)', fontsize=18)
    plt.savefig("test.pdf")


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

x_y_min_max = 0.5
Beam_Position = TH2D("Beam_Position", "X vs Y", bins, -x_y_min_max, x_y_min_max,
                     bins, -x_y_min_max, x_y_min_max,)
Beam_Position_X = TH1D("Beam_Position_X", "Beam_Position_X",
                       bins, -x_y_min_max, x_y_min_max)
Beam_Position_Y = TH1D("Beam_Position_Y", "Beam_Position_Y",
                       bins, -x_y_min_max, x_y_min_max)

Beam_Pos = {
    'Beam_Position': Beam_Position,
    'Beam_Position_X': Beam_Position_X,
    'Beam_Position_Y': Beam_Position_Y
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

theta_min = 10
theta_max = 60
phi_min = -360 / 2.0
phi_max = 360 / 2.0
sector_num = 6
min_phi = [0, 60, 120, -180, -120, -60]
max_phi = [60, 120, 180, -120, -60, 0]

fid_hist = TH2D("fid", "fid", bins, phi_min, phi_max, bins,
                theta_min, theta_max)

fid_sec_hist = []
for sec in range(sector_num):
    hname = "fid_sec%d" % (sec + 1)
    htitle = "fid_sec%d" % (sec + 1)

    fid_sec_hist.append(TH2D(hname, htitle, bins, min_phi[sec],
                             max_phi[sec], bins, theta_min, theta_max))

fid = {'fid_hist': fid_hist,
       'fid_sec_hist': fid_sec_hist
       }

histo = {}
for d in (W_Q2, P_B, MissMass, delta_t, ec, fid, Beam_Pos):
    histo.update(d)


def add_and_save(output, root_file):
    h10 = cppyy.gbl.H10()
    for _h in output:
        for key, value in histo.items():
            if key == 'fid_sec_hist':
                pass
                # for _i in xrange(len(6)):
                # print(type(_h['fid_sec_hist']))
                # print([type(_i) for _i in _h['fid_sec_hist'])
                # print([_i for _i in _h['fid_sec_hist'])
            else:
                value.Add(_h[key])

    histo['EC_sampling_fraction'].SetXTitle("Momentum (GeV)")
    histo['EC_sampling_fraction'].SetYTitle("Sampling Fraction")
    histo['EC_sampling_fraction'].Write()

    beam_folder = root_file.mkdir("beam")
    beam_folder.cd()
    histo['Beam_Position'].SetXTitle("X (cm)")
    histo['Beam_Position'].SetYTitle("Y (cm)")
    histo['Beam_Position'].Write()
    fitGaus(histo['Beam_Position_X'], -1.0, 1.0)
    fitGaus_py(histo['Beam_Position_X'], -1.0, 1.0)
    histo['Beam_Position_X'].SetXTitle("X (cm)")
    histo['Beam_Position_X'].Write()
    fitGaus(histo['Beam_Position_Y'], -1.0, 1.0)
    fitGaus_py(histo['Beam_Position_Y'], -1.0, 1.0)
    histo['Beam_Position_Y'].SetXTitle("Y (cm)")
    histo['Beam_Position_Y'].Write()

    WvsQ2_folder = root_file.mkdir("W vs Q2")
    WvsQ2_folder.cd()
    histo['WvsQ2_hist'].SetXTitle("W (GeV)")
    histo['WvsQ2_hist'].SetYTitle("Q^{2} (GeV^{2})")
    histo['WvsQ2_hist'].Write()
    histo['W_hist'].SetXTitle("W (GeV)")
    histo['W_hist'].Write()
    histo['Q2_hist'].SetXTitle("Q^{2} (GeV^{2})")
    histo['Q2_hist'].Write()
    histo['E_prime_hist'].SetXTitle("Energy (GeV)")
    histo['E_prime_hist'].Write()
    histo['Q2_vs_xb'].SetXTitle("X_b")
    histo['Q2_vs_xb'].SetYTitle("Q^{2} (GeV^{2})")
    histo['Q2_vs_xb'].Write()
    histo['WvsQ2_proton'].SetXTitle("W (GeV)")
    histo['WvsQ2_proton'].SetYTitle("Q^{2} (GeV^{2})")
    histo['WvsQ2_proton'].Write()
    histo['W_proton'].SetXTitle("W (GeV)")
    histo['W_proton'].Write()
    histo['Q2_proton'].SetXTitle("Q^{2} (GeV^{2})")
    histo['Q2_proton'].Write()
    histo['WvsQ2_pion'].SetXTitle("W (GeV)")
    histo['WvsQ2_pion'].SetYTitle("Q^{2} (GeV^{2})")
    histo['WvsQ2_pion'].Write()
    histo['W_pion'].SetXTitle("W (GeV)")
    histo['W_pion'].Write()
    histo['Q2_pion'].SetXTitle("Q^{2} (GeV^{2})")
    histo['Q2_pion'].Write()
    histo['WvsQ2_single_pi'].SetXTitle("W (GeV)")
    histo['WvsQ2_single_pi'].SetYTitle("Q^{2} (GeV^{2})")
    histo['WvsQ2_single_pi'].Write()
    histo['W_single_pi'].SetXTitle("W (GeV)")
    histo['W_single_pi'].Write()
    histo['Q2_single_pi'].SetXTitle("Q^{2} (GeV^{2})")
    histo['Q2_single_pi'].Write()
    histo['WvsQ2_single_proton'].SetXTitle("W (GeV)")
    histo['WvsQ2_single_proton'].SetYTitle("Q^{2} (GeV^{2})")
    histo['WvsQ2_single_proton'].Write()
    histo['W_single_proton'].SetXTitle("W (GeV)")
    histo['W_single_proton'].Write()
    histo['Q2_single_proton'].SetXTitle("Q^{2} (GeV^{2})")
    histo['Q2_single_proton'].Write()

    pb_dir = root_file.mkdir("Momentum vs Beta")
    pb_dir.cd()
    histo['MomVsBeta_hist'].SetXTitle("Momentum (GeV)")
    histo['MomVsBeta_hist'].SetYTitle("#beta")
    histo['MomVsBeta_hist'].Write()
    histo['MomVsBeta_hist_pos'].SetXTitle("Momentum (GeV)")
    histo['MomVsBeta_hist_pos'].SetYTitle("#beta")
    histo['MomVsBeta_hist_pos'].Write()
    histo['MomVsBeta_hist_neg'].SetXTitle("Momentum (GeV)")
    histo['MomVsBeta_hist_neg'].SetYTitle("#beta")
    histo['MomVsBeta_hist_neg'].Write()
    histo['MomVsBeta_proton_ID'].SetXTitle("Momentum (GeV)")
    histo['MomVsBeta_proton_ID'].SetYTitle("#beta")
    histo['MomVsBeta_proton_ID'].Write()
    histo['MomVsBeta_Pi_ID'].SetXTitle("Momentum (GeV)")
    histo['MomVsBeta_Pi_ID'].SetYTitle("#beta")
    histo['MomVsBeta_Pi_ID'].Write()
    histo['MomVsBeta_proton_Pi_ID'].SetXTitle("Momentum (GeV)")
    histo['MomVsBeta_proton_Pi_ID'].SetYTitle("#beta")
    histo['MomVsBeta_proton_Pi_ID'].Write()
    histo['Mom'].SetXTitle("Momentum (GeV)")
    histo['Mom'].Write()
    histo['Energy_hist'].SetXTitle("Energy (GeV)")
    histo['Energy_hist'].Write()

    mm_dir = root_file.mkdir("Missing Mass")
    mm_dir.cd()
    # h10.missing_mass_fit(histo['Missing_Mass'])
    histo['Mass'].SetXTitle("Mass (GeV)")
    histo['Mass'].Write()
    fitGaus(histo['Missing_Mass'], 0.88, 1.0)
    fitGaus_py(histo['Missing_Mass'], 0.88, 1.0)
    histo['Missing_Mass'].SetXTitle("Mass (GeV)")
    histo['Missing_Mass'].Write()

    fitGaus(histo['Missing_Mass_square'], 0.5, 1.1)
    histo['Missing_Mass_square'].SetXTitle("Mass^{2} (GeV^{2})")
    histo['Missing_Mass_square'].Write()

    DeltaT = root_file.mkdir("Delta_T")
    DeltaT.cd()
    h10.deltat_slicefit(delta_t['delta_t_mass_P'], delta_t['delta_t_mass_PIP'])
    for value in delta_t.values():
        value.SetXTitle("Momentum (GeV)")
        value.SetYTitle("#Deltat")
        value.Write()

    fid_dir = root_file.mkdir("Fid_cuts")
    fid_dir.cd()
    histo['fid_hist'].SetXTitle("#phi")
    histo['fid_hist'].SetYTitle("#theta")
    histo['fid_hist'].Write()

    for _h in output:
        for value in histo.values():
            del value
