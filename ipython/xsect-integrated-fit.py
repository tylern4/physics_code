#!/usr/bin/python2

"""Provides class, configuration functions, and main function meant to make
systematic studies of integrated background subtraction easier.  To add a new
fitting scheme, create a new FS* function that sets the adjustable parameters
of FitSchemeXsectInt and add corresponding scheme to 'fss' map in main().
"""

# TODO: simulation-inspired background function
# TODO: add radiative tail to signal function
# TODO: allow class to accumulate cross-sections
#       W, Q^2, cross-sections, error, and parameters see fit() of
#       xsect-integrated-fit.C. remember to propagate errors!
# README: background subtraction error propagation
#         average relative error of sideband yields
#         and apply to integral of background function under signal
# TODO: add background-subtracted histogram
#       for each scheme to canvas and legend
# TODO: add ability to change functional forms
#       number of parameters etc.
# TODO: pull reusables to module
#       e.g., MASS_P, goodcolors, common regexps, etc.
# TODO: ParticleConstants for python
# TODO: maybe add resolution convolution
# TODO: option to create histograms
#       rather than pulling from file (a la h3maker.h)
# TODO: disambiguate "sig"
#       which current can refer to sigma or signal

import sys
import math as m
import re

import ROOT as R
from ROOT import TFile
from ROOT import TF1
from ROOT import TCanvas
from ROOT import TLegend
from ROOT import TH1D

# from epxsectutils import vgflux

# R.gROOT.SetStyle("Simple")
R.gStyle.SetOptStat(0000)
R.gStyle.SetOptFit(0000)

MASS_P = 0.93827203
RE_FLOAT = '[-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?'
goodcolors = [R.kRed+1,
              R.kGreen+1,
              R.kBlue,
              R.kYellow+1,
              R.kMagenta+1,
              R.kCyan+1,
              9]


class FitSchemeXsectInt:
    """Primary functions include Setup(), Fit(), and Decorate().

        * Setup() determines parameters of main fit;
        * Fit(), well, fits; and
        * Decorate() adds background (bg), signal (sig), and combined (fn)
          functions, along with any other "decorations" or visual modification
          to the histogram.
    """
    def __init__(self, desc, name='', fnsetup=None):
        self.g = R.TMinuit()    # leave reference active for access to gMinuit
        self.name = name
        # ##########################################################
        # ######## default values of adjustable parameters #########
        # ##########################################################
        self.prefit = None
        self.wrange, self.q2range = [0, 0], [0, 0]
        self.drawrange = [0.4, 2]
        self.endAtEdge = False
        self.bgNparms = 7
        self.sigNparms = 5
        self.fnNparms = 10
        self.bg = TF1('fbg_%s' % self.name, self.d_pol4,
                      self.drawrange[0], self.drawrange[1],
                      self.bgNparms)
        self.sig = TF1('fsig_%s' % self.name, self.d_sig,
                       self.drawrange[0], self.drawrange[1],
                       self.sigNparms)
        self.fn = TF1('fbgsig_%s' % self.name, self.d_pol4gaus,
                      self.drawrange[0], self.drawrange[1],
                      self.fnNparms)
        # fit ranges
        self.bg.range = [0.4, 1.1]
        self.sig.range = [0.74, 0.85]
        self.fn.range = [0.4, 1.1]
        # initial parameters
        self.bg.setparms = {}
        self.sig.setparms = {1: 0.783, 2: 0.020}
        self.fn.setparms = {}
        # parameter limits
        self.bg.parlims = {}
        self.sig.parlims = {2: [0.01, 0.04]}
        self.fn.parlims = {7: [0.01, 0.04]}
        # fixed parameters
        self.bg.fixparms = {3: 0, 4: 0}
        self.fn.fixparms = {3: 0, 4: 0}
        # ranges to skip during fitting
        self.bg.skip = [[0.74, 0.85], [0.5, 0.6]]
        self.sig.skip = []
        self.fn.skip = [[0.5, 0.6]]
        # x values of modified step, i.e., x0/x1 of stepfactor()
        self.edgerange = [0.77, 0.9]
        # draw options during fitting
        self.doptions = 'QN'
        # ################## end default values ####################
        self.goptions = ''
        self.desc = desc
        self.converged = True
        # override default parameters
        if fnsetup is not None:
            fnsetup(self)

    # #######################################################################
    # #################### DEFAULT FIT FUNCTIONS ############################
    # #######################################################################
    def stepfactor(self, x, x0, x1):
        return 1 if x < x0 else (0 if x > x1 else 1-1/(x1-x0)*(x-x0))

    def d_sig(self, x, par):
        step = self.stepfactor(x[0], par[3], par[4])
        gaus = par[0]*m.exp(-0.5*pow(((x[0]-par[1])/par[2]), 2))
        return step*gaus

    def d_pol4(self, x, par):
        # reject points in ranges not modeled by bg function
        for skip in self.bg.skip:
            if x[0] > skip[0] and x[0] < skip[1]:
                TF1.RejectPoint()
                break
        step = self.stepfactor(x[0], par[5], par[6])
        retval = step*(par[0] + par[1]*x[0] + par[2]*x[0]**2
                       + par[3]*x[0]**3 + par[4]*x[0]**4)
        return 0 if retval < 0 else retval

    def d_pol4gaus(self, x, par):
        # reject points in ranges not modeled by function
        for skip in self.fn.skip:
            if x[0] > skip[0] and x[0] < skip[1]:
                TF1.RejectPoint()
                break
        step = self.stepfactor(x[0], par[8], par[9])
        return step*(par[0] + par[1]*x[0] + par[2]*x[0]**2
                     + par[3]*x[0]**3 + par[4]*x[0]**4
                     + par[5]*m.exp(-0.5*pow(((x[0]-par[6])/par[7]), 2)))
    # ################## END DEFAULT FIT FUNCTIONS ##########################

    def Setup(self, hist, wrange=None, q2range=None):
        # attempt to set modified step function (phase-space edge) parameters
        # and populate W and Q2 ranges
        # TODO: change re.search to re.findall
        #       and get wlow, whigh, q2low, q2high in one expression
        wlow, whigh = 0, 0
        if wrange is not None:
            wlow, whigh = wrange[0], wrange[1]
        else:
            # try to get wrange from histogram title
            # assuming histograms of h3maker
            # look for '(W = <number>' without returning '(W = '
            re_string = '(?<=W = \()'+RE_FLOAT
            wlow = float(re.search(re_string, hist.GetTitle()).group(0))
            # look for ',<number)' without returning ',' or ')'
            re_string = '(?<=,)'+RE_FLOAT+'(?<=\))?'
            whigh = float(re.search(re_string, hist.GetTitle()).group(0))
        # fix stepfactor parameters according to wrange
        if whigh > 0:
            x0 = m.sqrt(wlow**2+MASS_P**2-2*wlow*MASS_P)
            x1 = m.sqrt(whigh**2+MASS_P**2-2*whigh*MASS_P)
            self.edgerange = [x0, x1]
        self.wrange[0], self.wrange[1] = wlow, whigh

        q2low, q2high = 0, 0
        if q2range is not None:
            q2low, q2high = q2range[0], q2range[1]
        else:
            # try to get wrange from histogram title
            # assuming histograms of h3maker
            # look for '(W = <number>' without returning '(W = '
            re_string = '(?<=Q\^2 = \()'+RE_FLOAT
            q2low = float(re.search(re_string, hist.GetTitle()).group(0))
            # look for ',<number)' without returning ',' or ')'
            re_string = '(?<=,)'+RE_FLOAT+'(?<=\))?'
            q2high = float(re.findall(re_string, hist.GetTitle())[1][0])
        self.q2range[0], self.q2range[1] = q2low, q2high

        # estimate background subtraction
        for pidx, pval in self.bg.setparms.items():
            if hasattr(pval, '__call__'):
                pval = pval(self)
            self.bg.SetParameter(pidx, pval)
        for pidx, (plo, phi) in self.bg.parlims.items():
            self.bg.SetParLimits(pidx, plo, phi)
        for pidx, pval in self.bg.fixparms.items():
            if hasattr(pval, '__call__'):
                pval = pval(self)
            self.bg.FixParameter(pidx, pval)
        x1 = self.bg.range[1]
        if self.endAtEdge:
            x1 = x1 if x1 < self.edgerange[0] else self.edgerange[0]
        hist.Fit(self.bg, self.doptions, self.goptions,
                 self.bg.range[0], x1)
        bg = self.bg.Clone('bgtmp')
        for pidx, pval in self.bg.fixparms.items():
            if hasattr(pval, '__call__'):
                pval = pval(self)
            bg.FixParameter(pidx, pval)
        h = hist.Clone('%s_sig' % hist.GetName())
        h.Add(bg, -1)

        # estimate signal function parameters
        for pidx, pval in self.sig.setparms.items():
            if hasattr(pval, '__call__'):
                pval = pval(self)
            self.sig.SetParameter(pidx, pval)
        for pidx, (plo, phi) in self.sig.parlims.items():
            self.sig.SetParLimits(pidx, plo, phi)
        for pidx, pval in self.sig.fixparms.items():
            if hasattr(pval, '__call__'):
                pval = pval(self)
            self.sig.FixParameter(pidx, pval)
        x1 = self.sig.range[1]
        if self.endAtEdge:
            x1 = x1 if x1 < self.edgerange[0] else self.edgerange[0]
        h.Fit(self.sig, self.doptions, self.goptions,
              self.sig.range[0], x1)
        h.SetMarkerColor(R.kBlue)
        h.SetLineColor(R.kBlue)
        self.hist_sig = h

        # ################################################################
        # ################ bin-centric Xsect #############################
        # ################################################################
        sigma = self.sig.GetParameter(2)
        intrange = (h.FindBin(0.783-3*sigma), h.FindBin(0.783+3*sigma))
        self.N = 0
        self.Nerr = 0
        berr2 = 0
        for b in range(intrange[0], intrange[1]):
            byield = h.GetBinContent(b)
            berr2 += h.GetBinError(b)**2
            self.N += byield
        self.Nerr = m.sqrt(berr2)
        dw = self.wrange[1] - self.wrange[0]
        dq2 = self.q2range[1] - self.q2range[0]
        br = 0.891
        lum = 19.844
        self.xsect = self.N/(br*dw*dq2*lum*0.997)
        self.xsecterr = self.Nerr/(br*dw*dq2*lum*0.997)

        # TODO: !!! BEWARE of *Nparms confusion
        #           during param copying, code currently ASSUMES 2 extra
        #           "edge" parameters on bg and sig, so it throws off the
        #           bg/sig-to-fn parameter mapping
        # set initial parameters of bg+sig function
        [self.fn.SetParameter(i, v) for i, v
            in zip(range(0, self.bgNparms-2), bg.GetParameters())]
        [self.fn.SetParameter(i, v) for i, v
            in zip(range(self.bgNparms-2, self.fnNparms),
                   self.sig.GetParameters())]
        for pidx, (plo, phi) in self.fn.parlims.items():
            self.fn.SetParLimits(pidx, plo, phi)
        for pidx, pval in self.fn.fixparms.items():
            if hasattr(pval, '__call__'):
                pval = pval(self)
            self.fn.FixParameter(pidx, pval)
        return hist

    def Fit(self, hist):
        x1 = self.bg.range[1]
        if self.endAtEdge:
            x1 = x1 if x1 < self.edgerange[0] else self.edgerange[0]
        if self.prefit is not None:
            self.prefit(self)
        hist.Fit(self.fn, self.doptions, self.goptions,
                 self.fn.range[0], x1)
        self.converged = R.gMinuit.fCstatu.startswith('CONV')
        for pidx in range(0, self.bgNparms-2):
            self.bg.SetParameter(pidx, self.fn.GetParameter(pidx))
        for pidx in range(0, self.sigNparms-2):
            self.sig.SetParameter(pidx,
                                  self.fn.GetParameter(pidx+self.bgNparms-2))
        self.bg.SetParameter(self.bgNparms-2, self.fn.GetParameter(8))
        self.bg.SetParameter(self.bgNparms-1, self.fn.GetParameter(9))
        self.sig.SetParameter(self.sigNparms-2, self.fn.GetParameter(8))
        self.sig.SetParameter(self.sigNparms-1, self.fn.GetParameter(9))

        # ##################################################################
        # ############### calculate cross-section and error ################
        # ##################################################################
        self.sig.Nf, self.sig.Nferr = 0, 0
        Nferr2 = 0
        sigma = self.sig.GetParameter(2)
        intrange = (hist.FindBin(0.783-3*sigma), hist.FindBin(0.783+3*sigma))
        for b in range(intrange[0], intrange[1]):
            byield = hist.GetBinContent(b)
            berr2 = hist.GetBinError(b)**2
            bNf = self.sig.Eval(hist.GetXaxis().GetBinCenter(b))
            if byield > 0 and bNf > 0:
                self.sig.Nf += bNf
                Nferr2 += (self.sig.Nf/byield)*berr2
        self.sig.Nferr = m.sqrt(Nferr2)
        dw = self.wrange[1] - self.wrange[0]
        dq2 = self.q2range[1] - self.q2range[0]
        br = 0.891
        lum = 19.844
        self.sig.xsect = self.sig.Nf/(br*dw*dq2*lum)
        self.sig.xsecterr = self.sig.Nferr/(br*dw*dq2*lum)

    def Decorate(self, hist):
        fns = [self.bg.Clone('fbg'), self.sig.Clone('fsig'),
               self.fn.Clone('fbgsig')]
        for fn in fns:
            fn.SetRange(self.drawrange[0], self.drawrange[1])
            hist.GetListOfFunctions().Add(fn)


def sigmafix(fs):
    wval = (fs.wrange[0]+fs.wrange[1])/2
    return -0.0065+0.013*wval


def FSpol2trunc(fs):
    fs.endAtEdge = True
    fs.drawrange = [0.4, 2]
    fs.bg.nparms = 7
    fs.bg.range = [0.4, 1.1]
    fs.bg.skip = [[0.74, 0.85], [0.5, 0.6]]
    fs.bg.setparms = {}
    fs.bg.fixparms = {3: 0, 4: 0,
                      5: lambda s: s.edgerange[0],
                      6: lambda s: s.edgerange[1]}
    fs.sig.nparms = 5
    fs.sig.range = [0.74, 0.85]
    fs.sig.parlims = {2: [0.01, 0.04]}
    fs.sig.setparms = {1: 0.783, 2: 0.020}
    fs.sig.fixparms = {1: 0.783, 2: sigmafix,
                       3: lambda s: s.edgerange[0],
                       4: lambda s: s.edgerange[1]}
    fs.fn.nparms = 10
    fs.fn.range = [0.4, 1.1]
    fs.fn.skip = [[0.5, 0.6]]
    fs.fn.setparms = {}
    fs.fn.parlims = {7: [0.01, 0.04]}
    fs.fn.fixparms = {3: 0, 4: 0, 6: 0.783, 7: sigmafix,
                      8: lambda s: s.edgerange[0],
                      9: lambda s: s.edgerange[1]}


def FSpol4trunc(fs):
    # use 2nd order polynomial
    FSpol2trunc(fs)
    # remove restriction on 3rd and 4th order
    del fs.bg.fixparms[3]
    del fs.bg.fixparms[4]
    del fs.fn.fixparms[3]
    del fs.fn.fixparms[4]


def FSpol2trunc5(fs):
    FSpol2trunc(fs)
    fs.prefit = lambda s: s.fn.FixParameter(2, 1.05*s.fn.GetParameter(2))


def FSpol2truncM5(fs):
    FSpol2trunc(fs)
    fs.prefit = lambda s: s.fn.FixParameter(2, 0.95*s.fn.GetParameter(2))


def FSpol4full(fs):
    # use truncated 4th order polynomial
    FSpol4trunc(fs)
    #fs.doptions = ''
    # open range
    fs.bg.range = [0.4, 2]
    fs.fn.range = [0.4, 2]
    fs.endAtEdge = False


def gettext(x, y, s, f=12):
    t = R.TLatex(x, y, R.Form('#font[41]{}%s' % s))
    t.SetNDC()
    t.SetX(x)
    t.SetY(y)
    t.SetTextFont(f)
    t.SetTextAlign(12)
    t.SetTextSize(0.045)
    return t


def geterrrel(hist):
    h = TH1D('%s_relerr' % hist.GetName(), '%s_relerr' % hist.GetName(),
             hist.GetNbinsX(), hist.GetXaxis().GetBinLowEdge(1),
             hist.GetXaxis().GetBinLowEdge(hist.GetNbinsX()))
    for i in range(1, hist.GetNbinsX()):
        N = hist.GetBinContent(i)
        err = hist.GetBinError(i)
        relerr = err/N if N > 0 else 0
        h.SetBinContent(i, relerr)
    return h


def geterrsqrtN(hist):
    h = TH1D('%s_relerr' % hist.GetName(), '%s_relerr' % hist.GetName(),
             hist.GetNbinsX(), hist.GetXaxis().GetBinLowEdge(1),
             hist.GetXaxis().GetBinLowEdge(hist.GetNbinsX()))
    for i in range(1, hist.GetNbinsX()):
        N = hist.GetBinContent(i)
        err = hist.GetBinError(i)
        relerr = err/m.sqrt(N) if N > 0 else 0
        h.SetBinContent(i, relerr)
    return h


def main():
    R.gROOT.SetBatch(True)
    R.gErrorIgnoreLevel = R.kError
    fss = {'2^{nd}-order':
           FitSchemeXsectInt('2nd order polynomial background, Gauss signal.\
                             fit range = 0.4-1.1, skipping eta peak', 'p2t',
                             FSpol2trunc),
           '4^{th}-order':
           FitSchemeXsectInt('4th order polynomial background, Gauss signal.\
                             fit range = 0.4-2.0, skipping eta peak', 'p4t',
                             FSpol4full),
           '2^{nd}-order, par2 -5%':
           FitSchemeXsectInt('2th order polynomial background, Gauss signal.\
                             fit range = 0.4-1.1, skipping eta peak, -5\% p2',
                             'p2tm5', FSpol2truncM5),
           '2^{nd}-order, par2 +5%':
           FitSchemeXsectInt('2th order polynomial background, Gauss signal.\
                             fit range = 0.4-1.1, skipping eta peak, +5\% p2',
                             'p2t5', FSpol2trunc5)}
    orderedkeys = ['2^{nd}-order', '4^{th}-order',
                   '2^{nd}-order, par2 -5%', '2^{nd}-order, par2 +5%']
    fin = TFile('/home/ephelps/projects/phys-ana-omega/h3maker-hn.root')
    # fin = TFile('/home/ephelps/projects/phys-ana-omega/h3maker-hn-now8.root')
    # eff_cc = fin.Get('hq2w_eff_cc')
    # eff_acc = fin.Get('hq2w_eff_acc')
    # fin = TFile('/data/ephelps.bak/analysis/sandbox/h3maker-hn.root')
    nq2bins = 7        # TODO: make nq2bins dynamic
    outgroups = 0
    hss = [fin.Get('hs%d' % i).GetHists() for i in range(0, nq2bins)]
    outf = open('xsects.txt', 'w')
    # TODO: make number of xsects dynamic
    outf.write('W/D\tQ2\txs1\te1\txs2\te2\txs3\te3\txs4\te4')
    outf.write('\txs1b\te1b\txs2b\te2b\txs3b\te3b\txs4b\te4b\n')
    for h in [h for hs in hss for h in hs if h.Integral() > 10000]:
        c = TCanvas('cpreview', 'preview')
        c.cd()
        h.Draw()
        leg = TLegend(0.65, 0.67, 0.99, 0.91)
        wrange, q2range = [], []
        wval, q2val = 0, 0
        gsigs = []
        xss, xssb = [], []
        for i, k, v in [(i, orderedkeys[i], fss[orderedkeys[i]])
                        for i in range(0, len(orderedkeys))]:
            v.Setup(h)
            v.Fit(h)
            v.bg.SetLineColor(goodcolors[i])
            v.sig.SetLineColor(goodcolors[i])
            v.sig.SetLineStyle(2)
            v.fn.SetLineColor(goodcolors[i])
            v.fn.SetLineStyle(3)
            v.bg.SetNpx(500)
            v.sig.SetNpx(500)
            v.fn.SetNpx(500)
            h.GetListOfFunctions().Add(v.bg)
            h.GetListOfFunctions().Add(v.sig)
            h.GetListOfFunctions().Add(v.fn)
            # range only needs to be set once, but since range is accessible here,
            # h.GetXaxis().SetRangeUser(v.drawrange[0], v.edgerange[1]+0.1)
            lbl = k if v.converged else '%s (X)' % k
            gsig = v.fn.GetParameter(7)
            gsiginrange = gsig < 0.038 and gsig > 0.012
            if not v.converged:
                sys.stdout.write('C')
            elif not gsiginrange:
                sys.stdout.write('S')
            else:
                sys.stdout.write('.')
            lbl = lbl if gsiginrange else '%s (#sigma)' % lbl
            lbl += ' '
            leg.AddEntry(v.bg, lbl, 'l')
            wrange, q2range = v.wrange, v.q2range
            wval, q2val = (v.wrange[1]+v.wrange[0])/2, (v.q2range[1]+v.q2range[0])/2
            gsigs.append(gettext(0.15, 0.8-i*0.05, '#sigma = %.1f MeV' % (gsig*1000)))
            gsigs[i].SetTextColor(goodcolors[i])
            # gidx = eff_acc.FindBin(wval, q2val)
            # dw = wrange[1]-wrange[0]
            # dq2 = q2range[1]-q2range[0]
            # accf = eff_acc.GetBinContent(gidx)
            # ccf = eff_cc.GetBinContent(gidx)
            # flx = vgflux(wval, q2val)
            # brf = 0.891
            # lum = 19.844
            # corrinv = brf*lum*flx*accf*ccf*dw*dq2*1e6
            # corr = 0 if corrinv == 0 else 1/corrinv
            corr = 1/(1e6)
            xss.append((v.sig.xsect*corr, v.sig.xsecterr*corr))
            xssb.append((v.xsect*corr, v.xsecterr*corr))
        histsig = fss[orderedkeys[0]].hist_sig
        histsig.Draw('same')
        h.GetListOfFunctions().Add(leg)
        xstr = '%.3f\t%.3f' % (round(wval, 3), round(q2val, 3))
        for x, e in xss:
            # xstr += '\t%.0f\t%.0f' % (x, e)
            xstr += '\t%f\t%f' % (x, e)
        for x, e in xssb:
            xstr += '\t%f\t%f' % (x, e)
        outf.write(xstr + '\n')
        xhi = 0.933   # v.edgerange[1]+0.1
        xlo = 0.633   # v.drawrange[0]
        h.GetXaxis().SetTitle('M_{X} (GeV) in #gamma* p #rightarrow pX')
        h.GetXaxis().SetTitleSize(0.0375)
        h.GetXaxis().SetTitleOffset(1.2)
        h.GetXaxis().SetRangeUser(xlo, xhi)
        yo = 1.1*h.GetMaximum()
        h.GetYaxis().SetTitle('corrected yield')
        h.GetYaxis().SetTitleSize(0.0375)
        h.GetYaxis().SetRangeUser(0, yo)
        h.SetTitle('W = [%.3f, %.3f), Q^{2} = [%.3f, %.3f)'
                   % (wrange[0], wrange[1], q2range[0], q2range[1]))
        for gsig in gsigs:
            gsig.Draw('same')
        c.Modified()
        c.Update()
        c.SaveAs('xsect/%s.pdf' % h.GetName())
        herr = geterrsqrtN(h)
        herr.Draw()
        herr.GetXaxis().SetRangeUser(xlo, xhi)
        herr.SetMinimum(0)
        herr.SetTitle(h.GetTitle())
        c.Modified()
        c.Update()
        c.SaveAs('xsect/errnorm/%s.pdf' % herr.GetName())
        sys.stdout.flush()
        outgroups += 1
        if outgroups % 10 == 0:
            sys.stdout.write('\n')
        else:
            sys.stdout.write(' | ')
        c.Delete()
    print('\nFinished!')

if __name__ == '__main__':
    main()
