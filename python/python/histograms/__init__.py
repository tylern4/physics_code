# https://github.com/NichtJens/numpy-accumulative-histograms/blob/master/accuhist.py

import numpy as np


class Hist1D:
    def __init__(self, nbins, xlow, xhigh):
        self.nbins = nbins
        self.xlow = xlow
        self.xhigh = xhigh

        self.range = (xlow, xhigh)

        self.hist, edges = np.histogram([], bins=nbins, range=self.range)
        self.bins = (edges[:-1] + edges[1:]) / 2.

    def fill(self, arr):
        hist, _ = np.histogram(arr, bins=self.nbins, range=self.range)
        self.hist += hist

    @property
    def data(self):
        return self.bins, self.hist


class Hist2D:
    def __init__(self, nxbins, xlow, xhigh, nybins, ylow, yhigh):
        self.nxbins = nxbins
        self.xhigh = xhigh
        self.xlow = xlow

        self.nybins = nybins
        self.yhigh = yhigh
        self.ylow = ylow

        self.nbins = [nxbins, nybins]
        self.ranges = [[xlow, xhigh], [ylow, yhigh]]

        self.hist, xedges, yedges = np.histogram2d(
            [], [], bins=self.nbins, range=self.ranges)
        self.xbins = (xedges[:-1] + xedges[1:]) / 2.
        self.ybins = (yedges[:-1] + yedges[1:]) / 2.

    def fill(self, xarr, yarr):
        hist, _, _ = np.histogram2d(
            xarr, yarr, bins=self.nbins, range=self.ranges)
        self.hist += hist

    @property
    def data(self):
        return self.xbins, self.ybins, self.hist
