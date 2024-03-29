#!/usr/bin/env python

import ROOT
from ROOT import TLorentzVector
import uproot
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import sys
from multiprocessing import Pool
from glob import glob

ME = 0.00051099895
E0 = 4.81726


@np.vectorize
def q2_calc(px, py, pz):
    e_beam = np.array([0, 0, np.sqrt(E0 ** 2 - ME ** 2), E0])
    e_prime = np.array([px, py, pz, np.linalg.norm([px, py, pz, ME])])

    temp = e_beam - e_prime
    temp2 = (
        temp[0] * temp[0] + temp[1] * temp[1] + temp[2] * temp[2] - temp[3] * temp[3]
    )

    return temp2


def get_min_q2(file_name):
    mc_file = uproot.open(file_name)

    h10 = mc_file["h10"]
    px = h10.array("pxpart")
    py = h10.array("pypart")
    pz = h10.array("pzpart")
    min_q2 = 100
    for _px, _py, _pz in zip(px, py, pz):
        _q2 = q2_calc(_px[0], _py[0], _pz[0])
        min_q2 = np.minimum(_q2, min_q2)

    return f"{min_q2}: {file_name}"


if __name__ == "__main__":
    files = glob(sys.argv[2])

    with Pool(int(sys.argv[1])) as p:
        r = list(tqdm(p.imap(get_min_q2, files), total=len(files)))

    for x in r:
        print(x)

