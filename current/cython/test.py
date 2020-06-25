import multiprocessing as mp
import time
from multiprocessing import Pool, cpu_count

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

files = ["/Users/tylern/Data/e1d/h10_22904.root",
         "/Users/tylern/Data/e1d/h10_22905.root",
         "/Users/tylern/Data/e1d/h10_22906.root",
         "/Users/tylern/Data/e1d/h10_22907.root"]


def datahandler(file_name):
    from h10 import h10_data
    h10 = h10_data("h10", file_name)
    p_vs_b = np.asarray([[np.nan] * h10.num_entries,
                         [np.nan] * h10.num_entries])
    out = []
    for evt in h10:
        px = np.asarray([np.nan] * evt.gpart)
        py = np.asarray([np.nan] * evt.gpart)
        pz = np.asarray([np.nan] * evt.gpart)
        for i in range(0, evt.gpart):
            px[i] = evt.p[i] * evt.cx[i]
            py[i] = evt.p[i] * evt.cy[i]
            pz[i] = evt.p[i] * evt.cz[i]
            out.append({'p': evt.p[i], 'b': evt.b[i]})
    return pd.DataFrame(out)


if __name__ == '__main__':
    pool = Pool(processes=4)
    out = []
    for f in files:
        datahandler(f)
    #out = map(datahandler, files)
    pool.close()
    pool.join()

    print(out[0].head())
    print(out[0].describe())
