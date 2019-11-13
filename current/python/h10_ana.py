from matplotlib import pyplot as plt
from build.h10 import h10_data
from build.physics_vectors import LorentzVector
from python.reaction import reaction
from python.histograms import Hist1D, Hist2D

from tqdm import tqdm
import numpy as np


h1d = Hist1D(500, 0.5, 2)


data = h10_data()
data.add("~/Data/e1d/data/h10_r23501*.root")

for e in tqdm(range(0, data.num_entries)):
    data.get_entry(e)
    event = reaction(data)
    event.run()
    if event.PROT_EVENT:
        h1d.fill(event.W)

plt.step(*h1d.data)
plt.show()
