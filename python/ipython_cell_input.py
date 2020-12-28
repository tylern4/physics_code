import h10
import boost_histogram as bh
import numpy as np
from tqdm import tqdm

def good_event(event):
    if event.dc_sect.size == 0 or event.dc_vz.size == 0:
        return False
    if event.p.size == 0:
        return False
    
    return True
    

root_reader = h10.h10_data()
root_reader.add("/Users/tylern/Data/e1d/test/*.root")
root_reader.reset()
vz = dict()
vz_vs_phi = dict()

vz_vs_phi_all = bh.Histogram(bh.axis.Regular(100, -6, 10, metadata="Vertex_z"),
                                 bh.axis.Regular(
                                     100, -np.pi, np.pi, metadata="phi"),
                                 bh.axis.IntCategory([1,2,3,4,5,6], metadata="sector")
                                 )

for i in range(1, 7):
    vz[i] = bh.Histogram(bh.axis.Regular(100, -6, 10))
    vz_vs_phi[i] = bh.Histogram(bh.axis.Regular(100, -6, 10), 
                                bh.axis.Regular(100, -np.pi, np.pi))


with tqdm(total=root_reader.num_entries) as pbar:
    for event in root_reader:
        pbar.update(1)
        if not good_event(event):
            continue

        sector = event.dc_sect[0]
        
        vz[sector].fill(event.dc_vz[0])

        phi = np.arctan2(event.cx[0], event.cy[0])
        vz_vs_phi[sector].fill(event.dc_vz[0], phi)
        vz_vs_phi_all.fill(event.dc_vz[0], phi, int(sector))
