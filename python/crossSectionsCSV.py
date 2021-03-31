#!/usr/bin/env python
import warnings  # noqa
warnings.filterwarnings("ignore")  # noqa

from calc_xsections import *
from tqdm import tqdm
import argparse
import os
import time


def main(rec, mc_rec, mc_thrown, empty, binning, out_folder="plots", bins=12, radcorr=None, csvName="results"):
    results = []
    if not os.path.exists(f'{out_folder}/crossSections'):
        os.makedirs(f'{out_folder}/crossSections')
    # Make a set of values from 0 to 2Pi for plotting
    xs = np.linspace(0, 2 * np.pi, 250)
    if radcorr is not None:
        radcorr_df = pd.read_csv(radcorr)

    for w in tqdm(binning["wbins"]):
        for q2 in binning["q2bins"]:
            radcor_R = 1.0
            if radcorr is not None:
                cut = isclose(radcorr_df.w_left, w.left) & isclose(
                    radcorr_df.q2_left, q2.left)

                if(radcorr_df[cut].R.size == 0):
                    print(w.left, q2.left)
                    radcor_R = 1.0
                else:
                    radcor_R = np.array(radcorr_df[cut].R)

            for theta in binning["thetabins"]:
                # Cut data/mc for the w/q2/theta bin we're in
                data = make_cuts(rec, w, q2, theta)  # rec[].copy()
                data_e = make_cuts(empty, w, q2, theta)  # empty[].copy()
                data_mc = make_cuts(mc_rec, w, q2, theta)  # mc_rec[].copy()
                # mc_thrown[].copy()
                thrown = make_cuts(mc_thrown, w, q2, theta)
                num_good = np.sum(data.cut_fid)

                binCenter = binCetnerCorrection(w, q2, theta, num_bins=bins)

                cut_fids = {
                    "E1D Data": 0,
                    # "All Data": 2,
                    # "Fid cuts False": 1,
                }
                for name, cuts in cut_fids.items():
                    if cuts == 0:
                        # Histogram the data for plotting
                        _data_y, _x = hist_data(
                            data[data.cut_fid], density=False, bins=bins)
                        _empty_y, _ = hist_data(
                            data_e[data_e.cut_fid], density=False, bins=bins)
                        _mc_rec_y, _ = hist_data(
                            data_mc[data_mc.cut_fid], density=False, bins=bins)
                    elif cuts == 1:
                        # Histogram the data for plotting
                        _data_y, _x = hist_data(
                            data[~data.cut_fid], density=False, bins=bins)
                        _empty_y, _ = hist_data(
                            data_e[~data_e.cut_fid], density=False, bins=bins)
                        _mc_rec_y, _ = hist_data(
                            data_mc[~data_mc.cut_fid], density=False, bins=bins)
                    else:
                        _data_y, _x = hist_data(
                            data, density=False, bins=bins)
                        _empty_y, _ = hist_data(
                            data_e, density=False, bins=bins)
                        _mc_rec_y, _ = hist_data(
                            data_mc, density=False, bins=bins)

                    _thrown_y, _ = hist_data(thrown, density=False, bins=bins)

                    # Get numbers for error
                    N_y = _data_y
                    N_empty = _empty_y
                    # Empty target subtraction
                    _data_y = (_data_y/Q_FULL - _empty_y/Q_EMPTY)
                    # Remove points with 0 data count
                    cut = ~(_data_y == 0)
                    x = _x[cut]
                    mc_rec_y = _mc_rec_y[cut]
                    thrown_y = _thrown_y[cut]
                    data_y = _data_y[cut]
                    N_y = N_y[cut]
                    N_empty = N_empty[cut]

                    # Remove points with 0 mc_rec count and put 1???
                    mc_rec_y = np.where(mc_rec_y == 0, 1, mc_rec_y)
                    thrown_y = np.where(thrown_y == 0, 1, thrown_y)
                    acceptance = mc_rec_y / thrown_y

                    # Get bin widths
                    delta_W = (w.right-w.left)
                    delta_Q2 = (q2.right-q2.left)
                    delta_Theta = np.abs(theta.right-theta.left)
                    __phis = np.linspace(0, 2 * np.pi, bins)
                    delta_phi = __phis[1] - __phis[0]
                    kin_bin_width = delta_W * delta_Q2 * delta_Theta * delta_phi

                    # Calculate acceptance and correct data
                    flux = virtual_photon_flux(w.mid, q2.mid) * luminosity()
                    denom = kin_bin_width * flux * \
                        acceptance * radcor_R * binCenter(x)

                    stat_error = statistical(N_y, N_empty, denom)

                    # Normalize with bin widths
                    try:
                        y = data_y / denom
                    except ValueError:
                        continue

                    error_bar = get_error_bars(
                        y, mc_rec_y, thrown_y, stat_error)

                    for phi, cross, err in zip(x, y, error_bar):
                        results.append({"w_left": w.left,
                                        "w_right": w.right,
                                        "q2_left": q2.left,
                                        "q2_right": q2.right,
                                        "cos_theta": theta.left,
                                        "x": np.round(phi, 3),
                                        "y": np.round(cross, 5),
                                        "err": np.round(err, 5)})

    output = pd.DataFrame(results)
    print(output.head())
    output.to_csv(f"{out_folder}/crossSections/{csvName}.csv")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Make Cross Sections")
    parser.add_argument("--mc", dest="mc_data_file_path",
                        type=str, help="MC csv file", required=True)
    parser.add_argument("--data", dest="rec_data_file_path",
                        type=str, help="Data csv file", required=True)
    parser.add_argument("--empty", dest="empty_file_path",
                        type=str, help="Empty run csv file", required=True)
    parser.add_argument("--folder", dest="out_folder", type=str,
                        help="Folder for plots", required=False, default="plots")
    parser.add_argument("--radcorr", dest="radcorr", type=str, help="Location of radcorr data csv", required=False,
                        default=None)
    parser.add_argument("--highw", help="Use high W binning",
                        required=False, action='store_true', default=False)
    args = parser.parse_args()

    # Start to main

    # print("Start setup")
    start = time.time_ns()
    # Load mc file
    mc_rec, mc_thrown = read_csv(args.mc_data_file_path)
    # Load reconstructed file
    rec = read_csv(args.rec_data_file_path, True)
    empty_target = read_csv(args.empty_file_path, True)

    # Cut for missing mass
    rec, mc_rec, empty_target = cut_for_MM(rec, mc_rec, empty_target)
    end = time.time_ns()

    if args.highw:
        w_bins = w_bins_k
        q2_bins = q2_bins_k
        bins = 10
        csvName = "results_highw"
    else:
        w_bins = w_bins_e99
        q2_bins = q2_bins_e99
        bins = 10
        csvName = "full_results"

    # Make bins in the dataframes from the bins above
    _rec = prep_for_ana(rec, w_bins, q2_bins, theta_bins)
    _empty_target = prep_for_ana(empty_target, w_bins, q2_bins, theta_bins)
    _mc_rec = prep_for_ana(mc_rec, w_bins, q2_bins, theta_bins)
    _mc_thrown = prep_for_ana(mc_thrown, w_bins, q2_bins, theta_bins)

    # Create dict of sorted bins
    _binning = dict()
    _binning["wbins"] = pd.Index.sort_values(pd.unique(_rec.w_bin))
    _binning["q2bins"] = pd.Index.sort_values(pd.unique(_rec.q2_bin))
    _binning["thetabins"] = pd.Index.sort_values(pd.unique(_rec.theta_bin))
    print(f"Done setup: {(end-start)/1E9:0.2f}Sec")

    main(_rec, _mc_rec, _mc_thrown, _empty_target, _binning, bins=bins,
         out_folder=args.out_folder, radcorr=args.radcorr, csvName=csvName)

    del _rec, _mc_rec, _mc_thrown, _empty_target
