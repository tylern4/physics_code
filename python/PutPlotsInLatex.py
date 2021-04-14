
from glob import glob
import numpy as np

w_bins_e99 = np.array([1.1, 1.12, 1.14, 1.16, 1.18, 1.2, 1.22, 1.24, 1.26, 1.28, 1.3,
                       1.32, 1.34, 1.36, 1.38, 1.4, 1.42, 1.44, 1.46, 1.48, 1.5, 1.52,
                       1.54, 1.56, 1.58, 1.6, 1.62, 1.64, 1.66, 1.68, 1.7, 1.72, 1.74, 1.76, 1.78, 1.8, 1.82])

w_bins_k = np.array([1.605, 1.615, 1.625, 1.635, 1.645, 1.655, 1.665, 1.675, 1.685, 1.695,
                     1.705, 1.715, 1.725, 1.735, 1.745, 1.755, 1.765, 1.775, 1.785, 1.795, 1.805])

q2_bins_e99 = np.array([1.1, 1.33, 1.56, 1.87, 2.23, 2.66, 3.5])
q2_bins_k = np.array([1.8,  2.2,  2.6,  3.15, 4.0])

INSERT_FIG = """\\begin{{minipage}}{{.45\\textwidth}}
\\centering
\\includegraphics[width=\\linewidth]{{plots/CrossSections/{}}}
\\captionsetup{{width=0.8\\textwidth, font=scriptsize, list=no}}
\\captionof{{figure}}{{{}}}
\\label{{fig:xs{}}}
\\end{{minipage}}"""


# INSERT_FIG = "\\begin{{figure}}\n\\centering\n\\includegraphics[width=12cm]{{plots/CrossSections/{}}}\n\\caption{{ {} }}\n\\end{{figure}}\n"

CAPTION = "Cross-Sections $W~[{},{})~\\mathrm{{GeV}}~Q^2~[{},{})~\\mathrm{{GeV}}^2$."

files = glob("/Users/tylern/Desktop/show/plots/e1d/plots/crossSections/*.png")

for i, f in enumerate(np.sort(files)):
    x = f.split("/")
    name = x[-1]
    W = float(name[:-4].split("_")[2])
    Q2 = float(name[:-4].split("_")[4])
    try:
        W_min = w_bins_e99[np.where(w_bins_e99 == W)][0]
        W_max = w_bins_e99[np.where(w_bins_e99 == W)[0] + 1][0]

        Q2_min = q2_bins_e99[np.where(q2_bins_e99 == Q2)][0]
        Q2_max = q2_bins_e99[np.where(q2_bins_e99 == Q2)[0] + 1][0]

    except IndexError:
        W_min = w_bins_k[np.where(w_bins_k == W)][0]
        W_max = w_bins_k[np.where(w_bins_k == W)[0] + 1][0]

        Q2_min = q2_bins_k[np.where(q2_bins_k == Q2)][0]
        Q2_max = q2_bins_k[np.where(q2_bins_k == Q2)[0] + 1][0]

    caption = CAPTION.format(W_min, W_max, Q2_min, Q2_max)
    # if i > 11:
    #     continue

    print(INSERT_FIG.format(name, caption, i))
    if i % 2 != 0:
        print()
