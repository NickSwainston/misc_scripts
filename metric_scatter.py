#!/usr/bin/env python

# This script plots all ionospheric points with a median offset < 0.5 arcmin.

import os
import sys
import argparse
import itertools

import numpy as np
import pandas as pd
from astropy.time import Time

import matplotlib as mpl
mpl.rcParams["ps.useafm"] = True
mpl.rcParams["pdf.use14corefonts"] = True
mpl.rcParams["font.size"] = 10
mpl.rcParams["text.usetex"] = True
mpl.rcParams["font.serif"] = "Palatino"
from matplotlib.ticker import NullFormatter

import seaborn as sns
sns.set(color_codes=True)
sns.set_style({"font.family": [u"serif"]}, {"legend.frameon": True})
sns.set_context("paper", font_scale=1.5)


type_1_m_limits = (0, 0.14)
type_1_p_limits = (0.5, 0.63)

type_2_m_limits = (0.14, 0.5)
type_2_p_limits = (0.5, 0.7)

type_3_m_limits = (0, 0.14)
type_3_p_limits = (0.63, np.inf)

type_4_m_limits = (0.14, 0.5)
type_4_p_limits = (0.7, np.inf)

palette = itertools.cycle(sns.color_palette())
colours = []
# Preserving the colours used from Paper I.
colours.append(next(palette))
next(palette)
colours.append(next(palette))
colours.append(next(palette))
next(palette)
next(palette)
next(palette)
next(palette)
colours.append(next(palette))


def print_number_with_comma_and_curly_brackets(n):
    print("{", "{:,}".format(n), "}", sep="")


def do_the_plot(df, marker_size=2, plot=True, output=None):
    mag = df.iono_magnitude
    pca = df.iono_pca

    median = np.median(mag)
    p = np.percentile(mag, 95)
    print(f"median: {median}")
    print(f"90%: {p}")
    sys.exit(0)


    type_1 = df[(mag > type_1_m_limits[0]) & (mag <= type_1_m_limits[1]) &
                (pca > type_1_p_limits[0]) & (pca <= type_1_p_limits[1])]
    type_2 = df[(mag > type_2_m_limits[0]) & (mag <= type_2_m_limits[1]) &
                (pca > type_2_p_limits[0]) & (pca <= type_2_p_limits[1])]
    true_type_2 = df[(mag > type_2_m_limits[0]) &
                     (pca > type_2_p_limits[0]) & (pca <= type_2_p_limits[1])]
    type_3 = df[(mag > type_3_m_limits[0]) & (mag <= type_3_m_limits[1]) &
                (pca > type_3_p_limits[0]) & (pca <= type_3_p_limits[1])]
    type_4 = df[(mag > type_4_m_limits[0]) & (mag <= type_4_m_limits[1]) &
                (pca > type_4_p_limits[0]) & (pca <= type_4_p_limits[1])]
    true_type_4 = df[(mag > type_4_m_limits[0]) & (pca > type_4_p_limits[0])]
    total = len(type_1) + len(type_2) + len(type_3) + len(type_4)
    true_total = len(type_1) + len(true_type_2) + len(type_3) + len(true_type_4)

    num_nights = len(np.where(np.diff(df.obsid) > 30000)[0])

    # Stuff to copy+paste into the latex source.
    print("\n\\newcommand{\\noNights}", end="")
    print_number_with_comma_and_curly_brackets(num_nights)
    print("\\newcommand{\\noObs}", end="")
    print_number_with_comma_and_curly_brackets(true_total)
    print("\\newcommand{\\noObsInScatter}", end="")
    print_number_with_comma_and_curly_brackets(total)
    print("\\newcommand{\\noObsDiff}", end="")
    print_number_with_comma_and_curly_brackets(true_total - total)
    print("\\newcommand{\\firstObs}{%s}" % Time(df.iloc[0].obsid, scale="utc", format="gps").datetime.strftime("%B %Y"))
    print("\\newcommand{\\lastObs}{%s}" % Time(df.iloc[-1].obsid, scale="utc", format="gps").datetime.strftime("%B %Y"))
    for (t, l) in ([(len(type_1), "tone"),
                    (len(type_2), "ttwo"),
                    (len(type_3), "tthree"),
                    (len(true_type_4), "tfour"),
                    ]):
        print("\\newcommand{\\%s}" % l, end="")
        print_number_with_comma_and_curly_brackets(t)
    for (t, l) in ([(len(type_1), "tonef"),
                    (len(type_2), "ttwof"),
                    (len(type_3), "tthreef"),
                    (len(true_type_4), "tfourf"),
                    ]):
        print("\\newcommand{\\%s}{%.3g}" % (l, 100 * t / true_total))
    print("\n")

    # Plotting.
    pad = 0.1
    left_edge = 0.30
    bot_edge = 0.30
    width = 1 - left_edge - 0.01
    height = 1 - bot_edge - 0.01

    rect_scatter = (left_edge, bot_edge, width, height)
    rect_hist_x = (left_edge, pad, width, bot_edge - pad)
    rect_hist_y = (pad, bot_edge, left_edge - pad, height)

    fig = plt.figure(figsize=(12, 9))
    ax_scatter = plt.axes(rect_scatter)
    ax_hist_x = plt.axes(rect_hist_x)
    ax_hist_y = plt.axes(rect_hist_y)

    # Make all subplot edges thicker.
    for s in ["top", "bottom", "left", "right"]:
        for a in [ax_scatter, ax_hist_x, ax_hist_y]:
            a.spines[s].set_linewidth(3)
            a.spines[s].set_edgecolor("black")

    # Remove labels.
    nullfmt = NullFormatter()
    ax_scatter.xaxis.set_major_formatter(nullfmt)
    ax_scatter.yaxis.set_major_formatter(nullfmt)

    print(f"Total: {total}")

    mag_data = []
    pca_data = []
    for i, (t, l) in enumerate([(type_1, "Type 1"),
                                (type_2, "Type 2"),
                                (type_3, "Type 3"),
                                (type_4, "Type 4"),
                                ]):
        print(f"{l}: {len(t)} ({len(t)/total})")
        ax_scatter.scatter(t.iono_magnitude, 100*t.iono_pca,
                           edgecolors="face", s=marker_size,
                           c=(colours[i],), label=l, alpha=0.8)
        mag_data += t.iono_magnitude.to_numpy().tolist()
        pca_data += (100*t.iono_pca.to_numpy()).tolist()

    levels = [0.1, 0.4, 0.7]
    print("KDE levels: %s" % levels)
    plt.sca(ax_scatter)
    sns.kdeplot(mag_data, pca_data, levels=levels, linewidths=1, colors="black", cmap=None)

    # weights = [np.ones(x)/total for x in (len(type_1), len(type_2), len(type_3), len(type_4))]
    f = lambda d: d.iono_magnitude
    ax_hist_x.hist((f(type_1), f(type_2), f(type_3), f(type_4)), bins=np.arange(0, 0.5 + 0.01, 0.005),
                   stacked=True, density=True, color=colours)
    f = lambda d: d.iono_pca * 100
    ax_hist_y.hist((f(type_1), f(type_2), f(type_3), f(type_4)), bins=np.arange(50, 101, 1),
                   orientation="horizontal", stacked=True, density=True, color=colours)

    # Match the limits of the small subplots to the big one.
    ax_hist_x.set_xlim(ax_scatter.get_xlim())
    ax_hist_y.set_ylim(ax_scatter.get_ylim())

    # Make the extent of the small subplots a little bigger.
    ax_hist_x.set_ylim((0, ax_hist_x.get_ylim()[1] * 1.05))
    ax_hist_y.set_xlim((0, ax_hist_y.get_xlim()[1] * 1.05))

    leg = ax_scatter.legend(loc=4, markerscale=7, frameon=True)
    leg.get_frame().set_facecolor("#f8f8f8")
    ax_hist_x.set_xlabel("median(ionospheric offsets) [arcmin at 200 MHz]")
    ax_hist_x.set_ylabel("PDF of data")
    ax_hist_y.set_xlabel("PDF of data")
    ax_hist_y.set_ylabel("Value of the dominant eigenvalue\ndetermined by PCA [per cent]")
    plt.tight_layout()

    # Remove the last tick labels.

    # When using an Agg backend, this line will complain that it cannot show the
    # figure. BUT, not using this line will screw up the ticks :))))
    fig.show()
    fig.canvas.draw()
    ticks = ax_hist_x.yaxis.get_major_ticks()
    # I don't know why -2 is what -1 should be. Thanks, python! I can't wait for
    # you to die in a fire.
    # Commented out on 2021-01-12, because we can fit another tick in.
    # ticks[-2].label1.set_visible(False)

    ticks = ax_hist_y.xaxis.get_major_ticks()
    # ticks[-2].label1.set_visible(False)

    # Shift the "PDF of data" labels so they're not too close to each other.
    ax_hist_x.yaxis.set_label_coords(-0.06, 0.25)
    ax_hist_y.xaxis.set_label_coords(0.25, -0.05)

    if plot:
        plt.show()
    else:
        if output is None:
            output = args.output
        plt.savefig(output)
        print("Wrote:", output)
    plt.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--plot", action="store_true",
                        help="Plot interactively.")
    parser.add_argument("-o", "--output",
                        default="metric_scatter.pdf",
                        help="Output filename. Default: %(default)s")
    parser.add_argument("files", nargs='*', help="Files to be processed.")
    args = parser.parse_args()

    if not args.plot:
        mpl.use("Agg")
    import matplotlib.pyplot as plt

    df = pd.read_csv(os.path.join(os.path.dirname(__file__), "data/qadb.csv")).sort_values(by=["obsid"])
    do_the_plot(df, plot=args.plot, output=args.output)
