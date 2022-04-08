# -*- coding: utf-8 -*-
"""
Created on Sat Feb 19 17:37:14 2022

@author: fabian_balzer

Script for plotting the separations of the matched catalogue
"""
# %%

import os

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.offsetbox import AnchoredText
from scipy.stats import gaussian_kde, rayleigh

import configure_matplotlib as cm
import my_tools as mt

NBINS = 30


def arrange_subplots(axes):
    """Arranges the subplots of the given axes object
    and returns the important ones (main, xhist, yhist) for later use."""
    left, width, = 0.1, 0.6
    bottom, height = 0.1, 0.6
    spacing, sidesize = 0.000, 0.2
    xstart, ystart = bottom + height + spacing, left + width + spacing
    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, xstart, width, sidesize]
    rect_histy = [ystart, bottom, sidesize, height]

    main = axes[1][0]
    xhist = axes[0][0]
    yhist = axes[1][1]
    axes[0][1].axis("off")
    main.set_position(rect_scatter)
    xhist.set_position(rect_histx)
    yhist.set_position(rect_histy)
    main.tick_params(direction='in', top=True, right=True)
    xhist.tick_params(direction='in', labelbottom=False)
    yhist.tick_params(direction='in', labelleft=False)
    for ax in [main, xhist, yhist]:
        ax.grid(True)
        ax.minorticks_on()
    return main, xhist, yhist


def scatter_points(x, y, ax, lim, true_radius):
    """Produce the central scatter plot, including main lines and circles corresponding to the data used."""
    xy = np.vstack([x, y])
    z = gaussian_kde(xy)(xy)
    ax.scatter(x, y, label=f"{len(df)} sources", c=z,
               marker=".", cmap="viridis", s=1)
    ax.set_xlabel("$\Delta$ RA [arcsec]")
    ax.set_ylabel("$\Delta$ Dec [arcsec]", labelpad=0.1)
    ax.set_xlim((-lim, lim))
    ax.set_ylim((-lim, lim))
    ax.axvline(0, -1, 1, color="gray", linestyle="dashed", alpha=0.5)
    ax.axhline(0, -1, 1, color="gray", linestyle="dashed", alpha=0.5)
    # Encircle the data with the maximum matching radius:
    circle1 = plt.Circle((0, 0), lim, color='gray',
                         linestyle="dashed", alpha=0.5, fill=False, label="$r_{match}$")
    ax.add_patch(circle1)
    # Encircle the data that is actually used:
    xmed, ymed = x.median(), y.median()
    circle2 = plt.Circle((xmed, ymed), true_radius, color='r',
                         linestyle="dashed", alpha=0.85, fill=False, label="$r_{3\sigma}$")
    ax.add_patch(circle2)
    ax.legend(prop={"size": "x-small"})
    # fig.suptitle(f"Match separation of {len(df)} sources between {type_} and VHS data in eFEDS")
    # box = dict(boxstyle="round, pad=0.1", fc="gray", ec="k", lw=1, alpha=0.3)


def produce_histograms(x, y, x_ax, y_ax, lim):
    """Produces the histograms for the x and y axes including
    annotations for the medians and standard deviations."""
    binwidth = 2 * lim / NBINS
    bins = np.arange(-lim, lim + binwidth, binwidth)
    yvals1, _, _ = x_ax.hist(
        x, bins=bins, fc="white", ec="k", histtype="stepfilled")
    yvals2, _, _ = y_ax.hist(y, bins=bins, orientation='horizontal',
                             fc="white", ec="k", histtype="stepfilled")
    maximum = max(max(yvals1), max(yvals2)) * 1.05
    x_ax.set_ylim(0, maximum)
    y_ax.set_xlim(0, maximum)
    x_ax.set_xlim((-lim, lim))
    y_ax.set_ylim((-lim, lim))
    ticks = x_ax.get_yticks()[:-1]
    y_ax.set_xticks(ticks)
    x_ax.set_yticks(ticks)
    # Annotate the medians and standard deviations:
    xmed, ymed = x.median(), y.median()
    x_ax.axvline(xmed, color='k', linewidth=2)
    y_ax.axhline(ymed, color='k', linewidth=2)
    annotation = f"Median:\n  ${xmed:.3f}''$\n"
    annotation += "$\sigma=" + f"{x.std():.3f}''$"
    at = AnchoredText(annotation, prop=dict(size="x-small"),
                      frameon=True, loc="upper left")
    at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
    x_ax.add_artist(at)
    annotation = f"Median:\n  ${ymed:.3f}''$\n"
    annotation += "$\sigma=" + f"{y.std():.3f}''$"
    at = AnchoredText(annotation, prop=dict(size="x-small"),
                      frameon=True, loc="upper left")
    at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
    y_ax.add_artist(at)


def plot_separation_hist(ax, histdata, lim):
    binwidth = lim / NBINS
    bins = np.arange(0, lim + binwidth, binwidth)
    ax.hist(histdata, bins=bins, fc="white",
            ec="k", histtype="stepfilled", density=True)
    ax.set_xlim(0, lim)
    ax.set_xlabel("Corrected separation [arcsec]")
    ax.grid(True)
    ax.minorticks_on()
    med, std = histdata.median(), histdata.std()
    # annotation = f"Median:\n  ${med:.3f}''$\n"
    # annotation += "$\sigma=" + f"{std:.3f}''$"
    # at = AnchoredText(annotation, prop=dict(size="x-small"),
    #                   frameon=True, loc="upper right")
    # at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
    # ax.add_artist(at)
    ax.axvline(med, color='k', linewidth=2, label=f"Median ({med:.3f}'')")

    x = np.linspace(0, lim, 100)
    ax.plot(x, rayleigh.pdf(x, scale=std),
            label="Rayleigh dist.")
    ax.axvline(3 * std, color='r', linewidth=2, label=r"$3\sigma$ cut")
    ax.legend(prop={"size": "x-small"})


def plot_separation(df, type_, radius):
    fig, axes = plt.subplots(2, 2, figsize=cm.set_figsize(
        fraction=.5, subplots=(2, 2), aspect=True))
    fig.patch.set_facecolor('white')
    main, xhist, yhist = arrange_subplots(axes)
    df = df[df[f"delta_ra_{type_}"].notnull()]
    print(
        f"In total, there are {len(df)} matched sources within a matching radius of {radius}'' for {type_}.")
    x_data = df[f"delta_ra_{type_}"] * 3600
    y_data = df[f"delta_dec_{type_}"] * 3600
    true_sep = df[f"true_sep_{type_}"] * 3600
    true_radius = 3 * true_sep.std()
    # Calculate the point density:
    scatter_points(x_data, y_data, main, radius, true_radius)
    # Look at the histograms
    produce_histograms(x_data, y_data, xhist, yhist, radius)
    num_good = len(true_sep[true_sep < true_radius])
    print(true_radius / 3)
    print(
        f"The radius including sources with less than 3 \
\sigma separation from the center is {true_radius}'', \
corresponding to {num_good} ({num_good/len(true_sep)*100:.1f} %) sources.")
    cm.save_figure(fig, f"Data/separation_plot_{type_}", format_="pdf")

    # Produce a histogram with the given data
    fig, axes = plt.subplots(1, 1, figsize=cm.set_figsize(fraction=.5))
    plot_separation_hist(axes, true_sep, radius)
    cm.save_figure(fig, f"data/plots/separation_hist_{type_}", format_="pdf")
    plt.show()


if __name__ == "__main__":
    df1 = mt.read_fits_as_dataframe("data/input/plike_processed_input.fits")
    df2 = mt.read_fits_as_dataframe("data/input/ext_processed_input.fits")
    df = pd.concat([df1, df2])
    # plot_separation(df, "shu", 0.01)  # As we can see, no separation.
    # plot_separation(df, "vhs", 0.5)
    # plot_separation(df, "galex", 3.5)
    plot_separation(df, "eros", 3.5)
    # plot_separation(df, "hsc", 1)
    # plot_separation(df, "eros", 0.01, histplot=True)
