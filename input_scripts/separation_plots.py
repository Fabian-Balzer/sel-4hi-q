# -*- coding: utf-8 -*-
"""
Created on Sat Feb 19 17:37:14 2022

@author: fabian_balzer

Script for plotting the separations of the matched catalogue
"""
# %%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import util.configure_matplotlib as cm
import util.my_tools as mt
from matplotlib.offsetbox import AnchoredText
from scipy.stats import gaussian_kde, rayleigh

NBINS = 30


class ColumnMissingError(Exception):
    pass


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


def scatter_points(x, y, ax, lim, true_radius=None):
    """Produce the central scatter plot, including main lines and circles corresponding to the data used."""
    xy = np.vstack([x, y])
    z = gaussian_kde(xy)(xy)
    ax.scatter(x, y, label=f"{len(x)} sources", c=z,
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
    if true_radius is not None:
        circle2 = plt.Circle((xmed, ymed), true_radius, color='r',
                             linestyle="dashed", alpha=0.85, fill=False, label="$r_{3\sigma}$")
        ax.add_patch(circle2)
    ax.legend(prop={"size": "x-small"})
    # fig.suptitle(f"Match separation of {len(df)} sources between {ttype} and VHS data in eFEDS")
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


def plot_separation_hist(ax, histdata, lim, draw_cutoff):
    binwidth = lim / NBINS
    bins = np.arange(0, lim + binwidth, binwidth)
    ax.hist(histdata, bins=bins, fc="white",
            ec="k", histtype="stepfilled", density=True)
    ax.set_xlim(0, lim)
    ax.set_xlabel(f"{'Corrected s' if draw_cutoff else 'S'}eparation [arcsec]")
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
    if draw_cutoff:
        ax.axvline(3 * std, color='r', linewidth=2, label=r"$3\sigma$ cut")
    ax.legend(prop={"size": "x-small"})


def add_cols_if_necessary(df, ttype):
    """Checks wether the given columns are already present."""
    cols = df.columns
    if f"ra_{ttype}" not in cols:
        raise ColumnMissingError
    all_columns_found = f"true_sep_{ttype}" in cols
    for col in ["ra", "dec"]:
        colname = f"delta_{col}_{ttype}"
        if not colname in cols:
            df[colname] = df[f"{col}_{ttype}"] - df[col]
    colname = f"sep_to_{ttype}"
    if not colname in cols:
        df[colname] = (df[f"delta_ra_{ttype}"] **
                       2 + df[f"delta_dec_{ttype}"]**2)**0.5
    return df, all_columns_found


def plot_separation(df, ttype, radius, verbose=False):
    fig, axes = plt.subplots(2, 2, figsize=cm.set_figsize(
        fraction=.5, subplots=(2, 2), aspect=True))
    fig.patch.set_facecolor('white')
    main, xhist, yhist = arrange_subplots(axes)
    try:
        df, all_columns_found = add_cols_if_necessary(df, ttype)
    except ColumnMissingError:
        print(
            f"Could not find fitting columns to compare to for '{ttype}'. No plots are produced for it.")
        return
    df = df[df[f"delta_ra_{ttype}"].notnull()]
    if verbose:
        print(
            f"In total, there are {len(df)} matched sources within a matching radius of {radius}'' for {ttype}.")
    x_data = df[f"delta_ra_{ttype}"] * 3600
    y_data = df[f"delta_dec_{ttype}"] * 3600
    if all_columns_found:
        true_sep = df[f"true_sep_{ttype}"] * 3600
        true_radius = 3 * true_sep.std()
        # Calculate the point density:
        scatter_points(x_data, y_data, main, radius, true_radius)
        if verbose:
            num_good = len(true_sep[true_sep < true_radius])
            print(true_radius / 3)
            print(
                f"The radius including sources with less than 3 \
        \sigma separation from the center is {true_radius}'', \
        corresponding to {num_good} ({num_good/len(true_sep)*100:.1f} %) sources.")
    else:
        scatter_points(x_data, y_data, main, radius)
    # Look at the vertical and horizontal histograms
    produce_histograms(x_data, y_data, xhist, yhist, radius)
    cm.save_figure(
        fig, f"input_analysis/separation/separation_plot_{ttype}", format_="pdf")

    # Produce a histogram with the given data
    fig, axes = plt.subplots(1, 1, figsize=cm.set_figsize(fraction=.5))
    sep = df[f"sep_to_{ttype}"] if not all_columns_found else true_sep
    plot_separation_hist(axes, sep, radius, all_columns_found)
    cm.save_figure(
        fig, f"input_analysis/separation/separation_hist_{ttype}", format_="pdf")
    plt.show()


def plot_all_separations(df, verbose=False, names=None):
    """Produce separation plots for all columns"""
    radius_dict = {"opt_agn": 0.01, "vhs": 0.5,
                   "eros": 0.1, "hsc": 1, "galex": 3.5, "kids": 1.5}
    if names is None:
        names = (col for col in input_df.columns if col.startswith("ra"))
    for name in names:
        radius = radius_dict[names]
        plot_separation(df, name, radius, verbose=verbose)


if __name__ == "__main__":
    input_df = mt.read_plike_and_ext(prefix="matches/test2_",
                                     suffix="_processed_table.fits")
    input_df = mt.add_mag_columns(input_df)
    plot_all_separations(input_df, verbose=True)
