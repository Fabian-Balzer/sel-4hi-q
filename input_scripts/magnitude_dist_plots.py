# -*- coding: utf-8 -*-
"""
Created on Sat Feb 19 17:37:14 2022

@author: fabian_balzer

Script for plotting the spectroscopic availability of the matched catalogue
"""

# %%
from math import ceil, floor

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import util.configure_matplotlib as cm
import util.my_tools as mt
from matplotlib.patches import Patch


def plot_magnitude_dist(df, band_list=tuple(mt.BAND_LIST), split_ttype=False):
    num_plots = len(band_list)
    num_rows = 4
    num_cols = ceil(num_plots / num_rows)
    colors = plt.cm.viridis_r.colors[::int(256 / num_plots)]
    fig, axes = plt.subplots(
        num_rows, num_cols, figsize=cm.set_figsize(subplots=(num_rows, num_cols)), sharex=True, sharey=True)
    fig.suptitle(
        "Magnitude distribution for the different bands (bin width is 0.5 mag)")
    bin_heights = []
    plt.subplots_adjust(left=None, bottom=None, right=None,
                        top=None, wspace=0.1, hspace=0.1)
    plt.rcParams['xtick.labeltop'] = False
    for i, band in enumerate(band_list):
        row = floor(i / num_cols)
        col = i - row * num_cols
        ax = axes[row][col]
        df[f"mag_{band}"].hist(ax=ax, label=(band.replace("_", "-")),
                               color=colors[i], bins=np.arange(10, 25, 0.5))
        ax.legend(loc="upper left")
        if row == num_rows - 1:
            ax.set_xlabel("Mag")
        ax.set_xlim(10, 24)
        ax.set_xticks(range(10, 24 + 1, 2))
        for p in ax.patches:
            bin_heights.append(p.get_height())

    for ax in axes.reshape(-1):
        max_ = int(ceil(max(bin_heights) / 100.0)) * 100
        ax.set_ylim(0, max_)
        ax.set_yticks(range(0, max_ + 1, int(max_ / 4)))

    # Disable additional axes
    for i in range(num_plots, (num_rows * num_cols)):
        row = floor(i / num_cols)
        col = i - row * num_cols
        print(row, col)
        axes[row][col].axis('off')
    cm.save_figure(fig, "input_analysis/magnitude_distribution")



if __name__ == "__main__":
    df = mt.read_plike_and_ext(prefix="matches/test2_",
                               suffix="_processed_table.fits")
    df = mt.add_mag_columns(df)
    # plot_magnitude_dist(df)