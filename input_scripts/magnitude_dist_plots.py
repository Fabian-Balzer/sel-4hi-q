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


def plot_magnitude_dist(df: pd.DataFrame, band_list=tuple(mt.BAND_LIST), split_ttype=False):
    """Plot the distribution of sources for a given set of photometric bands.
        Parameters:
            [df]: Pandas DataFrame
                DataFrame containing the processed input data.
            [band_list]: tuple<str>
                Plots for these requested photometric bands will be produced.
            [context]: int = -1
                Only consider the bands fitting the context
    """
    num_plots = len(band_list)
    num_cols = min(floor(num_plots**0.5), 4)  # at most 4 cols
    num_rows = ceil(num_plots / num_cols)
    num_in_last_row = num_rows * num_cols - num_plots
    # colors = plt.cm.viridis_r.colors[::int(256 / num_plots)]
    fig, axes = plt.subplots(
        num_rows, num_cols, figsize=cm.set_figsize(subplots=(num_rows, num_cols)), sharex=True, sharey=True)
    fig.set_tight_layout(False)
    # fig.suptitle(
    #     "Relative magnitude distributions (bin width = 0.5 mag)")
    # bin_heights = []
    plt.subplots_adjust(left=None, bottom=None, right=None,
                        top=None, wspace=0.1, hspace=0.1)
    for i, band in enumerate(band_list):
        row = floor(i / num_cols)
        col = i - row * num_cols
        is_on_bottom = row == num_rows - 1 or row == num_rows - \
            2 and col >= num_cols - num_in_last_row
        ax = axes[row][col]
        if split_ttype:
            ext = df[df["Type"] == "extended"][f"mag_{band}"]
            plike = df[df["Type"] == "pointlike"][f"mag_{band}"]
            n_ext, n_plike = len(ext[ext > 0]), len(plike[plike > 0])
            ext.hist(ax=ax, label=f"Ext ({n_ext})",
                     bins=np.arange(10, 25, 0.5), density=True, histtype="stepfilled", alpha=0.5, ec="gray")
            plike.hist(ax=ax, label=f"Poi ({n_plike})",
                       bins=np.arange(10, 25, 0.5), density=True, histtype="stepfilled", alpha=0.5, ec="gray")
            title = mt.generate_pretty_band_name(band)
            surv_color = mt.give_survey_color(mt.give_survey_for_band(band))
            ax.text(0.1, 0.8, title, bbox=dict(
                facecolor=surv_color, alpha=0.3), transform=ax.transAxes)
            ax.legend(loc=(0.4, 0.6), prop={"size": "small"}, frameon=False)

        else:
            n_sources = mt.give_amount_of_good_photometry(df, band)
            label = f"{mt.generate_pretty_band_name(band)} ($N={n_sources}$)"
            surv_color = mt.give_survey_color(mt.give_survey_for_band(band))
            df[f"mag_{band}"].hist(ax=ax, label=label,
                                   color=surv_color, bins=np.arange(10, 25, 0.5), density=True, histtype="stepfilled", ec="black")
            ax.legend(loc="upper left")
        if is_on_bottom:
            ax.set_xlabel("Mag")
        ax.set_xlim(12, 24)
        ax.set_xticks(range(12, 24 + 1, 2))
        ax.set_xticklabels(range(12, 24 + 1, 2))
        ax.tick_params(axis="x", top=(row == 0),
                       bottom=is_on_bottom, label2On=(row == 0), label1On=is_on_bottom)
        # for patch in ax.patches:
        #     bin_heights.append(patch.get_height())

    for ax in axes.reshape(-1):
        max_ = 1  # int(ceil(max(bin_heights) / 100.0)) * 100
        ax.set_ylim(0, max_)
        # ax.set_yticks(range(0, max_ + 1, int(max_ / 4)))

    # Disable additional axes
    for i in range(num_plots, (num_rows * num_cols)):
        row = floor(i / num_cols)
        col = i - row * num_cols
        axes[row][col].axis('off')
        axes[row - 1][col].set_xticklabels(range(12, 24 + 1, 2))
    return fig
