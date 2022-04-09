# -*- coding: utf-8 -*-
"""
Created on Sat Feb 19 17:37:14 2022

@author: fabian_balzer

Script for plotting the spectroscopic availability of the matched catalogue
"""
# %%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import util.configure_matplotlib as cm
import util.my_tools as mt
from matplotlib.patches import Patch


def plot_r_band_magnitude(df):
    """Plots the number of sources for multiple r-band bins and shows how much spec-z is available for each."""
    fig, axes = plt.subplots(1, 1, figsize=cm.set_figsize(fraction=.5))
    with_z = df[df["ZSPEC"] > 0]
    without_z = df[~(df["ZSPEC"] > 0)]
    axes.hist([with_z["mag_r"], without_z["mag_r"]],
              bins=30, stacked=True, color=['orange', 'b'], label=["spec-z available", "no spec-z"], histtype="stepfilled")
    axes.set_xlabel("r band magnitude")
    axes.legend(loc="upper left")
    axes.set_title(
        f"r-band distribution ({len(df[df['mag_r']>0])} sources, eFEDS field)", size="small")
    cm.save_figure(fig, "input_analysis/r_band_hist")


def plot_input_distribution(df):
    fig, axes = plt.subplots(1, 1, figsize=cm.set_figsize(fraction=.8))
    total_count = len(df)
    # Count the number of available photometry for each band:
    counts = {}
    for ttype in ["extended", "pointlike"]:
        ttypedf = df[df["Type"] == ttype]
        total_length = len(ttypedf)
        ttypedf = ttypedf.rename(columns={"ZSPEC": "spec-z"})
        for band in mt.BAND_LIST + ["spec-z"]:
            refname = f"mag_{band}" if band != "spec-z" else band
            subset = ttypedf[ttypedf[refname] > 0]
            band = band.replace('_', '-')
            if band in counts:
                counts[band][ttype] = len(subset) / total_length
            else:
                counts[band] = {ttype: len(subset) / total_length}
    count_df = pd.DataFrame(counts).T
    colors = {"galex": "blue", "sweep": "green",
              "vhs": "red", "hsc": "lightgreen", "kids": "darkgreen"}
    color_dict = {band.replace("_", "-"): color for key, color in colors.items()
                  for band in mt.BAND_DICT[key]}
    color_dict["spec-z"] = "black"
    count_df["Color"] = pd.DataFrame.from_dict(color_dict, orient="index")
    axes.grid(True, axis="y")
    labels = count_df.index
    x = np.arange(len(labels))
    width = 0.4  # the width of the bars
    space = 0
    rects1 = axes.bar(
        x - width / 2 - space, count_df["extended"], width, label='Extended', color=count_df["Color"], edgecolor="k", alpha=0.7)
    rects2 = axes.bar(
        x + width / 2 + space, count_df["pointlike"], width, label='Pointlike', color=count_df["Color"], edgecolor="k")
    for bar in rects1:
        bar.set_linewidth(0.5)
        bar.set_linestyle("dashed")
    for bar in rects2:
        bar.set_linewidth(0.5)
        bar.set_hatch("///")

    # Add some text for labels, title and custom x-axis tick labels, etc.
    axes.set_ylabel('Relative availability')
    title = f"Photometry in the eFEDS field  ({total_count} in total)"
    axes.set_title(title, size="small")
    axes.set_xticks(x)
    axes.set_xticklabels(list(labels), minor=False, rotation=90)
    axes.axhline(1, color="k")
    legend_patches = [Patch(facecolor=color, edgecolor='k',
                            label=mt.SURVEY_NAME_DICT[key]) for key, color in colors.items()]
    ext_plike_patches = [Patch(facecolor="white", edgecolor='k',
                               label=label, linestyle=lstyle, hatch=hatch) for label, lstyle, hatch in [("Extended", "dashed", ""), ("Pointlike", "-", "///")]]
    legend_2 = plt.legend(handles=ext_plike_patches, prop={
        "size": "x-small"}, bbox_to_anchor=(1, 0.7), loc=2)
    fig.add_artist(legend_2)
    fig.legend(handles=legend_patches, prop={
        "size": "x-small"}, bbox_to_anchor=(0.9, 0.9), loc=2)
    cm.save_figure(fig, "input_analysis/input_distribution")
    return count_df


if __name__ == "__main__":
    df = mt.read_plike_and_ext(prefix="matches/test2_",
                               suffix="_processed_table.fits")
    df = mt.add_mag_columns(df)
    plot_r_band_magnitude(df)
    count_df1 = plot_input_distribution(df)
