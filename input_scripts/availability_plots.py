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
from matplotlib.patches import Patch

import configure_matplotlib as cm
import my_tools as mt


def plot_r_band_magnitude(df):
    fig, axes = plt.subplots(1, 1, figsize=cm.set_figsize(fraction=.5))
    with_z = df[df["ZSPEC"] > 0]
    without_z = df[~(df["ZSPEC"] > 0)]
    axes.hist([with_z["MAG_r"], without_z["MAG_r"]],
              bins=30, stacked=True, color=['orange', 'b'], label=["spec-z available", "no spec-z"], histtype="stepfilled")
    axes.set_xlabel("r band magnitude")
    axes.legend(loc="upper left")
    axes.set_title(
        f"r-band distribution ({len(df[df['MAG_r']>0])} sources, eFEDS field)", size="small")
    cm.save_figure(fig, "data/plots/r_band_hist", format_="pdf")


def plot_input_distribution(df):
    fig, axes = plt.subplots(1, 1, figsize=cm.set_figsize(fraction=.8))
    total_count = len(df)
    # Count the number of available photometry for each band:
    counts = {}
    for type_ in ["extended", "pointlike"]:
        type_df = df[df["Type"] == type_]
        total_length = len(type_df)
        type_df = type_df.rename(columns={"ZSPEC": "spec-z"})
        for band in mt.BAND_LIST + ["spec-z"]:
            subset = type_df[type_df[band] > 0]
            if band in counts:
                counts[f"{band}"][type_] = len(subset) / total_length
            else:
                counts[f"{band}"] = {type_: len(subset) / total_length}
    my_df = pd.DataFrame(counts).T
    colors = {"galex": "blue", "sweep": "green",
              "vhs": "red", "hsc": "lightgreen"}
    color_dict = {band: color for key, color in colors.items()
                  for band in mt.BAND_DICT[key]}
    color_dict["spec-z"] = "black"
    my_df["Color"] = pd.DataFrame.from_dict(color_dict, orient="index")
    axes.grid(True, axis="y")
    labels = my_df.index
    x = np.arange(len(labels))
    width = 0.4  # the width of the bars
    space = 0
    rects1 = axes.bar(
        x - width / 2 - space, my_df["extended"], width, label='Extended', color=my_df["Color"], edgecolor="k", alpha=0.7)
    rects2 = axes.bar(
        x + width / 2 + space, my_df["pointlike"], width, label='Pointlike', color=my_df["Color"], edgecolor="k")
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
        "size": "x-small"}, bbox_to_anchor=(1, 0.75), loc=2)
    fig.add_artist(legend_2)
    fig.legend(handles=legend_patches, prop={
        "size": "x-small"}, bbox_to_anchor=(0.9, 0.9), loc=2)
    cm.save_figure(fig, "data/plots/input_distribution")
    return my_df


if __name__ == "__main__":
    df = mt.read_plike_and_ext("data/input/plike_processed_input.fits",
                               "data/input/ext_processed_input.fits")
    weird_cols = [col for col in df.columns if col.startswith(
        "c_") and "flux" in col and "kids" not in col]
    renamed = {col: col[2:].split("_flux")[
        0] + "_err" if "err" in col else col[2:].split("_flux")[0] for col in weird_cols}
    df = df.rename(columns=renamed)
    df = mt.add_mag_columns(df)
    # plot_r_band_magnitude(df)
    my_df = plot_input_distribution(df)
