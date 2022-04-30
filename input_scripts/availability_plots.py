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


def plot_r_band_magnitude(df, stem=""):
    """Plots the number of sources for multiple r-band bins and shows how much spec-z is available for each."""
    fig, ax = plt.subplots(1, 1, figsize=cm.set_figsize(fraction=.5))
    with_z = df[df["ZSPEC"] > 0]
    without_z = df[~(df["ZSPEC"] > 0)]
    ax.grid(True, axis="y")
    ax.hist([with_z["mag_r"], without_z["mag_r"]],
            bins=30, stacked=True, color=['k', 'g'], alpha=0.7, label=["With spec-$z$", "Without spec-$z$"], histtype="stepfilled", edgecolor="k")
    ax.legend(loc="upper left")
    ax.set_xlabel(f"{mt.BAND_LABEL_DICT['r']}-band magnitude")
    ax.set_ylabel("Number of sources")
    ax.set_title(
        f"{mt.BAND_LABEL_DICT['r']}-band distribution ({len(df[df['mag_r']>0])} sources, eFEDS field)", size="small")
    cm.save_figure(fig, f"input_analysis/{stem}_r_band_hist")


def construct_band_availability_dataframe(df, desired_bands):
    """Count the number of available photometry for each band and assign colours to them depending on the surveys they stem from.
    returns a band_count_df with the bands as an index."""
    counts, source_nums = {}, {}
    for ttype in ["extended", "pointlike"]:
        ttypedf = df[df["Type"] == ttype]
        total_length = len(ttypedf)
        for band in desired_bands:
            refname = f"mag_{band}" if band != "ZSPEC" else band
            subset = ttypedf[ttypedf[refname] > 0]
            if band in counts:
                counts[band][ttype] = len(subset) / total_length
            else:
                counts[band] = {ttype: len(subset) / total_length}
        source_nums[ttype] = total_length
    band_count_df = pd.DataFrame(counts).T
    # Add colours:
    colours = {"galex": "blue", "sweep": "green",
               "vhs": "red", "hsc": "lightgreen", "kids": "darkgreen"}
    colour_dict = {band.replace("_", "-"): colour for key, colour in colours.items()
                   for band in mt.BAND_DICT[key]}
    colour_dict["ZSPEC"] = "black"
    band_count_df["Colour"] = pd.DataFrame.from_dict(
        colour_dict, orient="index")
    return band_count_df, colours, source_nums


def construct_num_band_dataframe(df, allowed_bands):
    """Count how many bands are available for each source and return a dataframe reflecting this distribution."""
    max_band_num = len(allowed_bands)
    count_dict = {}
    # Counts the number of not-nan-entries for each row
    df["num_bands"] = df[allowed_bands].count(axis=1)
    for i in range(max_band_num + 1):
        with_z = df[df["ZSPEC"] > 0]
        without_z = df[~(df["ZSPEC"] > 0)]
        count_dict[i] = {}
        count_dict[i]["with_z"] = len(with_z[with_z["num_bands"] == i])
        count_dict[i]["without_z"] = len(
            without_z[without_z["num_bands"] == i])
    band_count_df = pd.DataFrame(count_dict).T
    return band_count_df


def construct_availability_bar_plots(band_count_df, ax, title, source_nums):
    """Constructs the desired bar plots on the given ax."""
    ax.grid(True, axis="y")
    labels = [mt.BAND_LABEL_DICT[band] for band in band_count_df.index]
    x = np.arange(len(labels))
    width = 0.39  # the width of the bars
    space = 0
    rects1 = ax.bar(
        x - width / 2 - space, band_count_df["extended"], width, label='Extended', color=band_count_df["Colour"], edgecolor="k", alpha=0.7)
    rects2 = ax.bar(
        x + width / 2 + space, band_count_df["pointlike"], width, label='Pointlike', color=band_count_df["Colour"], edgecolor="k")
    for bar in rects1:
        bar.set_linewidth(0.5)
        bar.set_linestyle("dashed")
    for bar in rects2:
        bar.set_linewidth(0.5)
        bar.set_hatch("///")

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('Relative availability')
    ax.set_ylim(0, 1)
    if title:
        ax.set_title(title, size="small")
    ax.set_xticks(x)
    ax.set_xticklabels(list(labels), minor=False, rotation=90)
    ax.axhline(1, color="k", linewidth=0.6)
    ax.bar_label(rects1, labels=[
                 int(num * source_nums["extended"]) for num in band_count_df["extended"]], label_type='edge', rotation=90, padding=2.5, size="x-small")
    ax.bar_label(rects2, labels=[
                 int(num * source_nums["pointlike"]) for num in band_count_df["pointlike"]], label_type='edge', rotation=90, padding=2.5, size="x-small")


def plot_input_distribution(df, stem="", excluded_surveys=(), title=True):
    """Produce a bar plot of the relative input distribution for pointlike and extended.
    Parameters:
        df: Input or output Dataframe with magnitude columns in each band.
        excluded_surveys: i. e. hsc or kids not intended for plotting"""
    fig, axes = plt.subplots(1, 1, figsize=cm.set_figsize(fraction=.8))
    total_count = len(df)
    desired_bands = [
        band for band in mt.ORDERED_BANDS if not band in excluded_surveys] + ["ZSPEC"]
    band_count_df, colours, source_nums = construct_band_availability_dataframe(
        df, desired_bands)
    if title:
        title = f"Photometry in the eFEDS field  ({total_count} sources in total)"
    construct_availability_bar_plots(band_count_df, axes, title, source_nums)
    legend_patches = [Patch(facecolor=colour, edgecolor='k',
                            label=mt.SURVEY_NAME_DICT[key]) for key, colour in colours.items()]
    ext_plike_patches = [Patch(facecolor="white", edgecolor='k',
                               label=f"{label.capitalize()} ({source_nums[label]})", linestyle=lstyle, hatch=hatch) for label, lstyle, hatch in [("extended", "dashed", ""), ("pointlike", "-", "///")]]
    legend_2 = plt.legend(handles=ext_plike_patches, prop={
        "size": "x-small"}, bbox_to_anchor=(1, 0.7), loc=2)
    fig.add_artist(legend_2)
    fig.legend(handles=legend_patches, prop={
        "size": "x-small"}, bbox_to_anchor=(0.9, 0.9), loc=2)
    cm.save_figure(fig, "input_analysis/input_distribution")
    return band_count_df


def plot_band_number_distribution(df, stem="", excluded_surveys=("hsc", "kids"), suffix="", title=False):
    """Construct a bar plot with the number of possible bands for each of the sources."""
    total_count = len(df)
    allowed_bands = [f"mag_{band}" for band in mt.BAND_LIST if mt.give_survey_for_band(
        band) not in excluded_surveys]
    band_count_df = construct_num_band_dataframe(df, allowed_bands)
    fig, axes = plt.subplots(1, 1, figsize=cm.set_figsize(fraction=.8))
    ax = axes
    ax.grid(True, axis="y")
    labels = list(band_count_df.index)
    x = np.arange(len(labels))
    width = 0.6  # the width of the bars
    rects1 = ax.bar(
        x, band_count_df["without_z"], width, color="g", bottom=band_count_df["with_z"], label=f"Without spec-$z$ ({band_count_df['without_z'].sum()})", edgecolor="k", alpha=0.7)
    rects2 = ax.bar(
        x, band_count_df["with_z"], width, label=f"With spec-$z$ ({band_count_df['with_z'].sum()})", color="k", edgecolor="k", alpha=0.7)
    ax.legend()
    if title:
        title = "Number of bands available for the sources in the eFEDS field"
        ax.set_title(title)
    ax.set_xlabel("Number of bands with photometry available")
    ax.set_ylabel("Number of sources")
    ax.set_xticks(x)
    ax.set_xlim(4, max(x) + 1)
    cm.save_figure(fig, f"input_analysis/band_number_distribution{suffix}")


if __name__ == "__main__":
    df = mt.read_plike_and_ext(prefix="matches/test2_",
                               suffix="_processed_table.fits")
    df = mt.add_mag_columns(df)
    # df = df[df["Type"] == "extended"]
    plot_band_number_distribution(df, suffix="excluding_i_band", title=False)
    # plot_r_band_magnitude(df)
    # plot_input_distribution(df, title=False)
