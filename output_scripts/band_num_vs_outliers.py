
# %%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import util.configure_matplotlib as cm
import util.my_tools as mt


def plot_performance_against_band_num(df, ttype, stem_name):
    """Plot the photo-z performance against the number of bands used for the fit (according to LePhare)."""
    fig, ax = plt.subplots(
        1, 1, figsize=cm.set_figsize(fraction=.8))
    df["ZMeasure"] = (df["ZBEST"] - df["ZSPEC"]) / (1 + df["ZSPEC"])
    df = df[df["Type"] == ttype]
    df = df[df["ZSPEC"] > 0]
    sizes = df["ZSPEC"].apply(lambda x: x if x < 2 else 2)
    colors = df["PDZ_BEST"]
    df.plot.scatter("NBAND_USED", "ZMeasure", ax=ax, s=40 / sizes, marker="_",
                    c=colors, cmap=plt.cm.viridis)
    ax.set_xlim(6.5, 14.5)
    ax.set_xlabel("Number of bands used for the fit")
    ax.set_ylabel(r"$\Delta z/(1+z_{\rm spec})$")
    ax.axhline(0, color="k", lw=0.7)
    ax.axhline(0.05, color="r", lw=0.7)
    ax.axhline(-0.05, color="r", lw=0.7)
    ax.axhline(0.15, color="r", ls="--", lw=0.7)
    ax.axhline(-0.15, color="r", ls="--", lw=0.7)
    # colorax = fig.get_axes()[1]
    # colorax.set_ylim(0, 2)

    cm.save_figure(
        fig, f"output_analysis/{stem_name}_performance_against_band_num")


def plot_performance_for_each_survey(df, ttype, stem_name):
    """Plot the photo-z performance against the number of bands used for the fit (according to LePhare)."""
    fig, ax = plt.subplots(
        1, 1, figsize=cm.set_figsize(fraction=.8))
    df["ZMeasure"] = (df["ZBEST"] - df["ZSPEC"]) / (1 + df["ZSPEC"])
    df = df[df["Type"] == ttype]
    df = df[df["ZSPEC"] > 0]
    survey_dataframes = {}
    for i, survey in enumerate(mt.BAND_DICT.keys()):
        # any([x[f"mag_{band}"] > 0 for band in mt.BAND_DICT[survey]]))]
        subframe = df
        for band in mt.BAND_DICT[survey]:
            subframe = subframe[subframe[f"mag_{band}"] > 0]
        ax.scatter(np.ones(len(subframe)) * i,
                   subframe["ZMeasure"], s=1, label=f"{survey} ({len(subframe)})")
        survey_dataframes[survey] = subframe
    ax.legend()
    # sizes = df["ZSPEC"].apply(lambda x: x if x < 2 else 2)
    # colors = df["PDZ_BEST"]
    # df.plot.scatter("NBAND_USED", "ZMeasure", ax=ax, s=40 / sizes, marker="_",
    #                 c=colors, cmap=plt.cm.viridis)
    # ax.set_xlim(6.5, 14.5)
    ax.set_xlabel("Photometry used for the bands")
    ax.set_ylabel(r"$\Delta z/(1+z_{\rm spec})$")
    ax.axhline(0, color="k", lw=0.7)
    ax.axhline(0.05, color="r", lw=0.7)
    ax.axhline(-0.05, color="r", lw=0.7)
    ax.axhline(0.15, color="r", ls="--", lw=0.7)
    ax.axhline(-0.15, color="r", ls="--", lw=0.7)
    # colorax = fig.get_axes()[1]
    # colorax.set_ylim(0, 2)

    # cm.save_figure(
    #     fig, f"output_analysis/{stem_name}_performance_against_band_num")


if __name__ == "__main__":
    STEM = "new_test"
    output_df = mt.read_plike_and_ext(
        prefix=f"lephare_output/{STEM}_", suffix=".fits")
    output_df = mt.add_filter_columns(output_df)

    for ttype in ["pointlike", "extended"]:
        # plot_performance_against_band_num(output_df, ttype, STEM)
        plot_performance_for_each_survey(output_df, ttype, STEM)
