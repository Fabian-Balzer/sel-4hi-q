
# %%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import util.configure_matplotlib as cm
import util.my_tools as mt


def plot_performance_against_band_num(df, stem_name):
    """Plot the photo-z performance against the number of bands used for the fit (according to LePhare)."""
    fig, ax = plt.subplots(
        1, 1, figsize=cm.set_figsize(fraction=.8))
    df["ZMeasure"] = (df["ZBEST"] - df["ZSPEC"]) / (1 + df["ZSPEC"])
    df.plot.scatter("ZMeasure", "NBAND_USED", ax=ax, s=3,
                    c=df["PDZ_BEST"], cmap=plt.cm.viridis)
    ax.set_xlabel("$\Delta z/(1+z_{\rm spec})$")
    ax.set_ylabel("Number of bands used for the fit")

    cm.save_figure(fig, "output_analysis/performance_against_band_num")


if __name__ == "__main__":
    STEM = "new_test"
    output_df = mt.read_plike_and_ext(
        prefix=f"lephare_output/{STEM}_", suffix=".fits")
    output_df = mt.add_filter_columns(output_df)

    for ttype in ["pointlike", "extended"]:
        plot_performance_against_band_num(output_df, STEM)
