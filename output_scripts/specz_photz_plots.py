# %%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import util.configure_matplotlib as cm
import util.my_tools as mt


def make_scatter_plot(df, ax, label_prefix="", color="blue"):
    """Makes a scatter plot of the ZSPEC vs the ZBEST column on the given ax.
    Distinguishes the i and non-i-band."""
    bad, total = len(df[df['ISOUTLIER']]), len(df)
    perc = bad / total  # Outlier percentage
    label = label_prefix + f"{bad} of {total} ({perc*100:.1f} \%) outliers"
    df.plot.scatter("ZSPEC", "ZBEST", s=1,
                    label=label, color=color, ax=ax, marker=".")


def cut_i_band(df, ax):
    """Filter for i or i2 bands:
    # Without i or i2 band:"""
    with_i = df[df["filter_list"].apply(
        lambda x: 14 in x or 15 in x or 16 in x)]
    if len(with_i) == 0:
        label = "Sources without i band\n"
        make_scatter_plot(df, ax, label)
        return
    without_i = df[df["filter_list"].apply(
        lambda x: not(14 in x or 15 in x or 16 in x))]
    if len(without_i) == 0:
        label = "Sources with i band\n"
        make_scatter_plot(df, ax, label)
        return
    make_scatter_plot(without_i, ax, "Without i band\n")
    make_scatter_plot(with_i, ax, "Including i band\n", "black")


def plot_auxiliary_lines(ax, max_z=5, bounds=True):
    """Plots the necessary lines for outliers onto the given ax and provides the axis names."""
    z = np.arange(0, 10)
    ax.plot(z, mt.linearfunc(z), "r-", label="$z_{spec}$")
    ax.plot(z, mt.linearfunc(z, m=0.05), "r-",
            label="$z_{spec}\pm 0.05(1+z_{spec})$")
    ax.plot(z, mt.linearfunc(z, m=-0.05), "r-")
    ax.plot(z, mt.linearfunc(z, m=0.15), "r--",
            label="$z_{spec}\pm 0.15(1+z_{spec})$")
    ax.plot(z, mt.linearfunc(z, m=-0.15), "r--")
    ax.set_xlabel("Spectroscopic redshift")
    ax.set_ylabel("Best photometric redshift")
    ax.set_xlim(0, max_z)
    ax.set_ylim(0, max_z)
    if bounds:
        ax.axhspan(0, 0.5, xmin=0.5 / max_z, color="red", alpha=0.2,
                   label=r"$z_{phot} < 0.5$ (false negative)")
        ax.axvspan(0, 0.5, ymin=0.5 / max_z, color="green", alpha=0.2,
                   label=r"$z_{spec} < 0.5$ (false positive)")


def plot_photoz_vs_specz(df, ttype, stem_name="test", comparison=False, fnamesuffix=""):
    """Makes a scatter plot of photo z vs spectro-z."""
    fig, ax = plt.subplots(
        1, 1, figsize=cm.set_figsize(fraction=.8))
    ax.set_aspect(True)
    fig.set_facecolor("white")
    df = df[df["HASGOODZ"]]
    df = df[df["Type"] == ttype]
    if comparison:
        make_scatter_plot(df, ax)
    else:
        cut_i_band(df, ax)
    ax.set_title(f"{ttype.capitalize()} sources")
    max_z = 5 if ttype == "pointlike" else 2
    plot_auxiliary_lines(ax, max_z)
    # box = axes.get_position()
    # ax.set_position([box.x0, box.y0, box.width * 0.75, box.height * 0.75])
    fig.legend(loc=2, bbox_to_anchor=(0.77, 0.9), shadow=True,
               fancybox=True)
    ax.get_legend().remove()

    cm.save_figure(
        fig, f"output_analysis/{stem_name}_spec_z_phot_z_{ttype}{fnamesuffix}")


if __name__ == "__main__":
    output_df = mt.read_plike_and_ext(
        prefix="lephare_output/test2_", suffix=".fits")
    output_df = mt.add_filter_columns(output_df)

    for ttype in ["pointlike", "extended"]:
        plot_photoz_vs_specz(output_df, ttype)

    # mt.save_dataframe_as_fits(df["ra", "dec", "Type"], "Ra_dec_sources")

    # df = mt.read_plike_and_ext("data/input/plike_processed_input.fits",
    #                            "data/input/ext_processed_input.fits")
    # df = df.rename(columns={"PHOT_Z": "ZBEST"})
    # mt.add_outlier_information(df)
    # for ttype in ["pointlike", "extended"]:
    #     plot_photoz_vs_specz(df, ttype, comparison=True,
    #                          fnamesuffix="_shu_agn_")
