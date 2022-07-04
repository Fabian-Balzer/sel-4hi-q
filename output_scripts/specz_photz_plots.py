# %%
import matplotlib.pyplot as plt
import numpy as np
import util.configure_matplotlib as cm
import util.my_tools as mt
from matplotlib.offsetbox import AnchoredText


def make_scatter_plots(df, main_ax, secondary_ax, color="k"):
    """Makes a scatter plot of the ZSPEC vs the ZBEST column on the given ax.
    Distinguishes the i and non-i-band."""
    markerstyle = "."
    size = 1
    df[~df["IsOutlier"]].plot.scatter("ZSPEC", "ZBEST", s=size,  # label=label,
                                      color=color, ax=main_ax, marker=markerstyle)
    df[df["IsOutlier"]].plot.scatter(
        "ZSPEC", "ZBEST", s=size, color=color, ax=main_ax, marker=markerstyle, alpha=0.7)
    df.plot.scatter("ZSPEC", "ZMeasure", s=size,
                    color=color, ax=secondary_ax, marker=markerstyle)


def plot_auxiliary_lines(main_ax, secondary_ax, max_z=5, bounds=True, label=""):
    """Plots the necessary lines for outliers onto the given ax and provides the axis names."""
    lwidth = 0.7
    z = np.arange(0, 10)
    main_ax.plot(z, mt.linearfunc(z), "k-", label="$z_{spec}$", lw=lwidth)
    main_ax.plot(z, mt.linearfunc(z, m=0.05), "r-",
                 label="$z_{spec}\pm 0.05(1+z_{spec})$", lw=lwidth)
    main_ax.plot(z, mt.linearfunc(z, m=-0.05), "r-", lw=lwidth)
    main_ax.plot(z, mt.linearfunc(z, m=0.15), "r--",
                 label="$z_{spec}\pm 0.15(1+z_{spec})$", lw=lwidth)
    main_ax.plot(z, mt.linearfunc(z, m=-0.15), "r--", lw=lwidth)
    main_ax.set_ylabel("Best photometric redshift")
    main_ax.set_xlim(0, max_z)
    main_ax.set_ylim(0, max_z)
    main_ax.set_xticks([])
    if bounds:
        main_ax.axhspan(0, 0.5, xmin=0.5 / max_z, color="red", alpha=0.2,
                        label=r"False negative")
        main_ax.axvspan(0, 0.5, ymin=0.5 / max_z, color="green", alpha=0.2,
                        label=r"False positive")
    secondary_ax.set_ylabel("$\Delta z/(1+z)$")
    secondary_ax.set_xlabel("Spectroscopic redshift")
    secondary_ax.set_ylim(-0.3, 0.3)
    secondary_ax.set_xlim(0, max_z)
    secondary_ax.axhline(0, linestyle="-", color="k", lw=lwidth)
    secondary_ax.axhline(0.05, linestyle="-", color="r", lw=lwidth)
    secondary_ax.axhline(-0.05, linestyle="-", color="r", lw=lwidth)
    secondary_ax.axhline(0.15, linestyle="--", color="r", lw=lwidth)
    secondary_ax.axhline(-0.15, linestyle="--", color="r", lw=lwidth)
    at = AnchoredText(
        label, prop=dict(size="x-small"), frameon=True, loc='upper left')
    at.patch.set_boxstyle("round, pad=0.1, rounding_size=0.1")
    at.patch.set_alpha(0.8)
    main_ax.add_artist(at)
    main_ax.legend(loc="lower right", prop=dict(
        size="x-small"), frameon=True)
    legend = main_ax.get_legend()
    legend.legendPatch.set_boxstyle("round, pad=0.1, rounding_size=0.1")
    legend.legendPatch.set_edgecolor("k")
    legend.legendPatch.set_alpha(0.8)


def set_subplot_positions(main, secondary, is_left=False, is_right=False):
    """Sets the dimensions for the two subplots on top of each other"""
    # width and height have to be equal for the fig as well to preserve aspect ratio
    left, width, = 0.05, 0.8
    bottom, height = 0.23, 0.8
    if is_left:
        width = 0.4
        height = 0.8
    if is_right:
        left, width, = 0.55, 0.4
        height = 0.8
    spacing, secondheight = 0.000, 0.18
    rect_main = [left, bottom, width, height]
    rect_secondary = [left, bottom -
                      secondheight - spacing, width, secondheight]
    main.set_position(rect_main)
    secondary.set_position(rect_secondary)


def plot_ttype(df, ttype, main, secondary):
    """Wrapper function to filter the df for the ttype and plot it on the given axes"""
    df = df[df["HasGoodz"]]
    df = df[df["Type"] == ttype]
    make_scatter_plots(df, main, secondary)
    main.set_title(mt.give_plot_title(ttype, True))
    max_z = 5 if ttype == "pointlike" else 2
    plot_auxiliary_lines(main, secondary, max_z,
                         label=mt.give_photoz_performance_label(df))


def plot_photoz_vs_specz(df, ttype="both"):
    """Makes a scatter plot of photo z vs spectro-z."""
    df = df[df["HasGoodz"]]
    if ttype == "both":
        fig, axes = plt.subplots(2, 2, figsize=cm.set_figsize(
            width=cm.LATEXWIDTH, height=0.5 * cm.LATEXWIDTH, fraction=1))
        main_plike, secondary_plike = axes[0]
        main_ext, secondary_ext = axes[1]
        set_subplot_positions(main_plike, secondary_plike, is_left=True)
        set_subplot_positions(main_ext, secondary_ext, is_right=True)
        plot_ttype(df, "pointlike", main_plike, secondary_plike)
        plot_ttype(df, "extended", main_ext, secondary_ext)
    elif ttype is not None:
        fig, (main, secondary) = plt.subplots(
            2, 1, figsize=cm.set_figsize(width=0.5 * cm.LATEXWIDTH, height=0.5 * cm.LATEXWIDTH, fraction=1))
        set_subplot_positions(main, secondary)
        plot_ttype(df, ttype, main, secondary)
    else:
        mt.LOGGER.error(
            "Please specify a valid ttype [pointlike, extended or both] to plot the specz-photoz-plot.")

    # fig.legend(loc=2, bbox_to_anchor=(0.77, 0.9), shadow=False,
    #            fancybox=False)
    # main.get_legend().remove()
    stem = mt.CUR_CONFIG["LEPHARE"]["output_stem"]
    cm.save_figure(fig, f"{ttype}_spec_z_phot_z", "output_analysis", stem)


if __name__ == "__main__":
    STEM = "without_i"
    output_df = mt.read_plike_and_ext(
        prefix=f"lephare_output/{STEM}_", suffix=".fits")
    output_df = mt.add_filter_columns(output_df)
    context = mt.read_glb_context(f"lephare_output/{STEM}_pointlike.fits")
    # output_df = output_df[output_df["NBAND_USED"] > 10]
    for ttype in ["pointlike", "extended"]:
        plot_photoz_vs_specz(output_df, ttype, STEM)
