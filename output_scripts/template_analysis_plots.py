# %%
from math import ceil, floor

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import util.configure_matplotlib as cm
import util.my_tools as mt
from matplotlib.transforms import Bbox


def plot_color_versus_redshift(df: pd.DataFrame, ax, ttype: str,
                               temp_df: pd.DataFrame, c1: str, c2: str,
                               fitted_only: bool, plot_sources: bool):
    """Produce a color-redshift plot of colors c1-c2 for the sources and templates. """
    # Select a subset of the source dataframe with valid spec-z and actual data points in both of the requested bands
    df = df[df["Type"] == ttype]
    temp_df = temp_df[temp_df["Type"] == ttype]
    df = df[df["HasGoodz"]]
    for param in [f"mag_{c1}", f"mag_{c2}"]:
        df = df[(df[param] > 0) & (df[param] < 99)]
    colname = f"{mt.generate_pretty_band_name(c1)} - {mt.generate_pretty_band_name(c2)}"
    temp_nums_to_plot = give_count_df(
        df, ttype, 0.1).index if fitted_only else None
    temp_dict = mt.construct_template_dict(temp_df, temp_nums_to_plot)
    # Plot all templates:
    # info_dict = mt.provide_template_info(ttype)
    for temp_label in sorted(temp_dict):
        single_t_df = temp_dict[temp_label]
        y_data = single_t_df[f"mag_{c1}"] - single_t_df[f"mag_{c2}"]
        if not "E(B-V)" in temp_label:
            ax.plot(single_t_df["ZSPEC"], y_data,
                    "-", lw="0.5", label=temp_label)
        else:
            # Using this, the E(B-V)-lines will have the same color as their parents
            color = ax.get_lines()[-1].get_color()
            ax.plot(single_t_df["ZSPEC"], y_data, "--", lw="0.5", color=color)
        # print(temp_num)
    df[colname] = df[f"mag_{c1}"] - df[f"mag_{c2}"]
    zmax = 2.5 if ttype == "extended" else 5
    ax.set_xlim(0, zmax)
    ax.set_ylim(-1, 3)
    ax.set_xlabel("(spectroscopic) redshift $z$", labelpad=0.4)
    ax.set_ylabel(colname, labelpad=0.4)
    ax.grid(True)
    ax.set_title(mt.give_plot_title(ttype))
    if fitted_only:
        temp_num = len([label for label in temp_dict if not "E(B-V)" in label])
        ax.legend(ncol=ceil(temp_num / 3))
    if not plot_sources:  # Guard against plotting sources if unwanted
        return
    ymin, ymax = -1, 1
    bad, good = df[df["IsOutlier"]], df[~df["IsOutlier"]]
    for subset, label, style in [(bad, "Outliers", "rx"), (good, "Good sources", "kx")]:
        z, y = subset["ZSPEC"], subset[colname]
        if len(y) > 0:
            ymin = min(ymin, min(y))
            ymax = max(ymax, max(y))
            ax.plot(z, y, style, markersize=1,
                    label=f"{label} ({len(y)})", alpha=0.5)
    ax.set_ylim(floor(ymin * 2) / 2, ceil(ymax * 2) / 2)
    # for i, row in df.iterrows():
    #     if row["ZSPEC"] - row["ZBEST"] > 0.2:
    #         x = [row["ZSPEC"], row["ZBEST"]]
    #         y = [row[colname], row[colname]]
    #         ax.plot(x, y, "k-", linewidth=0.2)
    #         ax.plot(x[0], y[0], "rx", markersize=0.2)
    # for col, style in [("ZBEST", "rx"), ("ZSPEC", "kx")]:
    #     ax.plot(df[col], df[colname], style,
    #             markersize=0.5, label=col, alpha=0.5)


def plot_multiple_against_redshift(source_df, template_df, ttype, stem, bands=("g", "z"), templates_to_plot=None, onebigplot=False, joint_fig=False):
    """Creates a plot of color difference vs. redshift to analyse wether the
    templates cover the parameter space.
    Parameters:
        source_df: The pandas DataFrame containing output information of the sources.
        template_df: The pandas DataFrame hosting the template information
        ttype: pointlike, extended or
        bands: The bands of interest for the plots."""
    if ttype != "both":
        source_df = source_df[source_df["Type"] == ttype]
    template_df = template_df.sort_values(by=["ZSPEC", "model"])
    temp_dict = mt.construct_template_dict(template_df, templates_to_plot)
    if onebigplot:
        # Gaussian formula for the number of plots
        n = len(bands) - 1
        num_plots = n * (n + 1) / 2
        num_cols = 3
        num_rows = 2
        coords = [(i, j) for i in range(2) for j in range(3)]
        num_figs = ceil(num_plots / 6)
        fig_list = []
        for i in range(num_figs):
            fig, axes = plt.subplots(
                num_rows, num_cols, figsize=cm.set_figsize(subplots=(num_rows, num_cols)), sharex=True, sharey=False)
            fig_list.append((fig, axes))
        plt.subplots_adjust(left=None, bottom=None, right=None,
                            top=None, wspace=0.5, hspace=0.5)
    n = 0
    # Iterate over all possible combinations of the bands requested
    for k, c1 in enumerate(bands):
        for c2 in bands[k + 1:]:
            if not onebigplot:
                fig, ax = plt.subplots(
                    1, 1, figsize=cm.set_figsize(fraction=.5))
            else:
                fig, axes = fig_list[floor(n / 6)]
                i, j = coords[n % 6]
                ax = axes[i][j]
            plot_single_against_redshift(
                source_df, ttype, stem, ax, temp_dict, c1, c2)
            if templates_to_plot is not None:
                ax.legend(prop={"size": "x-small"})
            if not onebigplot:
                name = f"{stem}_{ttype}_template_plot_{c1}-{c2}"
                cm.save_figure(fig, name, "output_analysis/templates")
            n += 1
    if onebigplot:
        name = f"{ttype}_joint_template_plot"
        cm.save_current_figures(f"output_analysis/templates/{name}")
    elif joint_fig:
        filename = f"output_analysis/templates/{ttype}_all_temp_plots.pdf"
        cm.save_current_figures(filename)


def give_score_for_template(df, temp_num):
    """Calculates the score of a template with given template number."""
    subset = df[df["MOD_BEST"] == temp_num]
    if len(subset) == 0:
        return 0
    outlier_scores = np.exp(-abs(subset["ZMeasure"]) - subset["CHI_BEST"] / 2)
    return sum(outlier_scores) / len(subset)


def give_templates_to_keep(df, ttype, removal_frac=0.25):
    """Returns a list of the templates to keep via score selection + the a list of the discarded templates"""
    count_df = give_count_df(df, ttype)
    count_df.sort_values("Score", inplace=True, ascending=False)
    quartil = int(len(count_df) * (1 - removal_frac))
    temps_to_keep = pd.concat(
        [count_df.iloc[:quartil], count_df[count_df["Bad_frac"] == 0]]).drop_duplicates()
    temps_to_keep = temps_to_keep[temps_to_keep["Bad_frac"] < 1]
    return [temp_num for temp_num in temps_to_keep.index]


def give_templates_to_drop(df, ttype, removal_frac=0.25):
    """Returns the a list of the discarded templates"""
    count_df = give_count_df(df, ttype)
    temps_to_keep = give_templates_to_keep(df, ttype, removal_frac)
    # Find the dropped templates:
    ds1, ds2 = set(count_df.index), set(temps_to_keep)
    discarded_temps = ds1.difference(ds2)
    mt.LOGGER.info("Dropping %d of %d templates with numbers %s.",
                   len(discarded_temps), len(count_df), discarded_temps)
    return discarded_temps


def give_count_df(df, ttype: str, threshold_factor=0):
    """Analyzes an output dataframe for the models used for the best fits
    and returns a dataframe listing the counts of outliers/good fits and scores for
    model numbers (in index) that were matched more times than the threshold."""
    df = df[df["Type"] == ttype]
    df = df[df["HasGoodz"]]
    good = df[~df["IsOutlier"]]
    bad = df[df["IsOutlier"]]
    df1 = good["MOD_BEST"].value_counts().rename("Good")
    df2 = bad["MOD_BEST"].value_counts().rename("Bad")
    both = pd.concat([df1, df2], axis=1)
    if -99 in both.index:
        both = both.drop(labels=[-99])
    both["Total"] = both["Good"].fillna(0) + both["Bad"].fillna(0)
    both = both.sort_values(by="Total", ascending=False)
    both["Score"] = both.index.map(
        lambda temp_num: give_score_for_template(df, temp_num))
    both["Bad_frac"] = both["Bad"].fillna(0) / both["Total"]
    if threshold_factor == 0:
        return both
    # Otherwise, adopt a threshold and a column containing 'other' templates.
    threshold = max(1, max(both["Total"]) * threshold_factor)
    mt.LOGGER.info(
        "Adopting threshold (number of most occurences for models to be combined in 'other') of %.0f for %s.", threshold, ttype)
    other_dict = {"Good": {"Other": 0}, "Bad": {
        "Other": 0}, "Total": {"Other": 0}}
    for check in ["Good", "Bad"]:
        is_below_threshold = (both[check] <= threshold) & (
            both["Total"] <= threshold * 2)
        too_low = both[is_below_threshold]
        other_dict[check]["Other"] = len(too_low)
    is_below_threshold = (both["Bad"].fillna(0) <= threshold) & (
        both["Good"].fillna(0) <= threshold)
    both = both[~is_below_threshold]
    both = pd.concat([both, pd.DataFrame.from_dict(other_dict)])
    return both


def plot_bar_template_outliers(df, ax, ttype):
    """Produce a bar plot on the given ax for the templates of ttype."""
    df = df[df["Type"] == ttype]
    df = df[df["HasGoodz"]]
    both = give_count_df(df, ttype, threshold_factor=0.1)
    labels = both.index

    x = np.arange(len(labels))
    width = 0.5  # the width of the bars

    rects1 = ax.bar(x, both["Good"], width, bottom=both["Bad"].fillna(
        0), label='Good fits', align="center", color="green", edgecolor="k")
    rects2 = ax.bar(x, both["Bad"], width, label='Outliers',
                    color="firebrick", edgecolor="k")

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('Number of times used')
    # f"Models used for the best fits ({ttype})"
    ax.set_title(mt.give_plot_title(ttype, True))
    ax.set_xticks(x)
    ax.grid(True, axis="y")
    ax.set_xticklabels(list(labels), minor=False, rotation=90, )
    ax.legend()
    return both


def plot_bar_template_scores(df, ax, ttype):
    """Plot the scores of the templates in a bar plot."""
    count_df = give_count_df(df, ttype)
    count_df.sort_values("Score", inplace=True, ascending=False)
    labels = count_df.index

    x = np.arange(len(labels))
    width = 0.5  # the width of the bars
    sm = mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(
        vmin=0, vmax=1), cmap=mpl.cm.inferno_r)
    rects1 = ax.bar(x, count_df["Score"], width,
                    align="center", edgecolor="k", color=sm.to_rgba(count_df["Bad_frac"]))
    for i in x:
        num = int(count_df["Total"].iloc[i])
        offset = len(str(num)) * 0.1
        ax.text(i - offset, count_df["Score"].iloc[i] +
                count_df["Score"].max() * 0.012, s=num)
    quartil = int(len(labels) * (1 - 0.25))
    ax.axvline(x[quartil] - width, linestyle="--", color="red")
    cbar = plt.colorbar(sm, pad=0.02)
    cbar.set_label("Fraction of outliers", rotation=270, labelpad=10)

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('Score')
    # f"Models used for the best fits ({ttype})"
    title_text = mt.give_plot_title(
        ttype, True) + " (" + mt.give_temp_listname(ttype).split("/")[-1].replace("_", r"\_") + ")"
    ax.set_title(title_text)
    ax.set_xticks(x)
    ax.set_xlim(min(x) - width, max(x) + width)
    ax.grid(True, axis="y")
    ax.set_xticklabels(list(labels), minor=False, rotation=90, )
    dropping = [str(temp_num)
                for temp_num in sorted(give_templates_to_drop(df, ttype))]
    drop_text = "Dropping:\n"
    for i, temp_num in enumerate(dropping):
        if i % 5 == 0 and i != 0:
            drop_text += "\n"
        drop_text += temp_num + ", "
    ax.text(0.77, 0.8, drop_text.strip(", "), transform=ax.transAxes,
            bbox=dict(facecolor='white', alpha=1), fontsize="small")
    return count_df


if __name__ == "__main__":
    output_df = mt.read_saved_df(cat_type="out")

    # b = plot_bar_template_scores(output_df, ax, "extended")
    # print(give_templates_to_keep(output_df, "pointlike"))
    # for ttype, templates_to_plot in [("extended", [54, 13, 30, 15, 23, 16, 29, 85, 11, 28, 69, 68, "Other"]), ("pointlike", [27, 25, 28, 29, 18, 22, 23])]:
    #     # df2 = plot_problematic_templates(output_df, ttype)
    #     # templates_to_plot = list(df2.index)
    #     template_df = mt.read_template_library(f"{ttype}_mag_lib.dat")
    #     plot_multiple_against_redshift(
    #         output_df, template_df, ttype, bands=mt.ORDERED_BANDS, templates_to_plot=templates_to_plot, onebigplot=True, joint_fig=True)
    #     plot_multiple_against_redshift(
    #         output_df, template_df, ttype, bands=("W1", "W2"), onebigplot=False)
    #     plot_multiple_against_redshift(
    #         output_df, template_df, ttype, bands=("g", "z"), onebigplot=False)
    #     plot_multiple_against_redshift(
    #         output_df, template_df, ttype, bands=("i_kids", "z"), onebigplot=False)
    # filename = f"output_analysis/templates/selection_plots.pdf"
    # cm.save_current_figures(filename)
