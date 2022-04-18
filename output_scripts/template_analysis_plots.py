# %%
from math import ceil, floor

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import util.configure_matplotlib as cm
import util.my_tools as mt


def plot_single_against_redshift(df, ax, template_dict, c1, c2, ttype):
    # Select a subset of the source dataframe with valid spec-z and actual data points in both of the requested bands
    for param in ["ZSPEC", f"mag_{c1}", f"mag_{c2}"]:
        df = df[(df[param] > 0) & (df[param] < 99)]
    colname = f"{c1.replace('_', '-')} - {c2.replace('_', '-')}"

    # Plot all templates:
    info_dict = mt.provide_template_info(f"{ttype}_baseline_templates.list")
    for temp_num in template_dict:
        t_df = template_dict[temp_num]
        y_data = t_df[f"mag_{c1}"] - t_df[f"mag_{c2}"]
        label = ""  # str(info_dict.get(temp_num)).replace("_", "-")
        ax.plot(t_df["ZSPEC"], y_data, "-", lw="0.5", label=label)
    df[colname] = df[f"mag_{c1}"] - df[f"mag_{c2}"]
    zmax = 2.5 if ttype == "extended" else 4.5
    ymin, ymax = -1, 1
    # for i, row in df.iterrows():
    #     if row["ZSPEC"] - row["ZBEST"] > 0.2:
    #         x = [row["ZSPEC"], row["ZBEST"]]
    #         y = [row[colname], row[colname]]
    #         ax.plot(x, y, "k-", linewidth=0.2)
    #         ax.plot(x[0], y[0], "rx", markersize=0.2)
    # for col, style in [("ZBEST", "rx"), ("ZSPEC", "kx")]:
    #     ax.plot(df[col], df[colname], style,
    #             markersize=0.5, label=col, alpha=0.5)

    bad, good = df[df["ISOUTLIER"]], df[~df["ISOUTLIER"]]
    for subset, label, style in [(bad, "Outliers", "rx"), (good, "Good sources", "kx")]:
        z, y = subset["ZSPEC"], subset[colname]
        if len(y) > 0:
            ymin = min(ymin, min(y))
            ymax = max(ymax, max(y))
            ax.plot(z, y, style, markersize=1,
                    label=f"{label} ({len(y)})", alpha=0.5)
    ax.set_xlim(0, zmax)
    ax.set_ylim(floor(ymin), ceil(ymax))
    ax.set_xlabel("z", labelpad=0.4)
    ax.set_ylabel(colname, labelpad=0.4)


def plot_multiple_against_redshift(df, template_df, ttype, bands=("g", "z"), templates_to_plot=None, onebigplot=False, joint_fig=False):
    """Creates a plot of color difference vs. redshift to analyse wether the
    templates cover the parameter space.
    Parameters:
        source_df: The pandas DataFrame containing output information of the sources.
        template_df: The pandas DataFrame hosting the template information
        bands: The bands of interest for the plots."""
    df = df[df["Type"] == ttype]
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
            plot_single_against_redshift(df, ax, temp_dict, c1, c2, ttype)
            ax.legend(prop={"size": "x-small"})
            if not onebigplot:
                name = f"{ttype}_template_plot_{c1}-{c2}"
                cm.save_figure(fig, f"output_analysis/templates/{name}")
            n += 1
    if onebigplot:
        name = f"{ttype}_joint_template_plot"
        cm.save_current_figures(f"output_analysis/templates/{name}")
    elif joint_fig:
        filename = f"output_analysis/templates/{ttype}_all_temp_plots.pdf"
        cm.save_current_figures(filename)


def plot_problematic_templates(df, ttype):
    """Plots a histogram with the available templates and the percentage of usage for LePhare"""
    fig, axes = plt.subplots(1, 1, figsize=cm.set_figsize(fraction=.8))
    subset = df[df["Type"] == ttype]
    good = subset[~subset["ISOUTLIER"]]
    bad = subset[subset["ISOUTLIER"]]
    df1 = good["MOD_BEST"].value_counts().rename("Good")
    df2 = bad["MOD_BEST"].value_counts().rename("Bad")
    both = pd.concat([df1, df2], axis=1).drop(labels=[-99])
    both["Total"] = both["Good"].fillna(0) + both["Bad"].fillna(0)
    both = both.sort_values(by="Total", ascending=False)

    threshold = max(1, max(both["Total"]) * 0.07)
    print(
        f"Adopting threshold (number of most occurences for models to be combined in 'other') of {threshold:.0f} for {ttype}.")
    other_dict = {"Good": {"Other": 0}, "Bad": {
        "Other": 0}, "Total": {"Other": 0}}
    for check in ["Good", "Bad"]:
        is_single = (both[check] <= threshold) & (
            both["Total"] <= threshold * 2)
        singles = both[is_single]
        ind = [model for model in singles.index]
        # print(ind)
        other_dict[check]["Other"] = len(singles)
    is_single = (both["Bad"].fillna(0) <= threshold) & (
        both["Good"].fillna(0) <= threshold)
    both = both[~is_single]
    both = pd.concat([both, pd.DataFrame.from_dict(other_dict)])
    axes.grid(True, axis="y")
    labels = both.index

    x = np.arange(len(labels))
    width = 0.4  # the width of the bars

    rects1 = axes.bar(x, both["Good"], width, bottom=both["Bad"].fillna(
        0), label='Good fits', align="center")
    rects2 = axes.bar(x, both["Bad"], width, label='Bad fits', color="red")

    # Add some text for labels, title and custom x-axis tick labels, etc.
    axes.set_ylabel('Number of times used')
    title = f"Models used for the best fits ({ttype})"
    axes.set_title(title, size="small")
    axes.set_xticks(x)
    axes.set_xticklabels(list(labels), minor=False, rotation=90, )
    axes.legend(loc="upper right")
    fig.set_facecolor("white")
    cm.save_figure(fig, f"output_analysis/{ttype}_used_templates")
    return both


if __name__ == "__main__":
    output_df = mt.read_plike_and_ext(
        prefix="lephare_output/test2_", suffix=".fits")
    output_df = mt.add_filter_columns(output_df)

    for ttype, templates_to_plot in [("extended", [54, 13, 30, 15, 23, 16, 29, 85, 11, 28, 69, 68, "Other"]), ("pointlike", [27, 25, 28, 29, 18, 22, 23])]:
        # df2 = plot_problematic_templates(output_df, ttype)
        # templates_to_plot = list(df2.index)
        template_df = mt.read_template_library(f"{ttype}_mag_lib.dat")
        plot_multiple_against_redshift(
            output_df, template_df, ttype, bands=mt.ORDERED_BANDS, templates_to_plot=templates_to_plot, onebigplot=True, joint_fig=True)
    #     plot_multiple_against_redshift(
    #         output_df, template_df, ttype, bands=("W1", "W2"), onebigplot=False)
    #     plot_multiple_against_redshift(
    #         output_df, template_df, ttype, bands=("g", "z"), onebigplot=False)
    #     plot_multiple_against_redshift(
    #         output_df, template_df, ttype, bands=("i_kids", "z"), onebigplot=False)
    # filename = f"output_analysis/templates/selection_plots.pdf"
    # cm.save_current_figures(filename)
