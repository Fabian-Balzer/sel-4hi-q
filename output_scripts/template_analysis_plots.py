# %%
import util.configure_matplotlib as cm
import matplotlib.pyplot as plt
import util.my_tools as mt
import numpy as np
import pandas as pd


def produce_template_data(temp_df, band_1: str, band_2: str):
    """Read out a given lib template file and produce an x-array containing
    the redshift and an dict with x- and y-arrays and the template number as the keys."""
    temp_dict = {}
    for temp_number in set(temp_df["Model"]):
        subset = temp_df[temp_df["Model"] == temp_number]
        y_array = subset[band_1] - subset[band_2]
        x_array = subset["Redshift"]
        temp_dict[temp_number] = [x_array, y_array]
    return temp_dict


def plot_against_redshift(source_df, template_df, bands):
    """Creates a plot of color difference vs. redshift to analyse wether the
    templates cover the parameter space.
    Parameters:
        source_df: The pandas DataFrame containing output information of the sources.
        template_df: The pandas DataFrame hosting the template information
        bands: The bands of interest for the plots."""
    # Iterate over all possible combinations of the bands requested
    for i, lower_band in enumerate(bands):
        for higher_band in bands[i + 1:]:
            # Select a subset of the source dataframe with valid spec-z and actual data points in both of the requested bands
            for param in ["ZSPEC", lower_band, higher_band]:
                df = source_df[source_df[param] > 0]
                df = df[df[param] < 99]
    # TODO... !
    fig, axes = plt.subplots(figsize=(8, 8))
    temp_dict = produce_template_data(temp_df, band_1, band_2)
    for tempxy in temp_dict.values():
        axes.plot(tempxy[0], tempxy[1], "-")
    outliers = df[df["ISOUTLIER"]]
    good = df[~df["ISOUTLIER"]]
    xg = good["ZSPEC"]
    yg = good[band_1] - good[band_2]
    xb = outliers["ZSPEC"]
    yb = outliers[band_1] - outliers[band_2]
    axes.plot(xg, yg, "kx", markersize=2, label="Good sources")
    axes.plot(xb, yb, "rx", markersize=3, label="Outliers")
    zmin = 0
    zmax = max(xg)
    axes.set_xlim(zmin, zmax)
    axes.set_xlabel("Redshift z")
    axes.legend()
    band_1, band_2 = (band.split("_")[1] for band in (band_1, band_2))
    label = f"{band_1} - {band_2}"
    axes.set_ylabel(label)
    name = f"{sourcetype}_template_plot_{band_1}-{band_2}_{stem_name}"
    axes.set_title(name)
    name = name + ".png"
    save_to_stempath(fig, stem_name, name, sourcetype)


def plot_problematic_templates(df, ttype):
    """Plots a histogram with the available templates and the percentage of usage for LePhare"""
    fig, axes = plt.subplots(1, 1, figsize=cm.set_figsize(fraction=.8))
    subset = df[df["Type"] == ttype]
    good = subset[subset["ISOUTLIER"]]
    bad = subset[~subset["ISOUTLIER"]]
    df1 = good["MOD_BEST"].value_counts().rename("Good")
    df2 = bad["MOD_BEST"].value_counts().rename("Bad")
    both = pd.concat([df1, df2], axis=1)
    both["Total"] = both["Good"] + both["Bad"]
    both = both.sort_values(by="Total", ascending=False)

    axes.grid(True, axis="y")
    labels = both.index
    x = np.arange(len(labels))
    width = 0.4  # the width of the bars
    space = 0

    rects1 = axes.bar(
        x - width / 2 - space, both["Good"], width, label='Good fits')
    rects2 = axes.bar(
        x - width / 2 - space, both["Bad"], width, label='Bad fits')

    # Add some text for labels, title and custom x-axis tick labels, etc.
    axes.set_ylabel('Number of times used')
    title = f"Models used for the best fits ({ttype})"
    axes.set_title(title, size="small")
    axes.set_xticks(x)
    axes.set_xticklabels(list(labels), minor=False, rotation=90)
    axes.legend(loc="upper right")
    fig.set_facecolor("white")
    cm.save_figure(fig, f"data/plots/{ttype}_used_templates")
    return both


if __name__ == "__main__":
    df = mt.read_plike_and_ext("data/output/PLIKE_output_220126_test.fits",
                               "data/output/EXT_output_220126_test.fits")
    df = mt.read_plike_and_ext("data/output/test2_pointlike.fits",
                               "data/output/test2_extended.fits")
    df = mt.add_filter_columns(df)

    for ttype in ["pointlike", "extended"]:
        df2 = plot_problematic_templates(df, ttype)
