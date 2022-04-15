
# %%

import matplotlib.pyplot as plt
import numpy as np
import util.configure_matplotlib as cm
import util.my_tools as mt


def plot_color_difference(df, c1, c2, c3=None, c4=None, stem="", **kwargs):
    """Plots the color difference in two bands against the first color given.
    Parameters:
        df: DataFrame containing Type and magnitude columns.
        c1, c2: Colour on the x-axis (mag(c1) - mag(c2).
        c3, c4: If both are left out, mag(c2) is used for the y-axis.
            If only c3 is provided, mag(c2) - mag(c3) is used.
            If only c4 is provided, mag(c4) is used.
            If both are provided, mag(c3) - mag(c4) is used.
        stem: The stem name of the catalogue, considered when saving the figure.
    """
    fig, ax = plt.subplots(1, 1, figsize=cm.set_figsize(fraction=0.75))
    clabels = {colour: colour.replace(
        "_", "-") for colour in (c1, c2, c3, c4) if colour is not None}
    for ttype in ["pointlike", "extended"]:
        subset = df[df["Type"] == ttype]
        x_data = subset[f"mag_{c1}"] - subset[f"mag_{c2}"]
        if c4 is None:
            y_data = subset[f"mag_{c2}"] if c3 is None else subset[f"mag_{c2}"] - \
                subset[f"mag_{c3}"]
            ylabel = clabels[c2] if c3 is None else f"{clabels[c2]} - {clabels[c3]}"
        else:
            y_data = subset[f"mag_{c4}"] if c3 is None else subset[f"mag_{c3}"] - \
                subset[f"mag_{c4}"]
            ylabel = clabels[c4] if c3 is None else f"{clabels[c3]} - {clabels[c4]}"
        ax.scatter(x_data, y_data, s=0.9, label=ttype)
    xlabel = f'{clabels[c1]} - {clabels[c2]}'
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if "xlim" in kwargs:
        ax.set_xlim(kwargs["xlim"])
    if "ylim" in kwargs:
        ax.set_ylim(kwargs["ylim"])
    ax.legend()
    x = np.arange(-10, 10)
    y = 0.8 * x - 1.2
    ax.plot(x, y, "-", color="gray")
    ax.set_title("Colour-colour plot")
    cm.save_figure(
        fig, f"input_analysis/color_diff_{c1}-{ylabel.replace(' ', '')}_{stem}")


if __name__ == "__main__":
    df = mt.read_plike_and_ext(prefix="matches/test2_",
                               suffix="_processed_table.fits")
    df = mt.add_mag_columns(df)
    # plot_magnitude_dist(df)
    stem = "test2"
    plot_color_difference(df, "g", "r", "z", "W1", stem,
                          xlim=(-0.5, 2), ylim=(-3, 3))
    # plot_color_difference(df, "i_kids", "i_hsc", stem=stem)  # Why is this off?
