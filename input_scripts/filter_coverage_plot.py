# -*- coding: utf-8 -*-
"""
Created on Wed Feb 23 12:20:24 2022

@author: fabian_balzer
"""

# %%
import os

import matplotlib.pyplot as plt
import pandas as pd
import util.configure_matplotlib as cm
import util.my_tools as mt
from matplotlib.patches import Patch

import input_scripts.collect_filter_data as c_f


def produce_filter_plot(df):
    """Create a plot showing the normalised transmission curves of the filters."""
    fig, ax = plt.subplots(figsize=cm.set_figsize())
    surveys = []
    for band in df.keys().get_level_values(0).drop_duplicates().sort_values():
        x_data = df[band]["lambda"] / 10
        y_data = df[band]["trans"]
        lambda_mean, band, survey = band.split(">")
        surveys.append(survey)
        label = mt.generate_pretty_band_name(band)
        x, y = x_data[y_data.idxmax()], y_data.max()
        if band in ["i2_hsc", "g", "r"]:
            y -= 0.05
        ax.text(float(lambda_mean) * 0.9, y + 0.01, label)
        # y_data = y_data / max(y_data)
        ax.fill_between(x_data, y_data, alpha=0.5,
                        color=mt.give_survey_color_for_band(band))
        # ax.plot(x_data, y_data, label=label,
        #         color="k", linewidth=0.5)

    # ax.set_xscale("log")
    ax.set_xlim(1e2, 3e4)
    ax.set_xscale("log")
    ax.set_ylim(0, 1)
    ax.grid(True)
    ax.axhline(y=1, color="k", linewidth=0.5)
    ax.set_xlabel(r"$\lambda$ [nm]")
    ax.set_ylabel("Normalised filter transmission")
    survey_dict = {survey: mt.give_survey_color(
        survey) for survey in set(surveys)}
    legend_patches = [Patch(facecolor=color, edgecolor='k',
                            label=mt.give_survey_name(survey)) for survey, color in survey_dict.items()]
    fig.legend(handles=legend_patches, prop={
        "size": "x-small"}, bbox_to_anchor=(0.78, 0.94), loc=2)
    # ax.legend(ncol=6, prop={'size': "small"},
    #           bbox_to_anchor=(0, 1.2), loc="upper left")
    cm.save_figure(fig, "filter_coverage_plot", "input_analysis")


def latex_with_lines(df, *args, **kwargs):
    kwargs['column_format'] = ''.join([''] + ['l'] * df.index.nlevels
                                      + ['l'] * df.shape[1] + [''])
    res = df.to_latex(*args, **kwargs)
    return res


def read_filter_overview_file(stem, directory="lephare_files"):
    """Read the overview file and store the information in a dataframe. Format the column names to a LaTeX-readable format."""
    fname = mt.DATAPATH + directory + "/" + stem + "_overview.filt"
    with open(fname, "r") as f:
        for line in f.readlines():
            if line.startswith("# "):
                colnames = line.strip("# ").split()
    new_colnames = {"IDENT": "Context", "NAME": "Band", "Lbda_mean": r"$\Mean{\lambda}$ [$\mu$m]",
                    "Lbeff(Vega)": r"$\lambda_\tx{eff}^\tx{Vega}$ [$\mu$m]", "FWHM": r"FWHM [$\mu$m]", "AB-cor": r"$m_\tx{AB, c}$", "TG-cor": "TG-cor", "VEGA": r"$m_\tx{Vega}$", "M_sun(AB)": r"$M_{\odot,\tx{AB}}$", "CALIB": "Cal", "Lb_eff": r"$\lambda_0^{B_\nu}$ [$\mu$m]", "Fac_corr": r"$F_\tx{c}$"}
    colnames = [new_colnames[col] for col in colnames]
    df = pd.read_csv(fname, delim_whitespace=True,
                     skiprows=3, comment="#", names=colnames, skipfooter=1, engine='python')
    df["Context"] = df["Context"].apply(lambda x: int(2**(x - 1)))
    # for col in colnames:
    #     if "lambda" in col:
    # df[col] *= 1000  # In case we want the values to be in nm
    key_to_name_dict = {"FUV.pb": "FUV", "NUV.pb": "NUV", "r.pb": "r", "z.pb": "z", "newg.pb": "g", "wHSC_i.txt": "i-hsc",
                        "wHSC_i2.txt": "i2-hsc", "Y.lowres": "Y", "j.lowres": "J", "h.lowres": "H", "k.lowres": "Ks", "W1.res": "W1", "W2.res": "W2", "W3.res": "W3", "W4.res": "W4", "KiDSVIKING_aibn139_i.res": "i-kids"}
    key_to_name_dict = {"filters/" + key: val for key,
                        val in key_to_name_dict.items()}
    # for key, val in key_to_name_dict.items():
    #     vals = val.split("_")
    #     newval = vals[0] if len(vals) == 1 else f"{vals[0]}_\\tx{{{vals[1]}}}"
    #     key_to_name_dict[key] = f"{newval}"
    df = df.replace({"Band": key_to_name_dict})
    df["Survey"] = df["Band"].apply(
        lambda band: mt.SURVEY_NAME_DICT[mt.give_survey_for_band(band)])
    df["Band"] = df["Band"].apply(
        lambda band: mt.generate_pretty_band_name(band))
    return df


def save_filter_info(df, stem=""):
    """Saves the filter info overview in a LaTeX-readable document."""
    fname = mt.MYDIR + "other/" + stem + "_table_overview.tex"
    caption = r"""Overview on the important parameters for the table.\\
Here, the columns correspond to the mean wavelength $\Mean{\lambda}=\frac{\int R(\lambda)\lambda\dd \lambda}{\int R(\lambda)\dd\lambda}$, the Full Width at Half Maximum, the AB correction factor with $m_\tx{AB}=m_\tx{Vega} + m_\tx{AB, c}$, the Vega Magnitude $m_\tx{Vega}=-2.5\log_{10}\BracedIn{\frac{\int R(\lambda)\dd \lambda}{\int R(\lambda)F_\tx{Vega}(\lambda)\dd\lambda}}$, and the absolute magnitude of the sun in this filter.\\
In \LePhare{}, """
    cols = [col for col in df.columns if col not in [
        "Cal", "TG-cor", r"$\lambda_0^{B_\nu}$ [$\mu$m]", r"$F_\tx{c}$",
        r"$\lambda_\tx{eff}^\tx{Vega}$ [$\mu$m]"]]
    cols = cols[-1:] + cols[:-1]  # Move the 'Survey name' column to the front
    tabletext = latex_with_lines(df, na_rep="-", index=False, caption=caption,
                                 label="tab:filter_overview", position="htbp", float_format="{:0.2f}".format, columns=cols, escape=False)
    with open(fname, "w") as f:
        f.write(tabletext)
        LOGGER.info(
            f"The LaTeX input text for the filters has been written to {fname}")


if __name__ == "__main__":
    filter_df = c_f.read_filter_transmission_file(consider_context=True)
    produce_filter_plot(filter_df)

    # df = read_filter_overview_file("baseline_filters")
    # save_filter_info(df)
