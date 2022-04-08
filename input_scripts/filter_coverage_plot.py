# -*- coding: utf-8 -*-
"""
Created on Wed Feb 23 12:20:24 2022

@author: fabian_balzer
"""

# %%

import os
from glob import escape

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.offsetbox import AnchoredText
from scipy.stats import gaussian_kde

import configure_matplotlib as cm
import my_tools as mt


def read_filter_info_file(fname):
    """Reads and returns a .filt file from LePhare and converts it to a pandas
    dataframe."""
    filter_wls = {}
    filter_trans = {}
    with open(fname, "r") as f:
        for line in f.readlines():
            if line.startswith("# "):
                name = line.strip("# ")[:-1]
            else:
                if name in filter_wls:
                    filter_trans[name] = [
                        float(val) for val in line.strip("\n").split(", ")]
                else:
                    filter_wls[name] = [float(val)
                                        for val in line.strip("\n").split(", ")]
    df_dict = {}
    for key in filter_wls.keys():
        filter_df = pd.DataFrame(
            {"Wavelength": filter_wls[key], "Transmission": filter_trans[key]})
        df_dict[f"{filter_df['Wavelength'].mean():6.0f}>{key}"] = filter_df
    df = pd.concat(df_dict, axis=1)

    return df


def produce_filter_plot(df):
    """Create a plot showing the normalised transmission curves of the filters."""
    fig, ax = plt.subplots(figsize=cm.set_figsize(fraction=0.912))
    key_to_name_dict = {"FUV": "FUV_GALEX", "NUV": "NUV_GALEX", "r_decals": "r_DECam", "z_decals": "z_DECam", "decam_g": "g_DECam", "HSC_i2": "i_2, HSC",
                        "HSC-i": "i_HSC", "Y_VISTA": "Y_VISTA", "J_VISTA": "J_VISTA", "H_VISTA": "H_VISTA", "K_VISTA": "Ks_VISTA", "W1": "W_1", "W2": "W_2", "W3": "W_3", "W4": "W_4"}
    for key in df.keys().get_level_values(0).drop_duplicates().sort_values():
        x_data = df[key]["Wavelength"]
        y_data = df[key]["Transmission"]
        key = key.split(">")[1]
        key = key_to_name_dict[key]
        if "W_" in key:
            continue
        # Scale the y-data so the maximum is at 1
        while max(y_data) > 1:
            y_data = y_data / 10
        # y_data = y_data / max(y_data)
        keys = key.split("_")
        key = keys[0] if len(keys) == 1 else f"{keys[0]}_{{\\rm {keys[1]}}}"
        ax.plot(x_data, y_data, label=f"${key}$")
        ax.fill_between(x_data, y_data, alpha=0.5)
    # ax.set_xscale("log")
    ax.set_ylim(0, 1.3)
    ax.grid(True)
    ax.set_xlabel("$\lambda$ [$\AA$]")
    ax.set_ylabel("Normalised filter transmission")
    ax.legend(ncol=4, prop={'size': "x-small"})
    cm.save_figure(fig, "data/plots/filter_coverage_plot")


def latex_with_lines(df, *args, **kwargs):
    kwargs['column_format'] = ''.join([''] + ['l'] * df.index.nlevels
                                      + ['l'] * df.shape[1] + [''])
    res = df.to_latex(*args, **kwargs)
    return res


def read_filter_overview_file(fname):
    """Read the overview file and store the information in a dataframe. Format the column names to a LaTeX-readable format."""
    with open(fname, "r") as f:
        for line in f.readlines():
            if line.startswith("# "):
                colnames = line.strip("# ").split()
    new_colnames = {"IDENT": "Id", "NAME": "Band", "Lbda_mean": r"$\Mean{\lambda}$ [$\mu$m]",
                    "Lbeff(Vega)": r"$\lambda_\tx{eff}^\tx{Vega}$ [$\mu$m]", "FWHM": r"FWHM [$\mu$m]", "AB-cor": r"$m_\tx{AB, c}$", "TG-cor": "TG-cor", "VEGA": r"$m_\tx{Vega}$", "M_sun(AB)": r"$M_{\odot,\tx{AB}}$", "CALIB": "Cal", "Lb_eff": r"$\lambda_0^{B_\nu}$ [$\mu$m]", "Fac_corr": r"$F_\tx{c}$"}
    colnames = [new_colnames[col] for col in colnames]
    df = pd.read_csv(fname, delim_whitespace=True,
                     skiprows=3, comment="#", names=colnames, skipfooter=1, engine='python')
    # for col in colnames:
    #     if "lambda" in col:
    # df[col] *= 1000  # In case we want the values to be in nm
    key_to_name_dict = {"FUV.pb": "FUV_GALEX", "NUV.pb": "NUV_GALEX", "r.pb": "r_DECam", "z.pb": "z_DECam", "newg.pb": "g_DECam", "HSC_i2": "i_2, HSC",
                        "HSC-i": "i_HSC", "Y.lowres": "Y_VISTA", "j.lowres": "J_VISTA", "h.lowres": "H_VISTA", "k.lowres": "Ks_VISTA", "W1.res": "W_1", "W2.res": "W_2", "W3.res": "W_3", "W4.res": "W_4"}
    key_to_name_dict = {"filters/" + key: val for key,
                        val in key_to_name_dict.items()}
    for key, val in key_to_name_dict.items():
        vals = val.split("_")
        newval = vals[0] if len(vals) == 1 else f"{vals[0]}_\\tx{{{vals[1]}}}"
        key_to_name_dict[key] = f"${newval}$"
    df = df.replace({"Band": key_to_name_dict})
    return df


def save_filter_info(df, fname):
    """Saves the filter info overview in a LaTeX-readable document."""
    caption = r"""Overview on the important parameters for the table.\\
Here, the columns correspond to the mean wavelength $\Mean{\lambda}=\frac{\int R(\lambda)\lambda\dd \lambda}{\int R(\lambda)\dd\lambda}$, the effective wavelength w. r. t. Vega $\lambda_\tx{eff}=\frac{\int R(\lambda)\lambda F_\tx{Vega}(\lambda)\dd \lambda}{\int R(\lambda) F_\tx{Vega}(\lambda)\dd\lambda}$, the Full Width at Half Maximum, the AB correction factor with $m_\tx{AB}=m_\tx{Vega} + m_\tx{AB, c}$, the Vega Magnitude $m_\tx{Vega}=-2.5\log_{10}\BracedIn{\frac{\int R(\lambda)\dd \lambda}{\int R(\lambda)F_\tx{Vega}(\lambda)\dd\lambda}}$, the absolute magnitude of the sun in this filter, and the correction factor that is applied in \LePhare{}."""
    tabletext = latex_with_lines(df, na_rep="-", index=False, caption=caption,
                                 label="tab:filter_overview", position="htbp", float_format="{:0.2f}".format, columns=[col for col in df.columns if col not in ["Cal", "TG-cor", r"$\lambda_0^{B_\nu}$ [$\mu$m]"]], escape=False)
    with open(fname, "w") as f:
        f.write(tabletext)
        print(
            f"The LaTeX input text for the filters has been written to {fname}")


if __name__ == "__main__":
    df = read_filter_info_file(mt.DATAPATH + "/Input_data/filter_info.filt")
    produce_filter_plot(df)
    # df = read_filter_overview_file(
    #     mt.DATAPATH + "/Input_data/filter_overview.list")
    # save_filter_info(df, fname="Data/filter_overview.tex")
