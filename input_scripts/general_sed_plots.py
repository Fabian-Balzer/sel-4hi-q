# -*- coding: utf-8 -*-
"""
Created on Sat Feb 19 17:37:14 2022

@author: fabian_balzer

Script for plotting the SEDs of selected templates
"""

# %%
from math import ceil, floor

import matplotlib.pyplot as plt
import numpy as np
import output_scripts.template_analysis as t_a
import pandas as pd
import util.configure_matplotlib as cm
import util.my_tools as mt
from output_scripts.output_plot_container import OutputPlotContainer

import input_scripts.collect_filter_data as c_f


def get_sed_data(fname: str):
    df = pd.read_csv(fname, delim_whitespace=True, names=["lambda", "flux"])
    df["lambda"] /= 10  # Convert from Angstrom to nm
    return df


def plot_seds(ttype: str, templates: dict):
    # fname = mt.give_temp_listname(ttype, "baseline")
    # with open(fname, "r") as f:
    #     fnames = [line.replace("\n", "")
    #               for line in f.readlines() if not line.startswith("#")]
    fig, ax = plt.subplots(1, 1, figsize=cm.set_figsize())
    ax.set_xscale("log")
    # ax.set_yscale("log")
    # Plot two filters for comparison
    df = c_f.read_filter_transmission_file()
    for band in ["g", "r"]:
        key = [key for key in df.keys().get_level_values(0)
               if f">{band}>" in key][0]
        x, y = df[key]["lambda"] / 10, df[key]["trans"]
        ax.fill_between(x, y, alpha=0.5)
        label = mt.generate_pretty_band_name(band)
        ax.text(float(key.split(">")[0]) * 0.9, y.max() - 0.04, label)
    test_z = 0.5

    for fname, label in templates.items():
        try:
            fpath = "../" + fname
            df = get_sed_data(fpath)
            df = df[(df["lambda"] < 10000) & (df["lambda"] > 50)]
            df["flux"] /= (df["flux"].max())
            df["flux"] = df["flux"].replace({'0': np.nan, 0: np.nan})
            p = ax.plot(df["lambda"], df["flux"],
                        label=label, linewidth=0.6)[0]
            if label in ["QSO", "S0"]:
                lamb_redshift = df["lambda"] * (1 + test_z)
                color = p.get_color()
                ax.plot(lamb_redshift, df["flux"],
                        label=f"{label} at $z={test_z}$", linewidth=0.6, linestyle="--", color=color)
        except FileNotFoundError:
            mt.LOGGER.warning("Couldn't locate '%s'", fname)

    ax.set_xlim(70, 3e3)
    ax.set_ylim(0.01, 1.2)
    ax.set_xlabel("$\lambda$ [nm]")
    ax.set_ylabel("$F$ [arbitrary units]")
    ax.set_yticklabels([])
    lyman_x, balmer_x = 91.2, 364.6
    ax.axvline(lyman_x, linestyle="--", color="k")
    ax.text(lyman_x * 1.02, 1.13, "Lyman Break")
    ax.axvline(balmer_x, linestyle="--", color="k")
    ax.text(balmer_x * 1.02, 1.13, "Balmer Break")
    ax.axvline(balmer_x * (1 + test_z), linestyle="--",
               color="k", linewidth=0.5)
    ax.text(balmer_x * (1 + test_z) * 1.02, 1.05,
            f"Balmer Break at $z={test_z}$")
    ax.grid(True)
    ax.legend(loc="lower right")


if __name__ == "__main__":
    templates = {"MARA23032010/S0_template_norm.sed": "S0",
                 #   "MARA23032010/Sb_template_norm.sed": "SB",
                 #   "MARA23032010/Sey18_template_norm.sed": "Sey 1.8",
                 "MARA23032010/pl_QSOH_template_norm.sed": "QSO",
                 #  "MARA23032010/S0_50_QSO2_50.sed": "50/50 S0 to QSO",
                 #  "COSMOS_SED/Ell4_A_0.sed ": "Ell",
                 #  "MARA23032010/Spi4_template_norm.sed": "Spi4",
                 #  "BROWN_TONIMA/NGC4151_Central_temp_restframe.dat": "AGN (NGC4151 Central)"
                 "MARA23032010/M82_template_norm.sed": "M82",
                 "MARA23032010/I22491_template_norm.sed.save": "I22491"
                 }
    # templates = {f"MARA23032010/S0_{100-a}_QSO2_{a}.sed":
    #              f"{100-a} S0 to {a} QSO" for a in range(0, 100, 10)}
    # plot_seds("extended", templates)
    temp_num_dict = {mt.get_temp_num_for_name("pointlike",
                                              temp): label for temp, label in templates.items()}
    o_p_c = OutputPlotContainer()
    o_p_c.ext = False
    o_p_c.produce_both = False
    o_p_c.plot_color_vs_redshift(
        "g", "r", temp_nums=temp_num_dict, plot_sources=False, fitted_only=False)

# %%
