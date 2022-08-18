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
from cycler import cycler
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
    ymax = 3

    # Plot two filters for comparison
    df = c_f.read_filter_transmission_file()
    for band, color in zip(["g", "r"], ["k", "red"]):
        key = [key for key in df.keys().get_level_values(0)
               if f">{band}>" in key][0]
        x, y = df[key]["lambda"] / 10, df[key]["trans"] * ymax
        ax.fill_between(x, y, alpha=0.4, color=color)
        label = mt.generate_pretty_band_name(band)
        ax.text(float(key.split(">")[0]) * 0.83, -
                0.15, label + " band", fontweight="bold")
    test_z = 0.5
    test_z2 = 4
    custom_cycler = (cycler(color=['r', 'b', 'm', 'g', "k", "orange", "gray"]))

    # ax.set_prop_cycle(custom_cycler)

    i = 0
    for fname, label in templates.items():
        try:
            fpath = mt.GEN_CONFIG["PATHS"]["data"] + "seds/" + fname
            df = get_sed_data(fpath)
            df = df[(df["lambda"] < 10000) & (df["lambda"] > 0.01)]
            df["flux"] /= (df["flux"].max())
            df["flux"] = df["flux"].replace(
                {'0': np.nan, 0: np.nan}) + i * 0.4 + 0.1
            p = ax.plot(df["lambda"], df["flux"],
                        label=label, linewidth=0.6)[0]
            if label in ["S0", "Spiral bar", "Starburst (I22491)"]:
                lamb_redshift = df["lambda"] * (1 + test_z)
                color = p.get_color()
                ax.plot(lamb_redshift, df["flux"],
                        label=f"{label} at $z={test_z}$", linewidth=1, linestyle=":", color=color)
            if label in ["QSO", "Quasar"]:
                lamb_redshift = df["lambda"] * (1 + test_z2)
                color = p.get_color()
                ax.plot(lamb_redshift, df["flux"],
                        label=f"{label} at $z={test_z2}$", linewidth=1, linestyle=":", color=color)

        except FileNotFoundError:
            mt.LOGGER.warning("Couldn't locate '%s'", fname)
        i += 1
    ax.set_xlim(70, 3e3)
    ax.set_ylim(0.01, ymax)
    ax.set_xlabel(r"$\lambda$ [nm]")
    ax.set_ylabel("$F$ [arbitrary units]")
    ax.set_yticklabels([])
    lyman_x, balmer_x = 91.2, 364.6
    ax.axvline(lyman_x, linestyle="--", color="k")
    ax.text(lyman_x * 1.02, ymax * 0.95, "Lyman Break")
    ax.axvline(balmer_x, linestyle="--", color="k")
    ax.text(balmer_x * 1.02, ymax * 0.95, "Balmer Break")
    ax.axvline(balmer_x * (1 + test_z), linestyle="--",
               color="k", linewidth=0.5)

    ax.text(balmer_x * (1 + test_z) * 1.02, ymax * 0.53,
            f"Balmer Break at $z={test_z}$")

    ax.grid(True)
    ax.legend(loc="upper right")
    cm.save_figure(fig, "SED_comparison_plot_with_z",
                   "presentation", format_="png")


if __name__ == "__main__":
    templates = {
        'pl_QSO_DR2_029_t0.spec': "QSO",
        "S0_template_norm.sed": "S0",
        'Ell1_A_0.sed': "Ell",
        # 'M82_template_norm.sed': "M82",
        # 'NGC4151_Central_64.00_NGC4579.dat': "Composite Sey 1",
        'Sb_A_0.sed': "Sb",
        # 'Sc_A_0.sed': "Sc",
        'I22491_template_norm.sed.save': "Starburst (I22491)",
        'Sey2_template_norm.sed': "Sey 2",
        #  'Spi4_template_norm.sed': "Spiral",
    }

    templates = {
        'pl_QSO_DR2_029_t0.spec': "Quasar",
        # "S0_template_norm.sed": "S0",
        # 'Ell1_A_0.sed': "Ell",
        # 'M82_template_norm.sed': "M82",
        # 'NGC4151_Central_64.00_NGC4579.dat': "Composite Sey 1",
        'Sb_A_0.sed': "Spiral bar",
        # 'Sc_A_0.sed': "Sc",
        'I22491_template_norm.sed.save': "Starburst (I22491)",
        # 'Sey2_template_norm.sed': "Sey 2",
        #  'Spi4_template_norm.sed': "Spiral",
    }
    # templates = {f"MARA23032010/S0_{100-a}_QSO2_{a}.sed":
    #              f"{100-a} S0 to {a} QSO" for a in range(0, 100, 10)}
    plot_seds("pointlike", templates)

    # temp_num_dict = {mt.get_temp_num_for_name("pointlike",
    #                                           temp): label for temp, label in templates.items()}
    # o_p_c = OutputPlotContainer()
    # o_p_c.ext = False
    # o_p_c.produce_both = False
    # o_p_c.plot_color_vs_redshift(
    #     "g", "r", temp_nums=temp_num_dict, plot_sources=False, fitted_only=False)

# %%
