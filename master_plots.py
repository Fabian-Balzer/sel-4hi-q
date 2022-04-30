# -*- coding: utf-8 -*-
"""
Created on Fri Feb 18 13:40:36 2022

maybe needed:
%load_ext autoreload
%autoreload 2

@author: fabian_balzer
"""

# %%


import argparse
import os

import input_scripts.availability_plots as av
import input_scripts.filter_coverage_plot as fc
import input_scripts.separation_plots as sep
import output_scripts.specz_photz_plots as s_p
import output_scripts.template_analysis_plots as ta
import util.my_tools as mt


def read_args():
    """Reads out the arguments given by the user."""
    parser = argparse.ArgumentParser(
        description="Assess a LePhare output file.")
    stem_defaults = {"input": "baseline_input", "separation": "baseline_input",
                     "Filter": "compiled_filters", "output": "test", "template": "baseline_templates"}
    for argtype in ["input", "separation", "Filter", "output", "template"]:
        # Introduce a boolean call:
        short = f"-{argtype[0]}"  # i, s, F, o and t can be specified now.
        long = f"--produce_{argtype}_plots"
        helpstring = f"Specify wether to save the {argtype} plots"
        parser.add_argument(short, long, action="store_true",
                            help=helpstring, default=True)
        # Allow to specify the stem names:
        argname = f"--{argtype}_stem"
        default = stem_defaults[argtype]
        helpstring = f"Specify the stem name of the {argtype} data."
        parser.add_argument(argname, help=helpstring, default=default)
    parser.add_argument("--context", action="store_const",
                        help="The context for the LePhare run.", const=-1, type=int)
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Increase output verbosity.")
    args, _ = parser.parse_known_args()
    for argtype in [argtype.split('_')[1] for argtype, val in vars(args).items() if isinstance(val, bool) and val]:
        print(f"{argtype}_stem: {vars(args)[f'{argtype}_stem']}")
    return args


args = read_args()

# %%
# input dataframes:
if args.produce_input_plots:
    input_df = mt.read_plike_and_ext(prefix=f"matches/{args.input_stem}",
                                     suffix="_processed_table.fits")
    input_df = mt.add_mag_columns(input_df)
    av.plot_r_band_magnitude(input_df, args.input_stem)
    av.plot_input_distribution(input_df, args.input_stem)
    av.plot_band_number_distribution(input_df, args.input_stem)

# %% Filter analysis:
filter_df = fc.read_filter_info_file()
fc.produce_filter_plot(filter_df)
info_df = fc.read_filter_overview_file()
fc.save_filter_info(info_df)

# %% Separation plots:
# sep.plot_all_separations(input_df, verbose=True)


# %%

# output dataframes:
STEM_OUT = "without_i_suraj"
output_df = mt.read_plike_and_ext(
    prefix=f"lephare_output/{STEM_OUT}_", suffix=".fits")
output_df = mt.add_filter_columns(output_df)

for ttype in ["pointlike", "extended"]:
    # ta.plot_problematic_templates(output_df, ttype)
    s_p.plot_photoz_vs_specz(output_df, ttype, STEM_OUT, comparison=True)
