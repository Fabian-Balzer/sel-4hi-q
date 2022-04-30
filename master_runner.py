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
    parser.add_argument("-b", "--build_dirs", action="store_true",
                        help="Builds the directories necessary for running the code.")
    stem_defaults = {"input": "baseline_input", "separation": "baseline_input",
                     "Filter": "compiled_filters", "output": "test", "template": "baseline_templates"}
    for argtype in ["input", "Filter", "output", "template"]:
        # Introduce a boolean call on whether to produce any of the datasets:
        short = f"-{argtype[0]}"  # i, s, F, o and t can be specified now.
        long = f"--produce_{argtype}_data"
        helpstring = f"Specify whether to run scripts producing the {argtype} data"
        parser.add_argument(short, long, action="store_true",
                            help=helpstring, default=True)
        # Allow to specify the stem names:
        argname = f"--{argtype}_stem"
        default = stem_defaults[argtype]
        helpstring = f"Specify the stem name of the {argtype} data."
        parser.add_argument(argname, help=helpstring, default=default)
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Increase output verbosity.")
    args, _ = parser.parse_known_args()
    for argtype in [argtype.split('_')[1] for argtype, val in vars(args).items() if isinstance(val, bool) and val]:
        print(f"{argtype}_stem: {vars(args)[f'{argtype}_stem']}")
    return args


args = read_args()

# %%
if args.build_dirs:
    mt.init_plot_directory()

bands = [band for band in mt.BAND_LIST if band not in [
    "i_hsc", "i2_hsc", "i_kids"]]
context = mt.give_context(bands)
filters = mt.convert_context_to_filter_indices(context)
print(f"Use {context} for only using {bands} ({filters})")
# %%
# input dataframes:
input_df = mt.read_plike_and_ext(prefix="matches/test2_",
                                 suffix="_processed_table.fits")
input_df = mt.add_mag_columns(input_df)
# av.plot_r_band_magnitude(df)
av.plot_input_distribution(input_df)

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
