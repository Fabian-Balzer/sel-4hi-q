# -*- coding: utf-8 -*-
"""
Created on Fri Feb 18 13:40:36 2022

maybe needed:
%load_ext autoreload
%autoreload 2

@author: fabian_balzer
"""

# %%


import input_scripts.availability_plots as av
import input_scripts.filter_coverage_plot as fc
import input_scripts.separation_plots as sep
import output_scripts.specz_photz_plots as s_p
import output_scripts.template_analysis_plots as ta
import util.my_tools as mt

STEM_OUT = "without_i"
# %%
mt.init_plot_directory()

bands = [band for band in mt.BAND_LIST if band not in [
    "i-hsc", "i2-hsc", "i-kids"]]
context = mt.give_context(bands)
filters = mt.convert_context_to_band_indices(context)
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
