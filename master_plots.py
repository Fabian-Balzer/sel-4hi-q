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

# input dataframes:
df = mt.read_plike_and_ext(prefix="matches/test2_",
                           suffix="_processed_table.fits")
df = mt.add_mag_columns(df)
av.plot_r_band_magnitude(df)
av.plot_input_distribution(df)

# output dataframes:
# df = mt.read_plike_and_ext(prefix="lephare_output/test2_", suffix=".fits")
