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
import my_tools as mt


# input dataframes:
df = mt.read_plike_and_ext("data/output/test2_pointlike.fits",
                           "data/output/test2_extended.fits")


# output dataframes:
df = mt.read_plike_and_ext("data/output/test2_pointlike.fits",
                           "data/output/test2_extended.fits")
