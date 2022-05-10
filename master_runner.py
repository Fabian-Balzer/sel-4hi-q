# -*- coding: utf-8 -*-
"""
Created on Fri Feb 18 13:40:36 2022

maybe needed:
%load_ext autoreload
%autoreload 2

@author: fabian_balzer
"""

# %%


import shlex
import subprocess

import util.my_tools as mt
from util.assert_config import assert_all


def assemble_catalog():
    """Runs the jython script to match the input files."""
    run_jystilts = f"java -jar {mt.GEN_CONFIG['PATHS']['JYSTILTS']}"
    scriptpath = f"{mt.GEN_CONFIG['PATHS']['scripts']}jystilts_scripts/"
    match_table_string = f"{run_jystilts} '{scriptpath}match_tables.py'"
    subprocess.call(shlex.split(match_table_string))
    mt.LOGGER.info(
        "Successfully ran the jython program to match and write the tables.")
    return


if __name__ == "__main__":
    assert_all()

    # if mt.CUR_CONFIG["CAT_ASSEMBLY"].getboolean("assemble_cat"):
    #     assemble_catalog()


# %%
# Construct the input dataframe:
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
