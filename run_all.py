# -*- coding: utf-8 -*-
"""
Created on Fri Feb 18 13:40:36 2022

@author: fabian_balzer
"""
# %%

import util.my_tools as mt
import util.runner_commands as r_c
from input_scripts.input_plot_container import InputPlotContainer
from output_scripts.output_plot_container import OutputPlotContainer
from util.assert_config import assert_all

# %%


if __name__ == "__main__":
    mt.log_run_info()
    assert_all()

    if mt.CUR_CONFIG["CAT_ASSEMBLY"].getboolean("assemble_cat"):
        r_c.assemble_catalog()

    r_c.run_lephare_commands()

    # Input-related plots:
    i_p_c = InputPlotContainer(True)
    plot_con = mt.CUR_CONFIG["PLOTTING"]

    if plot_con.getboolean("input"):
        i_p_c.plot_input_dist(context=mt.CONTEXT)
        i_p_c.plot_magnitude_dist(context=mt.CONTEXT)

    if plot_con.getboolean("sep"):
        radius_dict = {"vhs": 0.5, "eros": 0.1, "hsc": 0.25,
                       "galex": 3.5, "kids": 1.5}
        i_p_c.plot_separation(radius_dict)

    if plot_con.getboolean("filters"):
        i_p_c.plot_filters()

    # Output-related plots:
    o_p_c = OutputPlotContainer(True)

    if plot_con.getboolean("output"):
        o_p_c.plot_specz_photo_z()
        o_p_c.plot_template_numbers()

    if plot_con.getboolean("template"):
        o_p_c.plot_color_vs_redshift("g", "r")
        o_p_c.plot_template_scores()

    # # %%
    # # Construct the input dataframe:
    # input_df = mt.read_plike_and_ext(prefix="matches/test2_",
    #                                  suffix="_processed_table.fits")
    # input_df = mt.add_mag_columns(input_df)
    # # av.plot_r_band_magnitude(df)
    # av.plot_input_distribution(input_df)

    # # %% Filter analysis:
    # filter_df = fc.read_filter_info_file()
    # fc.produce_filter_plot(filter_df)
    # info_df = fc.read_filter_overview_file()
    # fc.save_filter_info(info_df)
