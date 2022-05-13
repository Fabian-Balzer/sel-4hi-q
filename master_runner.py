# -*- coding: utf-8 -*-
"""
Created on Fri Feb 18 13:40:36 2022

maybe needed:
%load_ext autoreload
%autoreload 2

@author: fabian_balzer
"""

# %%
import util.my_tools as mt
from util.assert_config import assert_all
from input_scripts.filter_coverage_plot import read_filter_overview_file


def assemble_catalog():
    """Runs the jython script to match the input files."""
    mt.run_jystilts_program("match_tables.py")
    mt.LOGGER.debug(
        "Successfully ran the jython program to match and write the tables.")


def run_filters():
    """Runs the LePhare filter routine with the requested settings"""
    arg_dict = {"c": mt.give_parafile_fpath(),
                "FILTER_REP": f"{mt.GEN_CONFIG['PATHS']['params']}filters",
                "FILTER_FILE": mt.CUR_CONFIG["LEPHARE"]["filter_stem"]}
    additional = ">" + mt.give_filterfile_fpath()
    mt.run_lephare_command("filter", arg_dict, additional)
    info_df = fc.read_filter_overview_file()
    fc.save_filter_info(info_df)


def run_templates(ttype):
    """Runs the LePhare template routine with the requested settings"""
    arg_dict_sed = {"c": mt.give_parafile_fpath(),
                    "GAL_SED": mt.give_temp_listname(ttype),
                    "GAL_LIB": mt.give_temp_libname(ttype, "sed")}
    arg_dict_sed["t"] = "S" if ttype == "star" else "G"
    mt.run_lephare_command("sedtolib", arg_dict_sed)
    arg_dict_mag = {"c": mt.give_parafile_fpath(),
                    "GAL_LIB_IN": mt.give_temp_libname(ttype, "sed"),
                    "GAL_LIB_OUT": mt.give_temp_libname(ttype, "mag"),
                    "EM_LINES": "NO",
                    "LIB_ASCII": "YES"}
    arg_dict_mag["t"] = "S" if ttype == "star" else "G"
    if ttype == "pointlike":
        arg_dict_mag["EXTINC_LAW"] = "SMC_prevot.dat"
        arg_dict_mag["MOD_EXTINC"] = "11,23"
        arg_dict_mag["EB_V"] = "0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4"
    mt.run_lephare_command("mag_gal", arg_dict_sed)
    mt.run_jystilts_program("rewrite_fits_header", args=[ttype, "MAG"])


def run_zphota(ttype):
    """Runs the LePhare zphota routine with the requested settings"""
    arg_dict_sed = {"c": mt.give_parafile_fpath(),
                    "ZPHOTLIB": f"{mt.give_temp_libname(ttype, include_path=False)},{mt.give_temp_libname('star', include_path=False)}",
                    "CAT_IN": mt.give_lephare_filename(ttype),
                    "CAT_OUT": mt.give_lephare_filename(ttype, out=True),
                    "PARA_OUT": mt.give_parafile_fpath(out=True),
                    "GLB_CONTEXT": mt.CONTEXT,
                    "PDZ_OUT": mt.give_lephare_filename(ttype, suffix="")}
    if ttype == "pointlike":
        arg_dict_sed["MAG_REF"] = "7"
        arg_dict_sed["MAG_ABS"] = "-30,-20"
    mt.run_lephare_command("zphota", arg_dict_sed)
    mt.run_jystilts_program("rewrite_fits_header", args=[ttype, "OUT"])


if __name__ == "__main__":
    mt.LOGGER.info("Program started with the following requests:")
    cat_config = mt.CUR_CONFIG["CAT_ASSEMBLY"]
    if cat_config.getboolean("assemble_cat"):
        mt.LOGGER.info("Catalogue assembly with '%s' as a stem:",
                       cat_config["cat_stem"])
        for boolkey in ["use_matched", "use_processed", "reduce_to_specz", "write_lephare_input", "write_info_file"]:
            val = cat_config.getboolean(boolkey)
            mt.LOGGER.info("%s:\n\t\t%s", boolkey, str(val))
    lep_config = mt.CUR_CONFIG["LEPHARE"]
    if lep_config.getboolean("run_filters"):
        mt.LOGGER.info("LePhare filter run with '%s' as a stem.",
                       lep_config["filter_stem"])
    if lep_config.getboolean("run_templates"):
        mt.LOGGER.info("LePhare template run with '%s' as a stem.",
                       lep_config["template_stem"])
    if lep_config.getboolean("run_zphota"):
        mt.LOGGER.info("LePhare zphota run with '%s' as input and '%s' as output stem.",
                       lep_config["input_stem"], lep_config["output_stem"])
        mt.LOGGER.info("The provided global context is %s, corresponding to the following bands:\n%s",
                       mt.CONTEXT, mt.give_bands_for_context(mt.CONTEXT))
        if lep_config.getboolean("give_stats"):
            # TODO Stats file?
            mt.LOGGER.info(
                "Statistics about the LePhare run are going to be provided.")
    assert_all()

    if mt.CUR_CONFIG["CAT_ASSEMBLY"].getboolean("assemble_cat"):
        assemble_catalog()

    if mt.CUR_CONFIG["LEPHARE"].getboolean("run_filters"):
        run_filters()

    if mt.CUR_CONFIG["LEPHARE"].getboolean("run_templates"):
        if mt.CUR_CONFIG["GENERAL"].getboolean("use_pointlike"):
            run_templates("pointlike")
        if mt.CUR_CONFIG["GENERAL"].getboolean("use_extended"):
            run_templates("extended")
        run_templates("star")


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
