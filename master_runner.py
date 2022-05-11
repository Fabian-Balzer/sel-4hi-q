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
    mt.LOGGER.debug(
        "Successfully ran the jython program to match and write the tables.")


def run_lephare_command(command, arg_dict, additional=""):
    """Runs a given LePhare command in the LePhare source file."""
    main_command = f"{mt.GEN_CONFIG['PATHS']['lepharedir']}/source/" + command
    run_string = main_command + " " + \
        " ".join([f"-{arg} {val}" for arg, val in arg_dict.items()]
                 ) + " " + additional
    mt.LOGGER.debug("Running the following shell command:\n%s", run_string)
    try:
        subprocess.run(run_string, check=True, shell=True)
    except subprocess.CalledProcessError as err:
        mt.LOGGER.error(
            "The following error was thrown when running the last shell command:\n%s", err)


def run_filters():
    """Runs the LePhare filter routine with the requested settings"""
    arg_dict = {"c": mt.give_parafile_fpath(),
                "FILTER_REP": f"{mt.GEN_CONFIG['PATHS']['params']}filters",
                "FILTER_FILE": mt.CUR_CONFIG["LEPHARE"]["filter_stem"]}
    additional = ">" + mt.give_filterfile_fpath()
    run_lephare_command("filter", arg_dict, additional)


def run_templates(ttype):
    """Runs the LePhare template routine with the requested settings"""
    arg_dict_sed = {"c": mt.give_parafile_fpath(),
                    "GAL_SED": mt.give_temp_listname(ttype),
                    "GAL_LIB": mt.give_temp_libname(ttype, "sed")}
    arg_dict_sed["t"] = "S" if ttype == "star" else "G"
    run_lephare_command("sedtolib", arg_dict_sed)
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
    run_lephare_command("mag_gal", arg_dict_sed)


if __name__ == "__main__":

    [LEPHARE]
    glb_context = -1
    input_stem = baseline_input
    output_stem = baseline_output
    give_stats = True
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
