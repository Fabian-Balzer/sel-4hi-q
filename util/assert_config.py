"""Reads and asserts the config file.
Raises errors whenever a directory or file is not set up properly
(i. e. if a plot for filters is requested, but the filter file is not there."""
# %%


import os
from os.path import isfile
from sys import path

path.append(os.environ["LEPHARE"] + "/lephare_scripts")

import util.my_tools as mt

use_plike = mt.CUR_CONFIG["GENERAL"].getboolean("use_pointlike")
use_ext = mt.CUR_CONFIG["GENERAL"].getboolean("use_extended")

# %% General assertions


def assert_general():
    """General assertions (checking the fpath availability and context)"""
    fpath = mt.GEN_CONFIG["PATHS"]["config"] + \
        mt.GEN_CONFIG["PATHS"]["current_config"]
    mt.assert_file_exists(fpath, "configuration")
    assert (-1 <= mt.CONTEXT <= 2**(len(mt.BAND_LIST) + 1)
            ), f"Context {mt.CONTEXT} is not a viable context."


# %% Assert catalog availability:
def assert_catalog_assembly():
    """Assertions needed for running the catalog assembly"""
    cat_stem = mt.CUR_CONFIG["CAT_ASSEMBLY"]["cat_stem"]
    if mt.CUR_CONFIG["CAT_ASSEMBLY"].getboolean("assemble_cat"):
        assert mt.GEN_CONFIG["PATHS"]["cat"] != ""
        available_cats = os.listdir(mt.GEN_CONFIG["PATHS"]["cat"])
        for survey in mt.GEN_CONFIG["BAND_DICT"].keys():
            if survey != "galex":
                assert survey in available_cats

        if mt.CUR_CONFIG["CAT_ASSEMBLY"].getboolean("use_matched"):
            mt.assert_file_exists(mt.give_match_table_name(), "matched")

        if mt.CUR_CONFIG["CAT_ASSEMBLY"].getboolean("use_processed"):
            mt.assert_file_exists(mt.give_processed_table_name(
                "pointlike"), "processed pointlike")
            mt.assert_file_exists(mt.give_processed_table_name(
                "extended"), "processed extended")

        if mt.CUR_CONFIG["CAT_ASSEMBLY"].getboolean("write_lephare_input"):
            if use_plike:
                mt.assert_file_overwrite(mt.give_lephare_filename("pointlike"))
            if use_ext:
                mt.assert_file_overwrite(mt.give_lephare_filename("extended"))

# %% Assert LePhare stuff


def assert_lephare_assembly():
    """Assertions needed for running lephare"""
    filt = mt.CUR_CONFIG["LEPHARE"].getboolean("run_filters")
    temps = mt.CUR_CONFIG["LEPHARE"].getboolean("run_templates")
    zphot = mt.CUR_CONFIG["LEPHARE"].getboolean("run_zphota")
    if filt | temps | zphot:
        mt.assert_file_exists(mt.give_parafile_fpath(), "input parameter")

    if temps:
        if use_plike:
            mt.assert_file_exists(mt.give_temp_listname(
                "pointlike"), "pointlike template list")
        if use_ext:
            mt.assert_file_exists(mt.give_temp_listname(
                "extended"), "extended template list")
        mt.assert_file_exists(mt.give_temp_listname(
            "star"), "star template list")

    if zphot:
        mt.assert_file_exists(mt.give_parafile_fpath(
            out=True), "output parameter")
        # If no input cat is produced, check if it's available:
        if not mt.CUR_CONFIG["CAT_ASSEMBLY"].getboolean("write_lephare_input"):
            if use_plike:
                mt.assert_file_exists(mt.give_lephare_filename(
                    "pointlike"), "pointlike LePhare input")
            if use_ext:
                mt.assert_file_exists(mt.give_lephare_filename(
                    "extended"), "extended LePhare input")
        if not filt:
            # TODO: Look in LePhare instead
            mt.assert_file_exists(mt.give_filterfile_fpath(), "filter")
        # If no templates are generated, check if they are available:
        if not temps:
            if use_plike:
                mt.assert_file_exists(mt.give_temp_libname(
                    "pointlike", "mag", suffix=".doc", use_workpath=True), "pointlike template")
            if use_ext:
                mt.assert_file_exists(mt.give_temp_libname(
                    "extended", "mag", suffix=".doc", use_workpath=True), "extended template")
            mt.assert_file_exists(mt.give_temp_libname(
                "star", "mag", suffix=".doc", use_workpath=True), "star template")


# %% Assert plotting:
def assert_plots():
    """Assertions needed before running the plotting functions"""
    pass
# TODO!
# mt.CUR_CONFIG["PLOTTING"].getboolean("input")
# mt.CUR_CONFIG["PLOTTING"].getboolean("sep")
# mt.CUR_CONFIG["PLOTTING"].getboolean("filters")
# mt.CUR_CONFIG["PLOTTING"].getboolean("output")
# mt.CUR_CONFIG["PLOTTING"].getboolean("template")


def assert_all():
    """Goes through all config assertion functions."""
    assert_general()
    assert_catalog_assembly()
    assert_lephare_assembly()
    assert_plots()
    mt.LOGGER.info("All assertions have been successful.")
