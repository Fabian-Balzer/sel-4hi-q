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
    assert isfile(fpath), f"Unable to find config file at '{fpath}'"
    context = mt.CUR_CONFIG["LEPHARE"].getint("context")
    assert (-1 <= context <= (len(mt.BAND_LIST) + 1) **
            2), f"Context {context} is not a viable context."


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
            matchname = cat_stem + "_latest_match.fits"
            assert isfile(mt.GEN_CONFIG["PATHS"]["match"] + matchname)

        if mt.CUR_CONFIG["CAT_ASSEMBLY"].getboolean("use_processed"):
            procname = cat_stem + "pointlike_processed.fits"
            assert isfile(mt.GEN_CONFIG["PATHS"]["match"] + procname)
            procname = cat_stem + "pointlike_extended.fits"
            assert isfile(mt.GEN_CONFIG["PATHS"]["match"] + procname)

        input_stem = mt.CUR_CONFIG["LEPHARE"]["input_stem"]
        if mt.CUR_CONFIG["CAT_ASSEMBLY"].getboolean("write_lephare_input"):
            if use_plike:
                inputname = input_stem + "_pointlike.in"
                fpath = mt.GEN_CONFIG["PATHS"]["data"] + \
                    "lephare_input/" + inputname
                mt.assert_file_overwrite(fpath)
            if use_ext:
                inputname = input_stem + "_extended.in"
                fpath = mt.GEN_CONFIG["PATHS"]["data"] + \
                    "lephare_input/" + inputname
                mt.assert_file_overwrite(fpath)

# %% Assert LePhare stuff


def assert_lephare_assembly():
    """Assertions needed for running lephare"""
    filt = mt.CUR_CONFIG["LEPHARE"].getboolean("run_filters")
    temps = mt.CUR_CONFIG["LEPHARE"].getboolean("run_templates")
    zphot = mt.CUR_CONFIG["LEPHARE"].getboolean("run_zphota")
    if filt | temps | zphot:
        assert isfile(mt.give_parafile_fpath())

    if filt:
        mt.assert_file_overwrite(mt.give_filterfile_fpath())

    if temps:
        if use_plike:
            assert isfile(mt.give_temp_listname("pointlike"))
            mt.assert_file_overwrite(mt.give_temp_libname("pointlike", "mag", "fits"))
        if use_ext:
            assert isfile(mt.give_temp_listname("extended"))
            mt.assert_file_overwrite(mt.give_temp_libname("extended", "mag", "fits"))
        assert isfile(mt.give_temp_listname("star"))
        mt.assert_file_overwrite(mt.give_temp_libname("star", "mag", "fits"))

    if zphot:
        assert isfile(mt.give_parafile_fpath(out=True))
        # If no input cat is produced, check if it's available:
        if not mt.CUR_CONFIG["CAT_ASSEMBLY"].getboolean("write_lephare_input"):
            if use_plike:
                assert isfile(mt.give_lepharefile_fpath("pointlike"))
            if use_ext:
                assert isfile(mt.give_lepharefile_fpath("extended"))
        # If no templates are generated, check if they are available:
        if not temps:
            if use_plike:
                assert isfile(mt.give_lepharefile_fpath("pointlike"))
                assert isfile(mt.give_temp_libname("pointlike", "mag", "fits"))
            if use_ext:
                assert isfile(mt.give_lepharefile_fpath("extended"))
                assert isfile(mt.give_temp_libname("extended", "mag", "fits"))
            assert isfile(mt.give_temp_libname("star", "mag", "fits"))
        if use_plike:
            assert mt.assert_file_overwrite(
                mt.give_lepharefile_fpath("pointlike", out=True))
        if use_ext:
            assert mt.assert_file_overwrite(
                mt.give_lepharefile_fpath("extended", out=True))


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
