"""Reads and asserts the config file.
Raises errors whenever a directory or file is not set up properly
(i. e. if a plot for filters is requested, but the filter file is not there."""
# %%


import os
from os.path import isfile
from sys import path

path.append(os.environ["LEPHARE"] + "/lephare_scripts")

from util.my_tools import CUR_CONFIG, GEN_CONFIG, LOGGER, assert_file_overwrite

use_plike = CUR_CONFIG["GENERAL"].getboolean("use_pointlike")
use_ext = CUR_CONFIG["GENERAL"].getboolean("use_extended")

# %% General assertions
assert isfile(
    GEN_CONFIG["PATHS"]["config"] + GEN_CONFIG["PATHS"]["current_config"])


# %% Assert catalog availability:
cat_stem = CUR_CONFIG["CAT_ASSEMBLY"]["cat_stem"]
if CUR_CONFIG["CAT_ASSEMBLY"].getboolean("assemble_cat"):
    assert GEN_CONFIG["PATHS"]["cat"] != ""
    available_cats = os.listdir(GEN_CONFIG["PATHS"]["cat"])
    for survey in GEN_CONFIG["BAND_DICT"].keys():
        if survey != "galex":
            assert survey in available_cats

    if CUR_CONFIG["CAT_ASSEMBLY"].getboolean("use_matched"):
        matchname = cat_stem + "_latest_match.fits"
        assert isfile(GEN_CONFIG["PATHS"]["match"] + matchname)

    if CUR_CONFIG["CAT_ASSEMBLY"].getboolean("use_processed"):
        if use_plike:
            procname = cat_stem + "pointlike_processed.fits"
            assert isfile(GEN_CONFIG["PATHS"]["match"] + procname)
        if use_ext:
            procname = cat_stem + "pointlike_extended.fits"
            assert isfile(GEN_CONFIG["PATHS"]["match"] + procname)

    input_stem = CUR_CONFIG["LEPHARE"]["input_stem"]
    if CUR_CONFIG["CAT_ASSEMBLY"].getboolean("write_lephare_input"):
        if use_plike:
            inputname = input_stem + "_pointlike.in"
            fpath = GEN_CONFIG["PATHS"]["data"] + "lephare_input/" + inputname
            assert_file_overwrite(fpath)
        if use_ext:
            inputname = input_stem + "_extended.in"
            fpath = GEN_CONFIG["PATHS"]["data"] + "lephare_input/" + inputname
            assert_file_overwrite(fpath)

# %% Assert LePhare stuff
parampath = GEN_CONFIG["PATHS"]["params"]
filt = CUR_CONFIG["LEPHARE"].getboolean("run_filters")
temps = CUR_CONFIG["LEPHARE"].getboolean("run_templates")
zphot = CUR_CONFIG["LEPHARE"].getboolean("run_zphota")
if filt | temps | zphot:
    assert isfile(parampath + CUR_CONFIG["LEPHARE"]["para_stem"] + "_in.para")
    assert isfile(parampath + CUR_CONFIG["LEPHARE"]["para_stem"] + "_out.para")

if filt:
    fpath = GEN_CONFIG["PATHS"]["data"] + "lephare_files/" + \
        CUR_CONFIG["LEPHARE"]["filter_stem"] + ".filt"
    assert_file_overwrite(fpath)

if temps:
    if use_plike:
        fpath = GEN_CONFIG["PATHS"]["data"] + "lephare_files/templates/" + \
            CUR_CONFIG["LEPHARE"]["para_stem"] + "_pointlike_mag_lib.fits"
        assert_file_overwrite(fpath)
    if use_ext:
        fpath = GEN_CONFIG["PATHS"]["data"] + "lephare_files/templates/" + \
            CUR_CONFIG["LEPHARE"]["para_stem"] + "_extended_mag_lib.fits"
        assert_file_overwrite(fpath)
    fpath = GEN_CONFIG["PATHS"]["data"] + "lephare_files/templates/" + \
        CUR_CONFIG["LEPHARE"]["para_stem"] + "_star_mag_lib.fits"
    assert_file_overwrite(fpath)

if zphot:
    input_stem = GEN_CONFIG["PATHS"]["data"] + \
        "lephare_input/" + CUR_CONFIG["LEPHARE"]["input_stem"]
    if not CUR_CONFIG["CAT_ASSEMBLY"].getboolean("write_lephare_input"):
        if use_plike:
            fpath = input_stem + "_pointlike.in"
            assert isfile(fpath)
        if use_ext:
            fpath = input_stem + "_extended.in"
            assert isfile(fpath)
    output_stem = GEN_CONFIG["PATHS"]["data"] + \
        "lephare_output/" + CUR_CONFIG["LEPHARE"]["output_stem"]
    if use_plike:
        fpath = output_stem + "_pointlike.out"
        assert assert_file_overwrite(fpath)
    if use_ext:
        fpath = output_stem + "_extended.out"
        assert assert_file_overwrite(fpath)


assert -1 <= CUR_CONFIG["LEPHARE"].getint("context") <= (
    len(GEN_CONFIG["BANDS"]["listed"]) + 1)**2
# %% Assert plotting:
# TODO!
# CUR_CONFIG["PLOTTING"].getboolean("input")
# CUR_CONFIG["PLOTTING"].getboolean("sep")
# CUR_CONFIG["PLOTTING"].getboolean("filters")
# CUR_CONFIG["PLOTTING"].getboolean("output")
# CUR_CONFIG["PLOTTING"].getboolean("template")

LOGGER.info("All assertions have been successful.")
