# -*- coding: utf-8 -*-
"""
Script to be run to set a configfile (important for datapaths)

@author: fabian_balzer
"""

# %%


import os
from configparser import ConfigParser
from pathlib import Path

from util.my_logger import LOGGER

try:
    MYDIR = os.environ["LEPHARE"] + "/"
except KeyError as e:
    LOGGER.error(
        "Please initialize the 'LEPHARE' directory before running any further code.")
    raise Exception from e

try:
    CATPATH = os.environ["CATPATH"] + "/"
except KeyError:
    LOGGER.error("No catalog path could be found.")
    CATPATH = ""
try:
    LEPHAREDIR = os.environ["LEPHAREDIR"] + "/"
    LEPHAREWORK = os.environ["LEPHAREWORK"] + "/"
except KeyError:
    LOGGER.error("No path for a LePhare installation could be found.")
    LEPHAREDIR, LEPHAREWORK = "", ""

DATAPATH = MYDIR + "data/"
OUTPUTPATH = DATAPATH + "lephare_output/"
MATCHPATH = DATAPATH + "matches/"
PLOTPATH = MYDIR + "plots/"
OTHERPATH = MYDIR + "other/"
SCRIPTPATH = MYDIR + "lephare_scripts/"
PARAPATH = MYDIR + "lephare_scripts/lephare_parameters/"
CONFIGPATH = MYDIR + "lephare_scripts/config/"
JYPATH = f"{OTHERPATH}programs/jystilts.jar"


DEFAULTCONFIG = "my_config.ini"


# %% General configuration that is shared between all runs:
config = ConfigParser()


config["CONSTS"] = {
    "perc_efeds": 0.0145,
}

config["PATHS"] = {
    "data": DATAPATH,
    "match": MATCHPATH,
    "scripts": SCRIPTPATH,
    "plot": PLOTPATH,
    "other": OTHERPATH,
    "cat": CATPATH,
    "params": PARAPATH,
    "config": CONFIGPATH,
    "lepharedir": LEPHAREDIR,
    "lepharework": LEPHAREWORK,
    "jystilts": JYPATH,
    "current_config": DEFAULTCONFIG,
}
band_dict = {
    "galex": ["FUV", "NUV"],
    "sweep": ["g", "r", "z", "W1", "W2", "W3", "W4"],
    "vhs": ["Y", "J", "H", "Ks"],
    "hsc": ["i_hsc", "i2_hsc"],
    "kids": ["i_kids"],
    "ls10": ["i_ls10"],
}

config["BAND_DICT"] = band_dict


BAND_LIST = band_dict["galex"] + band_dict["sweep"][:3] + band_dict["vhs"] + \
    band_dict["sweep"][3:] + band_dict["hsc"] + \
    band_dict["kids"] + band_dict["ls10"]
ORDERED_LIST = band_dict["galex"] + band_dict["sweep"][:2] + band_dict["hsc"] + \
    band_dict["kids"] + band_dict["ls10"] + band_dict["sweep"][2:3] + \
    band_dict["vhs"] + band_dict["sweep"][3:]
GALEX_WL = [150, 220]
SWEEP_WL = [472, 641.5, 926, 3400, 4600, 12000, 22000]
VHS_WL = [1020, 1250, 1650, 2220]
HSC_WL = [806, 806]
KIDS_WL = [806]
LS10_WL = [806]
WL_LIST = GALEX_WL + SWEEP_WL[:3] + \
    VHS_WL + SWEEP_WL[3:] + HSC_WL + KIDS_WL + LS10_WL
FILTER_LIST = ['FUV.pb', 'NUV.pb',
               'newg.pb', 'r.pb', 'z.pb',
               'Y.lowres', 'j.lowres', 'h.lowres', 'k.lowres',
               'W1.res', 'W2.res', 'W3.res', 'W4.res',
               'wHSC_i.txt', 'wHSC_i2.txt',
               'KiDSVIKING_aibn139_i.res',
               'dr10i.pb']

config["BANDS"] = {
    "listed": BAND_LIST,
    "ordered": ORDERED_LIST,
}

config["WAVELENGTHS"] = {
    band: wl for band, wl in zip(BAND_LIST, WL_LIST)}
config["FILTERS"] = {
    band: wl for band, wl in zip(BAND_LIST, FILTER_LIST)}


fpath = CONFIGPATH + 'general.ini'
with open(fpath, 'w', encoding="utf8") as configfile:
    config.write(configfile)
    LOGGER.info("A general config file has been written to %s", fpath)


def init_plot_directory(ppath):
    """Constructs a plot directory with the necessary subfolders if it is missing."""
    for dirs in ["output_analysis/templates", "input_analysis/separation"]:
        path = ppath + dirs
        Path(path).mkdir(parents=True, exist_ok=True)
    LOGGER.info("Successfully initialized the path for plots at '%s'.", ppath)


def init_data_directory(dpath):
    """Constructs a data directory with the necessary subfolders if it is missing"""
    for dirs in ["lephare_files/templates", "lephare_input", "lephare_output", "matches", "raw_catalogues"]:
        path = dpath + dirs
        Path(path).mkdir(parents=True, exist_ok=True)
    LOGGER.info("Successfully initialized the path for data at '%s'.", dpath)


def init_other_directory(opath):
    """Constructs a data directory with the necessary subfolders if it is missing"""
    for dirs in ["latex", "programs"]:
        path = opath + dirs
        Path(path).mkdir(parents=True, exist_ok=True)
    jystiltsfpath = f"{opath}programs/jystilts.jar"
    if not os.path.isfile(jystiltsfpath):
        LOGGER.warning(
            "Jystilts couldn't be located. Please install it in '%s'", jystiltsfpath)
    LOGGER.info(
        "Successfully initialized the path for other stuff at '%s'.", opath)


init_plot_directory(PLOTPATH)
init_data_directory(DATAPATH)
init_other_directory(OTHERPATH)

# %% Individual configurations changing from run to run
config = ConfigParser()

config["GENERAL"] = {
    "logging_level": 20,  # 10 would be for DEBUG
    "use_pointlike": True,
    "use_extended": True,
    "ask_overwrite": True,
}


config["CAT_ASSEMBLY"] = {
    "assemble_cat": CATPATH != "",
    "cat_stem": "baseline_input",
    "use_matched": False,
    "use_processed": False,
    "reduce_to_specz": False,
    "write_lephare_input": True,
    "write_info_file": True,
}

config["LEPHARE"] = {
    "para_stem": "baseline",
    "run_filters": False,
    "filter_stem": "baseline_filters",
    "run_templates": False,
    "extinc_range_pointlike": "0,0,0,30",
    "template_stem": "baseline_templates",
    "run_zphota": LEPHAREDIR != "",
    "forbidden_bands": ['i_hsc', 'i2_hsc', 'i_kids', 'i_ls10'],
    "input_stem": "baseline_input",
    "output_stem": "baseline_output",
    "give_stats": True,
}

config["PLOTTING"] = {
    "input": False,
    "sep": False,
    "filters": False,
    "output": True,
    "template": False,
}


fpath = CONFIGPATH + DEFAULTCONFIG
with open(fpath, 'w', encoding="utf8") as configfile:
    config.write(configfile)
    LOGGER.info("A default config file has been written to %s", fpath)
