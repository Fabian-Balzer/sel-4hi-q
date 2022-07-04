# -*- coding: utf-8 -*-
"""
Created on Mon Jul 19 09:31:46 2021

@author: Fabian Balzer
"""


import os
import subprocess
from configparser import ConfigParser
from datetime import datetime

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.table import Table
from genericpath import isfile

from util.my_logger import LOGGER


def get_yes_no_input(question):
    """Tries to get user input for a yes/no question."""
    answer = input(question + "\n>>> ")
    while True:
        if answer.lower() in ["y", "yes", "yep"]:
            return True
        if answer.lower() in ["n", "no", "nope"]:
            return False
        answer = input("Please answer with 'yes' or 'no'\n>>> ")


def assert_file_overwrite(fpath):
    """Asks the user whether to really overwrite the given file."""
    if os.path.isfile(fpath):
        if CUR_CONFIG["GENERAL"].getboolean("ask_overwrite"):
            return get_yes_no_input(
                f"The file '{fpath}' already exists.\nContinue to overwrite it?")
        else:
            LOGGER.warning("Overwriting the file '%s'", fpath.split("/")[-1])
    return True


def assert_file_exists(fpath, ftype):
    """Checks whether the file at fpath exists."""
    fname = fpath.split("/")[-1]
    assert os.path.isfile(
        fpath), f"No {ftype} file with the name {fname} could be \
found, which is going to be needed in the process.\nFull path of \
the expected file:\n{fpath}"


CONFIGPATH = os.environ["LEPHARE"] + "/lephare_scripts/config/"
GEN_CONFIG = ConfigParser()
assert_file_exists(CONFIGPATH + "general.ini", "config")
GEN_CONFIG.read(CONFIGPATH + "general.ini")
CUR_CONFIG = ConfigParser()
CUR_CONFIG.read(CONFIGPATH + GEN_CONFIG["PATHS"]["current_config"])
LOGGER.setLevel(CUR_CONFIG.getint("GENERAL", "logging_level"))


def stringlist_to_list(stringlist):
    """Transform a string that is in the "['a', 'b', 'c']" pattern into a list of
    strings ["a", "b", "c"].
    Handy to read out lists from config files."""
    string = stringlist.strip("[] ")
    singles = string.split(", ")
    return [single.strip("'") for single in singles]


def read_list_from_config(section, key):
    """Returns a list of the values listed in the general config file."""
    stringlist = GEN_CONFIG.get(section, key)
    return stringlist_to_list(stringlist)


def give_context(bands, inverted=False):
    """Returns the context belonging to the set of bands provided."""
    unknown_bands = [band for band in bands if band not in BAND_LIST]
    if len(unknown_bands) > 0:
        LOGGER.warning(
            "You are trying to calculate the context of bands " +
            "that haven't been specified:\n%s", unknown_bands)
    if inverted:
        bands = [band for band in BAND_LIST if band not in bands]
        print(bands)
    return sum([2**i for i, band in enumerate(BAND_LIST) if band in bands])


def give_bands_for_context(context: int):
    """Returns the context belonging to the set of bands provided."""
    if context <= 0:
        return BAND_LIST
    band_indices = convert_context_to_band_indices(context)
    return [BAND_LIST[index - 1] for index in band_indices]


def give_survey_for_band(band: str) -> str:
    """Returns the survey that the band appeared in."""
    if band == "ZSPEC":
        return band
    surveys = [survey for survey, bands in BAND_DICT.items() if band.lower() in [
        band1.lower() for band1 in bands]]
    return surveys[0] if len(surveys) > 0 else "unknown"


def give_latex_band_name(band: str) -> str:
    """Returns the latex rep of the given band."""
    latex_map = {'fuv': "FUV", 'nuv': "NUV", 'g': "g", 'r': "r", 'z': "z", 'y': "Y",
                 'j': "J", 'h': "H", 'ks': "Ks", 'w1': "WOne", 'w2': "WTwo", 'w3': "WThree",
                      'w4': "WFour", 'i_hsc': "ihsc", 'i2_hsc': "iTwohsc", 'i_kids': "ikids", 'i_ls10': "ils"}
    return "\\" + latex_map[band.lower()] + r"band{}"


def give_latex_survey_name(survey: str, with_dr=False) -> str:
    """Returns the latex rep of the given survey. If [with_dr] it adds the DR rep."""
    latex_map = {"galex": "GALEX", "sweep": "LS", "hsc": "HSC",
                 "kids": "KIDS", "vhs": "VHS", "ls10": "LSTen"}
    suffix = r"DR{}" if with_dr else r"{}"
    return "\\" + latex_map[survey.lower()] + suffix


# Define the (hardcoded) path where the data sits in
CATPATH = GEN_CONFIG.get("PATHS", "cat")
# As we are adding these conversions in strings, they are stored as strings.
# Y, H, Ks Taken from Mara, J mag conversion from
# Blanton et al., Astronomical Journal 129, 2562 (2005), Eqs. (5) (2005AJ....129.2562B).
# OLD: {"Y": "0.938", "J": "0.91", "H": "1.379", "Ks": "1.85"}
VEGA_AB_DICT = {"Y": "0.60", "J": "0.92", "H": "1.37", "Ks": "1.83"}
GALEX_BANDS = read_list_from_config("BAND_DICT", "galex")
SWEEP_BANDS = read_list_from_config("BAND_DICT", "sweep")
VHS_BANDS = read_list_from_config("BAND_DICT", "vhs")
HSC_BANDS = read_list_from_config("BAND_DICT", "hsc")
KIDS_BANDS = read_list_from_config("BAND_DICT", "kids")
LS10_BANDS = read_list_from_config("BAND_DICT", "ls10")
BAND_LIST = read_list_from_config("BANDS", "listed")
SURVEYS = [pair[0] for pair in GEN_CONFIG.items("BAND_DICT")]
forbidden_bands = stringlist_to_list(
    CUR_CONFIG.get("LEPHARE", "forbidden_bands"))
CONTEXT = give_context(forbidden_bands, inverted=True)

BAND_DICT = {}

USED_TTYPES = []
if CUR_CONFIG["GENERAL"].getboolean("use_pointlike"):
    USED_TTYPES.append("pointlike")
if CUR_CONFIG["GENERAL"].getboolean("use_extended"):
    USED_TTYPES.append("extended")
USED_TTYPES = tuple(USED_TTYPES)  # Just to ensure they're not altered

for pair in GEN_CONFIG.items("BAND_DICT"):
    BAND_DICT[pair[0]] = stringlist_to_list(pair[1])


def give_nice_band_name(band, fluxtype="mag", err=False):
    """Helper function to unify band naming [SYNC with jy_tools!]"""
    errstring = "err_" if err else ""
    return fluxtype + "_" + errstring + band


def generate_pretty_band_name(band, in_math_environ=False):
    """Generates a band name that can be used in LaTeX."""
    if band == "ZSPEC":
        return r"{\rm spec-}z" if in_math_environ else "spec-$z$"
    band = band.replace("-", "_")
    suffix = r"}" if "_" in band else ""
    new = band.replace('_', r'_{\rm ') + suffix
    return f"${new}$" if not in_math_environ else new


def give_match_table_name():
    """Generates a uniform table name for the matched table.
    WARNING: Needs to be synced with jy_tools!"""
    path = GEN_CONFIG.get("PATHS", "match")
    stem = CUR_CONFIG.get("CAT_ASSEMBLY", "cat_stem")
    return path + stem + "_raw_match.fits"


def give_processed_table_name(ttype):
    """Generates a uniform table name for the matched table
    WARNING: Needs to be synced with jy_tools!"""
    path = GEN_CONFIG.get("PATHS", "match")
    stem = CUR_CONFIG.get("CAT_ASSEMBLY", "cat_stem")
    return path + stem + "_" + ttype + "_processed.fits"


def give_survey_name(survey_name: str):
    """Returns the survey name for the plots"""
    name_dict = {"galex": "GALEX (GR6+7)",
                 "vhs": "VHS (DR6)", "sweep": "LS (DR9)",
                 "hsc": "HSC (DR3)", "kids": "KiDS (DR4)",
                 "ls10": "LS (DR10)", "eros": "eFEDS ctps"}
    try:
        name = name_dict[survey_name]
    except KeyError:
        name = survey_name
        LOGGER.warning("Unidentified survey '%s' found.", name)
    return name


def give_survey_color(survey_name: str) -> str:
    """Returns a uniform color for each survey"""
    color_dict = {"galex": "blue", "sweep": "green",
                  "vhs": "red", "hsc": "lightgreen", "kids": "darkgreen",
                  "ls10": "darkolivegreen", "ZSPEC": "black"}
    try:
        color = color_dict[survey_name]
    except KeyError:
        color = "darkred"
        LOGGER.warning("Unidentified survey '%s' found.", color)
    return color


def give_survey_color_for_band(band: str) -> str:
    """Takes a photometric band, searches for the corresponding survey and returns the color associated with it."""
    survey = give_survey_for_band(band)
    return give_survey_color(survey)


def give_parafile_fpath(out=False):
    """Provides the name of the currently set LePhare parameter file.
    If out is True, the outputpara-name is used, else the inputparaname"""
    path = GEN_CONFIG.get('PATHS', 'params')
    suffix = "out" if out else "in"
    fname = CUR_CONFIG.get('LEPHARE', 'para_stem') + "_" + suffix + ".para"
    return path + fname


def give_filterfile_fpath(overview=True):
    """Provides the name of the requested filter file"""
    filtfilepath = GEN_CONFIG['PATHS']['params']
    filtstem = CUR_CONFIG["LEPHARE"]["filter_stem"]
    fname = filtstem + "_overview.filt" if overview else filtstem + "_transmissions.filt"
    return filtfilepath + fname


def give_lephare_filename(ttype, out=False, suffix: str = None, include_path=True) -> str:
    """Generates a uniform table name for the LePhare table
    WARNING: Needs to be synced with jy_tools!"""
    if out:
        path = GEN_CONFIG.get("PATHS", "data") + "lephare_output/"
        stem = CUR_CONFIG.get("LEPHARE", "output_stem")
        suffix = ".out" if suffix is None else suffix
    else:
        path = GEN_CONFIG.get("PATHS", "data") + "lephare_input/"
        stem = CUR_CONFIG.get("LEPHARE", "input_stem")
        suffix = ".in" if suffix is None else suffix
    path = path if include_path else ""
    return path + stem + "_" + ttype + suffix


def give_temp_listname(ttype: str, altstem: str = None, include_path=True):
    """Provides the name of the list file with the templates."""
    listpath = GEN_CONFIG["PATHS"]["params"] + "template_lists/"
    stem = CUR_CONFIG['LEPHARE']['template_stem'] if altstem is None else altstem
    fname = f"{stem}_{ttype}.list"
    return listpath + fname if include_path else fname


def give_statsfile_fname():
    """Provides the name of the stats file used to provide concise information on the LePhare run."""
    listpath = GEN_CONFIG["PATHS"]["params"]
    return listpath + "lephare_run_stats.txt"


def give_list_of_tempnames(ttype, altstem: str = None):
    """Reads the currently selected list of templates and returns a list of the active ones."""
    fname = give_temp_listname(ttype, altstem=altstem)
    with open(fname, "r") as f:
        templates = [line for line in f.readlines()
                     if not line.startswith("#")]
    templates = [temp.replace("\n", "").split()[0] for temp in templates]
    return templates


def give_list_of_tempnumbers(ttype, tempnames, altstem: str = None):
    """Reads the currently selected list of templates and returns a list of the active ones."""
    fname = give_temp_listname(ttype, altstem=altstem)
    with open(fname, "r") as f:
        templates = [line for line in f.readlines()
                     if not line.startswith("#")]
    templates = [temp.replace("\n", "").split()[0] for temp in templates]
    return templates


def get_temp_num_for_name(ttype: str, temp_name: str):
    """Scans the template lists for the template in question and returns its index.
    params:
        ttype: str
            'extended', 'pointlike'  Will check the requested list files in question
        temp_name: str
            name of the template
    returns:
        index: int || None
            The index of the template in question, or None if the template is not available
    """
    temps = give_list_of_tempnames(ttype)
    try:
        return temps.index(temp_name) + 1  # Template numbering starts at 1
    except ValueError:
        LOGGER.warning(
            "Couldn't find %s in the given %s template list.", temp_name, ttype)


def get_temp_name_for_num(ttype: str, temp_num: str):
    """Scans the template lists for the template in question and returns its index.
    params:
        ttype: str
            'extended', 'pointlike'  Will check the requested list files in question
        temp_name: str
            name of the template
    returns:
        index: int || None
            The index of the template in question, or None if the template is not available
    """
    try:
        # Template numbering starts at 1
        return give_list_of_tempnames(ttype)[temp_num - 1]
    except IndexError:
        LOGGER.warning(
            "There are less than %d templates in the given %s template list.", temp_num, ttype)


def give_temp_libname(ttype, libtype="mag", suffix="", include_path=True, use_workpath=False):
    """Provides the name of the compiled template file or the name
    of the mag_lib file.
    WARNING: Needs to be synced with jy_tools!"""
    temppath = GEN_CONFIG.get("PATHS", "data") + \
        "lephare_files/templates/"
    if use_workpath:
        temppath = GEN_CONFIG.get(
            "PATHS", "lepharework") + "lib_" + libtype + "/"
    temppath = temppath if include_path else ""
    fname = CUR_CONFIG.get('LEPHARE', 'template_stem') + \
        "_" + ttype + "_" + libtype + "_lib" + suffix
    return temppath + fname


def read_fits_as_dataframe(fname, saferead=False):
    """Read a given .fits file with the specified filename to return a pandas dataframe"""
    # This is the way to read the FITS data into a numpy structured array
    # (using astropy.io.fits.getdata didn't work out of the box
    # because it gives a FITSRec)
    if saferead:
        t = Table.read(fname, format="fits")
        safe_names = [name for name in t.colnames if len(
            t[name].shape) <= 1]
        data = t[safe_names].to_pandas()
    else:
        table = Table.read(fname, 1, format="fits")
        data = table.to_pandas()
    df = pd.DataFrame(data)
    return df


def read_ascii_as_df(fname, out=True):
    """Read a .out LePhare ASCII output file and return it as a dataframe."""
    with open(fname, "r", encoding="utf-8") as f:
        lines = f.readlines()
        index = lines.index("# Format topcat: \n")
        header = lines[index + 1]
    columns = header.split()[1:header.split().index("MAG_OBS0")]
    if out:
        additional1 = [give_nice_band_name(band) for band in BAND_LIST]
        additional2 = [give_nice_band_name(
            band, err=True) for band in BAND_LIST]
        columns = columns + additional1 + additional2
    else:
        raise NotImplementedError
    df = pd.read_csv(fname, comment="#", names=columns, delim_whitespace=True)
    return df


def save_dataframe_as_fits(df, filename, overwrite=False):
    """Store the given dataframe df as a fits file in 'filename'"""
    table = Table.from_pandas(df)
    fpath = GEN_CONFIG["PATHS"]["data"] + filename
    table.write(fpath, overwrite=overwrite)
    LOGGER.info("Successfully saved the dataframe at %s.", fpath)


def read_saved_df(cat_type="out"):
    """Reads the pointlike and extended fits input file."""
    df_list = []
    is_out = (cat_type == "out")
    for ttype in USED_TTYPES:
        if cat_type in ["in", "out"]:
            fpath = give_lephare_filename(ttype, out=is_out, suffix=".fits")
        if cat_type == "in_processed":
            fpath = give_processed_table_name(ttype)
        df = read_fits_as_dataframe(fpath)
        df["Type"] = ttype
        df_list.append(df)
    joined = pd.concat(df_list)
    joined = add_filter_columns(joined) if is_out else add_mag_columns(joined)
    return joined


def read_glb_context(fname):
    """Scans the given .fits file for the GLB_CONTEXT keyword"""
    fpath = GEN_CONFIG["PATHS"]["data"] + fname
    with open(fpath, "r") as f:
        while True:
            line = f.readline()
            if "GLB_CONTEXT" in line:
                break
    context = int(line.split(":")[1].strip())
    LOGGER.info(context)
    LOGGER.info(give_bands_for_context(context))


def read_ascii_as_dataframe(filename):
    """Read a .out ASCII-file as a dataframe"""
    table = Table.read(filename, format="ascii")
    data = table.to_pandas()
    df = pd.DataFrame(data)
    return df


def read_template_library():
    """Reads the template library .dat files"""
    df_list = []
    for ttype in USED_TTYPES:
        fpath = give_temp_libname(ttype, suffix=".dat")
        # The first line includes info about the columns
        with open(fpath, "r") as f:
            coltext = f.readline()
        cols = coltext.split()[1:]  # remove the hashtag
        cols = cols[:-3] + [f"mag_{band}"for band in BAND_LIST]
        cols = cols + [f"mag_err_{band}"for band in BAND_LIST]
        df = pd.read_csv(fpath, delim_whitespace=True,
                         header=None, skiprows=1, names=cols)
        df = df.rename(columns={"redshift": "ZSPEC"})
        df["Type"] = ttype
        df_list.append(df)
    joined = pd.concat(df_list)
    return joined


def construct_template_dict(temp_df, templates_to_plot=None):
    """Construct a dict with dataframes for each of the models
    requested, removing templates at E(B-V) != 0."""
    temp_dict = {}
    available_templates = set(temp_df["model"])
    templates_to_plot = available_templates if templates_to_plot is None \
        else available_templates.intersection(
            set(templates_to_plot))
    for temp_num in templates_to_plot:
        subset = temp_df[temp_df["model"] == temp_num]
        # Since there are usually different numbers of templates, iterate over all of them:
        for ebv in set(subset["E(B-V)"]):
            ebv_string = "" if ebv == 0 else f"; E(B-V)={ebv}"
            label = f"{temp_num}{ebv_string}"
            subsubset = subset[subset["E(B-V)"] == ebv]
            temp_dict[label] = subsubset.iloc[0:302]
            # .sort_values(by=["ZSPEC"])
            # subset = subset[(subset["ext_law"] == 0) & (subset["E(B-V)"] == 0)]
    return temp_dict


def provide_template_info(ttype):
    """Provides a dictionary linking the ID of a template to its name."""
    path = give_temp_listname(ttype)
    with open(path, "r") as f:
        number_to_name = {}
        for i, line in enumerate(f.readlines(), start=1):
            if len(line) > 4:
                if "/" in line:
                    line = line.split("/")[-1]
                number_to_name[i] = line.split(".")[0]
    return number_to_name


def add_mag_columns(df, verbose=False):
    """Adds magnitude columns to a dataframe with given fluxes for the
    columnlist"""
    for band in BAND_LIST:
        colname = f"c_flux_{band.replace('-', '_')}"
        errcolname = f"c_flux_err_{band.replace('-', '_')}"
        try:
            df.loc[df[colname] <= 0, colname] = None
            df.loc[df[errcolname] <= 0, errcolname] = None
            df["mag_" + band] = flux_to_AB(df[colname])
            df["mag_err_" + band] = flux_to_AB(df[errcolname])
        except KeyError:
            LOGGER.info(f"Could not find {colname} column in dataframe.")
            if verbose:
                LOGGER.info(f"Available columns: {' '.join(list(df.columns))}")
    return df


def convert_cols_to_deg(df):
    """Converts the ra and dec values for the VHS data to degrees as
    they are provided in radians."""
    for type_ in ["vhs"]:
        for axis in ["ra_", "dec_"]:
            df[axis + type_] = 180 / np.pi * df[axis + type_]
    return df


def prepare_dataframe(df, is_input=True):
    """Takes the dataframe and adds important columns."""
    if is_input:
        df = add_mag_columns(df)
    else:
        df = add_filter_columns(df)
    return df


def calculate_used_context(series, cols):
    """Calculates the context of bands actually used.
    TODO: Consider forbidden context"""
    context = 0
    for i, col in enumerate(cols):
        if series[col] > 0 and series[col] < 50:
            context += 2**i
    return context


def add_filter_columns(df):
    """Adds columns with lists of the used filters and the count of these
    filters.
    """
    mycols = [col for col in df.columns if "MAG" in col and "ERR" not in col]
    df["used_context"] = df[mycols].apply(
        lambda x: calculate_used_context(x, mycols), axis=1)
    df["filter_list"] = df.apply(convert_context_to_band_indices, axis=1)
    # this way, for each list the length is taken
    df["nfilters"] = df["filter_list"].str.len()
    # rename capitalized mag columns
    cols = [col for col in df.columns if "mag" in col]
    newnames = [col.replace("_", "-")
                .replace("mag-", "mag_")
                .replace("err-mag-", "mag_err_") for col in cols]
    df = df.rename(columns=dict(zip(cols, newnames)))
    cols = [col for col in df.columns if "mag" in col]
    for col in newnames:
        df.loc[(df[col] <= 0) | (df[col] > 98.0), col] = None
    return add_outlier_information(df)


def add_outlier_information(df):
    """Adds information on the outliers of the dataframe.
    Renames Z_BEST to ZBEST."""
    if "Z_BEST" in df.columns:
        df = df.rename(columns={"Z_BEST": "ZBEST"})
    df["ZMeasure"] = (df["ZBEST"] - df["ZSPEC"]) / (1 + df["ZSPEC"])
    df["IsOutlier"] = abs(df["ZMeasure"]) > 0.15
    df["HasGoodz"] = ((df["ZSPEC"] > 0) & (df["ZBEST"] > 0))
    df["IsFalsePositive"] = ((df["ZSPEC"] < 0.5) & (df["ZBEST"] > 0.5))
    df["IsFalseNegative"] = ((df["ZSPEC"] > 0.5) & (df["ZBEST"] < 0.5))
    df["TemplateScore"] = np.exp(-df["ZMeasure"] - df["CHI_BEST"] / 2)
    return df


def linearfunc(z, m=0):
    """Simple linear helper function with the z + m*(1+z) prescription"""
    return z + m * (1 + z)


def give_output_statistics(df, filters_used=False) -> dict:
    """Takes a LePhare output DataFrame, filters it for good specz and photo-z
    rows and calculates the statistics (outlier fraction eta, accuracy
    sig_NMAD, false positive fraction and false negative fraction) for them.
    returns:
        Dictionary with 'eta', 'sig_nmad', 'psi_pos' and 'psi_neg' as keys.
    """
    df = df[df["HasGoodz"]]
    source_num = len(df)
    outliers = df['IsOutlier'].sum()
    eta = outliers / source_num
    sig_nmad = 1.45 * (abs(df["ZMeasure"])).median()
    psi_pos = df['IsFalsePositive'].sum() / source_num
    psi_neg = df['IsFalseNegative'].sum() / source_num
    if filters_used:
        df_good = df[~df["IsOutlier"]]
        df_bad = df[df["IsOutlier"]]
        bads = [item for sublist in df_bad["filter_list"] for item in sublist]
        goods = [item for sublist in df_good["filter_list"]
                 for item in sublist]
        LOGGER.info("Filter\tgood photoz\tbad photoz")
        for n, filt in enumerate(BAND_LIST):
            LOGGER.info(f"{filt}:\t{goods.count(n+1)}\t{bads.count(n+1)}")
    stat_dict = {"eta": eta, "sig_nmad": sig_nmad,
                 "psi_pos": psi_pos, "psi_neg": psi_neg}
    return stat_dict


def give_row_statistics(df):
    """Takes a dataframe and computes the mean values for each band."""
    for column in BAND_LIST:
        try:
            LOGGER.info(f"{column}:\t{df[column].mean()*1e28:.4g}\t" + r"\pm"
                        f"{df[column].std()*1e28:.4g}\t 10**(-28) ergs/cm**2/Hz/s")
        except KeyError:
            LOGGER.info("Could not find the column %s in the dataframe."
                        "\nMoving on to the next one.", column)


def give_plot_title(ttype, with_info=False):
    """Provide a nice generic plot title that can optionally display the current context"""
    temp = CUR_CONFIG['LEPHARE']['template_stem']
    infostring = f" [$C={CONTEXT}$, {temp}]" if with_info else ""
    return f"{ttype.capitalize()} sources{infostring}"


def find_good_indices(df):
    """Returns the indices of the dataframe where photometry for all bands is
    available."""
    good_photo_indices = []
    total_count = len(df)
    for index in range(total_count):
        is_good = True
        for band in BAND_LIST:
            is_good = is_good and df[band][index] > 0
        if is_good:
            good_photo_indices.append(index)
    return good_photo_indices


def calculate_number_of_photometry(df):
    """Returns a dictionary with the number of photometric bands as keys and the
    number of sources that have this number available as values."""
    number_dict = {}
    for _, source in df.iterrows():
        i = 0
        for band in BAND_LIST:
            if source[band] > 0:
                i += 1
        if i in number_dict:
            number_dict[i] += 1
        else:
            number_dict[i] = 1
    return number_dict


def give_amount_of_good_photometry(df, band):
    """Returns the number of sources with good photometry for the given band in the dataframe provided."""
    return len(df[df[f"mag_{band}"] > 0])


def give_input_statistics(df, sourcetype):
    """Takes an input dataframe and constructs input statistics,
    including the number of sources and the the number this would
    correspond to on the whole sky.
    """
    source_num = len(df)
    total_source_num = source_num / GEN_CONFIG["CONSTS"].getfloat("perc_efeds")
    LOGGER.info("There are %d sources in the %s eFEDS dataset, corresponding \
to %d sources in the whole sky.", source_num, sourcetype, total_source_num)
    LOGGER.info(
        f"There are {len(find_good_indices(df))} sources with photometry in all bands.")
    bands = sorted(calculate_number_of_photometry(df))
    sources = [calculate_number_of_photometry(df)[band] for band in bands]
    tablestring = "_" * 80 + "\n"
    tablestring += "# bands:  " + \
        "|".join([f"{band:4}" for band in bands]) + "\n"
    tablestring += "# sources:" + \
        "|".join([f"{source:4}" for source in sources]) + "\n"
    tablestring += "_" * 80 + "\n"
    LOGGER.info(tablestring)
    for band_num, source_num in \
            sorted(calculate_number_of_photometry(df).items()):
        LOGGER.info(
            f"There are {source_num} sources with {band_num} \
bands of photometry available.")


def flux_to_AB(f):
    """Taxes a flux in Jy and converts it to AB magnitude."""
    return -2.5 * np.log10(f) - 48.60


def rad_to_deg(x):
    """Takes an angle in radians and converts it to degrees"""
    return 180 / np.pi * x


def find_lower_exponent(context):
    """Searches for the next lowest exponent of 2 and returns it"""
    i = 0
    while True:
        if 2**i > context:
            break
        i += 1
#    LOGGER.info(f"For a context of {context}, the lower exponent is {i-1}.")
    return i - 1


def convert_context_to_band_indices(context):
    """Decodes a given context and returns the used filter numbers.
    Works for a given pandas Series with a context column and also for a given
    number.
    Example:
        if context = 13, it will return [1, 3, 4], since 2^(1-1) + 2^2 + 2^3 = 13.
    """
    if isinstance(context, pd.Series):
        context = int(context["used_context"])
    if not isinstance(context, int):
        context = int(context)
        LOGGER.info(f"Forcing context to become {context}")
    if context <= 0:  # Return all bands if context is -1
        return list(range(1, len(BAND_LIST) + 1))
    filter_numbers = []
    while context > 0:
        last = find_lower_exponent(context)
        # +1 as the filter numbers start with 1
        filter_numbers.append(last + 1)
        context = context - 2**last  # We now want to work with the remains.
    return filter_numbers[::-1]  # lastly, we reverse the list


def run_jystilts_program(filename, *args, with_path=False):
    """Runs a .py file using the java jystilts implementation,
    assuming the file is located in the 'jystilts_scripts' directory."""
    run_jystilts = f"java -jar {GEN_CONFIG['PATHS']['JYSTILTS']}"
    scriptpath = "" if with_path else f"{GEN_CONFIG['PATHS']['scripts']}jystilts_scripts/"
    match_table_string = f"{run_jystilts} '{scriptpath}{filename}' {' '.join(args)}"
    if CUR_CONFIG["GENERAL"].getboolean("print_commands_only"):
        print(match_table_string)
        return
    try:
        subprocess.run(match_table_string, check=True, shell=True)
    except subprocess.CalledProcessError as err:
        LOGGER.error(
            "The following error was thrown when trying to run the jystilts code:\n%s", err)


def run_lephare_command(command, arg_dict, ttype, additional=""):
    """Runs a given LePhare command in the LePhare source file."""
    main_command = f"{GEN_CONFIG['PATHS']['lepharedir']}/source/" + command
    run_string = main_command + " " + \
        " ".join([f"-{arg} {val}" for arg, val in arg_dict.items()]
                 ) + " " + additional
    if CUR_CONFIG["GENERAL"].getboolean("print_commands_only"):
        print(run_string)
        return
    LOGGER.debug("Running the following shell command:\n%s", run_string)
    LOGGER.info("Running %s for %s. This could take a while...", command, ttype)
    try:
        subprocess.run(run_string, check=True, shell=True)
    except subprocess.CalledProcessError as err:
        LOGGER.error(
            "The following error was thrown when running the last shell command:\n%s", err)


def give_photoz_performance_label(df):
    """Produces a label that can be displayed in a spec-z-phot-z plot.
    Returns a text"""
    stat_dict = give_output_statistics(df)
    etalabel = r"$\eta_{\rm out} = " + f"{stat_dict['eta']:.3f}$\n"
    sig_nmadlabel = r"$\sigma_{\rm NMAD} = " + f"{stat_dict['sig_nmad']:.3f}$"
    fpos = df['IsFalsePositive'].sum()
    fposlabel = "\n" + r"$\psi_{\rm Pos} = " + \
        f"{stat_dict['psi_pos']:.3f}$ ({fpos})"
    fneg = df['IsFalseNegative'].sum()
    fneglabel = "\n" + r"$\psi_{\rm Neg} = " + \
        f"{stat_dict['psi_neg']:.3f}$ ({fneg})"
    label = f"{len(df)} sources\n{etalabel}{sig_nmadlabel}"
    return label + fposlabel + fneglabel


def log_run_info():
    """Function to log the input parameters"""
    LOGGER.info("Program started with the following requests:")
    cat_config = CUR_CONFIG["CAT_ASSEMBLY"]
    if cat_config.getboolean("assemble_cat"):
        LOGGER.info("Catalogue assembly with '%s' as a stem:",
                    cat_config["cat_stem"])
        for boolkey in ["use_matched", "use_processed", "reduce_to_specz",
                        "write_lephare_input", "write_info_file"]:
            val = cat_config.getboolean(boolkey)
            LOGGER.info("%s:\t%s", boolkey, str(val))
    lep_config = CUR_CONFIG["LEPHARE"]
    if lep_config.getboolean("run_filters"):
        LOGGER.info("LePhare filter run with '%s' as a stem.",
                    lep_config["filter_stem"])
    if lep_config.getboolean("run_templates"):
        LOGGER.info("LePhare template run with '%s' as a stem.",
                    lep_config["template_stem"])
    if lep_config.getboolean("run_zphota"):
        LOGGER.info("LePhare zphota run with '%s' as input and '%s' as output stem.",
                    lep_config["input_stem"], lep_config["output_stem"])
        LOGGER.info("The provided global context is %s, corresponding to the following bands:\n%s",
                    CONTEXT, give_bands_for_context(CONTEXT))
        if lep_config.getboolean("give_stats"):
            # TODO Stats file?
            LOGGER.info(
                "Statistics about the LePhare run are going to be provided.")


def assess_lephare_run(ttype, write=False):
    """Directly assess the quality of a photo-z run."""
    df = read_fits_as_dataframe(
        give_lephare_filename(ttype, out=True, suffix=".fits"))
    df = add_filter_columns(df)
    stat_dict = give_output_statistics(df)
    fpos = df['IsFalsePositive'].sum()
    fneg = df['IsFalseNegative'].sum()
    LOGGER.info("The outlier fraction for %s is eta = %.4f.",
                ttype, stat_dict['eta'])
    LOGGER.info("The accuracy for %s is sig_NMAD = %.4f",
                ttype, stat_dict['sig_nmad'])
    LOGGER.info(
        "The false pos fraction is psi_pos = %.4f (%d)", stat_dict['psi_pos'], fpos)
    LOGGER.info(
        "The false neg fraction is psi_neg = %.4f (%d)", stat_dict['psi_neg'], fneg)
    if not write:
        return
    # Read the template file and store the information in
    stat_dict["tempfile"] = "\\code{" + \
        give_temp_listname(ttype, include_path=False) + "}"
    with open(give_temp_listname(ttype, include_path=True), "r", encoding="utf-8") as f:
        lines = [line for line in f.readlines(
        ) if not line.startswith("#") and len(line) > 5]
    stat_dict["# templates"] = str(len(lines))
    stat_dict["time"] = datetime.now().strftime("%y-%m-%d %H:%M")
    stat_dict["context"] = str(CONTEXT)
    order = ["tempfile", "# templates", "eta", "sig_nmad",
             "psi_pos", "psi_neg", "context", "time"]
    text = " & ".join([f"${stat_dict[key]:.3f}$" if not isinstance(
        stat_dict[key], str) else stat_dict[key] for key in order]) + " \\\\\n"
    fname = give_statsfile_fname()
    if not isfile(fname):
        # Write a header if the stats file is not yet available:
        text = " & ".join(order) + " \\\\\n" + text
    with open(fname, "a", encoding="utf-8") as f:
        f.write(text)
    LOGGER.info("Wrote the results to the stats file at '%s'", fname)


def save_tex_file(fname, text):
    """Writes a file into the 'other' folder"""
    fpath = GEN_CONFIG["PATHS"]["other"] + "latex/" + fname
    with open(fpath, "w", encoding="utf8") as f:
        f.write(text)
        print(
            f"The LaTeX input text has been written to {fname}")
