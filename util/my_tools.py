# -*- coding: utf-8 -*-
"""
Created on Mon Jul 19 09:31:46 2021

@author: Fabian Balzer
"""


import logging
import os
from configparser import ConfigParser
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.table import Table

CONFIGPATH = os.environ["LEPHARE"] + "/lephare_scripts/config/"
GEN_CONFIG = ConfigParser()
GEN_CONFIG.read(CONFIGPATH + "general.ini")
CUR_CONFIG = ConfigParser()
CUR_CONFIG.read(CONFIGPATH + GEN_CONFIG["PATHS"]["current_config"])


def init_logger():
    """Initializes a logger."""
    # create logger
    logger = logging.getLogger('simple_logger')
    logger.setLevel(logging.DEBUG)
    # create console handler and set level to debug
    ch = logging.StreamHandler()
    ch.setLevel(CUR_CONFIG.getint("GENERAL", "logging_level"))
    # create formatter
    formatter = logging.Formatter('%(levelname)s: %(message)s')
    # add formatter to ch
    ch.setFormatter(formatter)
    # add ch to logger
    logger.addHandler(ch)
    return logger


LOGGER = init_logger()


def stringlist_to_list(stringlist):
    """Transform a string that is in the "['a', 'b', 'c']" pattern into a list of strings ["a", "b", "c"].
    Handy to read out lists from config files."""
    string = stringlist.strip("[]")
    singles = string.split(", ")
    return [single.strip("'") for single in singles]


def read_list_from_config(section, key):
    """Returns a list of the values listed in a config file."""
    stringlist = GEN_CONFIG.get(section, key)
    return stringlist_to_list(stringlist)


# Define the (hardcoded) path where the data sits in
CATPATH = GEN_CONFIG.get("PATHS", "cat")
# As we are adding these conversions in strings, they are stored as strings.
# Y, H, Ks Taken from Mara, J mag conversion from Blanton et al., Astronomical Journal 129, 2562 (2005), Eqs. (5) (2005AJ....129.2562B).
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

BAND_DICT = {}
for pair in GEN_CONFIG.items("BAND_DICT"):
    BAND_DICT[pair[0]] = stringlist_to_list(pair[1])


def generate_pretty_band_name(band, in_math_environ=False):
    """Generates a band name that can be used in LaTeX."""
    if band == "ZSPEC":
        return r"{\rm spec-}z" if in_math_environ else "spec-$z$"
    band = band.replace("-", "_")
    suffix = r"}" if "_" in band else ""
    new = band.replace('_', r'_{\rm ') + suffix
    return f"${new}$" if not in_math_environ else new


def give_survey_name(survey):
    """Returns the survey name for the plots"""
    name_dict = {"galex": "GALEX (GR6+7)",
                 "vhs": "VHS (DR6)", "sweep": "LS (DR9)", "hsc": "HSC (DR3)", "kids": "KiDS (DR4)", "ls10": "LS (DR10)"}
    return name_dict[survey]


def give_parafile_fpath(out=False):
    """Provides the name of the currently set LePhare parameter file.
    If out is True, the outputpara-name is used, else the inputparaname"""
    path = GEN_CONFIG['PATHS']['params']
    suffix = "out" if out else "in"
    fname = f"{CUR_CONFIG['LEPHARE']['para_stem']}_{suffix}.para"
    return path + fname


def give_filterfile_fpath():
    """Provides the name of the requested filter file"""
    return GEN_CONFIG["PATHS"]["data"] + "lephare_files/" + \
        CUR_CONFIG["LEPHARE"]["filter_stem"] + ".filt"


def give_lepharefile_fpath(ttype, out=False, suffix=None):
    """Provides the name of the requested lephare input- or output file.
    ttype: one of ["pointlike", "extended"]
    out: True or False (whether to consider in- or output)
    suffix: Custom name suffix
    """
    if out:
        stem = GEN_CONFIG["PATHS"]["data"] + \
            "lephare_input/" + CUR_CONFIG["LEPHARE"]["input_stem"]
        suffix = "in" if suffix is None else suffix
    else:
        stem = GEN_CONFIG["PATHS"]["data"] + \
            "lephare_output/" + CUR_CONFIG["LEPHARE"]["output_stem"]
        suffix = "out" if suffix is None else suffix
    return f"{stem}_{ttype}.{suffix}"


def give_temp_listname(ttype):
    """Provides the name of the list file with the templates."""
    listpath = GEN_CONFIG["PATHS"]["params"] + "template_lists/"
    fname = f"{CUR_CONFIG['LEPHARE']['template_stem']}_{ttype}.list"
    return listpath + fname


def give_temp_libname(ttype, libtype="mag", suffix=""):
    """Provides the name of the compiled template file or the name
    of the mag_lib file."""
    temppath = GEN_CONFIG["PATHS"]["data"] + "lephare_files/templates/"
    fname = f"{CUR_CONFIG['LEPHARE']['para_stem']}_{ttype}_{libtype}_lib{suffix}"
    return temppath + fname


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


def read_fits_as_dataframe(filename, saferead=False):
    """Read a given .fits file with the specified filename to return a pandas dataframe"""
    # This is the way to read the FITS data into a numpy structured array
    # (using astropy.io.fits.getdata didn't work out of the box
    # because it gives a FITSRec)
    if saferead:
        t = Table.read(filename, format="fits")
        safe_names = [name for name in t.colnames if len(
            t[name].shape) <= 1]
        data = t[safe_names].to_pandas()
    else:
        table = Table.read(filename, 1, format="fits")
        data = table.to_pandas()
    df = pd.DataFrame(data)
    return df


def save_dataframe_as_fits(df, filename, overwrite=False):
    """Store the given dataframe df as a fits file in 'filename'"""
    table = Table.from_pandas(df)
    fpath = DATAPATH + filename
    table.write(fpath, overwrite=overwrite)
    LOGGER.info(f"Successfully saved the dataframe at {fpath}.")


def read_plike_and_ext(prefix, suffix, fmt="fits"):
    """Reads a pointlike and extended fits file conforming to the prefix-ttype-suffix notation and merges them into a single pandas dataframe.
    Adds the 'Type' column set to either 'pointlike' or 'extended' and returns the df."""
    if fmt == "fits":
        df1 = read_fits_as_dataframe(f"{DATAPATH}{prefix}pointlike{suffix}")
        df2 = read_fits_as_dataframe(f"{DATAPATH}{prefix}extended{suffix}")
    else:
        df1 = read_ascii_as_dataframe(f"{DATAPATH}{prefix}pointlike{suffix}")
        df2 = read_ascii_as_dataframe(f"{DATAPATH}{prefix}extended{suffix}")
    # Add type columns for distinction
    df1["Type"] = "pointlike"
    df2["Type"] = "extended"
    df = pd.concat([df1, df2])
    return df


def read_glb_context(fname):
    """Scans the given .fits file for the GLB_CONTEXT keyword"""
    path = f"{DATAPATH}{fname}"
    with open(path, "r") as f:
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


def read_template_library(fname):
    """Reads a template library .dat file"""
    path = DATAPATH + "lephare_files/" + fname

    with open(path, "r") as f:
        coltext = f.readline()
    cols = coltext.split()[1:]
    cols = cols[:-3] + [f"mag_{band}"for band in BAND_LIST]
    cols = cols + [f"mag_err_{band}"for band in BAND_LIST]
    df = pd.read_csv(path, delim_whitespace=True,
                     header=None, skiprows=1, names=cols)
    df = df.rename(columns={"redshift": "ZSPEC"})
    return df


def construct_template_dict(temp_df, templates_to_plot):
    """Construct a dict with dataframes for each of the models requested, removing templates at E(B-V) != 0."""
    temp_dict = {}
    available_templates = set(temp_df["model"])
    templates_to_plot = available_templates if templates_to_plot is None else available_templates.intersection(
        set(templates_to_plot))
    for temp_number in templates_to_plot:
        subset = temp_df[temp_df["model"] == temp_number].sort_values(by=[
                                                                      "ZSPEC"])
        subset = subset[(subset["ext_law"] == 0) & (subset["E(B-V)"] == 0)]
        temp_dict[temp_number] = subset
    return temp_dict


def provide_template_info(fname):
    path = MYDIR + "lephare_scripts/lephare_parameters/template_lists/" + fname
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
    """Converts the ra and dec values for the VHS data to degrees as they are provided in radians."""
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
    cols = [col for col in df.columns if "MAG" in col]
    newnames = [col.replace("_", "-")
                .replace("MAG-", "mag_")
                .replace("ERR-mag-", "mag_err_") for col in cols]
    df = df.rename(columns=dict(zip(cols, newnames)))
    print(df["mag_i-ls10"].head())
    for col in newnames:
        df.loc[(df[col] <= 0) | (df[col] == 99.), col] = None
    print(df["mag_i-ls10"].head())
    return add_outlier_information(df)


def add_outlier_information(df):
    """Adds information on the outliers of the dataframe.
    Renames Z_BEST to ZBEST."""
    if "Z_BEST" in df.columns:
        df = df.rename(columns={"Z_BEST": "ZBEST"})
    df["ZMeasure"] = (df["ZBEST"] - df["ZSPEC"]) / (1 + df["ZSPEC"])
    df["ISOUTLIER"] = abs(df["ZMeasure"]) > 0.15
    df["HASGOODZ"] = ((df["ZSPEC"] > 0) & (df["ZBEST"] > 0))
    df["IsFalsePositive"] = ((df["ZSPEC"] > 0.5) & (df["ZBEST"] < 0.5))
    df["IsFalseNegative"] = ((df["ZSPEC"] < 0.5) & (df["ZBEST"] > 0.5))
    return df


def linearfunc(z, m=0):
    """Simple linear helper function with the z + m*(1+z) prescription"""
    return z + m * (1 + z)


def give_output_statistics(df, filters_used=False):
    """Takes a LePhare output DataFrame, filters it for good specz and photo-z
    rows and calculates the statistics (outlier fraction eta, accuracy
    sig_NMAD, false positive fraction and false negative fraction) for them.
    returns:
        Dictionary with 'eta', 'sig_nmad', 'psi_pos' and 'psi_neg' as keys.
    """
    df = df[df["HASGOODZ"]]
    source_num = len(df)
    outliers = df['ISOUTLIER'].sum()
    eta = outliers / source_num
    sig_nmad = 1.45 * (abs(df["ZMeasure"])).median()
    psi_pos = df['IsFalsePositive'].sum() / source_num
    psi_neg = df['IsFalseNegative'].sum() / source_num
    if filters_used:
        df_good = df[~df["ISOUTLIER"]]
        df_bad = df[df["ISOUTLIER"]]
        bads = [item for sublist in df_bad["filter_list"] for item in sublist]
        goods = [item for sublist in df_good["filter_list"]
                 for item in sublist]
        LOGGER.info("Filter\tgood photoz\tbad photoz")
        for n, filt in enumerate(BAND_LIST):
            LOGGER.info(f"{filt}:\t{goods.count(n+1)}\t{bads.count(n+1)}")
    return {"eta": eta, "sig_nmad": sig_nmad, "psi_pos": psi_pos, "psi_neg": psi_neg}


def give_photoz_performance_label(df):
    """Produces a label that can be displayed in a spec-z-phot-z plot."""
    stat_dict = give_output_statistics(df)
    etalabel = "$\eta = " + f"{stat_dict['eta']:.3f}$\n"
    sig_nmadlabel = r"$\sigma_{\rm NMAD} = " + f"{stat_dict['sig_nmad']:.3f}$"
    fpos = df['IsFalsePositive'].sum()
    fposlabel = "\n" + r"$\psi_{\rm Pos} = " + \
        f"{stat_dict['psi_pos']:.3f}$ ({fpos})"
    fneg = df['IsFalseNegative'].sum()
    fneglabel = "\n" + r"$\psi_{\rm Neg} = " + \
        f"{stat_dict['psi_neg']:.3f}$ ({fneg})"
    label = f"{len(df)} sources\n{etalabel}{sig_nmadlabel}"
    return label + fposlabel + fneglabel


def give_row_statistics(df):
    """Takes a dataframe and computes the mean values for each band."""
    for column in BAND_LIST:
        try:
            LOGGER.info(f"{column}:\t{df[column].mean()*1e28:.4g}\t\pm \
              {df[column].std()*1e28:.4g}\t 10**(-28) ergs/cm**2/Hz/s")
        except KeyError as k:
            LOGGER.info(f"Could not find the column {column} in the dataframe."
                        "\nMoving on to the next one.")


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


def give_input_statistics(df, sourcetype):
    """Takes an input dataframe and constructs input statistics, including the number of sources and the the number this would correspond to on the whole sky.
    """
    source_num = len(df)
    LOGGER.info(f"There are {source_num} sources in the {sourcetype} eFEDS dataset, \
        corresponding to {int(source_num/EFEDS_PERCENTAGE)} sources in the whole sky.")
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
            f"There are {source_num} sources with {band_num} bands of photometry available.")


def flux_to_AB(f):
    """Taxes a flux in Jy and converts it to AB magnitude."""
    return -2.5 * np.log10(f) - 48.60


def rad_to_deg(x):
    """Takes an angle in radians and converts it to degrees"""
    return 180 / np.pi * x


def find_lower_exponent(context):
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
    if type(context) is pd.Series:
        context = int(context["used_context"])
    if type(context) is not int:
        context = int(context)
        LOGGER.info(f"Forcing context to become {context}")
    if context == -1:  # Return all bands if context is -1
        return list(range(1, len(BAND_LIST) + 1))
    filter_numbers = []
    while context > 0:
        last = find_lower_exponent(context)
        # +1 as the filter numbers start with 1
        filter_numbers.append(last + 1)
        context = context - 2**last  # We now want to work with the remains.
    return filter_numbers[::-1]  # lastly, we reverse the list


def give_context(bands):
    """Returns the context belonging to the set of bands provided."""
    return sum([2**i for i, band in enumerate(BAND_LIST) if band in bands])


def give_bands_for_context(context):
    """Returns the context belonging to the set of bands provided."""
    band_indices = convert_context_to_band_indices(context)
    return [BAND_LIST[index - 1] for index in band_indices]


def give_survey_for_band(band):
    """Returns the survey that the band appeared in."""
    surveys = [key for key in BAND_DICT if band in BAND_DICT[key]]
    return surveys[0] if len(surveys) > 0 else "unknown"


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
        assert get_yes_no_input(
            f"The file '{fpath}' already exists.\nContinue to overwrite it?")
