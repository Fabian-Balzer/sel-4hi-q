# -*- coding: utf-8 -*-
"""
Created on Mon Jul 19 09:31:46 2021

@author: Fabian Balzer
"""


import os
import subprocess
from configparser import ConfigParser

import numpy as np
import pandas as pd
from astropy.table import Table

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
        return get_yes_no_input(
            f"The file '{fpath}' already exists.\nContinue to overwrite it?")
    return True


def assert_file_exists(fpath, ftype):
    """Asks the user whether to really overwrite the given file."""
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
    band_indices = convert_context_to_band_indices(context)
    return [BAND_LIST[index - 1] for index in band_indices]


def give_survey_for_band(band):
    """Returns the survey that the band appeared in."""
    surveys = [survey for survey, bands in BAND_DICT.items() if band in bands]
    return surveys[0] if len(surveys) > 0 else "unknown"


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


def give_survey_name(survey):
    """Returns the survey name for the plots"""
    name_dict = {"galex": "GALEX (GR6+7)",
                 "vhs": "VHS (DR6)", "sweep": "LS (DR9)",
                 "hsc": "HSC (DR3)", "kids": "KiDS (DR4)",
                 "ls10": "LS (DR10)"}
    return name_dict[survey]


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
    fname = filtstem + "_overview.filt" if overview else filtstem + "transmission.filt"
    return filtfilepath + fname


def give_lephare_filename(ttype, out=False, suffix=None, include_path=True):
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


def give_temp_listname(ttype, altstem=None):
    """Provides the name of the list file with the templates."""
    listpath = GEN_CONFIG["PATHS"]["params"] + "template_lists/"
    stem = CUR_CONFIG['LEPHARE']['template_stem'] if altstem is None else altstem
    fname = f"{stem}_{ttype}.list"
    return listpath + fname


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
    fname = CUR_CONFIG.get('LEPHARE', 'para_stem') + \
        "_" + ttype + "_" + libtype + "_lib" + suffix
    return temppath + fname


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
    fpath = GEN_CONFIG["PATHS"]["data"] + filename
    table.write(fpath, overwrite=overwrite)
    LOGGER.info("Successfully saved the dataframe at %s.", fpath)


def read_output_df():
    """Reads a pointlike and extended fits."""
    # Add type columns for distinction
    df_list = []
    for ttype in USED_TTYPES:
        fpath = give_lephare_filename(ttype, out=True, suffix=".fits")
        df = read_fits_as_dataframe(fpath)
        df["Type"] = ttype
        df_list.append(df)
    joined = pd.concat(df_list)
    joined = add_filter_columns(joined)
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


def read_template_library(fname):
    """Reads a template library .dat file"""
    fpath = GEN_CONFIG["PATHS"]["data"] + "lephare_files/" + fname

    with open(fpath, "r") as f:
        coltext = f.readline()
    cols = coltext.split()[1:]
    cols = cols[:-3] + [f"mag_{band}"for band in BAND_LIST]
    cols = cols + [f"mag_err_{band}"for band in BAND_LIST]
    df = pd.read_csv(fpath, delim_whitespace=True,
                     header=None, skiprows=1, names=cols)
    df = df.rename(columns={"redshift": "ZSPEC"})
    return df


def construct_template_dict(temp_df, templates_to_plot):
    """Construct a dict with dataframes for each of the models
    requested, removing templates at E(B-V) != 0."""
    temp_dict = {}
    available_templates = set(temp_df["model"])
    templates_to_plot = available_templates if templates_to_plot is None \
        else available_templates.intersection(
            set(templates_to_plot))
    for temp_number in templates_to_plot:
        subset = temp_df[temp_df["model"] == temp_number].sort_values(by=[
                                                                      "ZSPEC"])
        subset = subset[(subset["ext_law"] == 0) & (subset["E(B-V)"] == 0)]
        temp_dict[temp_number] = subset
    return temp_dict


def provide_template_info(fname):
    """Provides a dictionary linking the ID of a template to its name."""
    path = GEN_CONFIG["PATHS"]["params"] + "template_lists/" + fname
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
    return {"eta": eta, "sig_nmad": sig_nmad, "psi_pos": psi_pos, "psi_neg": psi_neg}


def give_row_statistics(df):
    """Takes a dataframe and computes the mean values for each band."""
    for column in BAND_LIST:
        try:
            LOGGER.info(f"{column}:\t{df[column].mean()*1e28:.4g}\t" + r"\pm"
                        f"{df[column].std()*1e28:.4g}\t 10**(-28) ergs/cm**2/Hz/s")
        except KeyError:
            LOGGER.info("Could not find the column %s in the dataframe."
                        "\nMoving on to the next one.", column)


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
    if context.isinstance(pd.Series):
        context = int(context["used_context"])
    if not context.isinstance(int):
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


def run_jystilts_program(filename, *args, with_path=False):
    """Runs a .py file using the java jystilts implementation,
    assuming the file is located in the 'jystilts_scripts' directory."""
    run_jystilts = f"java -jar {GEN_CONFIG['PATHS']['JYSTILTS']}"
    scriptpath = "" if with_path else f"{GEN_CONFIG['PATHS']['scripts']}jystilts_scripts/"
    match_table_string = f"{run_jystilts} '{scriptpath}{filename}' {' '.join(args)}"
    try:
        subprocess.run(match_table_string, check=True, shell=True)
    except subprocess.CalledProcessError as err:
        LOGGER.error(
            "The following error was thrown when trying to run the jystilts code:\n%s", err)


def run_lephare_command(command, arg_dict, additional=""):
    """Runs a given LePhare command in the LePhare source file."""
    main_command = f"{GEN_CONFIG['PATHS']['lepharedir']}/source/" + command
    run_string = main_command + " " + \
        " ".join([f"-{arg} {val}" for arg, val in arg_dict.items()]
                 ) + " " + additional
    LOGGER.debug("Running the following shell command:\n%s", run_string)
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


def assess_lephare_run(ttype):
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


def save_tex_file(fname, text):
    """Writes a file into the 'other' folder"""
    fpath = GEN_CONFIG["PATHS"]["other"] + "latex/" + fname
    with open(fpath, "w", encoding="utf8") as f:
        f.write(text)
        print(
            f"The LaTeX input text has been written to {fname}")
