# -*- coding: utf-8 -*-
"""
Created on Mon Jul 19 09:31:46 2021

@author: Fabian Balzer
"""


import os
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.table import Table

MYDIR = os.environ["LEPHARE"] + "/"
DATAPATH = MYDIR + "data/"
OUTPUTPATH = DATAPATH + "lephare_output/"
MATCHPATH = DATAPATH + "matches/"
PLOTPATH = MYDIR + "plots/"

VHS_BANDS = ["Y", "J", "H", "Ks"]
SWEEP_BANDS = ["g", "r", "z", "W1", "W2", "W3", "W4"]
GALEX_BANDS = ["FUV", "NUV"]
HSC_BANDS = ["i-hsc", "i2-hsc"]
KIDS_BANDS = ["i-kids"]


def generate_pretty_band_name(band, in_math_environ=False):
    """Generates a band name that can be used in LaTeX."""
    band = band.replace("-", "_")
    suffix = r"}" if "_" in band else ""
    new = band.replace('_', r'_{\rm ') + suffix
    return f"${new}$" if not in_math_environ else new


BAND_DICT = {"vhs": VHS_BANDS, "sweep": SWEEP_BANDS,
             "galex": GALEX_BANDS, "hsc": HSC_BANDS, "kids": KIDS_BANDS}
BAND_LIST = GALEX_BANDS + SWEEP_BANDS[:3] + \
    VHS_BANDS + SWEEP_BANDS[3:] + HSC_BANDS + KIDS_BANDS
# Used for LaTeX axis labels
BAND_LABEL_DICT = {band: generate_pretty_band_name(band) for band in BAND_LIST}
BAND_LABEL_DICT["ZSPEC"] = "spec-$z$"
ORDERED_BANDS = GALEX_BANDS + SWEEP_BANDS[:2] + HSC_BANDS + KIDS_BANDS\
    + SWEEP_BANDS[2:3] + VHS_BANDS + SWEEP_BANDS[3:]
VHS_WL = [1020, 1250, 1650, 2220]
SWEEP_WL = [472, 641.5, 926, 3400, 4600, 12000, 22000]
HSC_WL = [806, 806]
GALEX_WL = [150, 220]
KIDS_WL = [806]
BAND_DICT = {"galex": GALEX_BANDS, "vhs": VHS_BANDS,
             "sweep": SWEEP_BANDS, "hsc": HSC_BANDS,
             "kids": KIDS_BANDS}
SURVEY_NAME_DICT = {"galex": "GALEX (GR6+7)",
                    "vhs": "VHS (DR6)", "sweep": "LS (DR9)", "hsc": "HSC (DR3)", "kids": "KiDS (DR4)"}
WL_LIST = GALEX_WL + SWEEP_WL[:3] + \
    VHS_WL + SWEEP_WL[3:] + HSC_WL + KIDS_WL
EFEDS_PERCENTAGE = 0.0145
BAND_WL_DICT = {BAND_LIST[i]: WL_LIST[i]
                for i in range(len(BAND_LIST))}


def init_plot_directory():
    """Constructs a plot directory with the necessary subfolders if there is none already"""
    for dirs in ["output_analysis/templates", "input_analysis/separation"]:
        path = MYDIR + "plots/" + dirs
        Path(path).mkdir(parents=True, exist_ok=True)


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
    table.write(filename, overwrite=overwrite)
    print(f"Successfully saved the dataframe at {filename}.")


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
    print(context)
    print(give_bands_for_context(context))


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
    path = DATAPATH + "lephare_files/template_lists/" + fname
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
            print(f"Could not find {colname} column in dataframe.")
            if verbose:
                print(f"Available columns: {' '.join(list(df.columns))}")
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
    newnames = [col.replace("MAG", "mag").replace(
        "ERR_mag", "mag_err") for col in cols]
    df = df.rename(columns=dict(zip(cols, newnames)))
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


def give_output_statistics(df, sourcetype="", filters_used=False, verbose=True):
    """Takes a LePhare output DataFrame, filters it for good specz and photo-z
    rows and calculates the statistics (outlier fraction eta and accuracy
    sig_NMAD) for them.
    returns:
        eta, sig_NMAD
    """
    df = df[df["HASGOODZ"]]
    outliers = len(df[df["ISOUTLIER"]])
    eta = outliers / len(df)
    sig_nmad = 1.45 * (abs(df["ZMeasure"])).median()
    if verbose:
        print(f"The outlier fraction for {sourcetype} is eta = {eta:.4f}")
        print(f"The accuracy for {sourcetype} is sig_NMAD = {sig_nmad:.4f}")
    if filters_used:
        df_good = df[~df["ISOUTLIER"]]
        df_bad = df[df["ISOUTLIER"]]
        bads = [item for sublist in df_bad["filter_list"] for item in sublist]
        goods = [item for sublist in df_good["filter_list"]
                 for item in sublist]
        if verbose:
            print("Filter\tgood photoz\tbad photoz")
            for n, filt in enumerate(BAND_LIST):
                print(f"{filt}:\t{goods.count(n+1)}\t{bads.count(n+1)}")
    return eta, sig_nmad


def give_photoz_performance_label(df):
    """Produces a label that can be displayed in a spec-z-phot-z plot."""
    eta, sig_nmad = give_output_statistics(df, verbose=False)
    etalabel = "$\eta = " + f"{eta:.3f}$\n"
    sig_nmadlabel = r"$\sigma_{\rm NMAD} = " + f"{sig_nmad:.3f}$"
    length = len(df)
    fpos = df['IsFalsePositive'].sum()
    fposlabel = "\n" + r"$\psi_{\rm Pos} = " + f"{fpos/length:.3f}$ ({fpos})"
    fneg = df['IsFalseNegative'].sum()
    fneglabel = "\n" + r"$\psi_{\rm Neg} = " + f"{fneg/length:.3f}$ ({fneg})"
    label = f"{len(df)} sources\n{etalabel}{sig_nmadlabel}"
    return label + fposlabel + fneglabel


def give_row_statistics(df):
    """Takes a dataframe and computes the mean values for each band."""
    for column in BAND_LIST:
        try:
            print(f"{column}:\t{df[column].mean()*1e28:.4g}\t\pm \
              {df[column].std()*1e28:.4g}\t 10**(-28) ergs/cm**2/Hz/s")
        except KeyError as k:
            print(f"Could not find the column {column} in the dataframe."
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
    print(f"There are {source_num} sources in the {sourcetype} eFEDS dataset, \
        corresponding to {int(source_num/EFEDS_PERCENTAGE)} sources in the whole sky.")
    print(
        f"There are {len(find_good_indices(df))} sources with photometry in all bands.")
    bands = sorted(calculate_number_of_photometry(df))
    sources = [calculate_number_of_photometry(df)[band] for band in bands]
    tablestring = "_" * 80 + "\n"
    tablestring += "# bands:  " + \
        "|".join([f"{band:4}" for band in bands]) + "\n"
    tablestring += "# sources:" + \
        "|".join([f"{source:4}" for source in sources]) + "\n"
    tablestring += "_" * 80 + "\n"
    print(tablestring)
    for band_num, source_num in \
            sorted(calculate_number_of_photometry(df).items()):
        print(
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
#    print(f"For a context of {context}, the lower exponent is {i-1}.")
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
        print(f"Forcing context to become {context}")
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


# path = "../../data/"
# filename = path + "eFEDS_test_output_EXT.out"
#
# df = read_ascii_as_dataframe(filename)  # Takes some time


#
#    z = df["SPECZ_Redshift"]
#    data = plt.cm.jet(z)
#    im = axes.scatter(a, b, color=data, cmap="coolwarm")
#    fig.colorbar(im)
#    axes.grid(True)
#    axes.set_title(f"Available photometry in the efeds field for {name} sources")

# flatten filter lists
