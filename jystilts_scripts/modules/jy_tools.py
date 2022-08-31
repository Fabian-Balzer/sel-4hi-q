"""Helper module for jython scripts"""
import logging
import os

import stilts
from ConfigParser import ConfigParser

CONFIGPATH = os.environ["LEPHARE"] + "/sel-4hi-q/config/"
GEN_CONFIG = ConfigParser()
GEN_CONFIG.read(CONFIGPATH + "general.ini")
CUR_CONFIG = ConfigParser()
CUR_CONFIG.read(CONFIGPATH + GEN_CONFIG.get("PATHS", "current_config"))


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


def give_nice_band_name(band, fluxtype="mag", err=False):
    """Helper function to unify band naming [SYNC with my_tools!]"""
    errstring = "err_" if err else ""
    return fluxtype + "_" + errstring + band


def change_colnames(table, oldnames, newnames):
    """Changes all column names of a table from oldnames (list) to newnames (list)."""
    for oname, nname in zip(oldnames, newnames):
        table = table.cmd_colmeta("-name", nname, oname)
    return table


def keep_columns(table, columnlist):
    """Returns a table reduced to the columnames given via columnlist after checking whether they are available."""
    real_columns = [column.name.lower() for column in table.columns()]
    new_columns = [column.lower() for column in columnlist]
    # perform a check for all equal items:
    new_columns_as_set = set(new_columns)
    final_list = list(new_columns_as_set.intersection(real_columns))
    return table.cmd_keepcols(" ".join(final_list))


# %% Uniform name generation
def give_match_table_name():
    """Generates a uniform table name for the matched table"""
    path = GEN_CONFIG.get("PATHS", "match")
    stem = CUR_CONFIG.get("CAT_ASSEMBLY", "cat_stem")
    return path + stem + "_raw_match.fits"


def give_processed_table_name(ttype):
    """Generates a uniform table name for the matched table"""
    path = GEN_CONFIG.get("PATHS", "match")
    stem = CUR_CONFIG.get("CAT_ASSEMBLY", "cat_stem")
    return path + stem + "_" + ttype + "_processed.fits"


def give_lephare_filename(ttype, out=False, suffix=None):
    """Generates a uniform table name for the matched table
    WARNING: Needs to be synced with jystilts!"""
    if out:
        path = GEN_CONFIG.get("PATHS", "data") + "lephare_output/"
        stem = CUR_CONFIG.get("LEPHARE", "output_stem")
        suffix = "out" if suffix is None else suffix
    else:
        path = GEN_CONFIG.get("PATHS", "data") + "lephare_input/"
        stem = CUR_CONFIG.get("LEPHARE", "input_stem")
        suffix = "in" if suffix is None else suffix
    return path + stem + "_" + ttype + "." + suffix


def give_parafile_fpath(out=False):
    """Provides the name of the currently set LePhare parameter file.
    If out is True, the outputpara-name is used, else the inputparaname"""
    path = GEN_CONFIG.get('PATHS', 'params')
    suffix = "out" if out else "in"
    fname = CUR_CONFIG.get('LEPHARE', 'para_stem') + "_" + suffix + ".para"
    return path + fname


def give_temp_libname(ttype, libtype="mag", suffix="", include_path=True):
    """Provides the name of the compiled template file or the name
    of the mag_lib file.
    WARNING: Needs to be synced with jy_tools!"""
    temppath = GEN_CONFIG.get("PATHS", "data") + \
        "lephare_files/templates/" if include_path else ""
    fname = CUR_CONFIG.get('LEPHARE', 'para_stem') + \
        "_" + ttype + "_" + libtype + "_lib" + suffix
    return temppath + fname


# %% Reading and writing files
def write_match_as_backup(table):
    """Writes a table that includes matching to the given path."""
    path = give_match_table_name()
    table.write(path, fmt="fits")
    LOGGER.info("Successfully wrote a matched table to %s", path)


def read_match_from_backup():
    """Reads a matched table and returns it and the names of the tables that have been used for the matching process."""
    path = give_match_table_name()
    table = stilts.tread(path, fmt="fits")
    return table


def write_processed_table(table, ttype):
    """Writes a collated version of the matched table with the columns already processed."""
    path = give_processed_table_name(ttype)
    table.write(path, fmt="fits")
    LOGGER.info("Successfully wrote a matched and processed table to %s", path)


def read_processed_table(ttype):
    """Reads the processed versions of the table and returns the expected table."""
    path = give_processed_table_name(ttype)
    table = stilts.tread(path, fmt="fits")
    return table


def write_lephare_input(table, ttype):
    """Writes the table of type ttype to the test stempath"""
    path = give_lephare_filename(ttype)
    table.write(path, fmt="ASCII")
    LOGGER.debug(
        "Successfully wrote a matched and processed %s LePhare input table to '%s'.", ttype, path)
    LOGGER.info(
        "The %s LePhare input table contains %d sources.", ttype, table.count_rows())


def generate_info_text(pointlike, extended, stem):
    """Generates the text for the info file."""
    text = "Info on the matched and processed tables of '" + \
        stem + "'.\n" + "-" * 40 + "\n"
    text += "Catalogue:".ljust(15)
    info_dict = {}
    cats_used = [str(column.name.lower()[4:])
                 for column in pointlike.columns() if column.name.lower().startswith("dec_")]
    cats_used = [cat for cat in cats_used if "ivar" not in cat]
    for cat in cats_used:
        total = 0
        info_dict[cat] = {}
        for table, ttype in [(pointlike, "pointlike"), (extended, "extended")]:
            num = table.cmd_select("!NULL_dec_" + cat).count_rows()
            total += num
            info_dict[cat][ttype] = num
        info_dict[cat]["total"] = total
        text += " \t| " + cat.ljust(6)[:6]
    text += " \t| spec-z\n"
    total = 0
    info_dict["spec-z"] = {}
    for table, ttype in [(pointlike, "pointlike"), (extended, "extended")]:
        num = table.cmd_select("!NULL_ZSPEC").count_rows()
        total += num
        info_dict["spec-z"][ttype] = num
    info_dict["spec-z"]["total"] = total
    info_dict["total"] = {}
    for table, ttype in [(pointlike, "pointlike"), (extended, "extended")]:
        num = table.cmd_select("!NULL_dec").count_rows()
        total += num
        info_dict["total"][ttype] = num
    info_dict["total"]["total"] = total
    for ttype in ["pointlike", "extended", "total"]:
        row = ttype.capitalize().ljust(15) + " \t| " + \
            " \t| ".join([str(info_dict[cat][ttype]).ljust(6)
                         for cat in cats_used + ["spec-z", "total"]])
        text += row + "\n"
    return text


def write_info_file(pointlike, extended):
    """Writes an info file containing info about the the processed tables"""
    path = GEN_CONFIG.get("PATHS", "match")
    stem = CUR_CONFIG.get("CAT_ASSEMBLY", "cat_stem")
    if CUR_CONFIG.getboolean("CAT_ASSEMBLY", "write_info_file"):
        text = generate_info_text(pointlike, extended, stem)
        fpath = path + stem + "_info.txt"
        f = open(fpath, "w")
        f.write(text)
        f.close()
