# Jython code
import os
from sys import argv, path

WORKPATH = os.environ["LEPHARE"] + "/"
# This is needed to be able to import custom modules
path.append(WORKPATH + "sel-4hi-q/jystilts_scripts/modules")
path.append("/home/hslxrsrv3/p1wx150")
import jy_tools as jt
import stilts

try:
    TTYPE = argv[1]  # pointlike, extended, or star
    MAG_OR_OUT = argv[2]
except IndexError:
    raise IndexError(
        1, "Two arguments expected for rewriting the header: The ttype, and whether a mag lib file or an output file is expected to be rewritten.")


def rewrite_output_file():
    """Reads the LePhare output file and generates correct column names while converting it to a .fits file"""
    fpath = jt.give_lephare_filename(TTYPE, out=True)
    table = stilts.tread(fpath, fmt="ASCII")
    # cols = [str(col.name) for col in table.columns()]
    for i, band in enumerate(jt.BAND_LIST):
        table = table.cmd_colmeta(
            "-name", jt.give_nice_band_name(band), "MAG_OBS" + str(i))
        table = table.cmd_colmeta(
            "-name", jt.give_nice_band_name(band, err=True), "ERR_MAG_OBS" + str(i))
    new_fpath = jt.give_lephare_filename(TTYPE, out=True, suffix="fits")
    table.write(new_fpath, fmt="fits")


def rewrite_maglib_file():
    """Reads the .dat template file and generates correct column names while converting it to a .fits file"""
    fpath = jt.give_temp_libname(TTYPE, suffix=".dat")
    f = open(fpath, "r")
    columns = f.readline()[2:].split()
    entries = f.readline().split()
    f.close()
    vectors = [col for col in columns if "vector" in col]
    colnames = [col for col in columns if "vector" not in col]
    jt.LOGGER.debug(
        "Found the following bare columns in the %s template file: %s", TTYPE, colnames)
    table = stilts.tread(fpath, fmt="ASCII")
    # Just in case something doesn't match up:
    if len(colnames) + len(jt.BAND_LIST) * 2 != len(entries):
        jt.LOGGER.error(
            "Could not properly read the new column names of the mag_lib.dat file. The rewritten fits file is going to unmeaningful colnames.")
    else:
        for band in jt.BAND_LIST:
            colnames.append(jt.give_nice_band_name(band))
        for band in jt.BAND_LIST:
            colnames.append("kcor_" + band)
        for i, colname in enumerate(colnames):
            table = table.cmd_colmeta("-name", colname, "col" + str(i + 1))
    new_fpath = jt.give_temp_libname(TTYPE, suffix=".fits")
    table.write(new_fpath, fmt="fits")


if MAG_OR_OUT == "OUT":
    rewrite_output_file()
if MAG_OR_OUT == "MAG":
    rewrite_maglib_file()


# %% Outdated stuff:
def insert_filter_column(key, newname, colnames):
    """As some of the columns are only placeholders for the filters, check for them, insert
    a column with prefix 'newname' for each filter and return it.
    Pass if the given key is not found."""
    if not key in colnames:
        print(key + " is not present in the output para file.")
        return colnames
    i = colnames.index(key)
    return colnames[:i] + [newname + filt for filt in jt.BAND_LIST] + colnames[i + 1:]


def generate_names_from_parafile(parafile):
    """Reads the parafile to extract the correct column names"""
    colnames = []
    file = open(parafile, "r")  # 'with' context doesn't work in jython
    for line in file.readlines():
        # Empty strings are false by default
        if not line.startswith("#") and line.strip("\n"):
            colnames.append(line.strip("\n"))
    file.close()
    # Iterate over all sets where the filters are used.
    key_to_namedict = {
        "MAG_OBS()": "MAG_", "ERR_MAG_OBS()": "ERR_MAG_", "MAG_MOD()": "MOD_MAG_"}
    for key, name in key_to_namedict.items():
        colnames = insert_filter_column(key, name, colnames)
    print("Adopting the colnames from the file\n'%s'." % parafile)
    return colnames


def rewrite_colnames(filename, colnames):
    """Reads the file given in ASCII, replaces the colnames and rewrites it as a fits file."""
    table = stilts.tread(filename, fmt="ASCII")

    if len(colnames) == len(table.columns()):
        for i, colname in enumerate(colnames):
            table = table.cmd_colmeta("-name", colname, "col" + str(i + 1))
    else:
        print("Couldn't properly rename the column names. Keeping the old ones.")
        print(colnames)
        print(table.columns())
    pre, ext = os.path.splitext(filename)  # separate filename and extension
    newname = pre + ".fits"
    table.write(newname, fmt="fits")
    # os.remove(filename)
    print("Successfully rewrote the columnnames into the following fits file:")
    print(newname)


def new_way_to_rewrite(filename):
    """New way to rewrite colnames if stilts correctly reads out existing colnames"""
    table = stilts.tread(filename, fmt="ASCII")
    # cols = [str(col.name) for col in table.columns()]
    for i, band in enumerate(jt.BAND_LIST):
        table = table.cmd_colmeta("-name", "MAG_" + band, "MAG_OBS" + str(i))
        table = table.cmd_colmeta(
            "-name", "ERR_MAG_" + band, "ERR_MAG_OBS" + str(i))
    pre, ext = os.path.splitext(filename)  # separate filename and extension
    newname = pre + ".fits"
    table.write(newname, fmt="fits")


# colnames = generate_names_for_libfile(rewritefile) if parafile == "0" \
#     else generate_names_from_parafile(parafile)

# # print(rewritefile)
# # rewrite_colnames(rewritefile, colnames)
# new_way_to_rewrite(rewritefile)
