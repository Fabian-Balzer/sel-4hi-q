# Jython code
import os
from sys import argv, path

path.append("/home/hslxrsrv3/p1wx150/GitRepo/4hi-q_master/Skripte/Jython/modules")
path.append("/home/hslxrsrv3/p1wx150")
import stilts

WORKPATH = os.environ["LEPHARE"] + "/"

VHS_BANDS = ["Y", "J", "H", "Ks"]
SWEEP_BANDS = ["g", "r", "z", "W1", "W2", "W3", "W4"]
HSC_BANDS = ["i_hsc", "i2_hsc"]
GALEX_BANDS = ["FUV", "NUV"]
KIDS_BANDS = ["i_kids"]
FILTER_LIST = GALEX_BANDS + \
    SWEEP_BANDS[:3] + VHS_BANDS + SWEEP_BANDS[3:] + HSC_BANDS + KIDS_BANDS
try:
    parafile = argv[1]
    rewritefile = argv[2]
except IndexError:
    raise IndexError(
        1, "Please provide two arguments: a parafile (or 0 if you want to rewrite a lib file) and a filename of the file to rewrite.")


def insert_filter_column(key, newname, colnames):
    """As some of the columns are only placeholders for the filters, check for them, insert
    a column with prefix 'newname' for each filter and return it.
    Pass if the given key is not found."""
    if not key in colnames:
        print(key + " is not present in the output para file.")
        return colnames
    i = colnames.index(key)
    return colnames[:i] + [newname + filt for filt in FILTER_LIST] + colnames[i + 1:]


def generate_names_from_parafile(parafile):
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
    for key in key_to_namedict.keys():
        colnames = insert_filter_column(key, key_to_namedict[key], colnames)
    print("Adopting the colnames from the file\n'%s'." % parafile)
    return colnames


def generate_names_for_libfile(fname):
    """Reads the libfile to get the first few names"""
    f = open(fname, "r")
    columns = f.readline()[2:].split()
    entries = f.readline().split()
    f.close()
    vectors = [col for col in columns if "vector" in col]
    colnames = [col for col in columns if "vector" not in col]
    # Just in case something doesn't match up:
    has_weird_length = len(colnames) + len(FILTER_LIST) * 2 != len(entries)
    if has_weird_length:
        colnames = ["col " + str(i)
                    for i in range(len(entries) - len(FILTER_LIST))]
    for filt in FILTER_LIST:
        colnames.append("mag_" + filt)
    if not has_weird_length:
        for filt in FILTER_LIST:
            colnames.append("kcor_" + filt)
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


colnames = generate_names_for_libfile(rewritefile) if parafile == "0" \
    else generate_names_from_parafile(parafile)


rewrite_colnames(rewritefile, colnames)
