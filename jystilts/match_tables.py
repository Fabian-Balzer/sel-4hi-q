

import os
from sys import path

WORKPATH = os.environ["LEPHARE"] + "/"
# This is needed to be able to import custom modules
path.append(WORKPATH + "lephare_scripts/jystilts/modules")

path.append("/home/hslxrsrv3/p1wx150")

from sys import argv

# import select_region as s_region
import stilts
from Table_Container import TableContainer

args = ["stem", "lephare_stem", "use_matched", "use_processed", "test_only",
        "write_info", "omit_lephare_input"]
arg_dict = {}
for arg in args:
    arg_dict[arg] = False

try:
    arg_dict["stem"] = argv[1]
except IndexError:
    raise(RuntimeError("Please provide at least a stem name."))

if len(argv) - 1 == len(args):
    arg_dict["lephare_stem"] = argv[2]
    for i, arg in enumerate(args[2:]):
        arg_dict[arg] = argv[i + 3].lower() == "true"
else:
    arg_dict["lephare_stem"] = arg_dict["stem"]
    print("Defaulting to the following values:")
    print(arg_dict)


tables = TableContainer(arg_dict)

if not arg_dict["use_matched"] and not arg_dict["use_processed"]:
    tables.match_tables(cat_list=["vhs", "galex", "eros", "hsc", "kids"])
if not arg_dict["use_processed"]:
    tables.process_match()

if not arg_dict["omit_lephare_input"]:
    cats = ("galex", "sweep", "vhs", "hsc", "kids")
    excluded_extended = ()  # ("W3", "W4", "i_hsc", "i2_hsc")
    tables.process_for_lephare("extended", excluded_bands=excluded_extended,
                               used_cats=cats, make_fits=False, stem_alt=None)
    excluded_pointlike = ()  # ("i_hsc", "i2_hsc")
    tables.process_for_lephare("pointlike", excluded_bands=excluded_pointlike,
                               used_cats=cats, make_fits=False, stem_alt=None)


# region_boundaries = s_region.get_region(eFEDS=True)  # Get the region as a minmaxlist
# regions = s_region.write_region_list(region_boundaries)  # Separate it into brick strings


# shu_sweep = t_io.match_shu_sweep_by_id(shu, sweep)
# shu_sweep = stilts.tskymatch2(in1=shu, in2=sweep, error=1, find="best")
# print("In the skymatched catalogue are %i more objects." %(row_diff))
# t_io.write_matched_table(shu_sweep, "shu_sweep.fits")
# shu_vhs = stilts.tskymatch2(in1=shu_sweep, in2=vhs, error=10, join="1and2", find="best")  # Later research for a proper error radius needs to be conducted


# t_io.write_matched_table(shu_vhs, "eFEDS_match.fits")

# GALEX:
# EB-V column is needed for correction of UV flux
# stilts cdsskymatch cdstable=II/335/galex_ais
# find=each in=C001_BEST_CTPS_RADEC_11321.fits ra=BEST_LS8_RA dec=BEST_LS8_Dec radius=2 out=C001_GALEX_11321.fits

# LS DR9
# Divide by
