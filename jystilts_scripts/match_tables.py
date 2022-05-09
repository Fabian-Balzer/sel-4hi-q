

import os
from sys import path

from ConfigParser import ConfigParser

WORKPATH = os.environ["LEPHARE"] + "/"
# This is needed to be able to import custom modules
path.append(WORKPATH + "lephare_scripts/jystilts_scripts")

path.append("/home/hslxrsrv3/p1wx150")


from modules.Table_Container import TableContainer

tables = TableContainer()
tables.match_tables()
tables.process_match()
tables.process_for_lephare("extended", make_fits=False)
tables.process_for_lephare("pointlike", make_fits=False)


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
