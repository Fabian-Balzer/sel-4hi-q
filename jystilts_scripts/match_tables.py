

import os
from sys import path

WORKPATH = os.environ["LEPHARE"] + "/"
# This is needed to be able to import custom modules
path.append(WORKPATH + "sel-4hi-q/jystilts_scripts/modules")

import jy_tools as jt
import table_io as t_io


def load_and_clean_tables():
    """Read the tables."""
    # , "xray_agn"]  # "assef_agn"
    tables = {}
    for catname in CATS:
        table = t_io.read_table(catname)
        table = t_io.pre_clean_table(catname, table)
        tables[catname] = table
    jt.LOGGER.debug("Loaded tables from the following catalogs:\n%s", CATS)
    return tables


def match_tables(tables):
    """Matches all tables in the given table_list, or simply all tables provided."""
    if jt.CUR_CONFIG.getboolean("CAT_ASSEMBLY", "use_matched"):
        return jt.read_match_from_backup()
    else:
        _match = t_io.match_given_tables(tables)
        jt.write_match_as_backup(_match)
        return _match


def process_match(match):
    """Processes the matched table and writes pointlike and extended sub-tables."""
    tables = {}
    if jt.CUR_CONFIG.getboolean("CAT_ASSEMBLY", "use_processed"):
        for ttype in ["pointlike", "extended"]:
            tables[ttype] = jt.read_processed_table(ttype)
    else:
        with_sep = t_io.add_separation_columns(match)
        # Processing step to split the pointlike and extended part and convert all fluxes:
        tables["pointlike"], tables["extended"] = t_io.process_table(with_sep)
        for ttype in ["pointlike", "extended"]:
            # Get rid of possibly faultily-matched photometry:
            tables[ttype] = t_io.discard_problematic_matches(tables[ttype])
            jt.write_processed_table(tables[ttype], ttype)
    for ttype in ["pointlike", "extended"]:
        tables[ttype] = t_io.filter_for_testing(tables[ttype])
    jt.write_info_file(tables["pointlike"], tables["extended"])
    return tables


def process_for_lephare(tables):
    """Function to clean the tables and output them such that they can be used by LePhare."""
    for ttype in ["pointlike", "extended"]:
        if not ttype in tables.keys():
            jt.LOGGER.warning("You did not initialize the %s table.", ttype)
            return
        if jt.CUR_CONFIG.getboolean("CAT_ASSEMBLY", "write_lephare_input"):
            table = t_io.process_for_lephare(tables[ttype])
            jt.write_lephare_input(table, ttype)


if __name__ == "__main__":
    CATS = ["vhs", "sweep", "opt_agn", "eros", "hsc", "kids", "ls10"]
    STEM = jt.CUR_CONFIG.get("CAT_ASSEMBLY", "cat_stem")
    # Dictionary with the table names as keys and the tables as values
    TABLE_DICT = load_and_clean_tables()

    MATCH = match_tables(TABLE_DICT)
    PROCESSED = process_match(MATCH)
    process_for_lephare(PROCESSED)


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
