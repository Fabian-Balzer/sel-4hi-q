import os
from sys import path as p

WORKPATH = os.environ["LEPHARE"] + "/"
# This is needed to be able to import custom modules
p.append(WORKPATH + "lephare_scripts/jystilts_scripts/modules")
p.append("/home/hslxrsrv3/p1wx150")
p.append("C:/Program Files/Jystilts/stilts.jar")

import stilts

import table_io as t_io

# Define the (hardcoded) path where the data sits in
try:
    CATPATH = os.environ["CATPATH"] + "/"
except KeyError:
    CATPATH = "DUMMY"
VHS_BANDS = ["Y", "J", "H", "Ks"]
SWEEP_BANDS = ["g", "r", "z", "W1", "W2", "W3", "W4"]
GALEX_BANDS = ["FUV", "NUV"]
HSC_BANDS = ["i_hsc", "i2_hsc"]
KIDS_BANDS = ["i_kids"]
LS10_BANDS = ["i_ls10"]
BAND_DICT = {"vhs": VHS_BANDS, "sweep": SWEEP_BANDS,
             "galex": GALEX_BANDS, "hsc": HSC_BANDS, "kids": KIDS_BANDS, "ls10": LS10_BANDS}
BAND_LIST = GALEX_BANDS + SWEEP_BANDS[:3] + \
    VHS_BANDS + SWEEP_BANDS[3:] + HSC_BANDS + KIDS_BANDS + LS10_BANDS
# As we are adding these conversions in strings, they are stored as strings.
# Y, H, Ks Taken from Mara, J mag conversion from Blanton et al., Astronomical Journal 129, 2562 (2005), Eqs. (5) (2005AJ....129.2562B).
# OLD: {"Y": "0.938", "J": "0.91", "H": "1.379", "Ks": "1.85"}
FILTER_NAME_DICT = {"FUV": "filters/FUV.pb", "NUV": "filters/NUV.pb",
                    "g": "filters/newg.pb", "r": "filters/r.pb", "z": "filters/z.pb",
                    "Y": "filters/Y.lowres", "J": "filters/j.lowres", "H": "filters/h.lowres",
                    "Ks": "filters/k.lowres", "W1": "filters/W1.res", "W2": "filters/W2.res",
                    "W3": "filters/W3.res", "W4": "filters/W4.res", "i_hsc": "filters/wHSC_i.txt", "i2_hsc": "filters/wHSC_i2.txt", "i_kids": "kids/KiDSVIKING_aibn139_i.res", "i_dr10": "dr10i.pb"}
VEGA_AB_DICT = {"Y": "0.60", "J": "0.92", "H": "1.37", "Ks": "1.83"}

WORKPATH = os.environ["LEPHARE"] + "/"


def clean_ls10(table):
    """Rename some ls10 columns"""
    oldcols = ["ctp_ls8_ra", "ctp_ls8_dec", "CTP_LS8_Type"]
    newcols = ["ra", "dec", "Type"]
    for band in SWEEP_BANDS + ["i"]:
        oldcols.append("LU_flux_" + band)
        oldcols.append("LU_flux_" + band + "_err")
        newcols.append("c_flux_" + band)
        newcols.append("c_flux_err_" + band)
    table = t_io.change_colnames(table, oldcols, newcols)
    table = table.cmd_keepcols(" ".join(newcols))
    return table


agn = t_io.read_table("opt_agn")
ls10 = t_io.read_table("ls10")
vhs = t_io.read_table("vhs")
eros = t_io.read_table("eros")
agn = t_io.clean_opt_agn(agn)
ls10 = clean_ls10(ls10)
vhs = t_io.clean_vhs(vhs)
eros = t_io.clean_eros(eros)

match = t_io.match_opt_agn_sweep(agn, ls10, 0.1)
match = t_io.match_table_galex(match, 3.5)
match = t_io.match_table_vhs(match, vhs)
match = t_io.match_table_eros(match, eros)
match = t_io.correct_galex_fluxes(match)
pointlike, extended = t_io.split_by_type(match)
pointlike, extended = t_io.convert_vega_mag_to_flux(pointlike, extended)
