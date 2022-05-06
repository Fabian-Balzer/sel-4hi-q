import os
from sys import path as p

import stilts

p.append("C:/Program Files/Jystilts/stilts.jar")


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