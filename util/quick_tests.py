# %%
import shlex
import subprocess

import util.my_tools as mt

df = mt.read_fits_as_dataframe(mt.DATAPATH + "sweep_dr10_corr.fits")

df["ID"] = df.index
df["CONTEXT"] = -1
df["STRING"] = df[["CTP_LS8_RA", "CTP_LS8_DEC"]].astype(
    str).agg(' '.join, axis=1)
df["ZSPEC"] = df[df["SPECZ_NORMQ"] == 3]["SPECZ_REDSHIFT"]

fluxcols = [col for col in df.columns if "LU_" in col]

# Go through the band list and assign nice photometry to all columns desired
used_cols = []
for band in mt.BAND_LIST:
    avail_cols = [col for col in fluxcols if band.lower() in col.split("_")[2]]
    if len(avail_cols) == 2:
        print(avail_cols)
        df[band] = df[avail_cols[0]]
        df[band + "_err"] = df[avail_cols[1]]
    else:
        df[band] = -99.
        df[band + "_err"] = -99.
    used_cols += [band, band + "_err"]
# As the i band is not in mt.BAND_LIST, assign it manually:
df["i-lsdr"] = df["LU_flux_i"]
df["i-lsdr_err"] = df["LU_flux_i_err"]
used_cols += ["i-lsdr", "i-lsdr_err"]


# Prepare the dataframe for LePhare use, putting the columns into the right order:
needed_cols = ["ID"] + used_cols + \
    ["CONTEXT", "ZSPEC", "STRING", "CTP_LS8_TYPE"]
df = df[needed_cols]

# Split by ttype
pointlike = df[df["CTP_LS8_TYPE"].apply(lambda x: "PSF" in str(x))]
extended = df[df["CTP_LS8_TYPE"].apply(lambda x: "PSF" not in str(x))]

mt.save_dataframe_as_fits(
    pointlike, "lephare_input/dr10_test_pointlike.fits", True)
mt.save_dataframe_as_fits(
    extended, "lephare_input/dr10_test_extended.fits", True)

call = shlex.split(
    f"java -jar /home/hslxrsrv3/p1wx150/stilts.jar tcopy {mt.DATAPATH}lephare_input/dr10_test_extended.fits {mt.DATAPATH}lephare_input/dr10_test_extended.in ifmt=fits ofmt=ascii")
subprocess.call(call)
call = shlex.split(
    f"java -jar /home/hslxrsrv3/p1wx150/stilts.jar tcopy {mt.DATAPATH}lephare_input/dr10_test_pointlike.fits {mt.DATAPATH}lephare_input/dr10_test_pointlike.in ifmt=fits ofmt=ascii")
subprocess.call(call)
