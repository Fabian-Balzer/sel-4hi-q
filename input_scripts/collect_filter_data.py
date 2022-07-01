
# %%
import pandas as pd
import util.my_tools as mt


def read_filter_transmission_file(consider_context=False):
    """Reads and returns a .filt file from LePhare and converts it to a pandas
    dataframe."""
    fpath = mt.give_filterfile_fpath(overview=False)
    filter_wls = {}
    filter_trans = {}
    with open(fpath, "r") as f:
        for line in f.readlines():
            if line.startswith("# "):
                # We can be sure that this is always called first
                name = line.strip("# ")[:-1]
            else:
                if name in filter_wls:
                    filter_trans[name] = [
                        float(val) for val in line.strip("\n").split(", ")]
                else:
                    filter_wls[name] = [float(val)
                                        for val in line.strip("\n").split(", ")]
    if consider_context:
        bands = list(filter_wls.keys())
        for band in bands:
            if band not in mt.give_bands_for_context(mt.CONTEXT):
                del filter_wls[band]
                del filter_trans[band]
    df_dict = {}
    overview = read_filter_overview_file()
    for key, wl in filter_wls.items():
        band_info = overview[overview["Band_intern"] == key.lower()]
        filter_df = pd.DataFrame(
            {"lambda": wl, "trans": filter_trans[key]})
        lambda_mean = float(band_info['Lbda_mean']) * 1000  # to nm
        survey = band_info['Survey_intern'].iloc[0]
        df_dict[f"{lambda_mean:6.0f}>{key}>{survey}"] = filter_df
    df = pd.concat(df_dict, axis=1)
    return df


def latex_with_lines(df, *args, **kwargs):
    kwargs['column_format'] = ''.join([''] + ['l'] * df.index.nlevels
                                      + ['l'] * df.shape[1] + [''])
    res = df.to_latex(*args, **kwargs)
    return res


def read_filter_overview_file():
    """Read the overview file and store the information in a dataframe. Format the column names to a LaTeX-readable format."""
    fname = mt.give_filterfile_fpath()
    with open(fname, "r") as f:
        lines = f.readlines()
    for line in lines:
        if line.startswith("# NAME"):
            skipped = lines.index(line)
            break
    colnames = [name for name in line.strip("# ").split()]
    df = pd.read_csv(fname, delim_whitespace=True,
                     skiprows=skipped, comment="#", names=colnames, skipfooter=1, engine='python')
    fname_band_dict = {fname: band for band,
                       fname in mt.GEN_CONFIG["FILTERS"].items()}
    df["Band_intern"] = df["NAME"].apply(lambda name: fname_band_dict[name])
    df["Band"] = df["Band_intern"].apply(
        lambda name: mt.give_latex_band_name(name))
    df["Survey_intern"] = df["Band_intern"].apply(
        lambda name: mt.give_survey_for_band(name))
    df["Survey"] = df["Survey_intern"].apply(
        lambda name: mt.give_latex_survey_name(name, with_dr=True))
    return df


def save_filter_info(df):
    """Saves the filter info overview in a LaTeX-readable document."""
    caption = "Overview of the important parameters of the filters."
    desired_cols = ["IDENT", "Band", "Survey",
                    "Lbda_mean", "FWHM", "AB-cor", "VEGA", "M_sun(AB)"]
    tabletext = latex_with_lines(df[desired_cols], na_rep="-", index=False, caption=caption,
                                 label="tab:filter_overview", position="htbp", float_format="${:0.2f}$".format, columns=desired_cols, escape=False, )
    lines = tabletext.splitlines()
    header = r"""ID\hyperlink{filter_a}{$^a$} & Band &   Survey &          $\Mean{\lambda}$\hyperlink{filter_b}{$^b$} [$\mu$m] &  FWHM\hyperlink{filter_c}{$^c$} [$\mu$m] &  $m_\tx{AB, c}$\hyperlink{filter_d}{$^d$} &  $m_\tx{Vega}$\hyperlink{filter_e}{$^e$} &  $M_{\odot,\tx{AB}}$\hyperlink{filter_f}{$^f$} \\"""
    start = lines.index(r"\toprule")
    lines[start + 1] = header
    footer = r"""\\
{\raggedright\small
\hypertarget{filter_a}{$^a$} The context $C_i$ of the filter with ID $i$ can be calculated using $C_i=2^{i-1}$.\\
\hypertarget{filter_b}{$^b$} $\Mean{\lambda}=\frac{\int R(\lambda)\lambda\dd \lambda}{\int R(\lambda)\dd\lambda}$ is the mean wavelength of the filter.\\
\hypertarget{filter_c}{$^c$} Full Width at Half Maximum.\\
\hypertarget{filter_d}{$^d$} AB correction factor with $m_\tx{AB}=m_\tx{Vega} + m_\tx{AB, c}$.\\
\hypertarget{filter_e}{$^e$} Vega Magnitude $m_\tx{Vega}=-2.5\log_{10}\BracedIn{\frac{\int R(\lambda)\dd \lambda}{\int R(\lambda)F_\tx{Vega}(\lambda)\dd\lambda}}$.\\
\hypertarget{filter_f}{$^f$} Absolute magnitude of the sun in this filter.\\
}"""
    lines.insert(-1, footer)
    filetext = "\n".join(lines)
    fname = mt.GEN_CONFIG["PATHS"]["other"] + "latex/" + \
        mt.CUR_CONFIG["LEPHARE"]["filter_stem"] + "_table_overview.tex"
    with open(fname, "w") as f:
        f.write(filetext)
        mt.LOGGER.info(
            "The LaTeX input text for the filters has been written to '%s'", fname)


if __name__ == "__main__":
    df = read_filter_transmission_file()
    # print(df.columns)
    # overview_df = read_filter_overview_file()
    # save_filter_info(overview_df)
