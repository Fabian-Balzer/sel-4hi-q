import pandas as pd
import util.my_tools as mt


def read_filter_info_file():
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
    df_dict = {}
    for key, wl in filter_wls.items():
        filter_df = pd.DataFrame(
            {"Wavelength": wl, "Transmission": filter_trans[key]})
        df_dict[f"{filter_df['Wavelength'].mean():6.0f}>{key}"] = filter_df
    df = pd.concat(df_dict, axis=1)
    return df
