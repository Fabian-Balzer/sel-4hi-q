# %%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import util.my_tools as mt


class OutputPlotContainer:
    def __init__(self):
        df = mt.read_output_df()
        df_p = df[df["Type"] == "pointlike"]
        df_e = df[df["Type"] == "extended"]
        df_outliers = df[df["IsOutlier"]]
        df_p_outliers = df_p[df_p["IsOutlier"]]
        df_e_outliers = df_e[df_e["IsOutlier"]]
        df_good = df[~df["IsOutlier"]]
        df_p_good = df_p[~df_p["IsOutlier"]]
        df_e_good = df_e[~df_e["IsOutlier"]]


if __name__ == "__main__":
    o_p_c = OutputPlotContainer()
