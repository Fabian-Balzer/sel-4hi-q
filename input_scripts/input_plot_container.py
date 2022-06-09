# %%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import util.configure_matplotlib as cm
import util.my_tools as mt

import input_scripts.filter_coverage_plot as f_c
import input_scripts.magnitude_dist_plots as m_d
import input_scripts.separation_plots as sep_p


class InputPlotContainer:
    """Convenience class for containing and plotting the input data."""

    def __init__(self, save_figures=False):
        self.df = mt.read_saved_df(cat_type="in_processed")
        self.ext = "extended" in mt.USED_TTYPES
        self.plike = "pointlike" in mt.USED_TTYPES
        self.produce_both = self.ext & self.plike
        if not (self.plike | self.ext):
            mt.LOGGER.error(
                "Please specify a valid ttype [pointlike, extended or both] to plot the specz-photoz-plot.")
        self.stem = mt.CUR_CONFIG["LEPHARE"]["input_stem"]
        self.save_figures = save_figures

    def _quicksave_fig(self, fig, name, dir_=""):
        """Shortcut to save a given figure in the correct dir"""
        if self.save_figures:
            dir_ = "/" + dir_ if dir_ != "" else ""
            cm.save_figure(fig, name, f"input_analysis{dir_}", self.stem)

    def plot_separation(self, survey_name_to_radius: dict):
        """Generate separation scatter plots and histograms to each of the survey_names in [survey_name_to_radius].
        parameters:
            [survey_name_to_radius]: dict<str, float>
                Dictionary with the survey names as keys and the adopted radii for
                the initial matches as values.
        """
        for survey_name, radius in survey_name_to_radius.items():
            scatter_fig, hist_fig = sep_p.plot_separation(
                self.df, survey_name, radius)
            self._quicksave_fig(
                scatter_fig, f"{survey_name}_sep_histogram", "separation")
            self._quicksave_fig(
                hist_fig, f"{survey_name}_sep_scatter", "separation")


if __name__ == "__main__":
    i_p_c = InputPlotContainer(save_figures=False)
    radius_dict = {"vhs": 0.5,
                   "eros": 0.1, "hsc": 0.25, "galex": 3.5, "kids": 1.5}
    # radius_dict = {"vhs": 0.5}
    i_p_c.plot_separation(radius_dict)
