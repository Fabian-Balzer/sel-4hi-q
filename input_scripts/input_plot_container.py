# %%
import matplotlib.pyplot as plt
import numpy as np
import output_scripts.specz_photz_plots as s_p
import pandas as pd
import util.configure_matplotlib as cm
import util.my_tools as mt

import input_scripts.availability_plots as av
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

    def plot_magnitude_dist(self, band_list=tuple(mt.BAND_LIST), context=-1):
        """Plot the magnitude distribution of sources for a given set of photometric bands, each in their individual plot.
        Parameters:
            [bands]: tuple<str>
                Plots for these requested photometric bands will be produced.
            [context]: int = -1
                Only consider the bands fitting the context
        """
        band_list = tuple([
            band for band in band_list if band in mt.give_bands_for_context(context)])
        fig = m_d.plot_magnitude_dist(self.df, band_list, split_ttype=True)
        name = "magnitude_distribution"
        if context > 0:
            name += "_" + str(context)
        self._quicksave_fig(fig, name)

    def plot_input_dist(self, band_list=tuple(mt.BAND_LIST), context=-1, consider_specz=True):
        """Produce a bar plot of the distribution of sources for a given set of 
        photometric bands.
        Parameters:
            [bands]: tuple<str>
                Plots for these requested photometric bands will be produced.
            [context]: int = -1
                Only consider the bands fitting the context
            [consider_specz]: bool = True
                Whether to plot an additional bar showing the amount of spec-z
        """
        band_list = tuple([
            band for band in band_list if band in mt.give_bands_for_context(context)])
        fig = av.plot_input_distribution(self.df, band_list, consider_specz)
        name = "availability_distribution"
        if context > 0:
            name += "_" + str(context)
        self._quicksave_fig(fig, name)

    def plot_comparison_photz_vs_specz(self):
        """Produces a comparison plot to analyse the photo-z from the RF-classifier of the
        optically selected AGN."""
        df = self.df.rename(columns={"opt_agn_phot_z": "ZBEST"})
        df = mt.add_filter_columns(df)
        fig, axes = plt.subplots(2, 2, figsize=cm.set_figsize(
            width=cm.LATEXWIDTH, height=0.5 * cm.LATEXWIDTH, fraction=1))
        main_plike, secondary_plike = axes[0]
        main_ext, secondary_ext = axes[1]
        s_p.set_subplot_positions(
            main_plike, secondary_plike, is_left=True)
        s_p.set_subplot_positions(main_ext, secondary_ext, is_right=True)
        s_p.plot_ttype(df, "pointlike", main_plike, secondary_plike)
        s_p.plot_ttype(df, "extended", main_ext, secondary_ext)
        self._quicksave_fig(fig, "optically_selected_agn_spec_z_phot_z")


if __name__ == "__main__":
    i_p_c = InputPlotContainer(save_figures=True)
    radius_dict = {"vhs": 0.5,
                   "eros": 0.1, "hsc": 0.25, "galex": 3.5, "kids": 1.5}
    # i_p_c.plot_separation(radius_dict)
    # i_p_c.plot_comparison_photz_vs_specz()
    # i_p_c.plot_magnitude_dist(context=mt.CONTEXT)
    # i_p_c.plot_input_dist(context=mt.CONTEXT)
