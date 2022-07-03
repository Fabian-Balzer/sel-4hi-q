# %%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import util.configure_matplotlib as cm
import util.my_tools as mt

import output_scripts.specz_photz_plots as s_p
import output_scripts.template_analysis_plots as ta


class OutputPlotContainer:
    """Convenience class containing the output- and template data."""

    def __init__(self, save_figures=False):
        self.df = mt.read_saved_df(cat_type="out")
        # self.df = mt.read_saved_df(cat_type="in_processed")
        self.df_p = self.df[self.df["Type"] == "pointlike"]
        self.df_e = self.df[self.df["Type"] == "extended"]
        self.ext = "extended" in mt.USED_TTYPES
        self.plike = "pointlike" in mt.USED_TTYPES
        self.produce_both = self.ext & self.plike
        if not (self.plike | self.ext):
            mt.LOGGER.error(
                "Please specify a valid ttype [pointlike, extended or both] to plot the specz-photoz-plot.")
        self.stem = mt.CUR_CONFIG["LEPHARE"]["output_stem"]
        self.save_figures = save_figures
        try:
            self.template_df = mt.read_template_library()
        except FileNotFoundError:
            mt.LOGGER.error("Couldn't locate template files.")

    def _quicksave_fig(self, fig, name):
        """Shortcut to save a given figure in the correct dir"""
        if self.save_figures:
            cm.save_figure(fig, name, "output_analysis", self.stem)

    def plot_specz_photo_z(self):
        """Wrapper function for producing specz_photo_z scatter plots"""
        if self.produce_both:
            fig, axes = plt.subplots(2, 2, figsize=cm.set_figsize(
                width=cm.LATEXWIDTH, height=0.5 * cm.LATEXWIDTH, fraction=1))
            main_plike, secondary_plike = axes[0]
            main_ext, secondary_ext = axes[1]
            s_p.set_subplot_positions(
                main_plike, secondary_plike, is_left=True)
            s_p.set_subplot_positions(main_ext, secondary_ext, is_right=True)
            s_p.plot_ttype(self.df, "pointlike", main_plike, secondary_plike)
            s_p.plot_ttype(self.df, "extended", main_ext, secondary_ext)
            self._quicksave_fig(fig, "both_spec_z_phot_z")
            return
        fig, (main, secondary) = plt.subplots(
            2, 1, figsize=cm.set_figsize(width=0.5 * cm.LATEXWIDTH, height=0.5 * cm.LATEXWIDTH, fraction=1))
        s_p.set_subplot_positions(main, secondary)
        ttype = "extended" if self.ext else "pointlike"
        s_p.plot_ttype(self.df, ttype, main, secondary)
        self._quicksave_fig(fig, f"{ttype}_spec_z_phot_z")

    def plot_template_numbers(self):
        """Wrapper function to plot the template numbers"""
        if self.produce_both:
            fig, (ax_plike, ax_ext) = plt.subplots(1, 2, figsize=cm.set_figsize(
                width=cm.LATEXWIDTH, height=0.5 * cm.LATEXWIDTH, fraction=1))
            ta.plot_bar_template_outliers(self.df, ax_plike, "pointlike")
            ta.plot_bar_template_outliers(self.df, ax_ext, "extended")
            self._quicksave_fig(fig, "both_used_templates")
            return
        fig, ax = plt.subplots(1, 1, figsize=cm.set_figsize(
            width=0.5 * cm.LATEXWIDTH, height=0.5 * cm.LATEXWIDTH, fraction=1))
        ttype = "extended" if self.ext else "pointlike"
        ta.plot_bar_template_outliers(self.df, ax, ttype)
        self._quicksave_fig(fig, f"{ttype}_used_templates")

    def plot_template_scores(self):
        """Wrapper function to plot the template numbers"""
        ttype_list = ["extended"] if self.ext else ["pointlike"]
        if self.produce_both:
            ttype_list = ["extended", "pointlike"]
        for ttype in ttype_list:
            fig, ax = plt.subplots(1, 1, figsize=cm.set_figsize())
            ta.plot_bar_template_scores(self.df, ax, ttype)
            self._quicksave_fig(fig, f"{ttype}_used_templates")

    def plot_color_vs_redshift(self, c1: str, c2: str, fitted_only=True, plot_sources=True, temp_nums: dict = None):
        """Wrapper function to produce a color-redshift plot for the templates and sources.
        Parameters:
            c1: str
                Color 1
            c2: str
                Color 2 --> The color on the y-axis will be c1-c2.
            fitted_only: bool=True
                Only plot the templates that have been fit.
            plot_sources: bool=True
                Scatter all sources with spec-z.
            temp_nums: dict<Int, str>
        """
        temp_df = self.template_df[self.template_df["model"].apply(
            lambda num: num in temp_nums)] if temp_nums is not None else self.template_df
        if self.produce_both:
            fig, (ax_plike, ax_ext) = plt.subplots(1, 2, figsize=cm.set_figsize(
                width=cm.LATEXWIDTH, height=0.5 * cm.LATEXWIDTH, fraction=1))
            fig.set_tight_layout(True)
            ta.plot_color_versus_redshift(
                self.df, ax_plike, "pointlike", temp_df, c1, c2, fitted_only, plot_sources)
            ta.plot_color_versus_redshift(
                self.df, ax_ext, "extended", temp_df, c1, c2, fitted_only, plot_sources)
            self._quicksave_fig(fig, f"both_template_plot_{c1}-{c2}")
            return
        fig, ax = plt.subplots(
            1, 1, figsize=cm.set_figsize(width=0.5 * cm.LATEXWIDTH, height=0.5 * cm.LATEXWIDTH, fraction=1))
        ttype = "extended" if self.ext else "pointlike"
        ta.plot_color_versus_redshift(
            self.df, ax, ttype, temp_df, c1, c2, fitted_only, plot_sources)
        self._quicksave_fig(fig, f"{ttype}_template_plot_{c1}-{c2}")


if __name__ == "__main__":
    o_p_c = OutputPlotContainer(save_figures=True)
    # o_p_c.plot_specz_photo_z()
    # o_p_c.plot_template_numbers()
    # o_p_c.plot_color_vs_redshift(
    #     "W1", "W2", fitted_only=True, plot_sources=True)
    # o_p_c.plot_color_vs_redshift(
    # "g", "r", fitted_only=True, plot_sources=True)
    # o_p_c.plot_color_vs_redshift(
    #     "FUV", "r", fitted_only=True, plot_sources=True)
    # o_p_c.plot_template_numbers()
    o_p_c.plot_template_scores()

# %%
