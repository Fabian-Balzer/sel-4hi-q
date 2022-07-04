#!/usr/bin/python
"""Helper file containing tools and the Spectrum class for the
spec.py program to analyse LePhare .spec files."""

import argparse
import logging
import os
import sys
from enum import Enum

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import FuncFormatter

custom_params = {"axes.titlesize": 30,
                 "axes.labelsize": 24,
                 "lines.linewidth": 2,
                 "lines.markersize": 10,
                 "xtick.labelsize": 16,
                 "ytick.labelsize": 16,
                 "font.size": 20}
mpl.rcParams.update(custom_params)


def init_logger():
    """Initializes a logger."""
    # create logger
    logger = logging.getLogger('simple_logger')
    logger.setLevel(logging.DEBUG)
    # create console handler and set level to debug
    ch = logging.StreamHandler()
    # create formatter
    formatter = logging.Formatter('%(levelname)s: %(message)s')
    # add formatter to ch
    ch.setFormatter(formatter)
    # add ch to logger
    while logger.hasHandlers():
        logger.removeHandler(logger.handlers[0])
    logger.addHandler(ch)
    return logger


LOGGER = init_logger()


class Device(Enum):
    """Possible outputs."""
    multi = 'multi'
    pdf = 'pdf'
    png = 'png'
    eps = 'eps'
    ps = 'ps'
    screen = 'screen'

    def __str__(self):
        return self.value


class YAxisUnit(Enum):
    """Possible y axis units."""
    flux = 'flux'
    mag = 'mag'

    def __str__(self):
        return self.value


def find_lower_exponent(context):
    """Searches for the next lowest exponent of 2 and returns it"""
    i = 0
    while True:
        if 2**i > context:
            break
        i += 1
    return i - 1


def convert_context_to_band_indices(context):
    """Decodes a given context and returns the used filter numbers.
    Works for a given pandas Series with a context column and also for a given
    number.
    Example:
        if context = 13, it will return [1, 3, 4], since 2^(1-1) + 2^2 + 2^3 = 13.
    """
    if context in [-1, 0]:  # Return all bands if context is -1
        return -1
    filter_numbers = []
    remaining_context = context
    while remaining_context > 0:
        last = find_lower_exponent(remaining_context)
        # +1 as the filter numbers start with 1
        filter_numbers.append(last + 1)
        # We now want to work with the remains.
        remaining_context = remaining_context - 2**last
    # lastly, we reverse the list to start with lowest nums
    ordered_nums = filter_numbers[::-1]
    LOGGER.info("Only plotting the following filters due to the provided context of %d:\n %s",
                context, str(ordered_nums))
    return ordered_nums


def read_args():
    """Reads out the arguments given by the user."""
    parser = argparse.ArgumentParser(
        description="Script to analyze LePhare's .spec output files."
        "SYNTAX: >>> python spec.py file[s].spec [OPTIONS].")
    parser.add_argument("-d", "--device", type=Device, choices=list(Device),
                        help="""Select the output device.
If 'multi' is chosen, all the plots are collected in a multi pdf file.
If the option is not set, then print(on screen (with
a limit of 1 object/window).""", default=Device.screen)
    parser.add_argument("-i", "--input_file", type=str, default="",
                        help="Specify an input file with a list of the desired .spec files.")
    parser.add_argument("-o", "--output", type=str, default="",
                        help="""If --device = 'multi', this option specifies the name of
the pdf file. With any other --device values, a string STR is appended at
the end of the filename, just before the extension (e.g. *STR.png).
Does nothing when printing on screen.""")
    parser.add_argument("-c", "--context", type=int, default=-1,
                        help="Specify the CONTEXT of the filters to be plotted.")
    parser.add_argument("-y", "--ytype", type=Device, choices=list(YAxisUnit),
                        help="""Select the y axis unit for the plot""", default=YAxisUnit.mag)
    parser.add_argument("-l", "--log_level", type=int, default=20,
                        help="""Set the logging level (conforming with the logging module)
CRITICAL: 50; ERROR: 40; WARNING: 30; INFO: 20; DEBUG: 10""")
    parser.add_argument("-f", "--find_input", action="store_true",
                        help="""Plot spectra for all .spec files
found in the current directory. This works in addition to specifying
.spec files directly.""")
    args, _ = parser.parse_known_args()
    loglevel = args.log_level if 10 <= args.log_level <= 50 else 20
    LOGGER.setLevel(loglevel)
    if args.log_level < 10 or args.log_level > 50:
        LOGGER.warning("Unknown log level set. Reverting to 'INFO'.")
    if args.device not in list(Device):
        LOGGER.error("Unknown device provided!")
        LOGGER.info("Please choose one of the following devices for output:")
        LOGGER.info(list(Device))
        sys.exit()
    else:
        LOGGER.info("Output device will be in ***%s*** format.", args.device)
        args.device = f".{args.device}"  # Convert from enum to string.
    args.allowed_filters = convert_context_to_band_indices(args.context)
    return args


def convert_mag_to_flux(mag_list):
    """Takes the magnitudes of a given np mag list and returns a list in fluxes"""
    ibad = np.where((mag_list <= 0) | (mag_list > 35))
    flux_list = -0.4 * (mag_list - 23.91)  # uJy
    flux_list[ibad] = -99.
    return flux_list


class Spectrum:
    """Convenience class containing the spectrum, responsible for
    producing the plots."""

    def __init__(self, fname: str, args: argparse.Namespace):
        with open(fname, "r") as f:
            lines = f.readlines()
        self.y_is_mag = args.ytype == YAxisUnit.mag  # flux or mag for y axis
        self.args = args
        self.fontsize = 20
        self._read_general_info(lines)
        self._read_model_params(lines)
        self._read_filter_data(lines)
        self._read_pdf_data(lines)
        self._read_model_data(lines)
        self._init_plot()
        self._plot_data()

    def _read_general_info(self, lines):
        self.id = lines[1].split()[0]
        self.zspec = float(lines[1].split()[1])
        self.zphot = float(lines[1].split()[2])
        self.num_filt = int(lines[3].split()[1])
        if 2**self.num_filt < self.args.context:
            LOGGER.warning(
                "The provided context exceeds the number of filters available.")
        self.num_pdf_entries = int(lines[5].split()[1])

    def _read_model_params(self, lines):
        model_keys = [name for name in lines[6].split() if "#" not in name]
        # Read out information on each of the models and store them in dicts:
        model_dict = {}
        model_names = ["gal_1", "gal_2", "gal_fir", "gal_stoch", "qso", "star"]
        model_colors = ['b', 'g', 'k', 'm', 'y', 'gray']
        start_index = 7
        for i, name in enumerate(model_names):
            values = lines[start_index + i].split()
            color = model_colors[i]
            model_dict[name] = self._read_single_model_info(
                values, model_keys, color)
        self.model_dict = model_dict

    def _read_single_model_info(self, values, keys, color):
        """Helper method to store the model information in a dictionary"""
        model_specs = dict(zip(keys, values))
        for key, val in model_specs.items():
            if key != "Type":
                if "." in val:
                    model_specs[key] = float(val)
                else:
                    model_specs[key] = int(val)
        model_specs["color"] = color
        return model_specs

    def _read_filter_data(self, lines):
        start_index = 13
        filter_keys = ["mag", "mag_err", "lambda_eff",
                       "filter_width", "_", "_", "_", "fcorr", "model_mag"]
        filter_dict = {}
        for key in filter_keys:
            filter_dict[key] = np.zeros(self.num_filt)
        # Iterate over all lines with filter information
        for i in range(self.num_filt):
            values = lines[start_index + i].split()
            for j, key in enumerate(filter_keys):
                filter_dict[key][i] = float(values[j])
        # convert mag(AB syst.) in log(flux)
        filter_dict["flux"] = convert_mag_to_flux(filter_dict["mag"])
        filter_dict["flux_err"] = convert_mag_to_flux(filter_dict["mag_err"])
        df = pd.DataFrame(filter_dict)  # Use pandas dataframe
        df = df.replace(99.0, np.NaN)
        df = df.replace(-99.0, np.NaN)
        df["lambda_eff"] /= 10000
        df["filter_width"] /= 10000
        df["context"] = 2**df.index
        if isinstance(self.args.allowed_filters, list):
            df["include_in_plot"] = pd.Series(df.index).apply(
                lambda filt_index: filt_index in self.args.allowed_filters)
        else:
            df["include_in_plot"] = True
        self.filter_df = df

    def _read_pdf_data(self, lines):
        start_index = 13 + self.num_filt
        zpdf = np.zeros([3, self.num_pdf_entries])
        for i in range(self.num_pdf_entries):
            line_vals = lines[start_index + i].split()
            zpdf[:, i] = np.array(line_vals)
        val_dict = {"z": zpdf[0, :], "p_1": zpdf[1, :], "p_2": zpdf[2, :]}
        df = pd.DataFrame(val_dict)
        # Norm the probabilities:
        df["p_1"] = df["p_1"] / max(df["p_1"])
        df["p_2"] = df["p_2"] / max(df["p_2"])
        self.pdf_df = df

    def _read_model_data(self, lines):
        start_index = 13 + self.num_filt + self.num_pdf_entries
        for model in self.model_dict.values():
            line_num = model["Nline"]
            data = np.zeros([2, line_num])
            if line_num > 0:
                for i in range(line_num):
                    values = lines[start_index + i].split()
                    data[:, i] = np.array(values)
            val_dict = {"wavelength": data[0, :] / 10000.,
                        "mag": data[1, :], "flux": convert_mag_to_flux(data[1, :])}
            df = pd.DataFrame(val_dict)
            model["data"] = df
            start_index += line_num  # start at the next line number

    def _init_plot(self):
        fig = plt.figure(1, figsize=(15, 10))
        fig.subplots_adjust(
            left=0.10, right=0.94, top=0.94, bottom=0.10, wspace=0.05, hspace=0.05)
        ax1 = fig.add_subplot(111)
        ylabel = r'$magnitude$' if self.y_is_mag else r'$F_{\nu}/\mu Jy$'
        ax1.set_xlabel(r'$\lambda/\mu m$')
        ax1.set_ylabel(ylabel)
        ax1.yaxis.set_label_coords(-0.07, 0.5)
        ax1.set_xscale("log")
        ax1.xaxis.set_major_formatter(
            FuncFormatter(lambda x, pos: str(int(round(x)))))
        ax1.set_xticks([1, 10, 100, 1000])
        ax1.set_title(
            f"ID: {self.id}, zspec: {self.zspec:.3f}, zphot: {self.zphot:.3f}", color='black')
        # PDZ plot
        pdf_ax = fig.add_axes([0.65, 0.15, 0.25, 0.20])
        pdf_ax.set_ylim([0.0, 1.1])
        pdf_ax.set_yticks([0.25, 0.5, 0.75])
        pdf_ax.set_xlabel("zphot", size=14)
        pdf_ax.set_ylabel("p(z)", size=14)
        pdf_ax.grid(True)
        pdf_ax.yaxis.tick_right()
        pdf_ax.tick_params(axis='both', which='major', labelsize=10)
        pdf_ax.set_title("Probability distribution", size=20)
        for tick in pdf_ax.yaxis.get_major_ticks():
            tick.label1On = False
            tick.label2On = True
        for tick in pdf_ax.xaxis.get_major_ticks():
            tick.label1On = False
            tick.label2On = True
        self.pdf_ax = pdf_ax
        self.main = ax1
        self.fig = fig

    def _plot_data(self):
        df = self.filter_df.sort_values(by="lambda_eff")
        df = df[df["include_in_plot"]]
        xmax = df.iloc[-1]["lambda_eff"] + 2 * df.iloc[-1]["filter_width"]
        xmin = df.iloc[0]["lambda_eff"] - 2 * df.iloc[0]["filter_width"]
        xmin = xmin if xmin > 0 else 1
        ymax = self.filter_df["mag"].min() - 1.5  # Inverted for magnitudes!
        ymin = self.filter_df["mag"].max() + 1.5
        info_cols = {"Type": "Type", "Model": "Model", "Zphot": "z",
                     "Chi2": r"$\chi^2$", "Extlaw": "Ext_law", "EB-V": r"$E_{B-V}$"}
        for i, val in enumerate(info_cols.values()):
            self.main.text(0.05 + 0.1 * i, 0.95, val,
                           transform=self.main.transAxes)
        # Plot the lines for the models:
        counter = 1
        for model in self.model_dict.values():
            if model["Library"] < 0:
                continue  # plot only models really used
            x, y = model["data"]["wavelength"], model["data"]["mag"],
            color = model["color"]
            self.main.plot(x, y, "-", color=color, label=model["Type"])
            for i, key in enumerate(info_cols.keys()):
                text = str(model[key]) if not isinstance(
                    model[key], float) else f"{model[key]:.3f}"
                self.main.text(0.05 + 0.1 * i, 0.95 - 0.05 * counter, text,
                               transform=self.main.transAxes, color=color)
            counter += 1
        good_vals = self.filter_df[self.filter_df["mag"].notna(
        ) & self.filter_df["mag_err"].notna() & (self.filter_df["mag_err"] < 2.)]
        x, y = good_vals["lambda_eff"], good_vals["mag"]
        xerr, yerr = good_vals["filter_width"] / 2, good_vals["mag_err"]
        self.main.errorbar(x, y, xerr=xerr, yerr=yerr, color="black",
                           fmt="o")
        # Plot values with low S/N separately in gray
        bad_vals = self.filter_df[self.filter_df["mag"].notna(
        ) & self.filter_df["mag_err"].notna() & (self.filter_df["mag_err"] >= 2.)]
        x, y = bad_vals["lambda_eff"], bad_vals["mag"]
        xerr, yerr = bad_vals["filter_width"] / 2, bad_vals["mag_err"]
        self.main.errorbar(x, y, xerr=xerr, yerr=yerr, color="gray",
                           fmt="o")
        self.main.set_xlim(xmin, xmax)
        self.main.set_ylim(ymin, ymax)

        # Plot the pdz:
        self.pdf_ax.plot(self.pdf_df["z"],
                         self.pdf_df["p_1"], 'b-', label="Best")
        self.pdf_ax.plot(self.pdf_df["z"],
                         self.pdf_df["p_2"], 'r-', label="Second")
        self.pdf_ax.axvline(self.zphot, linestyle="--",
                            color="gray", label=f"zphot-best ({self.zphot:.3f})")
        self.pdf_ax.axvline(self.zspec, linestyle="--",
                            color="black", label=f"zspec-best ({self.zspec:.3f})")
        self.pdf_ax.set_xlim(self.pdf_df["z"].min(), self.pdf_df["z"].max())
        self.pdf_ax.legend(prop={"size": 10})

    def save_single_plot(self, fname):
        """Save the plot of the spectrum to the given fname."""
        self.fig.savefig(fname, format=self.args.device.strip("."),
                         dpi=300, bbox_inches='tight')

    def save_plot_to_multi(self, multi_doc):
        """Append the plot on the multi-page plot"""
        self.fig.savefig(multi_doc, format='pdf',
                         dpi=300, bbox_inches='tight')


def get_specfiles_from_working_dir():
    """Finds and returns a list of files ending with .spec in the
    current working directory."""
    return [fname for fname in os.listdir() if fname.endswith(".spec")]


def get_specfiles_from_input_file(input_fname):
    """Finds and returns a list of files ending with .spec in the
    current working directory."""
    with open(input_fname, "r") as file:
        lines = file.readlines()
    return [fname for fname in lines if fname.endswith(".spec")]


def read_specnames(args: argparse.Namespace):
    """Reads the .spec filenames specified in the various possible ways.
    Already performs a check whether the files exist on the system.

    Returns:
        Filtered list of available .spec filenames
    """
    fnames = [fname for fname in sys.argv if fname.endswith(".spec")]

    if args.find_input:
        extra_fnames = [fname for fname in get_specfiles_from_working_dir()
                        if fname not in fnames]
        fnames += extra_fnames
        LOGGER.info(
            "Added %d .spec files from the working dir.", len(extra_fnames))

    if args.input_file != "":
        try:
            extra_fnames = [fname for fname in get_specfiles_from_input_file(args.input_file)
                            if fname not in fnames]
            fnames += extra_fnames
            LOGGER.info(
                "Added %d .spec files in the specified input file '%s'.", len(extra_fnames), args.input_file)
        except FileNotFoundError:
            LOGGER.error(
                "Unable to find the file '%s' containing a list of .spec files", args.input_file)
    available = [fname for fname in fnames if os.path.isfile(fname)]
    not_available = [fname for fname in fnames if fname not in available]
    if len(not_available) > 0:
        LOGGER.warning(
            "You provided %d files that could not be \
located and are skipped over.", len(not_available))
        LOGGER.info("These unavailable files are '%s'.",
                    "'\n'".join(not_available))
    return available
