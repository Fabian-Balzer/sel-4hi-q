#!/usr/bin/python

""" Plot (file or window) the observed magnitudes
of the objects, along with best-fit templates,
reading the info form .zsp and .pdz output
files of LePhare.
author: F. Balzer (adopted from scripts by M. Salvato and O. Ilbert)
"""


import sys

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from util.spec_helper import LOGGER, Spectrum, read_args, read_specnames

######## GET FILENAMES AND OPTIONS ########


ARGS = read_args()

SPEC_FNAMES = read_specnames(ARGS)

NUM_SPEC = len(SPEC_FNAMES)  # N of .spec files

LOGGER.debug(SPEC_FNAMES)

if NUM_SPEC == 0:
    LOGGER.error('Please specify at least one valid .spec file')
    LOGGER.error('Try -h or --help options to get help.')
    sys.exit()


if ARGS.device == ".screen":
    if NUM_SPEC > 1:
        LOGGER.warning("""Since you have chosen 'screen'
as output device, only the first spectrum will
be plotted although you have provided more.""")
    LOGGER.info("Producing plot for single spectrum, showing it on screen.")
    spec = Spectrum(SPEC_FNAMES[0], ARGS)
    plt.ioff()
    plt.show()
    sys.exit()

LOGGER.info("\n%s\n  Creating %d %s files...   \n%s\n",
            40 * "*", NUM_SPEC, ARGS.device, 40 * "*")

####### LOOP OVER .SPEC FILES ########

if ARGS.device == ".multi":
    multipage_name = f"{ARGS.output}.pdf" if ARGS.output != "" else "MULTISPEC.pdf"
    LOGGER.info("All objects will be collected in a single pdf file named:")
    LOGGER.info("--> %s", multipage_name)
    PDF_MULTIPAGE = PdfPages(multipage_name)
    for spec_fname in SPEC_FNAMES:
        spec = Spectrum(spec_fname, ARGS)
        spec.save_plot_to_multi(PDF_MULTIPAGE)
    PDF_MULTIPAGE.close()
    LOGGER.info("Successfully created plots for %d .spec files, saving them in %s.",
                NUM_SPEC, multipage_name)
    sys.exit()


# Plot SEDs
for spec_fname in SPEC_FNAMES:
    spec = Spectrum(spec_fname, ARGS)
    savename = spec_fname.split(".")[0] + ARGS.output + ARGS.device
    LOGGER.debug(savename)
    spec.save_single_plot(savename)
LOGGER.info("Successfully created plots for %d .spec files.", NUM_SPEC)
LOGGER.info("These are in the pattern *old_name*%s%s",
            ARGS.output, ARGS.device)
