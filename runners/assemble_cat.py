"""Script acting as a decoy pipe to the jython action going on in the background."""

import argparse
import os
import shlex
import subprocess

import util.my_tools as mt
from util.my_logger import logger


def read_args():
    """Reads out the arguments given by the user."""
    parser = argparse.ArgumentParser(
        description="Collect data for a LePhare run and run it.")
    parser.add_argument(
        "stem_name", help="The name of the directory to be created for the run.")
    parser.add_argument("-m", "--use_matched", action="store_true",
                        help="Specify wether to use already matched data with the given stem name")
    parser.add_argument("-p", "--use_processed", action="store_true",
                        help="Specify wether to use already processed data with the given stem name")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Increase output verbosity.")
    parser.add_argument("-z", "--reduce_to_specz",
                        help="Specify on wether to discard all sources without spec-z", action="store_true")
    parser.add_argument("-l", "--omit_lephare_input",
                        help="Specify on wether to write the LePhare input file", action="store_true")
    parser.add_argument("-w", "--write_info_file",
                        help="Write a short info file about the matches", action="store_true")

    # parser.add_argument("-r", "--rerun_filters", action="store_true", help="Use LePhare to rerun the filter files")
    parser.add_argument(
        "--lephare_stem", help="Set the unique stem for LePhare if required (defaults to the other given stem if not provided)", default="other_stem")
    args = parser.parse_args()
    args.lephare_stem = args.stem_name if args.lephare_stem == "other_stem" else args.lephare_stem
    print_verbose(args.verbose, "This script is run with the following arguments:\n" +
                  "\n".join([f"{key}:\t{val}" for key, val in vars(args).items()]))
    return args


def get_input_assembly_specifications(input_stem, lephare_stem):
    """Prompts the user for several steps about using default or testing data."""
    answer_dict = {"use_matched": False,
                   "use_processed": False,
                   "reduce_to_specz": False,
                   "omit_lephare_input": False,
                   "write_output_file": True}
    matched_fpath = mt.DATAPATH + \
        f"matches/{input_stem}_latest_match_backup.fits"
    if os.path.isfile(matched_fpath):
        answer_dict["use_matched"] = mt.get_yes_no_input(
            f"Would you like to use the catalogue that has already been matched for {input_stem}?")
    processed_fpath = mt.DATAPATH + \
        f"matches/{input_stem}_pointlike_processed_table.fits"
    if os.path.isfile(processed_fpath):
        answer_dict["use_processed"] = mt.get_yes_no_input(
            f"Would you like to use already processed data for {input_stem}?")
    answer_dict["reduce_to_specz"] = mt.get_yes_no_input(
        "Would you like your input catalogue to only contain sources with spec-z available?")
    answer_dict["reduce_to_specz"] = mt.get_yes_no_input(
        "Would you like to omit writing the LePhare input files?")


def arrange_input_files(args):
    """Runs the jython script to match the input files."""
    match_table_string = f"{JYSTILTS} '{WORKPATH}lephare_scripts/jystilts/match_tables.py' {args.stem_name} {args.lephare_stem} {args.use_matched} {args.use_processed} {args.reduce_to_specz} {args.write_info_file} {args.omit_lephare_input}"
    subprocess.call(shlex.split(match_table_string))
    print_verbose(
        args.verbose, "Successfully ran the jython program to match and write the tables.")
    return


if __name__ == "__main__":
    user_args = read_args()
    arrange_input_files(user_args)
