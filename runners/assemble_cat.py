"""Script acting as a decoy pipe to the jython action going on in the background."""

import argparse
import os
import shlex
import subprocess

WORKPATH = os.environ["LEPHARE"] + "/"  # Path for everything
LEPHAREDIR = os.environ["LEPHAREDIR"] + "/"
JYSTILTS = f"java -jar '{WORKPATH}static_files/programs/jystilts.jar'"


def print_verbose(verbose, statement):
    """Helper function to print only if verbose and separate the text from the rest"""
    if verbose:
        print("-*" * 40)
        print(statement)
        print("-*" * 40)


def get_yes_no_input(question):
    """Tries to get user input for a yes/no question."""
    answer = input(question + "\n>>> ")
    while True:
        if answer.lower() in ["y", "yes", "yep"]:
            return True
        if answer.lower() in ["n", "no", "nope"]:
            return False
        answer = input("Please answer with 'yes' or 'no'\n>>> ")


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
