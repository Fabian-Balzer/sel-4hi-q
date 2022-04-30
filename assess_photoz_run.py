"""Script to directly assess the quality of a photo-z run."""

import argparse
import os

import util.my_tools as mt



def read_args():
    """Reads out the arguments given by the user."""
    parser = argparse.ArgumentParser(
        description="Assess a LePhare output file.")
    parser.add_argument(
        "output_fpath", help="The file path of the output file.")
    parser.add_argument(
        "ttype", help="The table type of the output file (extended or pointlike).")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Increase output verbosity.")
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = read_args()
    df = mt.read_fits_as_dataframe(args.output_fpath)
    df = mt.add_filter_columns(df)
    stat_dict = mt.give_output_statistics(df)
    etalabel = "$\eta = " + f"{stat_dict['eta']:.3f}$\n"
    sig_nmadlabel = r"$\sigma_{\rm NMAD} = " + f"{stat_dict['sig_nmad']:.3f}$"
    fpos = df['IsFalsePositive'].sum()
    fposlabel = "\n" + r"$\psi_{\rm Pos} = " + \
        f"{stat_dict['psi_pos']:.3f}$ ({fpos})"
    fneg = df['IsFalseNegative'].sum()
    fneglabel = "\n" + r"$\psi_{\rm Neg} = " + \
        f"{stat_dict['psi_neg']:.3f}$ ({fneg})"
    label = f"{len(df)} sources\n{etalabel}{sig_nmadlabel}"
    if args.verbose:
        print(
            f"The outlier fraction for {args.ttype} is eta = {stat_dict['eta']:.4f}")
        print(
            f"The accuracy for {args.ttype} is sig_NMAD = {stat_dict['sig_nmad']:.4f}")
        print(
            f"The false pos fraction is psi_pos = {stat_dict['psi_pos']:.4f}")
        print(
            f"The false neg fraction is psi_neg = {stat_dict['psi_neg']:.4f}")
    else:
        print(
            f"Eta = {stat_dict['eta']:.4f},\tSig_NMAD = {stat_dict['sig_nmad']:.4f}")
