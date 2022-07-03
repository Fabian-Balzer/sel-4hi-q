"""Script to read a zfix run, plot the used templates,
and write a new template list file based on the old one."""
# %%
import output_scripts.template_analysis_plots as t_a

import util.my_tools as mt


def write_zfix_improved_temp_list(ttype):
    """Writes a new template list containing only the templates that were fit
    to any sources of the current run."""
    templist_fpath = mt.give_temp_listname(ttype)
    # Open the old template file
    with open(templist_fpath, "r", encoding="utf-8") as f:
        old_temps = [line for line in f.readlines()
                     if not line.startswith("#")]
    temp_dict = {num: temp for num, temp in enumerate(old_temps, start=1)}
    output_df = mt.read_saved_df()
    models_used = set(output_df[output_df["Type"] == ttype]["MOD_BEST"])

    correct_temps = [temp_dict[num] for num in models_used if num != -99]
    mt.LOGGER.info("Only %d templates of initially %d templates remaining for %s (%s)",
                   len(correct_temps), len(old_temps), ttype, templist_fpath.split("/")[-1])
    new_templist_fpath = mt.give_temp_listname(
        ttype, mt.CUR_CONFIG['LEPHARE']['template_stem'] + "_reduced")

    with open(new_templist_fpath, "w", encoding="utf-8") as f:
        f.writelines(correct_temps)


def write_combined_temp_list(ttype, stem1: str, stem2: str):
    """Writes a new template list containing only three quarters of the templates with the highest scores."""
    templist_fpath_1 = mt.give_temp_listname(ttype, altstem=stem1)
    templist_fpath_2 = mt.give_temp_listname(ttype, altstem=stem2)
    # Open both template files
    with open(templist_fpath_1, "r", encoding="utf-8") as f:
        temps_1 = [line for line in f.readlines()
                   if not line.startswith("#")]
    with open(templist_fpath_2, "r", encoding="utf-8") as f:
        temps_2 = [line for line in f.readlines()
                   if not line.startswith("#")]
    combined_temps = []
    for temp in temps_1 + temps_2:
        if temp not in combined_temps:
            combined_temps.append(temp)
    new_templist_fpath = mt.give_temp_listname(ttype, altstem="combined")
    with open(new_templist_fpath, "w", encoding="utf-8") as f:
        f.writelines(combined_temps)
    mt.LOGGER.info(
        "Combined templates have been written to '%s'", new_templist_fpath)


def write_score_improved_temp_list(ttype):
    """Writes a new template list containing only three quarters of the templates with the highest scores."""
    df = mt.read_saved_df()
    good_score_temps = t_a.give_templates_to_keep(df, ttype)
    # Convert the numbers into names
    temp_list = [mt.get_temp_name_for_num(
        ttype, num) for num in good_score_temps]
    new_templist_fpath = mt.give_temp_listname(
        ttype, altstem=mt.CUR_CONFIG["LEPHARE"]["template_stem"] + "_score_reduced")
    with open(new_templist_fpath, "w", encoding="utf-8") as f:
        f.write('\n'.join(temp_list))
    mt.LOGGER.info(
        "Combined templates have been written to '%s'", new_templist_fpath)


if __name__ == "__main__":
    for TTYPE in ["pointlike", "extended"]:
        write_combined_temp_list(
            TTYPE, "baseline_score_reduced", "ananna_score_reduced")
        write_score_improved_temp_list(TTYPE)
    #     write_zfix_improved_temp_list(TTYPE)
