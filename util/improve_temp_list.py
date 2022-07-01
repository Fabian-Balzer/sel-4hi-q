"""Script to read a zfix run, plot the used templates,
and write a new template list file based on the old one."""
# %%
import util.my_tools as mt

ttype = "extended"
templist_fpath = mt.give_temp_listname(ttype)

with open(templist_fpath, "r", encoding="utf-8") as f:
    old_temps = [line for line in f.readlines() if not line.startswith("#")]
temp_dict = {num: temp for num, temp in enumerate(old_temps, start=1)}

zfix_fpath = mt.give_lephare_filename(ttype, out=True, suffix=".fits")
df = mt.read_fits_as_dataframe(zfix_fpath)
df = df["MOD_BEST"]
mods_used = set(df)

correct_temps = [temp_dict[num] for num in mods_used if num != -99]
mt.LOGGER.info("%d templates of initially %d templates remaining for %s", len(
    correct_temps), len(old_temps), ttype)
new_templist_fpath = mt.give_temp_listname(ttype, "improved_templates")

with open(new_templist_fpath, "w", encoding="utf-8") as f:
    f.writelines(correct_temps)
