# %%
import pandas as pd
import util.my_tools as mt


def generate_template_list_info(ttype, stem=""):
    """Read the template list file and store the information in an output file."""
    fname = mt.DATAPATH + \
        f"lephare_files/templates/lists/{ttype}_{stem}.list"
    templates = {}
    with open(fname, "r") as f:
        for line in f.readlines():
            if not line.startswith("#"):
                tempname = line.split("/")[-1]  # Remove filepath
                directory = line.split("/")[0]  # Additionally save directory
                tempname = tempname.split(".")[0]  # Remove ending
                templates[tempname] = directory
    name_dict = {"Sb": "Spiral bREF1", "Spi4": "Spiral cREF1",
                 "M82": "StarburstREF1", "Sey18": "Seyfert 1.8REF1", "Sey2": "Seyfert 2REF1", "pl_QSOH": "High lum. QSOREF4", "pl_QSOREF4": "Low lum. QSOREF4", "pl_TQSO1": "High IR lum. QSOREF4", "BC03": "Blue starformingREF2", "S0": "S0REF1"}
    type_dict = {"Name": {}, "Type": {}}
    for i, temp in enumerate(templates, start=1):
        if temp.endswith("_template_norm"):
            temp = temp[:-len("_template_norm")]
            temptype = temp
        if temp.endswith("_temp_restframe"):
            temp = temp[:-len("_temp_restframe")]
            temptype = temp
        temptype = name_dict[temp] if temp in name_dict else "TODO"
        if ("S0" in temp and "QSO2" in temp):
            temptype = "HybridREF3"
        if temptype == "TODO" and ttype == "star":
            temptype = templates[temp].lower().capitalize()
        type_dict["Name"][i] = temp
        type_dict["Type"][i] = temptype
    df = pd.DataFrame(type_dict)
    df = df.rename(columns={"Name": "NameREF0"})
    caption = f"Template Library used for the {ttype} sources"
    label = f"tab:{ttype}_temp_lib"
    table_text = df.to_latex(
        caption=caption, label=label, position="htbp")
    table_text = table_text.replace(r"{}", "ID").replace("TODO", "\TODO{}")
    footnotetexts = ["Names as in the original libraries",
                     r"\cite{2009IlbertRedshiftPrecision}", r"\cite{2007Polletta}", r"\cite{2009Salvato}", r"Template from \url{http://classic.sdss.org/dr5/algorithms/spectemplates/}"]
    # Still need two closing brackets at the end
    notes = r"\multicolumn{3}{l}{\makecell[l]{\textbf{Notes.}\\"
    for i, letter in enumerate("abcde"):
        # Then create references to that footnote
        if f"REF{i}" in table_text:
            notes += f"$^{letter}$ {footnotetexts[i]}." + r"\\"
        table_text = table_text.replace(
            f"REF{i}", f"$^{letter}$")
    table_text = table_text.replace(
        r"\bottomrule", r"\bottomrule" + notes + r"}}")
    fname = mt.MYDIR + f"other/{ttype}_templates.tex"
    with open(fname, "w", encoding="utf8") as f:
        f.write(table_text)
        print(
            f"The LaTeX input text for the {ttype} template library has been written to {fname}")


if __name__ == "__main__":
    generate_template_list_info("pointlike", "baseline_templates")
    generate_template_list_info("extended", "baseline_templates")
    generate_template_list_info("star", "baseline_templates")
