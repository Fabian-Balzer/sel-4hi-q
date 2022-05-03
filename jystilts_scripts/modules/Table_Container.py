import table_io as t_io


class TableContainer(dict):
    """Custom table container with each of the tables stored in a dictionary."""

    def __init__(self, arg_dict):
        self.arg_dict = arg_dict
        self.stem = arg_dict["stem"]
        self.match = None
        self.required_tables = ["vhs", "sweep", "opt_agn",
                                "eros", "hsc", "kids", "ls10"]  # , "xray_agn"]  # "assef_agn"
        if not arg_dict["use_matched"]:
            self.load_tables()
            self.pre_clean_tables()
        else:
            self.match, self.cats_used_for_matching = t_io.read_match_from_backup(
                self.stem, self.required_tables + ["galex"])
        if arg_dict["use_processed"]:
            print("Loading the processed tables.")
            for ttype in ["pointlike", "extended"]:
                self[ttype] = t_io.read_processed_input_for_analysis(
                    self.stem, ttype)
                if arg_dict["test_only"]:  # In this case, select only a subset with SPECZ available
                    self[ttype] = t_io.filter_for_testing(self[ttype])
            if arg_dict["write_info"]:
                t_io.write_info_file(
                    self["pointlike"], self["extended"], self.stem, self.cats_used_for_matching)
        super(TableContainer, self).__init__()

    def __str__(self):
        text = ""
        for table in self:
            text = text + "There are %s objects in the joint %s table.\n" % (
                self[table].count_rows(), table)
        return text

    def pre_clean_tables(self):
        """Perform the necessary cleaning steps for the tables even before the matching process."""
        for tablename in self:
            self[tablename] = t_io.pre_clean_table(tablename, self[tablename])

    def load_tables(self):
        """Read the tables."""
        for tablename in self.required_tables:
            self[tablename] = t_io.read_table(tablename)

    def match_tables(self, cat_list=None):
        """Matches all tables in the given table_list, or simply all tables provided."""
        cat_list = self.keys() if cat_list is None else cat_list
        self.match = t_io.match_given_tables(self, cat_list)
        t_io.write_match_as_backup(self.match, self.stem)
        self.cats_used_for_matching = cat_list

    def process_match(self):
        """Processes the matched table and writes pointlike and extended sub-tables."""
        self.match = t_io.add_separation_columns(
            self.match, self.cats_used_for_matching)
        # Processing step to split the pointlike and extended part and convert all fluxes:
        self["pointlike"], self["extended"] = t_io.process_table(
            self.match, self.cats_used_for_matching)
        for ttype in ["pointlike", "extended"]:
            # Get rid of possibly faultily-matched photometry:
            self[ttype] = t_io.discard_problematic_matches(self[ttype])
            t_io.write_processed_input_for_analysis(
                self[ttype], self.stem, ttype)
            # In this case, select only a subset with SPECZ available
            if self.arg_dict["test_only"]:
                self[ttype] = t_io.filter_for_testing(self[ttype])
        if self.arg_dict["write_info"]:
            t_io.write_info_file(
                self["pointlike"], self["extended"], self.stem, self.cats_used_for_matching)

    def process_for_lephare(self, ttype, excluded_bands=(), used_cats=("sweep", "vhs", "galex"), make_fits=False, stem_alt=None):
        """Function to clean the tables and output them such that they can be used by LePhare.
        Parameters:
            ttype: pointlike or extended (the subset of the general catalogue
            context: the lephare context provided

            used_cats: Tuple of catalogue names that you want to be included in LePhare (order doesn't matter as the order is forced anyways)

            make_fits: Bool to determine wether to additionally write a fits output file

            stem_alt: Alternative string name for the stem (handy if doing testing with different contexts)"""
        if not ttype in self.keys():
            print("You did not initialize the " + ttype +
                  " table. Maybe you'll need to call process_match() first.")
            return
        self[ttype + "_lephare"] = t_io.process_for_lephare(
            self[ttype], excluded_bands, used_cats)
        stem = stem_alt if stem_alt is not None else self.arg_dict["lephare_stem"]
        t_io.write_lephare_input(self[ttype + "_lephare"], ttype, stem)
        if make_fits:
            t_io.write_lephare_input(
                self[ttype + "_lephare"], ttype, stem, ofmt="fits")

    # def write_ra_and_dec(self):
    #     ra_dec = self.tables["matches"].cmd_keepcols("ra dec")
    #     t_io.write_matched_table(ra_dec, "ra_and_dec_sources.fits")
