import jy_tools as jt
import table_io as t_io


class TableContainer(dict):
    """Custom table container with each of the tables stored in a dictionary."""

    def __init__(self):
        super(TableContainer, self).__init__()
        self.match = None
        self.stem = jt.CUR_CONFIG.get("CAT_ASSEMBLY", "cat_stem")
        self.load_tables()
        self.pre_clean_tables()
        self.match_tables()
        self.process_match()
        for ttype in ["pointlike", "extended"]:
            self.process_for_lephare(ttype)

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
        # , "xray_agn"]  # "assef_agn"
        cats = ["vhs", "sweep", "opt_agn", "eros", "hsc", "kids", "ls10"]
        for tablename in cats:
            self[tablename] = t_io.read_table(tablename)
        jt.LOGGER.debug(
            "Loaded tables from the following catalogs:\n%s", cats)

    def match_tables(self):
        """Matches all tables in the given table_list, or simply all tables provided."""
        if jt.CUR_CONFIG.getboolean("CAT_ASSEMBLY", "use_matched"):
            self.match = jt.read_match_from_backup()
        else:
            self.match = t_io.match_given_tables(self)
            jt.write_match_as_backup(self.match)

    def process_match(self):
        """Processes the matched table and writes pointlike and extended sub-tables."""
        if jt.CUR_CONFIG.getboolean("CAT_ASSEMBLY", "use_processed"):
            for ttype in ["pointlike", "extended"]:
                self[ttype] = jt.read_processed_table(ttype)
        else:
            with_sep = t_io.add_separation_columns(self.match)
            # Processing step to split the pointlike and extended part and convert all fluxes:
            self["pointlike"], self["extended"] = t_io.process_table(with_sep)
            for ttype in ["pointlike", "extended"]:
                # Get rid of possibly faultily-matched photometry:
                self[ttype] = t_io.discard_problematic_matches(self[ttype])
                jt.write_processed_table(self[ttype], ttype)
        for ttype in ["pointlike", "extended"]:
            self[ttype] = t_io.filter_for_testing(self[ttype])
        jt.write_info_file(self["pointlike"], self["extended"])

    def process_for_lephare(self, ttype):
        """Function to clean the tables and output them such that they can be used by LePhare.
        Parameters:
            ttype: pointlike or extended"""
        if not ttype in self.keys():
            jt.LOGGER.warning("You did not initialize the %s table.", ttype)
            return
        if jt.CUR_CONFIG.getboolean("CAT_ASSEMBLY", "write_lephare_input"):
            self[ttype + "_lephare"] = t_io.process_for_lephare(self[ttype])
            jt.write_lephare_input(self[ttype + "_lephare"], ttype)

    # def write_ra_and_dec(self):
    #     ra_dec = self.tables["matches"].cmd_keepcols("ra dec")
    #     t_io.write_matched_table(ra_dec, "ra_and_dec_sources.fits")
