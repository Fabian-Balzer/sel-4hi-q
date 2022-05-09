"""Helper module for jython scripts"""
import logging
import os

from ConfigParser import ConfigParser

CONFIGPATH = os.environ["LEPHARE"] + "/lephare_scripts/config/"
GEN_CONFIG = ConfigParser()
GEN_CONFIG.read(CONFIGPATH + "general.ini")
CUR_CONFIG = ConfigParser()
CUR_CONFIG.read(CONFIGPATH + GEN_CONFIG.get("PATHS", "current_config"))


def init_logger():
    """Initializes a logger."""
    # create logger
    logger = logging.getLogger('simple_logger')
    logger.setLevel(logging.DEBUG)
    # create console handler and set level to debug
    ch = logging.StreamHandler()
    ch.setLevel(CUR_CONFIG.get("GENERAL", "logging_level"))
    # create formatter
    formatter = logging.Formatter('%(levelname)s: %(message)s')
    # add formatter to ch
    ch.setFormatter(formatter)
    # add ch to logger
    logger.addHandler(ch)
    return logger


LOGGER = init_logger()
