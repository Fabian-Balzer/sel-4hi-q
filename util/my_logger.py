"""Short script to set up a reliable logger."""
import logging


def init_logger():
    """Initializes a logger."""
    # create logger
    logger = logging.getLogger('simple_logger')
    logger.setLevel(logging.DEBUG)
    # create console handler and set level to debug
    ch = logging.StreamHandler()
    # create formatter
    formatter = logging.Formatter('%(levelname)s: %(message)s')
    # add formatter to ch
    ch.setFormatter(formatter)
    # add ch to logger
    while logger.hasHandlers():
        logger.removeHandler(logger.handlers[0])
    logger.addHandler(ch)
    return logger


LOGGER = init_logger()
