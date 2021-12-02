import logging
import sys

STREAM_FORMAT = logging.Formatter('%(asctime)s- %(levelname)s - %(message)s')
FILE_FORMAT = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')


def get_console_handler():
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(STREAM_FORMAT)
    console_handler.setLevel(logging.INFO)
    return console_handler


def get_file_handler(log_file):
    file_handler = logging.FileHandler(log_file, mode='w')
    file_handler.setFormatter(FILE_FORMAT)
    file_handler.setLevel(logging.INFO)
    return file_handler


def get_logger(name, log_file):
    #  Create a custom logger
    logger = logging.getLogger(name)
    logger.propagate = False
    logger.setLevel(logging.INFO)
    # Create handlers
    logger.addHandler(get_console_handler())
    logger.addHandler(get_file_handler(log_file))
    return logger
