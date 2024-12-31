
import logging
import sys

def setup_logger(name="dedup_logger", log_file="dedup.log", verbosity=logging.DEBUG):
    """
    Set up a logger config
    
    Args:
        verbosity (int): The verbosity level for the logger. Defaults to logging.DEBUG.

    Returns:
        logging.Logger: The configured logger instance.
    """    
    # Create formater
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')

    # Create file handler
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(verbosity)
    file_handler.setFormatter(formatter)

    # Create stream handler for writing to console
    console_handler = logging.StreamHandler(sys.stdout)  
    console_handler.setLevel(verbosity)
    console_handler.setFormatter(formatter)

    # Add the file handler to the logger
    logger = logging.getLogger(name)
    logger.setLevel(verbosity)
    logger.addHandler(file_handler)
    logger.addHandler(console_handler)
    logger.propagate = False  # Prevent duplicate logs

    return logger