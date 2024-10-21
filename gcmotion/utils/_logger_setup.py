"""
Logger setup
------------

Sets up the global logger of the system.
"""

from loguru import logger

# Setup logger
logger.remove()

# Format templates
# fmt = "{time:HH:mm:ss:SSS} | {function: <20} |  {level: ^7} | {message}"
# fmt = "{time:HH:mm:ss:SSS} | {name: <18} |  {level: >6} | {message}"

fmt = "{time:HH:mm:ss:SSS} |  {level: ^7} | {message}"

level = "DEBUG"

logger.add("log.txt", delay=True, level=level, format=fmt, mode="w")


def add_logger():  # For re-enabling logger
    logger.add("log.txt", delay=True, level=level, format=fmt, mode="a")


add_logger()

logger.info(f"Logger added on {level} level.")
logger.info("---------------------------------")
