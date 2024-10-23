"""
Logger setup
------------

Sets up the global logger of the system.

Here is how to supress and re-enable the logger:

.. code-block:: python

    >>> from gcmotion.utils._logger_setup import logger, add_logger
    >>> ...
    >>> logger.disable("gcmotion")
    >>> function() # function that logs stuff but we want to disable
    >>> logger.enable()
    >>> ...

"""

from loguru import logger


# Setup logger
logger.remove()  # Removes the default one which prints on sys.stderr


# Format templates
# fmt = "{time:HH:mm:ss:SSS} | {function: ^25} |  {level: ^7} | {message}"
# fmt = "{time:HH:mm:ss:SSS} | {name: <18} |  {level: >6} | {message}"
fmt = "{time:HH:mm:ss:SSS} |  {level} | {message}"

level = "DEBUG"

logger.add("log.txt", level=level, mode="w")

logger.info(f"Logger added on {level} level.")
logger.info("---------------------------------")
