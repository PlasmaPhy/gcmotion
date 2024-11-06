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
    >>> logger.enable("gcmotion")
    >>> ...

"""

from loguru import logger


# Setup logger
logger.remove()  # Removes the default one which prints on sys.stderr


# Format template
fmt = "{time:HH:mm:ss:SSS} |{level: ^7}| {message}"

level = "DEBUG"

logger.add("log.log", format=fmt, level=level, mode="w")

logger.info(f"Logger added on {level} level.")
