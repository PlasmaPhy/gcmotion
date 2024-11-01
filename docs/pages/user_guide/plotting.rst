Plotting
========

Single-particle plotting functions must always take the :code:`cwp` as their first arguement. The rest of the arguements
are optional and have default values, so every plotting function can be called as such:

.. code-block:: python

    >>> gcm.time_evolution(cwp)

or

.. code-block:: python

    >>> gcm.time_evolution(cwp, percentage=10, units="NU")

If the code is run in a .py file, the plots will open in a new window in interactive mode.

If the code is run in a Jupyter notebook, then the result will be a static image. Its dpi can be set by importing 
the :code:`matplotlib` package and setting the default dpi as such:

.. code-block:: python

    >>> import matplotlib.pyplot as plt
    >>> plt.rcParams["figure.dpi"] = 300

Jupyter Notebooks can also open interactive windows, though its not recommended and the :code:`pyqt6` package is
required. This is done as such:

.. code-block::

    >>> %matplotlib qt
    >>> gcm.time_evolution(cwp)

The :code:`%matplotlib qt` command tells the notebook to open every plot as an interactive window. To change back to
static images, use :code:`%matplotlib inline`

