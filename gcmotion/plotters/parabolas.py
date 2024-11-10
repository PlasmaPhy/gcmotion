r"""
In the absence of electric field and in LAR configuration, it plots
the orbit type parabolas as described by White.

Can also plot the cwp's orbit type point.

Example
-------

.. code-block:: python

    gcm.parabolas(plot_point=True)

.. rubric:: Function:
    :heading-level: 4

"""

import matplotlib.pyplot as plt

from gcmotion.utils.logger_setup import logger

from gcmotion.classes.parabolas import Construct

from gcmotion.configuration.plot_parameters import orbit_point as config


def parabolas(cwp, plot_point: bool = True, plot_label: bool = True, **params):
    """Constructs and plots the orbit type parabolas.

    Returns early if there is no Electric field.

    :meta public:

    Parameters
    ----------
    cwp : :py:class:`~gcmotion.classes.particle.Particle`
        The current working particle.
    plot_point : bool
        Whether or not to plot cwp's orbit point. Defaults to True.
    plot_label : bool
        Whether or not to plot a label above particle.
    params : dict, optional
        Extra plotting parameters:

            * labels : bool
                Whether or not to plot a label above particle.
            * autoscale : bool
                Whether to autoscale the x axis in order to get
                all the points int view, or limit the plot to the
                parabolas' area.
    """
    # Unpack params
    _internal_call = params.pop("_internal_call", False)  # POP
    canvas = params.pop("canvas", None)  # POP!
    autoscale = params.get("autoscale", False)

    logger.info("Plotting orbit type Parabolas:")

    if cwp.has_efield or not cwp.Bfield.is_lar:
        string = "Electric field is present, or Magnetic field is not LAR. Orbit type parabolas do not work."
        print(string)
        logger.info("\t" + string)
        return

    logger.debug("Calling 'Construct' class")
    obj = Construct(cwp, canvas=canvas)
    canvas = obj.get_canvas()
    logger.info("--> Parabolas and Boundary plotted successfully.\n")

    fig, ax = canvas

    if plot_point:
        orbit_point(cwp, canvas=canvas)

    if autoscale:
        plt.autoscale(axis="x")
    # if not _internal_call:
    #     plt.ion()
    #     plt.show(block=True)


# ----------------------------------------------------------------


def orbit_point(cwp, **params):
    r"""Plots the particle point on the :math:`\mu-P_\zeta` (normalized) plane.

    :meta private:

    Parameters
    ----------

    canvas : tuple
        A 2-tuple containing the (fig, ax) object upon which to plot
        the point.
    plot_label : bool
        Whether or not to plot a label above particle.
    """
    if cwp.has_efield or not cwp.Bfield.is_lar:
        string = "Electric field is present, or Magnetic field is not LAR. Orbit type point cannote be plotted."
        print(string)
        logger.info("\t" + string)
        return

    # Unpack params
    _internal_call = params.pop("_internal_call", False)  # POP!
    canvas = params.pop("canvas", None)  # POP!
    plot_label = params.get("plot_label", True)

    if _internal_call:
        logger.disable("gcmotion")
    else:
        logger.enable("gcmotion")

    logger.info("Plotting orbit type point on parabolas plot...")

    # Get all needed attributes first
    orbit_x = cwp.orbit_x
    orbit_y = cwp.orbit_y
    t_or_p = cwp.t_or_p
    l_or_c = cwp.l_or_c

    fig, ax = canvas

    ax.scatter(orbit_x, orbit_y, **config["orbit_point_kw"])

    if plot_label:
        label = "  Particle " + f"({t_or_p[0]}-{l_or_c[0]})"
        ax.annotate(label, (orbit_x, orbit_y), color="b")
        logger.debug("\tPlotting particle's labels.")

    ax.set_xlabel(r"$P_\zeta/\psi_p$")
    logger.info("--> Plotted orbit type point successfully.\n")
