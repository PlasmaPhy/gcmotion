"""
Runs a 3d animation of the particle's orbit inside the tokamak.

.. caution::
    This feature is experimental.

Animation parameters include:

1] percentage (int)
    The percentage of the orbit to be animated, by default 100.

2] truescale (bool)
    Whether or not to plot the orbit in True scale, by default True.

3] min_step (float)
    The minimum stepsize to take per animation step, with respect to R.
    The higher the number the better the performance, but the more jagged
    the particle's trail. By default 0.01.

4] seconds (int)
    The total amount of seconds for the animation to take place. The animation
    rate is automatically adjusted internally to match the duration.

Example
-------

.. code-block:: python

    params = {"percentage": 40, "truescale": True, "min_step": 0.02, "seconds": 300}
    gcm.animation.run(cwp, params)

.. rubric:: Functions
    :heading-level: 4

"""

import multiprocessing as mp
import vpython as vp
import numpy as np
from scipy.signal import argrelextrema as ex
from tqdm import tqdm

from gcmotion.configuration.animation_parameters import animation_kw as config
from gcmotion.utils.canonical_to_toroidal import canonical_to_toroidal


spawn_ctx = mp.get_context("spawn")


def run(cwp, params: dict = {}):
    r"""
    Runs the script in a separate process. Necessary when running interactively.

    Parameters
    ----------

    cwp : particle object.
        The Current Working Particle.
    params : dict, optional
        The animation parameters, by default {}

    """
    proc = spawn_ctx.Process(target=animate, args=(cwp, params), daemon=True)
    proc.start()
    try:
        proc.join()
    finally:
        proc.kill()


def animate(cwp, params: dict = {}):
    r"""Sets up animation parameters and runs it.

    Parameters
    ----------
    cwp : particle object
        The Current Working Particle.
    params : dict, optional
        The animation parameter, by default {}

    """

    for key in config.keys():
        if key not in params:
            params[key] = config[key]

    percentage = params.get("percentage", 100)
    truescale = params.get("truescale", True)
    min_step = params.get("min_step", 0.01)
    seconds = params.get("seconds", 60)

    R, a, r_torus, theta_torus, z_torus = canonical_to_toroidal(
        cwp, percentage=percentage, truescale=truescale
    )
    # Cartesian (y and z are switched in vpython!)
    x = (R + r_torus * np.cos(theta_torus)) * np.cos(z_torus)
    z = (R + r_torus * np.cos(theta_torus)) * np.sin(z_torus)
    y = r_torus * np.sin(theta_torus)
    rate = 30
    running = 0

    print(f"Minimum step size:\t{np.around(min_step,5)}.")
    print(f"Animation duration:\t{seconds} seconds.")

    def dist_compress(min_step: None):

        n = len(x)  # number of steps

        step = np.sqrt((x[1:] - x[:-1]) ** 2 + (y[1:] - y[:-1]) ** 2 + (z[1:] - z[:-1]) ** 2)

        # Create arrays used for plotting
        plotx, ploty, plotz = [np.zeros(n) for _ in range(3)]
        plotx[0], ploty[0], plotz[0] = x[0], y[0], z[0]

        step_buffer = 0
        skipped = 0
        j = 0
        for i in range(n - 2):
            if step_buffer > min_step:
                plotx[j], ploty[j], plotz[j] = x[i], y[i], z[i]
                j += 1
                step_buffer = 0
                skipped = 0
            else:
                step_buffer += step[i]
                skipped += 1

        plotx = np.trim_zeros(plotx, trim="b")
        ploty = np.trim_zeros(ploty, trim="b")
        plotz = np.trim_zeros(plotz, trim="b")

        compression = int((x.shape[0] - plotx.shape[0]) / x.shape[0] * 100)
        print(f"\nCompression level: -{compression}%\n")

        return plotx, ploty, plotz

    def adjust_rate(x, seconds=seconds):
        rate = int(x.shape[0] / seconds)
        print(f"Rate = {rate}/s")
        return rate

    def _setup_scene(truescale: bool = True):

        # Canvas
        width = 1920
        height = 850
        scene = vp.canvas(width=width, height=height, background=vp.color.white)
        scene.userpan = False
        scene.center = vp.vector(0, -0.55 * R, 0)

        # Vertical axis
        height = 1 * (2 * R)
        pos = [vp.vector(0, -height / 2, 0), +vp.vector(0, height / 2, 0)]
        vaxis = vp.curve(pos=pos, color=eval("vp.color." + config["vaxis_color"]), radius=0.004 * R)

        # Torus walls
        shape = vp.shapes.circle(radius=float(a), np=60)
        path = vp.paths.circle(radius=float(R), np=60)
        torus = vp.extrusion(
            pos=vp.vector(0, 0, 0),
            shape=shape,
            path=path,
            color=eval("vp.color." + config["torus_color"]),
            opacity=0.4,
        )

        return scene, vaxis, torus

    def _flux_surface():

        # Get theta of 1 period
        if cwp.t_or_p == "Trapped":
            span = ex(cwp.theta, np.greater)[0][:2]
            theta = cwp.theta[span[0] : span[1]]

        if cwp.t_or_p == "Passing":
            condition = (np.abs(cwp.theta) > cwp.theta0) & (
                np.abs(cwp.theta) < cwp.theta0 + 2 * np.pi
            )
            theta = cwp.theta[condition]

        r = r_torus[: theta.shape[0]]
        xflux = r * np.cos(theta)
        zflux = r * np.sin(theta)
        zero = np.zeros(theta.shape[0])
        points = np.vstack((xflux, zflux, zero)).T
        vectors = []

        for i in range(len(points)):
            vectors.append(vp.vector(points[i, 0], points[i, 1], points[i, 2]))

        shape = np.vstack((xflux, zflux)).T.tolist()
        shape.append(shape[0])
        path = vp.paths.circle(radius=float(R), np=60)
        flux_surface = vp.extrusion(
            pos=vp.vector(0, 0, 0),
            shape=shape,
            path=path,
            color=eval("vp.color." + config["flux_surface_color"]),
            opacity=config["flux_surface_opacity"],
        )

        return flux_surface

    def _setup_particle():

        # Particle
        pos = vp.vector(x[0], y[0], z[0])
        p = vp.sphere(
            pos=pos,
            radius=a / 20,
            color=eval("vp.color." + config["particle_color"]),
            make_trail=True,
            trail_radius=R / 1000,
            interval=1,
        )

        return p

    class buttons:
        def runAnimation(foo):
            nonlocal running
            running = not running
            if running:
                foo.text = "Pause"
            else:
                foo.text = "Run"

        def ptrail(foo):
            nonlocal p
            p.make_trail = foo.checked
            if not foo.checked:
                p.clear_trail()

        def rate_slider(foo):
            nonlocal rate
            rate = foo.value

        def restart(foo):
            nonlocal i, p, pbar
            i = 0
            p.clear_trail()
            p.visible = False
            p = _setup_particle()
            pbar.reset()

        def setup():
            vp.radio(bind=buttons.runAnimation, text="Run")
            vp.checkbox(bind=buttons.ptrail, text="Show Trail", checked=True)
            vp.slider(bind=buttons.rate_slider, text="Rate", min=30, max=10000, step=5, value=rate)
            vp.button(bind=buttons.restart, text="Restart")

    x, y, z = dist_compress(min_step)
    rate = adjust_rate(x, seconds=seconds)
    scene, vaxis, torus = _setup_scene()
    # flux = _flux_surface()
    p = _setup_particle()
    buttons.setup()

    i = 0
    with tqdm(total=len(x), leave=True) as pbar:
        while True:
            vp.rate(rate)

            if running:
                p.pos = vp.vector(x[i], y[i], z[i])
                pbar.update(1)
                i += 1

            if i == len(x):
                pbar.refresh()
                scene.waitfor("click")
                buttons.restart("foo")

        pass
