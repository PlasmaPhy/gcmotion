===================
Creating a Particle
===================

Setting up tokamak configuration
--------------------------------

First we need to configure our tokamak:

.. code-block:: python
    
    >>> R, a = 6, 2  # Major/Minor Radius in [m]
    >>> q = gcm.qfactor.Hypergeometric(R, a)
    >>> Bfield = gcm.bfield.LAR(i=0, g=1, B0=5)
    >>> Efield = gcm.efield.Radial(R, a, q, Ea=75000, minimum=0.9, waist_width=50)


#. In the first line, we set the tokamak's major and minor radii in [m]. Note that 
   first thing in particle initialization is normalizing the lengths to the major 
   radius R. Therefore, during calculations, the quantity :math:`r_{wall} = \alpha/R` 
   is used to derive quantities like :math:`\psi_{wall}`, :math:`$B_{wall}`, etc.

#. In the second line, we set up the q factor profile. In this example we use a 
   hypergeometric profile, but other configurations are :ref:`availiable <qfactors>`. To choose a 
   different profile, we simple set :code:`q = gcm.qfactor.<q factor>(parameters)`. 
   Different q profiles require different parameters. We can also easily create 
   new profiles with however many parameters we want. For more info, see the 
   :ref:`documentation <about_qfactors>`.

#. In the third line we set up the Magnetic field inside the tokamak. Here we use 
   the class Large Aspect Ratio configuration, which is the only :ref:`availiable <bfields>` at the 
   moment. We can easily create new ones though, see :ref:`documentation <about_bfields>`. 
   Note that during calculations, all values are normalized to the Magnetic field 
   axis strength $B_0$. $B_0$ itself never shows up in calculations, but only when 
   we want to convert the results to lab units.

#. In the forth line, we set up the Electric field. Similarly to the q factor, 
   we can choose from the :ref:`availiable <efields>` configurations, or create new ones, see 
   :ref:`documentation <about_efields>`.

.. note::

    The `q`, `Bfield` and `Efield` are objects created from their respective classes, 
    and the solver doesn't care about what they do and how, as long as they support 
    the required querry methods described in the documentation. For example, whenever 
    the solver wants to get the value of the magnetic field stength, it just calls
    :code: `Bfield.B(psi, theta)`, which returns the numerical vaule of the field 
    strength. That way, we never have to tinker with the solver itself whenever 
    we change configuration, and we can also use whichever method we want for 
    returning the required values (e.g. analytical calculation, lookup tables etc.)

Setting up particle's initial conditions
----------------------------------------

.. code-block:: python
    
    >>> species = "p"
    >>> mu = 1e-5
    >>> theta0 = 0
    >>> psi0 = 0.5  # times psi_wall
    >>> zeta0 = np.pi
    >>> Pzeta0 = -0.027
    >>> t_eval = np.linspace(0, 1000000, 100000)  # t0, tf, steps


#. First we must setup the particle's species. This is needed since all calculations 
   are normalized to the protons mass, so if we want to study other particles, we 
   need to account for that. The availiable particle species are "p" (proton), "e" (electron),
   "D" (deuterium), "T" (Tritium), "He3" and "He4".

#. The rest are pretty straight forward, with the exception of :math:`\psi_0`. Since 
   :math:`\psi_{wall}` is directly dependant on the tokamak's :math:`R` and :math:`\alpha` 
   through :math:`\psi_{wall} = \dfrac{r_{wall}^2}{2} = \dfrac{1}{2}\bigg(\dfrac{\alpha}{R}\bigg)^2`, 
   it is usually a weird number (in our example, it is 0.05). Therefore, we set up $\psi_0$ **with respect to $\psi_{wall}$**, and the actual $\psi_0$ is calculated upon particle initialization.
  
#. :code:`teval` contains the times for which we want to now the orbit. As with all solvers, 
   the stepsize of :code:`teval` is not the actual step size of the solver, since it uses 
   an adaptive step size.

Initializing particle
---------------------

.. code-block:: python

   >>> tokamak = {"R": R, "a": a, "q": q, "Bfield": Bfield, "Efield": Efield}
   >>> init_cond = {"theta0": theta0, "psi0": psi0, "zeta0": zeta0, "Pzeta0": Pzeta0}

   >>> particle1 = gcm.Particle(tokamak, t_eval, init_cond, mu, species)
   >>> cwp = particle1

Here we instantiate a single particle from the :py:class:`~gcmotion.classes.particle.Particle` 
class. We then set it as :code:`cwp` (current working particle), to avoid confusion.

As we can see, we simply create it by passing the parameters we created above. We can 
create multiple particles that way, which are all stored in memory until we restart 
Python, and we can switch between them by simply setting them as :code:`cwp`.

Now :code:`cwp` represents a fully-fledged particle. It stores the parameters we gave 
it, as well as all later derived and calculated quantities **inside** it 
(called its *attributes*). We can easily access them and print them as such:

.. code-block:: python

   >>> print(cwp.theta0)
   0
   >>> print(cwp.t_eval)
   [0.00000e+00 1.00001e+01 2.00002e+01 ... 9.99980e+05 9.99990e+05
   1.00000e+06]
   >>> print(cwp.Bfield.B(0.5, np.pi)) # B field stength at r=0.5, θ=π
   1.5

Calculating particle's obrit
----------------------------

One of the first things we wanna do is calculate the particles full obrit.

We simply run:

.. code-block:: python

   >>> cwp.run()

The :code:`run()` method goes through some specific steps in a specific order. 

#. It calculates the convertion factors from [NU] (Normalized Units) to lab units. This 
   is needed first for converting the Electric field strength from [V/m] to [NU].

#. It calculates the particle's energy. The energy can be directly calculated from the 
   initial conditions alone, since it is a constant of motion.

#. If and only if we use LAR Magnetic field **and** no Electric field, it calculates 
   the particle's orbit type (Trapped/Passing and Lost/Confined) as described by White.

#. Calculates the actual orbit and stores the results inside itself.

#. Prints the particle's calculated attributes.

The :code:`run()` method also accepts 3 optional parameters:

#. :code:`orbit: bool`: Whether or not to actually calculate the orbit or skip it. This 
   is used in e.g. studying the results that are calculated from the initial conditions 
   only, especially in the many-particles case. Defaults to True.

#. :code:`info:bool`:  Whether or not to print the information message after the orbit 
   is calculated. Defaults to True.

#. :code:`events:list`: a list with events to provide to the solver's *event locator*. 
   Defaults to an empty list.

For example, if we want to calculate the particles orbit up until :math:`\theta` returns
to its initial value 10 times, but not print the info message, we can run:

.. code-block:: python
   
   >>> events = [gcm.events.when_theta(theta0, 10)]
   >>> cwp.run(events=events, info=False)s

.. note::
   
   The solver solves the differential equations with respect to the dynamical variables 
   :math:`\theta, \psi, \zeta` and :math:`\rho_{||}`. The quantities 
   :math:`\psi_p, P_\theta` and :math:`P_\zeta`` are calculated afterwards.
