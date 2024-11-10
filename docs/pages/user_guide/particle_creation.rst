.. _user_guide_particle_creation:

===================
Creating a Particle
===================

Setting up the Quantity Constructor
-----------------------------------

The quantity constructor is an easy-to-use but powerful tool provided by the 
`pint <https://pint.readthedocs.io/en/stable/>`_ package. It offers the ability
create our own units and unit systems, and makes unit conversions a breeze, thus
saving us the headache of manually converting them every time it is needed. The units
are defined once in the :py:mod:`gcmotion.utils.setup_pint` module and are henceforth
globally available. A quick read-through the 
`pint user guide <https://pint.readthedocs.io/en/stable/user/defining-quantities.html>`_
is recommended.

From our theory we know that the normalization depends on the Magnetic field strength on
the magnetic axis :math:`B_0`, the tokamak's major radius :math:`R` and the particle's 
species, through its atomic and mass number. To make our lives easier when defining 
initial conditions for :math:`\psi_0`, we also define :math:`[\psi_{wall}]` as a 
unit of magnetic flux itself, so we can define :math:`\psi_0` relative to that.
Thus, the tokamak's minor radius :math:`\alpha` is also needed in defining our new
units system:

.. code-block:: python

   >>> Rnum = 1.65
   >>> anum = 0.5
   >>> B0num = 1
   >>> species = "p"
   >>> ureg, Q = gcm.setup_pint(R=Rnum, a=anum, B0=B0num, species=species)

The values passed to ``setup_pint()`` must be purely numerical, and **in [m] and [T]**
accordingly. The object we are interested in is Q. Through Q we can define any *Quantity* 
we want, in both S.I. and normalized units. ``ureg`` (UnitRegistry) is just a (very big)
dictionary containing all the defined units pint offers. We don't really needed but it is nice
to have around.

.. tip::
   *Quantity* objects are floats (or numpy arrays) with their units attached to them. All operations
   between Quantities work the same, as long as their dimensionality allows it. For example, pint won't
   allow us to convert from [Joule] to [kilograms] or add 1[meter] and 5[Tesla] (why would we anyway?).

Configuring the tokamak profile
--------------------------------------

First we need to configure our tokamak. We create all the Quantities needed for the q-factor,
magnetic and electric field:

.. code-block:: python

   >>> R = Q(Rnum, "meters")
   >>> a = Q(anum, "meters")
   >>> B0 = Q(B0num, "Tesla")
   >>> i = Q(0, "NUPlasma_current")
   >>> g = Q(1, "NUPlasma_current")
   >>> Ea = Q(73500, "Volts/meter")

The prefix "NU" means that the quantity is defined in normalized units.

.. tip::
   We are not restricted to "meters" or "Tesla". We might as well define :math:`R` in kilometers
   and :math:`B0` in Webers per square mile. The correct SI values are calculated internally.
   This applies to every Quantity used to create a particle.

Now we are ready to configure the tokamak profile. This is done by creating a dictionary containing
all the necessary parameters and objects:

.. code-block:: python

   >>> tokamak = {
   >>>     "R": R,
   >>>     "a": a,
   >>>     "B0": B0,
   >>>     "qfactor": gcm.qfactor.Hypergeometric(a, B0, q0=1.1, q_wall=3.8, n=2),
   >>>     "bfield": gcm.bfield.LAR(B0=B0, i=i, g=g),
   >>>     "efield": gcm.efield.Radial(a, Ea, B0, peak=0.98, rw=1 / 50),
   >>> } 

#. The first 3 entries are simply the Quantities we defined above.

#. In the 4th entry, we set up the q factor profile. In this example we use a 
   hypergeometric profile, but other configurations are :ref:`available <qfactors>`. To choose a 
   different profile, we simple set :code:`qfactor = gcm.qfactor.<qfactor>(parameters)`. 
   Different q profiles require different parameters. We can also easily create 
   new profiles with however many parameters we want. For more info, see the 
   :ref:`documentation <about_qfactors>`.

#. In the 5th entry we set up the Magnetic field inside the tokamak. Here we use 
   the Large Aspect Ratio configuration, which is the only :ref:`available <bfields>` at the 
   moment. We can easily create new ones though, see :ref:`documentation <about_bfields>`. 

#. In the forth line, we set up the Electric field. Similarly to the q factor, 
   we can choose from the :ref:`available <efields>` configurations, or create new ones, see 
   :ref:`documentation <about_efields>`.

.. note::

    The `qfactor`, `bfield` and `efield` are objects created from their respective classes, 
    and the solver doesn't care about what they do and how, as long as they support 
    the required querry methods described in the documentation. For example, whenever 
    the solver wants to get the value of the magnetic field stength it just calls
    :code: `Bfield.bigNU(psi, theta)`, which returns the numerical vaule of the field 
    strength. That way, we never have to tinker with the solver itself whenever 
    we change configuration, and we can also use whichever method we want for 
    returning the required values (e.g. analytical calculation, lookup tables etc.)

Setting up particle's initial conditions and parameters
-------------------------------------------------------

Particle initial and conditions are Quantities as well:

.. code-block:: python
    
   >>> parameters = {
   >>>     "species": species,
   >>>     "mu/muB": Q(5, "keV"),
   >>>     "theta0": 0,
   >>>     "zeta0": 0,
   >>>     "psi0": Q(0.78, "psi_wall"),
   >>>     "Pzeta0": Q(-0.0272, "NUMagnetic_flux"),
   >>>     "t_eval": Q(np.linspace(0, 1e-3, 100000), "seconds"),
   >>> }


#. The first entry is  the particle's species. The available particle species are 
   "p" (proton), "e" (electron), "D" (deuterium), "T" (Tritium), "He3" and "He4". 
   The species string is case-insensitive.

#. The ``"mu/muB"`` entry is interesting. Sometimes we need to define the particle's
   magnetic moment :math:`\mu`, or the product :math:`\mu B_{initial}`. If the Quantity
   has *dimensionality* of [magnetic moment] = :math:`[current]\cdot[length]^2`, the parameter is parsed
   as the particle's :math:`\mu`. If it has *dimensionality* of [energy] = :math:`[mass]\cdot[length]^2/[time]^2`, 
   the parameter is parsed as the particle's :math:`\mu B_{initial}`.

#. The ``"theta0"`` and ``"zeta0"`` entries can be either purely numerical, or Quantities with
   units of [radians].

#. ``"psi0"`` can be defined in many ways. It can either be defined in units of "Magnetic_flux" 
   (:math:`[Tesla\cdot meters^2]`), "NUMagnetic_flux" (:math:`[NUTesla\cdot NUmeters^2]`),
   "psi_wall" (:math:`B_0 a^2/2` [Magnetic_flux]), or "NUpsi_wall" 
   (:math:`(a/R)^2/2` [NUMagnetic_flux]).
  
#. :code:`teval` contains the times for which we want to now the orbit. We can define it in any unit of time
   we want, such as "miliseconds", or "NUseconds" for normalized time units. As with all solvers, 
   the stepsize of :code:`teval` is not the actual step size of the solver, since it uses 
   an adaptive step size.

Initializing particle
---------------------

Now, particle creation is as simple as:

.. code-block:: python

   >>> cwp = gcm.Particle(tokamak, parameters)

It is handy to symbolize the particle as ``cwp`` (current working particle). We can create as 
many different particles as we want, and they will stay in memory until we close the session.

Now :code:`cwp` represents a fully-fledged particle. It stores the parameters we gave 
it, as well as all later derived and calculated quantities **inside** it 
(called its *attributes*). We can easily access them and print them as such:

.. code-block:: python

   >>> cwp["psi0"]
   0.0975 Tm^2
   >>> cwp["psi0NU"] # psi0 in [NUMagnetic_flux]
   0.0358127 NUmf
   >>> cwp["muB"], cwp["mu"], cwp["muNU"]
   5 keV
   6.82714 keV / T
   2.61793e-05 NUmu
   >>> cwp["t_eval"]
   [0 1.00001e-08 2.00002e-08 ... 0.000777408 0.000777418 0.000777428] s
   >>> cwp["bfield"]
   LAR: B0=1 T, I=0 Tm, g=1.65 Tm.

Calculating particle's obrit
----------------------------

One of the first things we wanna do is calculate the particles full obrit.

We simply run:

.. code-block:: python

   >>> cwp.run()

Once the calculation is complete, we should see an info message in our screen.

The :code:`run()` method also accepts 3 optional parameters:

#. :code:`orbit: bool`: Whether or not to actually calculate the orbit or skip it. This 
   is used in e.g. studying the results that are calculated from the initial conditions 
   only, especially in the many-particles case. Defaults to True.

#. :code:`info:bool`:  Whether or not to print the information message after the orbit 
   is calculated. Defaults to True.

#. :code:`events:list`: a list with :py:mod:`~gcmotion.scripts.events` to provide to the solver's *event locator*. 
   Defaults to an empty list.

For example, if we want to calculate the particles orbit up until :math:`\theta` returns
to its initial value 10 times, but not print the info message, we can run:

.. code-block:: python
   
   >>> events = [gcm.events.when_theta(parameters["theta0"], 10)]
   >>> cwp.run(events=events)

.. note::
   
   The solver solves the differential equations with respect to the dynamical variables 
   :math:`\theta, \psi, \zeta` and :math:`\rho_{||}`. The quantities 
   :math:`\psi_p, P_\theta` and :math:`P_\zeta`` are calculated afterwards.
