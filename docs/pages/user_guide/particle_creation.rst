Configuring the tokamak and creating a particle
===============================================

Tokamak configuration
---------------------

We first configure our tokamak by defining its major and minor radii (R, :math:`\alpha`),
its q-factor, magnetic and electric fields:

.. code-block:: python

    >>> R, a    = 12, 2
    >>> q       = gcm.qfactor.Hypergeometric(R, a)
    >>> Bfield  = gcm.bfield.LAR(i=0, g=1, B0=5)
    >>> Efield  = gcm.efield.Radial(R, a, q, Ea=75000, minimum=0.9, waist_width=50)

Now all the needed configuration objects are ready to be passed to whatever
particle we want to create!

Setting up initial conditions
-----------------------------

The necessary