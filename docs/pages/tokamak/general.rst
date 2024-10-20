Tokamak Configuration
=====================

How q-factors and magnetic/electric
fields are implemented in GCMotion. 

.. toctree::
   :maxdepth: 1

   module_qfactor
   module_bfield
   module_efield

To use a q-factor/magnetic/electric field, simply construct an object:

.. code-block:: python

   R = 12
   a = 2
   q = gcm.qfactor.Hypergeometric(R, a)
   Bfield = gcm.bfield.LAR(i=0, g=1, B0=5)
   Efield = gcm.efield.Radial(R, a, q, Ea=75000, minimum=0.9, waist_width=50)

and pass it to the particle you want to create.


Currently supported configurations:
-----------------------------------

.. toctree::
   :maxdepth: 2

   qfactors
   bfields
   efields