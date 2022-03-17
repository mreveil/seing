====================
Code Documentation
====================


Basic Classes
---------------

.. doxygenclass:: AtomicSystem
   :project: SEING
   :members:


.. doxygenclass:: Atom
   :project: SEING
   :members:


Neighbor Searching
--------------------

.. doxygenclass:: NeighborList
   :project: SEING
   :members:


Input Parsing
-----------------

.. doxygenstruct:: fingerprintProperties
   :project: SEING
   :members:


Fingerprint Generation Utilities
---------------------------------

.. doxygenclass:: FingerprintGenerator
   :project: SEING
   :members:

.. doxygenclass:: GenericLocalCalculator
   :project: SEING
   :members:


Fingerprint Calculators
------------------------

All fingerprint calculator classes have the same constructor signature and expose
a get_size and calculate_fingerprint functions to allow for easy switching between
different types.
The AGNI fingerprint calculator will be used here as an example.

.. doxygenclass:: AGNICalculator
   :project: SEING
   :members:


Utility Functions
--------------------

.. doxygenfunction:: cutoff_func
   :project: SEING
