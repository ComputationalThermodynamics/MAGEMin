.. MAGEMin documentation

.. image:: /figs/eq.png
   :width: 256
   :align: right

Analytical formulation		  
======================

|

.. image:: /figs/simplex.png
   :width: 256
   :align: right

Generate pseudocompounds
========================

To further improve performance, pseudocompounds of the solution phases are pregenerated using a *python* script. The goal here, is to get a set of points uniformely covering (and restricted) to the feasible solution phase space.

This is achieved in two stage. First, the solution phase space is explored within allowed compositional variable range to generate potential pseudocompound candidate. Then, only pseudocompounds satisfying the site fraction constraint ( > 0.0) are saved to a list.


