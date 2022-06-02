.. MAGEMin documentation

Julia interface to MAGEMin
===========================

.. image:: /figs/julia.png
   :width: 120
   :align: right

The julia package :guilabel:`MAGEMin_C.jll` is a Julia wrapper for :guilabel:`MAGEMin`.
The interface is well-fitted for geodynamic coupling as the ``C`` functions can be directly called from Julia. 
Moreover, the inferface allows to get rid of file I/O (i.e., no slow disk access).

Installation
============

.. toctree::
   :maxdepth: 2
   
   julia/Installation


Use MAGEMin_C
=============

.. toctree::
   :maxdepth: 2
   
   julia/Run
 
  