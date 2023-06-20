.. MAGEMin documentation

.. image:: /figs/julia.png
   :width: 100
   :align: right

Julia interface to MAGEMin
===========================

The julia package :guilabel:`MAGEMin_C.jll` is a Julia wrapper for :guilabel:`MAGEMin`.
The interface is well-fitted for geodynamic coupling as the ``C`` functions can be directly called from Julia. 
Moreover, the inferface allows to get rid of file I/O (i.e., no slow disk access).


.. image:: /figs/install.png
   :width: 80
   :align: right
  
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
 
  