.. MAGEMin documentation

Introduction
============

Finding the most stable phase assemblage boils down to a constrained optimization problem that has to ensure three main points:

- The objective function to minimize is the weighted sum of the Gibbs free energy of the solution and pure phases. 
- The equality constraint is the weighted sum of the phase composition that must be equal to the prescribed bulk rock composition.
- Inequality constraints are expressed in complex solution phases as site fractions and represent "forbidden" compositional space.

Gibbs free energy
=================

The Gibbs energy of a multi-component multiphase system is given by the weighted summation of the chemical potentials of each phase:

.. math::
   G_{sys} = \sum_{\lambda=1}^{\Lambda} n_{\lambda} \sum_{i=1}^{N_{\lambda}} \mu_{i}(\lambda) \chi_{i}(\lambda) + \sum_{\omega=1}^{\Omega} \mu_{\omega} n_{\omega}

where  

.. math:: 
   n_{\lambda} 
   
is are the molar fraction of the solution phase.

.. math:: 
   n_{\omega} 
   
is are the molar fraction of the pure phase.

The Gibbs energy of the system can alternatively be represented with respect to the chemical potentials [J.mol-1] and the number of moles [mol] of the system for which there are C elements in the system.

.. math:: 
   G_{sys} = \sum_{j=1}^{C} \Gamma_{j} b_{j}

The chemical potential of a phase is either a constant for a condensed (pure) phase:

.. math:: 
   \mu_{i} = G_{i}^{0}

or a function for a phase within a solution:

.. math:: 
   \mu_{i} = G_{i}^{0} + RTln(a_{i}) + G_{i}^{ex}

where

.. math::
   a_{i}

is the thermochemical activity is related to the mole fraction and the activity coefficient by

.. math::
   a_{i} = x_{i} \gamma_{i}


For the case of ideal mixing, the activity coefficient is unity. The mixing of a species dissolved in a condensed phase, however rarely behaves ideally
and is typically a function of both temperature and composition.


Global minimum and G-hyperplane
===============================

.. image:: /figs/Levelling_sketch_Manual.png
   :align: center
   
