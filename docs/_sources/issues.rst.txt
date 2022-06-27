.. MAGEMin documentation

.. image:: /figs/tips.png
   :width: 80
   :align: right

Quick tips!
===========

Minimum pressure
****************

* Minimum pressure for phase equilibrium calculation should be set to 0.01 kbar. A choosen pressure of 0.0 will lead to failed minimization in some cases. 

Reduced chemical system
***********************

* While :literal:`H2O` can be set to 0.0, other system components should be set to low value to ensure stable and efficient computation. 

* We recommend a minimum value of 0.01 :literal:`[mol%]` (0.0001 :literal:`[mol fraction]`)


|

.. image:: /figs/wip.png
   :width: 80
   :align: right

Known problems
==============

Reduced chemical systems
************************

* While :literal:`H2O` can be prescribed to 0.0 without impacting the efficiency of the minimization we found out that :literal:`Cr2O3`, :literal:`TiO2`, :literal:`O` etc. should not be set to 0.0 without strongly degrading the stability of the minimization. 

* Instead setting a small amount of 0.01 mol\% stabilizes the minimization and yields the expected results. 

* This also applies when an oxide has to be turned off e.g., during the computing liquid line of descent

* The current implementation of the solution phase models (EOS) from Holland et al. (2018) do not allow removing several oxides such as :literal:`MgO` and :literal:`CaO`. This would lead to ill-defined solution phases models such as for olivine and plagioclase.


Low temperature stabilitzation
******************************

* In a number of cases, mainly water under-saturated cases, and at temperature below 600-650°C, the current version of the solver fails to provide consistent minimizations. We are working on an algorithm upgrade in order to fix this issue. 

* While a working alternative approach is already being tested, we are further improving the performances before releasing it.


HP-LT melt prediction
*********************

* For water-saturated composition, we found that at pressure > 20 kbar and temperature < 650°C melt is predicted to be stable. 

* This is a known problem from the THERMOCALC developer community, and a fix should be provided in the following versions of the igenous thermodynamic dataset.