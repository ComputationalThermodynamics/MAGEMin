.. MAGEMin documentation

.. image:: /figs/tips.png
   :width: 80
   :align: right

Quick tips!
===========

Compiler (Linux)
****************

* When using the C version of MAGEMin on linux, the :literal:`gcc` compiler yields slightly better performances (5-10% faster) over the :literal:`clang`.

Minimum pressure
****************

* Minimum pressure for phase equilibrium calculation should be set to 0.01 kbar. A choosen pressure of 0.0 will lead to failed minimization in some cases. 

Reduced chemical system
***********************

* While :literal:`H2O` can be set to 0.0, other system components should be set to low value to ensure stable and efficient computation. 

* We recommend a minimum value of 0.01 :literal:`[mol%]` (0.0001 :literal:`[mol fraction]`)


Set Oxygen content
******************

* Such like in THERMOCALC, the oxygen content :literal:`O` is set as :literal:`O [mol%]` = :literal:`Fe2O3 [mol%]`. 

* If the :literal:`Fe2O3`-content is not available you can turn off :literal:`O` by setting it to 0.01 :literal:`[mol%]`.

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


Low temperature stabilization
******************************

* In a number of cases using the Igneous database, mainly water under-saturated cases, and at temperature below 600-650°C, the PGE solver fails to provide consistent minimizations. 

* In nearly all of these cases, the Gibbs-Duhem constraint cannot be enforced by the PGE (which relies on being able to enforce it), and an alternative "legacy" has now been added. The "legacy" solver uses the approach presented by de Capitani & Brown (1987). This approach enforces the Gibbs-Duhem constraint only on the effective composition of the solution phase and not on each constitutive endmember of the solution phase such as in the PGE. When the PGE method fails the "legacy" solution can be seen as a relaxed Gibbs-Duhem constraint solution. 

* Note that in the cases of non PGE failure, the PGE and the "legacy" solvers yield identical results, which points out to thermodynamic database limitation outside calibration range.


HP-LT melt prediction
*********************

* For water-saturated composition, we found that at pressure > 20 kbar and temperature < 650°C melt is predicted to be stable. 

* This is a known problem from the THERMOCALC developer community, and a fix should be provided in a future version of the igneous thermodynamic dataset.