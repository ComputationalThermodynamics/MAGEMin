.. MAGEMin documentation


Single point calculation Tutorial
=================================

First open a terminal and execute Julia, then type :literal:`using MAGEMin_C` to load the package.

|

Initialize database 
*******************
.. code-block:: shell

   gv, z_b, DB, splx_data      = init_MAGEMin();

This initiatizes the global variables and the Database.

Set P-T-(pressure temperature)
**********************************************************
.. code-block:: shell

   P           = 8.
   T           = 800.
   sys_in      = "mol"     # wt or mol, default is mol
   test        = 0;
   gv          = use_predefined_bulk_rock(gv, test)

:literal:`get_bulk_rock` retrieves the saved bulk-rock composition 0, which corresponds to KLB-1 peridotite. 

|

Set bulk-rock composition
**********************************************************

Use a pre-defined bulk-rock "test" composition

.. code-block:: shell

   test        = 0;
   gv          = use_predefined_bulk_rock(gv, test)

|

Set the level of verbose :literal:`[-1,0,1]`
********************************************
.. code-block:: shell   

   gv.verbose  = -1    # switch off any verbose

:literal:`-1`, none; :literal:`0`, stable phase assemblage; :literal:`1`, full verbose. By default :literal:`gv.verbose` = 0.

|

Call optimization routine for given P-T-X
*****************************************
.. code-block:: shell   

   gv.verbose  = -1    # switch off any verbose
   out         = point_wise_minimization(P,T, gv, z_b, DB, splx_data, sys_in)

|

Display minimized point
************************

.. code-block:: shell   

   @show out

The command :literal:`@show` allows you to display the minimized point:

.. image:: /figs/julia_out.png
   :width: 420
   :align: center

|

Access output structure
************************

In Julia all the informations stored in the output structure ``stb_systems`` can be listed by typing ``out.`` and hitting the tab key twice :kbd:`tab` + :kbd:`tab`:

.. code-block:: shell   

   out.

which displays the content of structure ``out``:

.. image:: /figs/julia_out_struct.png
   :width: 640
   :align: center

|

The displayed informations are part of the ``C`` output structure ``stb_systems``, and can be accessed individually (e.g., ``out.Gamma``) or displayed all at once using 

.. code-block:: shell   

   print_info(out)

The full description of what contains the output structure is given in the CookBook: :doc:`/ckbk/out_struct`.

|

Provide custom bulk rock composition
====================================

Bulk-rock compositions can be provided either by directly using the system unit of MAGEMin or by using other system unit and converting them.

MAGEMin system unit
*******************

To define and use your own bulk rock composition you can provide it directly in the format used in MAGEMin: :literal:`[mol,wt]` fraction using the following chemical system i.e.

+-------+--------+-------+-------+-------+------+------+------+------+-------+-------+
| SiO2  |  Al2O3 |  CaO  | MgO   | FeOt  |  K2O | Na2O | TiO2 |   O  | Cr2O3 |  H2O  |
+-------+--------+-------+-------+-------+------+------+------+------+-------+-------+

For instance:

.. code-block:: shell

   bulk_rock   = [0.38451319035870185, 0.017740308257833806, 0.028208688355924924, 0.5050993397328966, 0.0587947378409965, 9.988912307338855e-5, 0.0024972280768347137, 0.0009988912307338856, 0.0009589355815045301, 0.0010887914414999351, 0.0]
   gv          = define_bulk_rock(gv, bulk_rock)


Variable system unit
********************

However, measured bulk-rock compositions usually provide :literal:`FeO` and :literal:`Fe2O3` instead of :literal:`FeOt` and :literal:`O`. Therefore, we provide an additional routine that convert bulk-rock composition into the right MAGEMin format:

.. code-block:: shell

   bulk_in_ox = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "Fe2O3"; "K2O"; "Na2O"; "TiO2"; "Cr2O3"; "H2O"]
   bulk_in    = [48.43; 15.19; 11.57; 10.13; 6.65; 1.64; 0.59; 1.87; 0.68; 0.0; 3.0]
   sys_in     = "wt"

   bulk_rock   = convertBulk4MAGEMin(bulk_in,bulk_in_ox,sys_in)
   gv          = define_bulk_rock(gv, bulk_rock)

where :literal:`bulk_in_ox` is a :literal:`Vector(String)` containing the oxide names and :literal:`bulk_in` is a :literal:`Vector(Float)` of the :literal:`[mol,wt]` fraction of the bulk-rock composition. The function 

.. code-block:: shell

   bulk_rock   = convertBulk4MAGEMin(bulk_in,bulk_in_ox,sys_in)` 

converts

 :literal:`SiO2`, ..., :literal:`FeO` and :literal:`Fe2O3` in system unit :literal:`[mol,wt]`
 
to:
 
 :literal:`SiO2`, ..., :literal:`FeOt` and :literal:`O` in system unit :literal:`[mol]`.

Note that if the provided bulk-rock composition includes more oxides than supported, they will be ignored and the composition will be renormalized accordingly. Moreover, if both :literal:`Fe2O3` and :literal:`O` are provided, :literal:`O` will be recalculated as function of :literal:`Fe2O3`. Thus, if you want to prescribe a different :literal:`O` content, do not define :literal:`Fe2O3`!

A full Julia script demonstrating how to use this function is provided below:

.. code-block:: shell

   #load MAGEMin
   using MAGEMin_C 

   #initialize
   gv, z_b, DB, splx_data      = init_MAGEMin();         

   # provide bulk-rock composition
   bulk_in_ox = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "Fe2O3"; "K2O"; "Na2O"; "TiO2"; "Cr2O3"; "H2O"]
   bulk_in    = [48.43; 15.19; 11.57; 10.13; 6.65; 1.64; 0.59; 1.87; 0.68; 0.0; 3.0]
   sys_in     = "wt"

   # convert bulk rock
   bulk_rock  = convertBulk4MAGEMin(bulk_in,bulk_in_ox,sys_in)

   # send bulk-rock to MAGEMin
   gv         = define_bulk_rock(gv, bulk_rock)

   # provide pressure and temperature conditions
   P          = 10.0
   T          = 1100.0

   # switch off any verbose
   gv.verbose = -1    

   # perform minimization    
   out        = point_wise_minimization(P,T, gv, z_b, DB, splx_data, sys_in)

   # print output
   print_info(out)

   # free memory
   finalize_MAGEMin(gv,DB)
