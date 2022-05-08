.. MAGEMin documentation

Single point calculation
=========================

First open a terminal and execute Julia, then type :literal:`using MAGEMin_C` to load the package.


Initialize database 
*******************
.. code-block:: shell

   gv, DB = init_MAGEMin();

This initiatize the global variables and the Database.

Set P-T-X (pressure temperature and bulk rock composition)
**********************************************************
.. code-block:: shell

   P           = 8.
   T           = 800.
   test        = 0;
   bulk_rock   = get_bulk_rock(gv, test)

:literal:`get_bulk_rock` retrieve the saved bulk-rock composition 0, which corresponds to KLB-1 peridotite. 
This is where you can pass your own bulk-rock and P-T conditions.

Set the level of verbose :literal:`[-1,0,1]`
********************************************
.. code-block:: shell   

   gv.verbose  = -1    # switch off any verbose

:literal:`-1`, none; :literal:`0`, stable phase assemblage; :literal:`1`, full verbose. By default :literal:`gv.verbose` = 0.


Call optimization routine for given P-T-X
*****************************************
.. code-block:: shell   

   gv.verbose  = -1    # switch off any verbose
   out         = point_wise_minimization(P,T, bulk_rock, gv, DB);

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

The displayed informations are part of the ``C`` output structure ``stb_systems``, and can be access individually (e.g., ``out.Gamma``) or displayed all at once using 

.. code-block:: shell   

   print_info(out)

The full description of what contains the output structure is given in the CookBook: :doc:`/ckbk/out_struct`.
