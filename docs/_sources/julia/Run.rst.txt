.. MAGEMin documentation


Serial point calculation Tutorial
=================================

First open a terminal and execute Julia, then type :literal:`using MAGEMin_C` to load the package.

|

Initialize database 
*******************
.. code-block:: shell

   # Initialize database  - new way
   data        =   Initialize_MAGEMin("ig", verbose=true); # database: ig, igneous (Holland et al., 2018); mp, metapelite (White et al 2014b); mb, metabasite (Green et al.,2016); um, ultramafic (Evans & Frost 2021)

This initiatizes the global variables and the Database.

|

Set P-T-(pressure temperature)
**********************************************************
.. code-block:: shell

   P           =   8.0
   T           =   800.0

:literal:`use_predefined_bulk_rock` retrieves the saved bulk-rock composition 0 from database ig, which corresponds to KLB-1 peridotite. 

|

Set bulk-rock composition
**********************************************************

Use a pre-defined bulk-rock "test" composition

.. code-block:: shell
   test        =   0         #KLB1
   data        =   use_predefined_bulk_rock(data, test);

or a custom bulk-rock composition:

.. code-block:: shell

    Xoxides = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "Fe2O3"; "K2O"; "Na2O"; "TiO2"; "Cr2O3"; "H2O"];
    X       = [48.43; 15.19; 11.57; 10.13; 6.65; 1.64; 0.59; 1.87; 0.68; 0.0; 3.0];
    sys_in  = "wt"

Note that the system unit :literal:`[mol,wt]` has to be provided here.

where :literal:`Xoxides` is a :literal:`Vector(String)` containing the oxide names and :literal:`X` is a :literal:`Vector(Float)` of the :literal:`[mol,wt]` fraction of the bulk-rock composition.
The function converts

 :literal:`SiO2`, ..., :literal:`FeO` and :literal:`Fe2O3` in system unit :literal:`[mol,wt]`
 
to:
 
 :literal:`SiO2`, ..., :literal:`FeOt` and :literal:`O` in system unit :literal:`[mol]`.

Note that if the provided bulk-rock composition includes more oxides than supported, they will be ignored and the composition will be renormalized accordingly. Moreover, if both :literal:`Fe2O3` and :literal:`O` are provided, :literal:`O` will be recalculated as function of :literal:`Fe2O3`. Thus, if you want to prescribe a different :literal:`O` content, do not define :literal:`Fe2O3`!

|

Call optimization routine for given P-T-X
*****************************************
.. code-block:: shell   

   out         =   point_wise_minimization(P,T, data);

if a predefined test is used (see :doc:`/ckbk/predef`) or:

.. code-block:: shell   

   out     = single_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in=sys_in)

if a custom bulk-rock composition is provided.

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

Examples of serial point calculation
************************************


.. code-block:: shell

   #load MAGEMin
   using MAGEMin_C 

   data    = Initialize_MAGEMin("ig", verbose=false);

   # One bulk rock for all points
   P,T     = 10.0, 1100.0
   Xoxides = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "Fe2O3"; "K2O"; "Na2O"; "TiO2"; "Cr2O3"; "H2O"];
   X       = [48.43; 15.19; 11.57; 10.13; 6.65; 1.64; 0.59; 1.87; 0.68; 0.0; 3.0];
   sys_in  = "wt"    
   out     = single_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in=sys_in)
   Finalize_MAGEMin(data)


for the metapelite database:

.. code-block:: shell

   #load MAGEMin
   using MAGEMin_C 

   #initialize
   data    = Initialize_MAGEMin("mp", verbose=false);

   # provide bulk-rock composition
   P,T      = 2.0, 650.0
   Xoxides  = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "Fe2O3"; "K2O"; "Na2O"; "TiO2"; "MnO"; "H2O"]
   X        = [69.64; 13.76; 1.77; 1.73; 4.32; 0.4; 2.61; 2.41; 0.80; 0.07; 10.0]
   sys_in   = "wt"    
   out      = single_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in=sys_in)
   Finalize_MAGEMin(data)

for the ultramafic database:

.. code-block:: shell

   #load MAGEMin
   using MAGEMin_C 

   #initialize
   data    = Initialize_MAGEMin("um", verbose=false);

   # provide bulk-rock composition
   P,T      = 2.0, 650.0
   out      = single_point_minimization(P, T, data, test=0)
   Finalize_MAGEMin(data)



Parallel point calculation Tutorial
===================================

To compute a list of single point calculation in parallel your can use the native Julia multi-threading. To activate multi-threading simply launch the Julia terminal as:

.. code-block:: shell

   julia -t 4

or 

.. code-block:: shell

   julia --threads 4

where the number of threads depends on your system, generally twice the number of cores. 

|

Examples of serial point calculation
************************************

To run :literal:`n` points, using database :literal:`ig` and :literal:`test 0` (see :doc:`/ckbk/predef`):

.. code-block:: shell

   #load MAGEMin
   using MAGEMin_C 

   #initialize
   data    =   Initialize_MAGEMin("ig", verbose=false);
   n       =   100;
   P       =   fill(8.0,n)
   T       =   fill(800.0,n)
   out     =   multi_point_minimization(P, T, data, test=1);
   Finalize_MAGEMin(data)

Here the results are stored in :literal:`out` as :literal:`out[1:end]`. Various bulk-rock compositions can be prescribed as:

.. code-block:: shell

   #load MAGEMin
   using MAGEMin_C 

   #initialize
   data    = Initialize_MAGEMin("ig", verbose=false);

   #set P-T-X conditions
   P       = [10.0, 10.0];
   T       = [1100.0, 1100.0];
   Xoxides = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "Fe2O3"; "K2O"; "Na2O"; "TiO2"; "Cr2O3"; "H2O"];
   X1      = [48.43; 15.19; 11.57; 10.13; 6.65; 1.64; 0.59; 1.87; 0.68; 0.0; 3.0];
   X2      = [49.43; 14.19; 11.57; 10.13; 6.65; 1.64; 0.59; 1.87; 0.68; 0.0; 3.0];
   X       = [X1,X2];
   sys_in  = "wt"    
   out     = multi_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in=sys_in)


Other examples
==============

Several additional tests are provided in :literal:`./test/runtests.jl`.
