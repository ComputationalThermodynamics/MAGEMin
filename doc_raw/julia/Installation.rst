.. MAGEMin documentation


MAGEMin_C.jll installation		  
==========================

.. image:: /figs/julia.png
   :width: 120
   :align: right

:guilabel:`MAGEMin_C.jll` like any other Julia package is easy to install and use.

Install Julia
=============

If you don't have Julia already you can download it here:

https://julialang.org/downloads/

Then, extract the archive and create a symbolic link. On linux this can be achieved as:

.. code-block:: shell

	sudo ln -s /path_to_julia/bin/julia /usr/local/bin/julia

On linux, the symbolic link makes Julia accessible from the terminal (:kbd:`ctrl` + :kbd:`alt` + :kbd:`t`) by simply typing in terminal:

.. code-block:: shell

	julia

Running Julia from the command line should open the next window:

.. image:: /figs/julia_term.png
   :width: 640
   :align: center

|

Install MAGEMin_C
==================

Open the package manager by typing :literal:`]`:

.. image:: /figs/julia_pack.png
   :width: 640
   :align: center

|

To test the package type :literal:`test MAGEMin_C`. You should get the test summary:

.. image:: /figs/julia_test.png
   :width: 360
   :align: center

|

The tests run a serie of minimizations and compare them with the expected results. If everything went fine all tests should pass.
