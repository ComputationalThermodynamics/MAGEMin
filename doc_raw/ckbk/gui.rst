.. MAGEMin documentation

.. image:: /figs/GUI_label.png
	:width: 640
	:align: center
   
|

Usage		  
=====

In order to generate pseudosection (including isochemical phase diagrams) we provide a in-house GUI developped using MATLAB 2020b and tested on MATLAB 2019b.

**Specifications**

- The GUI sends a list of pressure-temperature points to MAGEMin for a specified bulk-rock composition and receives back the stable phase mineral assemblage.
- The GUI uses adaptive mesh refinement to refine the reactions lines.
- Possibility to save and load pseudosection
- Various display options including, stable phases, Gibbs energy, phase fractions, phase compositions (isopleth), densities etc...

KLB-1 Pseudosection tutorial
============================

As a simple tutorial we present the necessary steps needed to produce a pseudosection for KLB-1 peridotite in the KNCFMASTOCr chemical system and from 0.0 to 50 kbar and 800 to 2000°C.

Sarting GUI
***********

Using MATLAB (2019+) launch the GUI within MAGEMin working directory

.. code-block:: shell
	
	cd /your_path_to_MAGEMin/MAGEMin_1.x/
	run PlotPseudosection.mlapp
	
This opens the following window:
	
.. image:: /figs/GUI_1.png
   :width: 480
   :align: center
   
----

- On the displayed window (simulation) you can choose and modifiy several parameters/options. First, you can change the name of the simulation (or save/load a pre-existing simulation)

.. image:: /figs/GUI_infos.png
   :width: 240
   :align: center
   
|

#. Change the filename from "test" to "KLB-1_pseudosection"
   
Set P-T conditions
******************

- Below you can accept or modify the suggested pressure and temperature range. The step parameter defines the number of points along the pressure and temperature axes.

.. image:: /figs/GUI_PT.png
   :width: 240
   :align: center
   
|

#. Use the default pressure and temperature range and change the pressure step to 4 kbar and the temperature to 100 °C

Mesh refinement and computation options
***************************************

- Under the P-T range selection, you can access the refinement options and the calculation parameters. By default the refinement is set on 2 levels and the parallel calculation is deactivated.

.. image:: /figs/GUI_refine.png
   :width: 240
   :align: center
   
|
 
#. Change the default refinement value to 4
#. Use default value for verbose
#. Activate parallel computation by ticking the box "Compute all points at once"

----

- Ticking "Compute all points at once" releases the following options: "# of MPI ranks", "mpiexec path" and "Perform computations on remote server"

.. image:: /figs/GUI_mpi.png
   :width: 240
   :align: center
   
|

#. Change "# of MPI ranks" according to the number of core available on your machine
#. Change "mpiexec path" to your own mpi path e.g. commonly "/usr/bin/" for Linux-based system
#. Leave "Perform computations on remote server" unticked

Choose bulk-rock composition
****************************

- In the middle "Composition" panel you can either specify your own bulk-rock composition or access the in-built list of bulk-rock compositions

.. image:: /figs/GUI_bulk.png
   :width: 240
   :align: center
   
|
 
#. Use default option "pre-defined test 0" (KLB-1 peridotite)

Compute pseudosection
*********************

- The set of parameters should look like this.

.. image:: /figs/GUI_ready.png
   :width: 480
   :align: center
   
|

#. Click on "Start new computation" at the bottom of the "Thermodynamics" panel to launch the computation. If everything is correctly setup it should opens the folloying window

.. image:: /figs/GUI_running.png
   :width: 240
   :align: center
   
|

- At the end of each refinement step, a set of 3 figures are displayed and updated. They show the "melt fraction", the "adaptive mesh refinement" and "the number of stable phases"

.. image:: /figs/GUI_figs.png
   :width: 800
   :align: center

----

- Once the computation is over, the GUI displays the pseudosection. The number of minimization points is shown at the bottom right, for this example "3493 points". Note that by default the number of stable phases is displayed together wit the mesh refinement grid. Note that, you are in the "Pseudosection" window.

.. image:: /figs/GUI_pseudo.png
   :width: 480
   :align: center
   
|

Refine pseudosection
********************
   
- From there you can decide to refine the pseudosection as the resolution is rather low.

.. image:: /figs/GUI_PS_refine.png
   :width: 480
   :align: center
   
|

#. click on "Refine all phase boundaries" in the bottom right corner, and then "Start computation"

Change displayed field
**********************

- Upon refinement you can choose to change the field that is displayed

.. image:: /figs/GUI_PS_refine.png
   :width: 480
   :align: center
   
|

1. To display the variance click on the "Plot color" menu and select "Variance"

.. image:: /figs/GUI_pcolor.png
   :width: 240
   :align: center
   
|
 
2. To hide the mesh refinement grid click, set "Show edges" to "Off", which gives

.. image:: /figs/GUI_variance.png
   :width: 480
   :align: center
   
|

Display stable phases labels
****************************
 
- In the bottom right menu, you can now access the plotting options

.. image:: /figs/GUI_plot_menu.png
	:width: 240
	:align: center
   
|

#. Tick "Add phase labels" and then click on the different fields to display the stable mineral assemblage directly on the grid such as

.. image:: /figs/GUI_label.png
	:width: 480
	:align: center
   
|
	
Note that the label box position, with respect to the field you click on, remains the same e.g. in a top-right position. If you miss-clicked and two boxes overlap, don't panick, you can always click on "remove last label", or even "remove all labels" and redo the process


Save pseudosection data
************************

- You have the option to save the pseudosection data in a MATLAB "mat" file.

.. image:: /figs/GUI_save.png
	:width: 240
	:align: center
   
|

1. To save the pseudosection data come back to the simulation window, by clicking on the top left "simulation" tab. Then click on "save". This created the following file

.. code-block:: shell
	
	KLB-1_pseudosection.mat
 
2. If you close and relaunch the GUI, you can load the saved pseudosection by entering its name in the "Filename" window and click on "load". Note that loading can potentially take a long time if the pseudosection you saved contains a large number of points


