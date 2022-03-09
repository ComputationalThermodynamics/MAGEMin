.. MAGEMin documentation

List of arguments		  
=================

MAGEMin is run by using command line arguments when executing the binary file.

The list of valid command line arguments is the following

+---------------+-----------------------------------------------+
|  Arguments    |                Application                    | 
+===============+===============================================+
| -\-Verb=x     | Verbose option, 0. inactive, 1. active        |
+---------------+-----------------------------------------------+
| -\-File=path  | Given file for multiple point calculation     |
+---------------+-----------------------------------------------+
| -\-n_points=x | Number of points when using *File* argument   |
+---------------+-----------------------------------------------+
| -\-test=x     | Run calculation on included test compositions |
+---------------+-----------------------------------------------+
| -\-Pres=y     | Pressure in kilobar                           |
+---------------+-----------------------------------------------+
| -\-Temp=y     | Temperature in celsus                         |
+---------------+-----------------------------------------------+
| -\-Bulk=[y]   | Bulk rock composition in molar amount         |
+---------------+-----------------------------------------------+
| -\-Gam=[y]    | Gamma, when a guess of gamma is known         |
+---------------+-----------------------------------------------+

where *x* is an integer, *y* a float/double and *[]* an array of size *number of oxydes*.


Single point calculation		  
========================

Using previously defined arguments, a valid command to run a single point calculation with MAGEMin is for instance:

.. code-block:: shell

	./MAGEMin --Verb=1 --Temp=718.750 --Pres=30.5000 --test=0 >&log.txt

Here, the verbose is active and the bulk rock composition of *test 0* is selected. The output of the verbose is saved as a log file *log.txt*.

Multiple points	calculation	  
===========================

When running multiple points using the command line argument *-\-file=path_to_file* several parameters can be passed to MAGEMin.

One line per point has to be provided as follows

+------------+----------------+----------------+---------+---------+-------+---------+
|  Mode(0-1) | Pressure(kbar) | Temperature(C) | Gamma_1 | Gamma_2 | .\.\. | Gamma_n |
+------------+----------------+----------------+---------+---------+-------+---------+

Note that mode = 0 for global minimization


