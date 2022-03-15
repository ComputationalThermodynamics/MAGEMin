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

To run multiple points at once you can pass an input file containing the list of points such as

.. code-block:: shell

	./MAGEMin --Verb=1 --File=path_to_file --n_points=x
	
where "path_to_file" is the location of the file and "x" is an integer corresponding to the total number of points contained in the file. The file must have one point per line using the following structure

+------------+----------------+----------------+---------+---------+-----+---------+
|  Mode(0-1) | Pressure(kbar) | Temperature(C) | Gamma_1 | Gamma_2 | ... | Gamma_n |
+------------+----------------+----------------+---------+---------+-----+---------+

- mode = 0 for global minimization
- "Gamma_n" is the array of oxide chemical potential (defining the Gibbs-hyperplane). This functionality is in development and will be used in a future update to provide an initial guess of the Gibbs-hyperplane to improve performances.

A valid list of points is for instance:

.. code-block:: shell

	MAGEMin_input.dat:
	
	0 0.0 800.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
	0 4.0 800.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
	0 8.0 800.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
	0 8.0 800.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
	(...)


Multiple points in parallel 
===========================

To run a list of points in parallel, you simply need to call "MPI" before MAGEMin and give the number of cores you want to use. Valid calls using previously defined input file are for instance:

.. code-block:: shell

	mpirun -np 8 ./MAGEMin --Verb=1 --File=/path_to_file/MAGEMin_input.dat --n_points=4
	mpiexec -n 8 ./MAGEMin --Verb=1 --File=/path_to_file/MAGEMin_input.dat --n_points=4

where "y" is the desired number of cores.
