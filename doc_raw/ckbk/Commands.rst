.. MAGEMin documentation

List of arguments		  
=================

MAGEMin is run by using command line arguments when executing the binary file.

The list of valid command line arguments is the following

+---------------+-----------------------------------------------+
|  Arguments    |                Application                    | 
+===============+===============================================+
| -\-version    | Displays MAGEMin version                      |
+---------------+-----------------------------------------------+
| -\-help       | Displays help                                 |
+---------------+-----------------------------------------------+
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
| -\-Temp=y     | Temperature in Celsius                        |
+---------------+-----------------------------------------------+
| -\-Bulk=[y]   | Bulk rock composition in molar amount         |
+---------------+-----------------------------------------------+
| -\-Gam=[y]    | Gamma, when a guess of gamma is known         |
+---------------+-----------------------------------------------+
| -\-sys_in=""  | system comp "mol" or "wt", default is "mol"   |
+---------------+-----------------------------------------------+

where *x* is an ``integer``, *y* a ``float``/``double`` and *[]* a comma-separated list of size *number of oxydes*. 

The list of oxides must be given in the following order: 

+------+-------+-----+-----+------+-----+------+------+---+-------+-----+
| SiO2 | Al2O3 | CaO | MgO | FeOt | K2O | Na2O | TiO2 | O | Cr2O3 | H2O |
+------+-------+-----+-----+------+-----+------+------+---+-------+-----+

Note that FeOt is the total iron!

|

Single point calculation		  
========================

Using previously defined arguments, a valid command to run a single point calculation with MAGEMin is for instance:

.. code-block:: shell

	./MAGEMin --Verb=1 --Temp=718.750 --Pres=30.5000 --test=0 >&log.txt

Here, the verbose is active and the bulk rock composition of *test 0* is selected. The output of the verbose is saved as a log file *log.txt*.

If you want to do a computation using a different bulk rock composition you can pass the custom bulk such as:

.. code-block:: shell

	./MAGEMin --Verb=1 --Temp=488.750 --Pres=3.5000 --Bulk=41.49,1.57,3.824,50.56,5.88,0.01,0.25,0.10,0.1,0.0 --sys_in=mol
	
|

Multiple points	calculation	  
===========================

To run multiple points at once you can pass an input file containing the list of points such as

.. code-block:: shell

	./MAGEMin --Verb=1 --File=path_to_file --n_points=x
	
where "path_to_file" is the location of the file and "x" is an integer corresponding to the total number of points contained in the file. The file must have one point per line using the following structure

+------------+----------------+----------------+---------+---------+-----+---------+
|  Mode(0-1) | Pressure(kbar) | Temperature(C) | Bulk_1  | Bulk_2  | ... | Bulk_n  |
+------------+----------------+----------------+---------+---------+-----+---------+

- mode = 0 for global minimization
- "Bulk_n" is the bulk rock composition in oxides (mol or wt fraction). 

A valid list of points for an in-built test is given by is for instance:

.. code-block:: shell

	MAGEMin_input.dat:
	
	0 0.0 800.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
	0 4.0 800.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
	0 8.0 800.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
	0 8.0 800.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
	(...)

You can compute the list of points using:

.. code-block:: shell

	./MAGEMin --Verb=1 --File=/path_to_file/MAGEMin_input.dat --n_points=4 --test=0

To compute a custom liast of bulk-rock compositions you have to provide the oxide composition and replace the "0.0" such as:

.. code-block:: shell

	MAGEMin_input.dat:
	
	0 0.0 800.0 41.49 1.57 4.824 52.56 5.88 0.01 0.25 0.10 0.1 0.0
	0 4.0 800.0 44.49 1.56 3.24 48.56 5.2 0.01 0.25 0.10 0.1 0.0
	0 8.0 800.0 42.49 1.27 3.84 51.56 4.28 0.01 0.25 0.10 0.1 0.0
	0 8.0 800.0 40.49 1.87 1.824 50.56 6.08 0.01 0.25 0.10 0.1 0.0
	(...)

Then compute the list of points while indicating the system composition unit (:literal:`mol` or :literal:`wt` fraction):

.. code-block:: shell

	./MAGEMin --Verb=1 --File=/path_to_file/MAGEMin_input.dat --n_points=4 --sys_in=mol


|

Multiple points in parallel 
===========================

To run a list of points in parallel, you simply need to call "MPI" before MAGEMin and give the number of cores you want to use. Valid calls using previously defined input file are for instance:

.. code-block:: shell

	mpirun -np 8 ./MAGEMin --Verb=1 --File=/path_to_file/MAGEMin_input.dat --n_points=4 --test=0
	mpiexec -n 8 ./MAGEMin --Verb=1 --File=/path_to_file/MAGEMin_input.dat --n_points=4 --test=0

where 8 is the desired number of cores.
