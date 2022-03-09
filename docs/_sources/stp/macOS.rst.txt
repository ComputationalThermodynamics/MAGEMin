.. MAGEMin documentation

External libraries		  
==================

**1. C and fortran compilers**

Using either *gcc* or *clang* to compile *MAGEMin* is up to you as the runtime performances are similar. However, in the event you want to modify *MAGEMin* for your own use, I would advice that you compile *MAGEMin* with *clang* as the error handling system is more strict which will save you from having unexpected segfaults errors...

.. code-block:: shell

	sudo apt-get install gcc
	sudo apt-get install clang
	sudo apt-get install gfortran

**2. Open MPI** (Message Passing Interface)

.. code-block:: shell

	sudo apt-get install openmpi-bin

Note that *mpich* can equally be used.

**3. lapacke** (C version of the fotran lapack library)

.. code-block:: shell

	sudo apt-get install liblapacke-dev

**4. NLopt** (Non Linear optimization library)

First *cmake* must be installed on your machine

.. code-block:: shell
	
	sudo apt-get install cmake

Then NLopt can be downloaded and installed

.. code-block:: shell

	git clone https://github.com/stevengj/nlopt.git
	cd nlopt
	mkdir build
	cd build
	cmake ..
	make
	sudo make install


MAGEMin
=======
	
Choose the C compiler in the first line of the *Makefile* by commenting out one or another

.. code-block:: shell

	#CC=gcc
	CC=clang

Make sure the open MPI paths for libraries and include directory in *Makefile* are correct. By default the paths to openmpi are the following:

.. code-block:: shell

	LIBS   += (...) -L/usr/lib/x86_64-linux-gnu/openmpi/lib -lmpi
	INC     = (...) -I/usr/lib/x86_64-linux-gnu/openmpi/include/

Depending on the machine on which you want to install MAGEMin, you might need to manually give the paths to NLopt libraries and include directory too

.. code-block:: shell

	LIBS   += (...) --L/local/home/kwak/nlopt_install/install/lib -lnlopt 
	INC     = (...) -I/local/home/kwak/nlopt_install/install/include


Then compile MAGEMin:

.. code-block:: shell
	
	make clean; make all;

Note that by default the optimization flag *-O3* and debugging flag *-g* are used.
