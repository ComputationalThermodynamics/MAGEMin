
.. MAGEMin documentation

.. image:: /figs/macOS.png
   :width: 64
   :align: right
   
External libraries		  
==================

The installation details for Mas OS X use ``Homebrew``. However, the libraries can also be installed using ``MacPorts``.

**1. C and fortran compilers**

Using either ``gcc`` or ``clang`` to compile MAGEMin is up to you as the runtime performances are similar. However, in the event you want to modify MAGEMin for your own use, I would advice that you compile MAGEMin with ``clang`` as the error handling system is more strict which will save you from having unexpected segfaults errors...

.. code-block:: shell

	brew install llvm
	brew install gcc

Note that the ``gcc`` package comes with ``gcc``, ``g++`` and ``gfortran``

**2. MPICH** (Message Passing Interface)

.. code-block:: shell

	brew install mpich

Note that ``openmpi`` can equally be used.

**3. lapacke** (``C`` version of the fortran ``lapack`` library, should now be included in the ``lapack`` libraries)

.. code-block:: shell

	brew install lapack
	
If the ``lapacke`` libraries are not included you can download the ``lapack`` package from netlib that includes it.

**4. NLopt** (Non Linear optimization library)

``NLopt`` can be installed using

.. code-block:: shell

	brew install nlopt

.. image:: /figs/macOS.png
   :width: 64
   :align: right
   

MAGEMin
=======
	
Choose the ``C`` compiler in the first line of the ``Makefile`` by commenting out one

.. code-block:: shell

	#CC=gcc
	CC=clang

Make sure the ``MPICH`` paths for libraries and include directory in the ``Makefile`` are correct for instance:

.. code-block:: shell

	LIBS   += (...) /opt/homebrew/lib/libmpich.dylib
	INC     = (...) -I/opt/homebrew/include 

Do the same for ``lapacke``:

.. code-block:: shell

	LIBS   += (...) /opt/homebrew/opt/lapack/lib/liblapacke.dylib
	INC    += (...) -I/opt/homebrew/opt/lapack/include
	
And ``NLopt``:

.. code-block:: shell

	LIBS   += (...) /opt/homebrew/lib/libnlopt.dylib
	
	
Which should give:

.. code-block:: shell

	LIBS    = -lm -framework Accelerate /opt/homebrew/opt/lapack/lib/liblapacke.dylib /opt/homebrew/lib/libnlopt.dylib /opt/homebrew/lib/libmpich.dylib
	INC     = -I/opt/homebrew/opt/lapack/include -I/opt/homebrew/include 
	
Note that this setup is provided by default in the ``Makefile`` for Mac OS X. As long as you installed every package using ``Homebrew`` you should be able to install MAGEMin without modifying these entries.

If you decided to use ``openmpi`` instead of ``mpich`` your ``Makefile`` should look like:

.. code-block:: shell

	LIBS    = -lm -framework Accelerate /opt/homebrew/opt/lapack/lib/liblapacke.dylib /opt/homebrew/opt/nlopt/lib/libnlopt.dylib /opt/homebrew/opt/openmpi/lib/libmpi.dylib  
	INC     = -I/opt/homebrew/opt/openmpi/include/ -I/opt/homebrew/opt/lapack/include -I/usr/local/include -I/opt/homebrew/opt/nlopt/include/

Then simply enter MAGEMin directory and compile MAGEMin as:

.. code-block:: shell
	
	make clean; make all;

Note that by default the optimization flag ``-O3`` and debugging flag ``-g`` are used.

To test if MAGEMin compilation was successful you can for instance check the version of MAGEMin by running:

.. code-block:: shell
	
	./MAGEMin --version

