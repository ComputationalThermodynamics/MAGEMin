.. MAGEMin documentation

|

This section reports encountered issues with MAGEMin and provide a way to fix these problems. If your problem is not listed here, please open an issue on:
https://github.com/ComputationalThermodynamics/MAGEMin

Firewall warning (Apple M1)
===========================

When using the cross-platform Binarybuilder users have reported the following message when running MAGEMin on Apple M1:

.. code-block:: shell

	Do you want the application “MAGEMin” to accept incoming network connections?

First of all, MAGEMin **does not** require an internet connection to be executed and **no data is ever exchanged** with a distant server.

We provide a fix in ``/julia/firewall_macos.jl`` which essentially blocks incoming traffic for the MAGEMin binary. It can be run from the terminal (upon having ``Julia`` installed on your device) with

.. code-block:: shell

	$julia firewall_macos.jl

Note that you do need to have the sudo password for your machine. If that is not the case, you will have you ask your system administrator for help.

You will need to run this once for every version of MAGEMin (if you update at some stage in the future, this will likely have to be repeated).
