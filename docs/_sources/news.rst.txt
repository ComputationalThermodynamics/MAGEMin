.. MAGEMin documentation
  
Seismic velocity calculation
============================

We are working on adding the calculation of seismic velocities using the method described in Connolly and Kerrick (2002) such as

.. math:: 
	v_p = \sqrt{ \frac{K_s + 4 \mu / 3}{\rho}}

and

.. math:: 
	v_s = \sqrt{ \frac{\mu}{\rho}}


where :math:`\rho` is the density, :math:`K_s` is the adiabatic shear modulus and :math:`\mu` is the shear modulus.
