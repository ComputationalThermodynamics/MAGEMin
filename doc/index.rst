.. MAGEMin documentation

.. image:: /figs/opt.png
   :width: 400
   :align: center

|

Mineral Assemblage Gibbs Energy Minimization		  
============================================

**Contributors:**

- MAGEMin developpers: Riel N., Kaus. B.
- Database translation and debugging: Green E., Berlie N., and Rummel L. 
     
**Contacts: nriel@uni-mainz.de, kaus@uni-mainz.de**	 
                                           
C routine to compute a stable mineralogical assemblage using Gibbs Free Energy Minimization. We employ a combination of linear programming, partitioning Gibbs free energy and non-linear constrained minimization to solve for phase equilibria at constant pressure, temperature and bulk-rock composition.     
                              
- Igneous thermodynamic database (Holland et al., 2008)
- KNCFMASHTOCr chemical space                           
- Solid solutions (biotite, clinopyroxene, cordierite, epidote, fluids, garnet, hornblende, ilmenite, K-felspar, liquid, muscovite, olivine, orthopyroxene, plagioclase, spinel) 

                                                    
Imported libraries                       
==================

- LAPACKE (C version of LAPACK)                         
- NLopt  (https://nlopt.readthedocs.io/)                
- uthash (https://troydhanson.github.io/uthash/)        
- ketopt (https://github.com/attractivechaos/klib/blob/master/ketopt.h) 
   
.. toctree::
   :maxdepth: 2
   :caption: Installation:
   
   stp
                                                   
.. toctree::
   :maxdepth: 2
   :caption: Methods and examples:
   
   man
   
.. toctree::
   :maxdepth: 2
   :caption: Cookbook:
   
   ckbk

.. toctree::
   :maxdepth: 2
   :caption: Developers section:
   
   dev

.. toctree::
   :maxdepth: 2
   :caption: Source code:

   src
