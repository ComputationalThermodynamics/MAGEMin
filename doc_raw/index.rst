.. MAGEMin documentation

.. image:: /figs/MAGMA_Logo2.png
   :width: 420
   :align: right

.. image:: /figs/GUI.png
   :width: 720
   :align: center

|

**MAGEMin v1.7.3**
==================
                                          
MAGEMin (Mineral Assemblage Gibbs Energy Minimization) is a parallel C library callable from any petrological/geodynamic tool. For a given set of pressure, temperature and bulk-rock composition MAGEMin uses a combination of linear programming, extended Partitioning Gibbs free Energy and gradient-based local minimization to compute the most stable mineral assemblage     
 
A full description of the minimization approach used in MAGEMin is given in:

Riel, N., Kaus, B. J. P., Green, E. C. R., & Berlie, N. (2022). MAGEMin, an efficient Gibbs energy minimizer: Application to igneous systems. Geochemistry, Geophysics, Geosystems, 23, e2022GC010427. https://doi.org/10.1029/2022GC010427 

|
     
Available thermodynamic dataset                       
================================

Igneous thermodynamic dataset
*****************************
                    
- Holland et al., 2018 (see http://hpxeosandthermocalc.org)
- K2O-Na2O-CaO-FeO-MgO-Al2O3-SiO2-H2O-TiO2-O-Cr2O3 chemical system
- Equations of state for
	- Pure stoichiometric phases quartz (q), cristobalite (crst), tridymite (trd), coesite (coe), stishovite (stv), kyanite (ky), sillimanite (sill), andalusite (and), rutile (ru) and sphene (sph). 
	- Solution phases spinel (spl), biotite (bi), cordierite (cd), clinopyroxene (cpx), orthopyroxene (opx), epidote (ep), garnet (g), hornblende (amp), ilmenite (ilm), silicate melt (liq), muscovite (mu), olivine (ol), ternary feldspar (pl4T), and aqueous fluid (fl).
- added May 2022


Ultramafic thermodynamic dataset
********************************
- Evans & Frost, 2021 (see http://hpxeosandthermocalc.org)
- SiO2-Al2O3MgO-FeO-O-H2O-S chemical system
- Equations of state for
	- Pure stoichiometric phases quartz (q), cristobalite (crst), tridymite (trd), coesite (coe), stishovite (stv), kyanite (ky), sillimanite (sill), pyrite (pyr)
	- Solution phases fluid (fluid), brucite (br), antigorite (atg), garnet (g), talc (t), chlorite (chl), spinel (spi), orthopyroxene (opx), pyrrhotite (po) and anthophylite (anth)
- added May 2023

Metapelite thermodynamic dataset
********************************
                    
- White et al., 2014a, 2014b (see http://hpxeosandthermocalc.org)
- K2O-Na2O-CaO-FeO-MgO-Al2O3-SiO2-H2O-TiO2-O-MnO chemical system
- Equations of state for
	- Pure stoichiometric phases quartz (q), cristobalite (crst), tridymite (trd), coesite (coe), stishovite (stv), kyanite (ky), sillimanite (sill), andalusite (and), rutile (ru) and sphene (sph). 
	- Solution phases spinel (spl), biotite (bi), cordierite (cd), orthopyroxene (opx), epidote (ep), garnet (g),  ilmenite (ilm), silicate melt (liq), muscovite (mu),  ternary feldspar (pl4T), sapphirine (sa), staurolite (st), magnetite (mt), chlorite (chl), chloritoid (ctd) and margarite (ma).
- added March 2023

Metabasite thermodynamic dataset
********************************
                    
- Green et al., 2016 (see http://hpxeosandthermocalc.org)
- K2O-Na2O-CaO-FeO-MgO-Al2O3-SiO2-H2O-TiO2-O chemical system
- Equations of state for
	- Pure stoichiometric phases quartz (q), cristobalite (crst), tridymite (trd), coesite (coe), stishovite (stv), kyanite (ky), sillimanite (sill), andalusite (and), rutile (ru) and sphene (sph). 
	- Solution phases spinel (sp), biotite (bi), orthopyroxene (opx), epidote (ep), garnet (g), ilmenite (ilm), silicate melt (liq), muscovite (mu),  ternary feldspar (pl4T), chlorite (chl), Omphacite(omph) and Augite(aug).
- added October 2023
	

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
   :caption: Julia interface:
   
   julia
   

.. toctree::
   :maxdepth: 2
   :caption: Known issues & Tips!
   
   issues
   
.. toctree::
   :maxdepth: 2
   :caption: FAQs:
   
   faq
   
.. toctree::
   :maxdepth: 2
   :caption: Developers section:
   
   dev

.. toctree::
   :maxdepth: 2
   :caption: Source code:

   src

|

Contributors
============

- MAGEMin developpers: Riel N., Kaus. B.
- Database translation and debugging: Green E., Berlie N., and Rummel L. 
     
*Contacts: nriel@uni-mainz.de, kaus@uni-mainz.de*

.. image:: /figs/Magma_alone_Logo.png
   :width: 240
   :align: right
   
This open source project was funded by the European Research Council through the `MAGMA`_ project, ERC Consolidator Grant #771143

.. _MAGMA: https://magma.uni-mainz.de
