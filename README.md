
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://computationalthermodynamics.github.io/MAGEMin/index.html)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10573416.svg)](https://doi.org/10.5281/zenodo.10573416)
<img src="./pics/GUI.png" alt="drawing" width="640" alt="centered image"/>

# Mineral Assemblage Gibbs Energy Minimization (MAGEMin)
`MAGEMin` is a Gibbs energy minimization solver package, which computes the thermodynamically most stable assemblage for a given bulk rock composition and pressure/temperature condition. It also returns parameters such as melt fraction or density, which can be combined with geodynamic/petrological tools to simulate, for example, the evolving chemistry of a crystallising melt.

`MAGEMin` is written as a parallel C library and uses a combination of linear programming, the extended Partitioning Gibbs free Energy approach and gradient-based local minimization to compute the most stable mineral assemblage. In this, it differs from existing approaches which makes it particularly suitable to utilize modern multicore processors.

We also provide a MATLAB-based graphical user interface to help computing pseudosections for given bulk rock composition and a [julia interface](https://github.com/ComputationalThermodynamics/MAGEMin_C.jl).

## Documentation
Full support to install and use MAGEMin is available here: https://computationalthermodynamics.github.io/MAGEMin/index.html

## Installing MAGEMin

**Quick start**

**Julia interface**

To make it easier to interface `MAGEMin` with other (geodynamic) codes, we provide a [julia interface](https://github.com/ComputationalThermodynamics/MAGEMin_C.jl) to the `MAGEMin` C library, with which you can perform pointwise calculations. 

**with Matlab using Julia**

You can use `MAGEMin` is through the MATLAB graphical user interface, which has an installation script to download the correct parallel binaries for your system (created using [BinaryBuilder](https://binarybuilder.org) & [julia](https://julialang.org)).

Follow these steps:
1) Download a zip file with the most recent release of `MAGEMin` (click on the green `Code` button @ the top of this page) and unzip it on your machine.
2) Open the `PlotPseudosection` graphical user interface from MATLAB (2020+). 
3) Follow the binary installation instructions (which requires you to install a recent [julia](https://www.julialang.org) version).
4) After this you are ready to get started, for example by pushing the `Start new computation` button. 

Note that the Matlab GUI is not maintained anymore as the primary phase diagram generator is now the Julia app!

**Manual compilation**

if you wish, you can also compile MAGEMin yourself, which requires you to install these packages as well:
- MPICH (to allow parallel computations)
- LAPACKE (C version of LAPACK)
- NLopt (https://nlopt.readthedocs.io/)

Details and guidelines are given in the extended documentation: https://computationalthermodynamics.github.io/MAGEMin/index.html  

In addition, we make use of [uthash](https://troydhanson.github.io/uthash/) and [ketopt](https://github.com/attractivechaos/klib/blob/master/ketopt.h).


## Available thermodynamic datasets
The MAGEMin algorithm is general and can be used with any thermodynamic database that are hardcoded for speed reasons. Presently the igneous (Holland et al., 2018) and the metapelite (White et al., 2014) database are available.


**Igneous database**

The hydrous mafic melting model of Holland et al. 2018 can be used to simulate the fractional crystallisation from a hydrous basalt to a felsic melt. 

The details of this thermodynamic solid solution and endmember database are:
- Holland et al., 2018 (see http://hpxeosandthermocalc.org)
- K2O-Na2O-CaO-FeO-MgO-Al2O3-SiO2-H2O-TiO2-O-Cr2O3 chemical system
- Equations of state for
	- Pure stoichiometric phases quartz (q), cristobalite (crst), tridymite (trd), coesite (coe), stishovite (stv), kyanite (ky), sillimanite (sill), andalusite (and), rutile (ru) and sphene (sph). 
	- Solution phases spinel (spn), biotite (bi), cordierite (cd), clinopyroxene (cpx), orthopyroxene (opx), epidote (ep), garnet (g), hornblende (hb), ilmenite (ilm), silicate melt (liq), muscovite (mu), olivine (ol), ternary feldspar (pl4T), and aqueous fluid (fl).


**Metapelite database**

The metapelitic model (extended with MnO, White et al., 2014) allows to compute the mineral assemblage from low temperature to supra-solidus conditions.

- Added March 2023, `MAGEMin v1.3.0` 
- White et al., 2014a, 2014b (see http://hpxeosandthermocalc.org)
- K2O-Na2O-CaO-FeO-MgO-Al2O3-SiO2-H2O-TiO2-O-MnO chemical system
- Equations of state for
	- Pure stoichiometric phases quartz (q), cristobalite (crst), tridymite (trd), coesite (coe), stishovite (stv), kyanite (ky), sillimanite (sill), andalusite (and), rutile (ru) and sphene (sph). 
	- Solution phases spinel (spn), biotite (bi), cordierite (cd), orthopyroxene (opx), epidote (ep), garnet (g),  ilmenite (ilm), silicate melt (liq), muscovite (mu),  ternary feldspar (pl4T), sapphirine (sa), staurolite (st), magnetite (mt), chlorite (chl), chloritoid (ctd) and margarite (ma).

**Ultramafic thermodynamic dataset**
- aadded May 2023, `MAGEMin v1.3.2` 
- Evans & Frost, 2021 (see http://hpxeosandthermocalc.org)
- SiO2-Al2O3MgO-FeO-O-H2O-S chemical system
- Equations of state for
	- Pure stoichiometric phases quartz (q), cristobalite (crst), tridymite (trd), coesite (coe), stishovite (stv), kyanite (ky), sillimanite (sill), pyrite (pyr)
	- Solution phases fluid (fluid), brucite (br), antigorite (atg), garnet (g), talc (t), chlorite (chl), spinel (spi), orthopyroxene (opx), pyrrhotite (po) and anthophylite (anth)
- added May 2023

**Metabasite thermodynamic dataset**

- Added October 2023, `MAGEMin v1.3.5`
- Green et al., 2016 (see http://hpxeosandthermocalc.org)
- K2O-Na2O-CaO-FeO-MgO-Al2O3-SiO2-H2O-TiO2-O chemical system
- Equations of state for
	- Pure stoichiometric phases quartz (q), cristobalite (crst), tridymite (trd), coesite (coe), stishovite (stv), kyanite (ky), sillimanite (sill), andalusite (and), rutile (ru) and sphene (sph). 
	- Solution phases spinel (sp), biotite (bi), orthopyroxene (opx), epidote (ep), garnet (g), ilmenite (ilm), silicate melt (liq), muscovite (mu),  ternary feldspar (pl4T), chlorite (chl), Omphacite(omph) and Augite(aug).

	


Please keep in mind that the datasets are only calibrated for a limited range of `P`,`T` and `bulk rock` conditions. If you go too far outside those ranges, `MAGEMin` (or most other thermodynamic software packages for that matter) may not converge or give bogus results. 
Developing new, more widely applicable, thermodynamic datasets is a huge research topic, which will require funding to develop the models themselves, as well as to perform targeted experiments to calibrate those models.


## Citation
An open-acces paper describing the methodology is:

- Riel N., Kaus B.J.P., Green E.C.R., Berlie N., (2022) MAGEMin, an Efficient Gibbs Energy Minimizer: Application to Igneous Systems. *Geochemistry, Geophysics, Geosystems* 23, e2022GC010427 [https://doi.org/10.1029/2022GC010427](https://doi.org/10.1029/2022GC010427)

If you use the software, we really appreciate if you cite this study. We also appreciate stars (see the top of this page). 

## Contributing
You are very welcome to request new features and point out bugs by opening an issue (top of this page). You can also help by adding features and creating a pull request. 

## Funding
Development of this software package was funded by the European Research Council under grant ERC CoG #771143 - [MAGMA](https://magma.uni-mainz.de).

<img src="./pics/MAGMA_Logo.png" alt="drawing" width="480" alt="centered image"/>
