<img src="./pics/MAGMA_Logo.png" alt="drawing" width="480" alt="centered image"/>

<img src="./pics/GUI.png" alt="drawing" width="640" alt="centered image"/>

# Mineral Assemblage Gibbs Energy Minimization (MAGEMin)
MAGEMin is a Gibbs energy minimization solver package, which computes the thermodynamically most stable assemblage for a given bulk rock composition and pressure/temperature condition. It also returns parameters such as melt fraction or density, which can be combined with geodynamic/petrological tools to simulate, for example, the evolving chemistry of a crystallising melt.

MAGEMin is written as a parallel C library and uses a combination of linear programming, extended Partitioning Gibbs free Energy and gradient-based local minimization to compute the most stable mineral assemblage. In this, it differs from exisisting approaches which makes it particularly suitable to utilize modern multicore processors.

We also provide a MATLAB-based graphical user interface to help computing pseudosections for given bulk rock composition.
  
 
## Available thermodynamic dataset
The MAGEMin algorithm is general and can in principle be used with any thermodynamic database. Yet, for speed reasons, the current implementation hardcoded the hydrous mafic melting model of Holland et al. 2018 which is callibrated for hydrous mafic melts and can be used to simulate the fractional crystallisation from a hydrous basalt to a felsic melt. 

The details of this thermodynamic solid solution and endmember database are:
- Holland et al., 2018 (see http://hpxeosandthermocalc.org)
- K2O-Na2O-CaO-FeO-MgO-Al2O3-SiO2-H2O-TiO2-O-Cr2O3 chemical system
- Equations of state for
	- Pure stoichiometric phases quartz (q), cristobalite (crst), tridymite (trd), coesite (coe), stishovite (stv), kyanite (ky), sillimanite (sill), andalusite (and), rutile (ru) and sphene (sph). 
	- Solution phases spinel (spn), biotite (bi), cordierite (cd), clinopyroxene (cpx), orthopyroxene (opx), epidote (ep), garnet (g), hornblende (hb), ilmenite (ilm), silicate melt (liq), muscovite (mu), olivine (ol), ternary feldspar (pl4T), and aqueous fluid (fl).
     
## Installing MAGEMin

**Quick start**
The easiest way to use `MAGEMin` is through the MATLAB graphical user interface, which has an installation script to download the correct parallel binaries for your system (creauing using [BinaryBuilder](https://binarybuilder.org) & julia).

Do these 3 steps:
1) Download a zip file with the most recent release of `MAGEMin` (click on the green Code button @ the top of this page) and unzip it on your machine .
2) Open the `PlotPseudosection` graphical user interface from matlab. 
3) Follow the binary installation instructions (which requires you to install a recent [julia](https://www.julialang.org) version).

**Manual compilation**
if you wish, you can also compile MAGEMin yourself, which requires you to install these packages as well:
- MPICH (to allow parallel computations)
- LAPACKE (C version of LAPACK)
- NLopt (https://nlopt.readthedocs.io/)
  
In addition, we make use of [uthash](https://troydhanson.github.io/uthash/) and [ketopt](https://github.com/attractivechaos/klib/blob/master/ketopt.h).

## Documentation
Full support to install and use MAGEMin is available here: https://computationalthermodynamics.github.io/MAGEMin/index.html

## Contributing
You are very welcome to request new features and point out bugs by opening an issue. You can also help by adding features and creating a pull request.

## Funding
Development of this software package was funded by the European Research Council under grant ERC CoG #771143 - [MAGMA](https://magma.uni-mainz.de).
