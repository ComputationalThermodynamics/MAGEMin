<img src="./pics/GUI.png" alt="drawing" width="640" alt="centered image"/>

# Mineral Assemblage Gibbs Energy Minimization (MAGEMin)

MAGEMin is written as a parallel C library callable from any petrological/geodynamic tool. For a given set of pressure, temperature and bulk-rock composition MAGEMin uses a combination of linear programming, extended Partitioning Gibbs free Energy and gradient-based local minimization to compute the most stable mineral assemblage     
  
## Available thermodynamic dataset
- Holland et al., 2018 (see http://hpxeosandthermocalc.org)
- K2O-Na2O-CaO-FeO-MgO-Al2O3-SiO2-H2O-TiO2-O-Cr2O3 chemical system
- Equations of state for
	- Pure stoichiometric phases quartz (q), cristobalite (crst), tridymite (trd), coesite (coe), stishovite (stv), kyanite (ky), sillimanite (sill), andalusite (and), rutile (ru) and sphene (sph). 
	- Solution phases spinel (spn), biotite (bi), cordierite (cd), clinopyroxene (cpx), orthopyroxene (opx), epidote (ep), garnet (g), hornblende (hb), ilmenite (ilm), silicate melt (liq), muscovite (mu), olivine (ol), ternary feldspar (pl4T), and aqueous fluid (fl).
     
## Imported libraries
- LAPACKE (C version of LAPACK)
- NLopt (https://nlopt.readthedocs.io/)
- uthash (https://troydhanson.github.io/uthash/)
- ketopt (https://github.com/attractivechaos/klib/blob/master/ketopt.h)

## Usage 

## Installation 

## Contributing
You are very welcome to request new features and point out bugs by opening an issue. You can also help by adding features and creating a pull request.

## Development roadmap

## Funding
Development of this software package was funded by the European Research Council under grant ERC CoG #771143 - [MAGMA](https://magma.uni-mainz.de).
