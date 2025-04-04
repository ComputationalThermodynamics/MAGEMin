
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://computationalthermodynamics.github.io/MAGEMin_C.jl/dev/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10573416.svg)](https://doi.org/10.5281/zenodo.10573416)

<img src="https://raw.githubusercontent.com/ComputationalThermodynamics/repositories_pictures/main/MAGEMinApp/readme_pd.png?raw=true" alt="drawing" width="480" alt="centered image"/>

*Example of auto-labelled isochemical phase diagram for KLB-1 peridotite computed using MAGEMinApp*

# Mineral Assemblage Gibbs Energy Minimization (MAGEMin)
`MAGEMin` is a Gibbs energy minimization solver package, which computes the thermodynamically most stable assemblage for a given bulk rock composition and pressure/temperature condition. It also returns parameters such as melt fraction or density, which can be combined with geodynamic/petrological tools to simulate, for example, the evolving chemistry of a crystallising melt.

`MAGEMin` is written as a parallel C library and uses a combination of linear programming, the extended Partitioning Gibbs free Energy approach and gradient-based local minimization to compute the most stable mineral assemblage. In this, it differs from existing approaches which makes it particularly suitable to utilize modern multicore processors.

While `MAGEMin` is the engine for the prediction of the stable phases, using it is more convenient through the [julia interface](https://github.com/ComputationalThermodynamics/MAGEMin_C.jl) `MAGEMin_C` and/or the [web-browser julia app](https://github.com/ComputationalThermodynamics/MAGEMinApp.jl) `MAGEMinApp`. 

## Available thermodynamic database
 **Mantle** (Holland et al., 2013), **Metapelite** (White et al., 2014), **Metabasite** (Green et al., 2016), **Igneous** (Holland et al., 2018) and **Ultramafic** (Evans & Frost, 2021).

### MAGEMin_C

[MAGEMin_C](https://github.com/ComputationalThermodynamics/MAGEMin_C.jl) allow to quickly and easily compute single point minimization (serial and parallel) using `Julia` and retrieve the results in a structure. This gives flexibility to the user on how to treat the data. Some programming experience in `Julia` are necessary.

### MAGEMinApp
[MAGEMinApp](https://github.com/ComputationalThermodynamics/MAGEMinApp.jl) is the graphic user interface developped to compute various type of phase diagrams (P-T, T-X, P-X, PT-X) and PTX paths. MAGEMinApp offers the options to display isocontours, select solution phase for the calculation, automatic labelling of the phase fields, saving the diagrams as vector graphic files and export the data as tables and csv files...

## Installing and using MAGEMinApp

`MAGEMinApp` is a web-browser application ([see repository](https://github.com/ComputationalThermodynamics/MAGEMinApp.jl)) developped in Julia using `Dash.jl` that relies on `MAGEMin_C` (which relies `MAGEMin`) to compute phase diagrams (PT, TX and PX) but also fractional melting and crystallization paths.

To install `MAGEMinApp` simply do:
```julia
julia> ]
pkg> add MAGEMinApp
```

To run the app:

```julia
julia -t 6 			# here 6 is the number of used threads. You can adjust the value to your machine to compute the diagrams faster!
julia> using MAGEMinApp
julia> App()
[ Info: Listening on: 127.0.0.1:8050, thread id: 2
```

Then copy and paste the address `127.0.0.1:8050` in your web-browser 

<img src="https://raw.githubusercontent.com/ComputationalThermodynamics/repositories_pictures/main/MAGEMinApp/MAGEMin_app.png?raw=true" alt="drawing" width="640" alt="centered image"/>

\
## With Matlab using Julia (not maintained anymore)

You can use `MAGEMin` is through the MATLAB graphical user interface, which has an installation script to download the correct parallel binaries for your system (created using [BinaryBuilder](https://binarybuilder.org) & [julia](https://julialang.org)).

Follow these steps:
1) Download a zip file with the most recent release of `MAGEMin` (click on the green `Code` button @ the top of this page) and unzip it on your machine.
2) Open the `PlotPseudosection` graphical user interface from MATLAB (2020+). 
3) Follow the binary installation instructions (which requires you to install a recent [julia](https://www.julialang.org) version).
4) After this you are ready to get started, for example by pushing the `Start new computation` button. 

Note that the Matlab GUI is not maintained anymore as the primary phase diagram generator is now the Julia app!

## Manual compilation

if you wish, you can also compile MAGEMin yourself, which requires you to install these packages as well:
- MPICH (to allow parallel computations)
- LAPACKE (C version of LAPACK)
- NLopt (https://nlopt.readthedocs.io/)

Details and guidelines are given in the extended documentation: https://computationalthermodynamics.github.io/MAGEMin_C.jl/dev/

In addition, we make use of [uthash](https://troydhanson.github.io/uthash/) and [ketopt](https://github.com/attractivechaos/klib/blob/master/ketopt.h).


## Available thermodynamic datasets

The MAGEMin algorithm is general and can be used with any thermodynamic database that are hardcoded for speed reasons.

**Igneous database**

The hydrous mafic melting model of Holland et al. 2018 can be used to simulate the fractional crystallisation from a hydrous basalt to a felsic melt. 

The details of this thermodynamic solid solution and endmember database are:
- Holland et al., 2018 (see http://hpxeosandthermocalc.org)
- K2O-Na2O-CaO-FeO-MgO-Al2O3-SiO2-H2O-TiO2-O-Cr2O3 chemical system
- Equations of state for
	- Pure stoichiometric phases quartz (q), cristobalite (crst), tridymite (trd), coesite (coe), stishovite (stv), kyanite (ky), sillimanite (sill), andalusite (and), rutile (ru) and sphene (sph). 
	- Solution phases spinel (spl), biotite (bi), cordierite (cd), clinopyroxene (cpx), orthopyroxene (opx), epidote (ep), garnet (g), clino-amphibole  (amp), ilmenite (ilm), silicate melt (liq), muscovite (mu), olivine (ol), ternary feldspar (pl4T), and aqueous fluid (fl).

**Metapelite database**

The metapelitic model (extended with MnO, White et al., 2014) allows to compute the mineral assemblage from low temperature to supra-solidus conditions.

- Added March 2023, `MAGEMin v1.3.0` 
- White et al., 2014a, 2014b (see http://hpxeosandthermocalc.org)
- K2O-Na2O-CaO-FeO-MgO-Al2O3-SiO2-H2O-TiO2-O-MnO chemical system
- Equations of state for
	- Pure stoichiometric phases quartz (q), cristobalite (crst), tridymite (trd), coesite (coe), stishovite (stv), kyanite (ky), sillimanite (sill), andalusite (and), rutile (ru) and sphene (sph). 
	- Solution phases spinel (spl), biotite (bi), cordierite (cd), orthopyroxene (opx), epidote (ep), garnet (g),  ilmenite (ilm), silicate melt (liq), muscovite (mu),  ternary feldspar (pl4T), sapphirine (sa), staurolite (st), magnetite (mt), chlorite (chl), chloritoid (ctd) and margarite (ma).

**Ultramafic thermodynamic dataset**
- Added May 2023, `MAGEMin v1.3.2` 
- Evans & Frost, 2021 (see http://hpxeosandthermocalc.org)
- SiO2-Al2O3-MgO-FeO-O-H2O-S chemical system
- Equations of state for
	- Pure stoichiometric phases quartz (q), cristobalite (crst), tridymite (trd), coesite (coe), stishovite (stv), kyanite (ky), sillimanite (sill), pyrite (pyr)
	- Solution phases fluid (fluid), olivine (ol), brucite (br), antigorite (atg), garnet (g), talc (t), chlorite (chl), spinel (spi), orthopyroxene (opx), pyrrhotite (po) and anthophylite (anth)

**Metabasite thermodynamic dataset**

- Added October 2023, `MAGEMin v1.3.5`
- Green et al., 2016 (see http://hpxeosandthermocalc.org)
- K2O-Na2O-CaO-FeO-MgO-Al2O3-SiO2-H2O-TiO2-O chemical system
- Equations of state for
	- Pure stoichiometric phases quartz (q), cristobalite (crst), tridymite (trd), coesite (coe), stishovite (stv), kyanite (ky), sillimanite (sill), andalusite (and), rutile (ru) and sphene (sph). 
	- Solution phases spinel (sp), biotite (bi), orthopyroxene (opx), epidote (ep), garnet (g), ilmenite (ilm), silicate melt (liq), muscovite (mu),  ternary feldspar (pl4T), chlorite (chl), Omphacite(omph), Augite(aug) and clino-amphibole (amp).

**Extended Ultramafic thermodynamic dataset**
- Added September 2024, `MAGEMin v1.5.3`
- Na2O-CaO-FeO-MgO-Al2O3-SiO2-H2O-S-O chemical system
- Equations of state for
	- Pure stoichiometric phases quartz (q), cristobalite (crst), tridymite (trd), coesite (coe), stishovite (stv), kyanite (ky), sillimanite (sill), pyrite (pyr)
	- Solution phases from Evans & Frost, 2021: fluid (fluid), olivine (ol), brucite (br), antigorite (atg), garnet (g), talc (t), chlorite (chl), spinel (spi), orthopyroxene (opx), pyrrhotite (po) and anthophylite (anth)
	- Solution phases Green et al., 2016: ternary feldspar (pl4T), Augite (aug) and clino-amphibole (amp).

**Mantle dataset (Transition Zone into the Uppermost Lower Mantle)**
- Added October 2024,`MAGEMin v1.5.5`
- Holland et al., 2013 (see https://academic.oup.com/petrology/article/54/9/1901/1514886)
- Na2O–CaO–FeO–MgO–Al2O3–SiO2 (NCFMAS) system
- Equations of state for
	- Pure stoichiometric phases quartz (q), cristobalite (crst), tridymite (trd), coesite (coe), stishovite (stv), kyanite (ky), sillimanite (sill) and andalusite (and). 
	- Solution phases garnet (g), clinopyroxene (cpx), orthopyroxene (opx) and its high-P polymorph (hpx), olivine (ol), wadsleyite (wad), ringwoodite (ring), akimotoite (ak), MgSi-perovskite (mpv), CaSi–perovskite (cpv), cf, nal, corundum (cor) and ferropericlase (fp)

Please keep in mind that the datasets are only calibrated for a limited range of `P`,`T` and `bulk rock` conditions. If you go too far outside those ranges, `MAGEMin` (or most other thermodynamic software packages for that matter) may not converge or give bogus results. 
Developing new, more widely applicable, thermodynamic datasets is a huge research topic, which will require funding to develop the models themselves, as well as to perform targeted experiments to calibrate those models.

**Igneous dataset (update and correction)**
- Added December 2024,`MAGEMin v1.6.2`
- Green et al., 2025, corrected from Holland et al., 2018 (see http://hpxeosandthermocalc.org)
- K2O-Na2O-CaO-FeO-MgO-Al2O3-SiO2-H2O-TiO2-O-Cr2O3 chemical system
- Equations of state for
	- Pure stoichiometric phases quartz (q), cristobalite (crst), tridymite (trd), coesite (coe), stishovite (stv), kyanite (ky), sillimanite (sill), andalusite (and), rutile (ru) and sphene (sph). 
	- Solution phases spinel (spl), biotite (bi), cordierite (cd), clinopyroxene (cpx), orthopyroxene (opx), epidote (ep), garnet (g), clino-amphibole  (amp), ilmenite (ilm), silicate melt (liq), muscovite (mu), olivine (ol), ternary feldspar (pl4T), and aqueous fluid (fl).

**Igneous alkaline dry dataset**
- Added December 2024,`MAGEMin v1.6.2`
- Weller et l., 2024 (see doi:10.1093/petrology/egae098)
- K2O-Na2O-CaO-FeO-MgO-Al2O3-SiO2-TiO2-O-Cr2O3 chemical system
- Equations of state for
	- Pure stoichiometric phases quartz (q), cristobalite (crst), tridymite (trd), coesite (coe), stishovite (stv), kyanite (ky), sillimanite (sill), andalusite (and), rutile (ru) and sphene (sph). 
	- Solution phases spinel (spl), clinopyroxene (cpx), orthopyroxene (opx), garnet (g), ilmenite (ilm), silicate melt (liq), olivine (ol), ternary feldspar (pl4T), Nepheline (ness), Kalsilite (kals), Leucite (lct) and Melilite (mel).


## Documentation
Full support to install and use MAGEMin is available [here](https://computationalthermodynamics.github.io/MAGEMin_C.jl/dev/).

## Citation
An open-acces paper describing the methodology is:

- Riel N., Kaus B.J.P., Green E.C.R., Berlie N., (2022) MAGEMin, an Efficient Gibbs Energy Minimizer: Application to Igneous Systems. *Geochemistry, Geophysics, Geosystems* 23, e2022GC010427 [https://doi.org/10.1029/2022GC010427](https://doi.org/10.1029/2022GC010427)

If you use the software, we really appreciate if you cite this study. We also appreciate stars (see the top of this page). 

## Contributing
You are very welcome to request new features and point out bugs by opening an issue (top of this page). You can also help by adding features and creating a pull request. 

## Funding
Development of this software package was funded by the European Research Council under grant ERC CoG #771143 - [MAGMA](https://magma.uni-mainz.de).

<img src="./pics/MAGMA_Logo.png" alt="drawing" width="480" alt="centered image"/>

## References

- Green, ECR, Holland, TJB, Powell, R, Weller, OM, & Riel, N (2025).
XXXXXX Journal of Petrology, doi: XXXXXX

- Weller, OM, Holland, TJB, Soderman, CR, Green, ECR, Powell, R, 
Beard, CD & Riel, N (2024). New Thermodynamic Models for Anhydrous
Alkaline-Silicate Magmatic Systems. Journal of Petrology, 65,
doi: 10.1093/petrology/egae098

- Holland, TJB, Green, ECR & Powell, R (2022). A thermodynamic model
for feldspars in KAlSi3O8-NaAlSi3O8-CaAl2Si2O8 for mineral 
equilibrium calculations. Journal of Metamorphic Geology, 40, 587-600, 
doi: 10.1111/jmg.12639

- Tomlinson, EL & Holland, TJB (2021). A Thermodynamic Model for the
Subsolidus Evolution and Melting of Peridotite. Journal of Petrology,
62, doi: 10.1093/petrology/egab012

- Holland, TJB, Green, ECR & Powell, R (2018). Melting of Peridotites
through to Granites: A Simple Thermodynamic Model in the System
KNCFMASHTOCr. Journal of Petrology, 59, 881-900, 
doi: 10.1093/petrology/egy048

- Green, ECR, White, RW, Diener, JFA, Powell, R, Holland, TJB & 
Palin, RM (2016). Activity-composition relations for the calculation
of partial melting equilibria in metabasic rocks. Journal of 
Metamorphic Geology, 34, 845-869, doi: 10.1111/jmg.12211

- White, RW, Powell, R, Holland, TJB, Johnson, TE & Green, ECR (2014). 
New mineral activity-composition relations for thermodynamic calculations 
in metapelitic systems. Journal of Metamorphic Geology, 32, 261-286,
doi: 10.1111/jmg.12071

- Holland, TJB & Powell, RW (2011). An improved and extended internally
consistent thermodynamic dataset for phases of petrological interest,
involving a new equation of state for solids. Journal of Metamorphic 
Geology, 29, 333-383, doi: 10.1111/j.1525-1314.2010.00923.x


<!-- ## Extended informations

G0HSC(T,P)  =  Tr·Selements(Tr,Pr) +  G0SUPCRT(T,P)
Tr = 298.15K
HSC -> thermocalc
SUPCRT -> DEW

Dean, John A. Lange’s Handbook of Chemistry, 12th ed.; McGraw-Hill: New York, New York, 1979; p 9-4–9-94. -->
