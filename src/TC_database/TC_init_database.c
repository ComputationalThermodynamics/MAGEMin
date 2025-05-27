/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **   Project      : MAGEMin
 **   License      : GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007
 **   Developers   : Nicolas Riel, Boris Kaus
 **   Contributors : Dominguez, H., Assunção J., Green E., Berlie N., and Rummel L.
 **   Organization : Institute of Geosciences, Johannes-Gutenberg University, Mainz
 **   Contact      : nriel[at]uni-mainz.de, kaus[at]uni-mainz.de
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <complex.h> 

#include "../uthash.h"
#include "../MAGEMin.h"
#include "../initialize.h"
#include "../all_endmembers.h"
#include "TC_init_database.h"


/* 	Note that "ecp" is the electrochemical potential "component" used for partitioning gibbs energy when using aqueous fluid model 	*/
oxide_data oxide_info = {
	15,						/* number of endmembers */
	{"SiO2"	,"Al2O3","CaO"	,"MgO"	,"FeO"	,"K2O"	,"Na2O"	,"TiO2"	,"O"	,"MnO"	,"Cr2O3","H2O"	,"CO2"	,"S"	,"Cl"	,"ecp"		},
	{"Si"	,"Al"	,"Ca"	,"Mg"	,"Fe"	,"K"	,"Na"	,"Ti"	,"O"	,"Mn"	,"Cr"	,"H"	,"C"	,"S"	,"Cl"	,"ecp"		},
	{60.08  ,101.96 ,56.08  ,40.30  ,71.85  ,94.2   ,61.98  ,79.88  ,16.0   ,70.94	,151.99 ,18.015	,44.01	, 32.06	,35.453	,0.0		},
	{3.0	,5.0	,2.0	,2.0	,2.0	,3.0	,3.0	,3.0	,1.0	,2.0 	,5.0	,3.0	,3.0	, 1.0	,1.0	,0.0		},
	{66.7736,108.653,42.9947,40.3262,38.7162,69.1514,61.1729,70.3246,30.5827,40.1891,106.9795,69.5449,62.8768,9.5557,33.2556,0.0		},
	{2,3,1,1,1,1,1,2,1,1,3,1,2,0,0,0}, //opo
	// {3,5,2,2,2,3,3,3,0,2,5,3,3,1,1,0},  //apo
	{1,2,1,1,1,2,2,1,1,1,2,2,1,1,1,0}	//ape
	// for the standard molar entropy the values are already normalized by the reference temperature = 298.15K (25°C) and expressed in kJ
};


metapelite_aq_dataset metapelite_aq_db = {
	62,							/* Endmember default dataset number */
	11,							/* number of oxides */			
	26,							/* number of pure phases */
	18,							/* number of solution phases */
	{"SiO2"	,"Al2O3","CaO"	,"MgO"	,"FeO"	,"K2O"	,"Na2O"	,"TiO2"	,"O"	,"MnO"	,"H2O"											},
	{"q"	,"crst"	,"trd"	,"coe"	,"stv"	,"ky"	,"sill"	,"and"	,"ru"	,"sph"	,"O2"  ,"H2O"	,"zo"	,
	"qfm"	,"mw"	,"qif"	,"nno"	,"hm"	,"iw"	,"cco"	,"aH2O"	, "aO2"	,"aMgO"	,"aFeO"	,"aAl2O3"		,"aTiO2"								},
	{"liq"	,"fsp"	,"bi"	,"g"	,"ep"	,"ma"	,"mu"	,"opx"	,"sa"	,"cd"	,"st"	,"chl"	,"ctd"	,"sp"  ,"mt"  ,"ilm"  ,"ilmm"  ,"aq17"},
	
	{1		,1		,1		,1		,1		,1		,1		,1		,1 		,1 		,1 		,1 		,1 		,1 		,1 		,1 		,1 		,1		},  // allow solvus?
	{3150	,231 	,981	,757	,110 	,1875	,1877	,1277	,230	,343	,540	,2270	,216	,405 	,87 	,130 	,1430 	,1		},  // # of pseudocompound
	{0.15	,0.049	,0.19	,0.19	,0.049	,0.19	,0.19	,0.249	,0.19	,0.145	,0.19	,0.249	,0.19	,0.124 	,0.099 	,0.09 	,0.099 	,1.0	},  // discretization step

	6.0, 						/* max dG under which a phase is considered to be reintroduced  					*/
	473.15,						/* max temperature above which PGE solver is active 								*/
	773.15,						/** minimum temperature above which melt is considered 								*/

	4,							/** number of inner PGE iterations, this has to be made mass or dG dependent 		*/
	0.025,						/** maximum mol% phase change during one PGE iteration in wt% 						*/
	2.5,						/** maximum delta_G of reference change during PGE 									*/
	1.0,						/** maximum update factor during PGE under-relax < 0.0, over-relax > 0.0 	 		*/

	1e-1,						/** merge instances of solution phase if norm < val 								*/
	1e-4,						/** fraction of solution phase when re-introduced 									*/
	1e-6						/** objective function tolerance 				 									*/
};



metapelite_dataset metapelite_db = {
	62,							/* Endmember default dataset number */
	11,							/* number of oxides */			
	26,							/* number of pure phases */
	17,							/* number of solution phases */
	{"SiO2"	,"Al2O3","CaO"	,"MgO"	,"FeO"	,"K2O"	,"Na2O"	,"TiO2"	,"O"	,"MnO"	,"H2O"											},
	{"q"	,"crst"	,"trd"	,"coe"	,"stv"	,"ky"	,"sill"	,"and"	,"ru"	,"sph"	,"O2"  ,"H2O"	,"zo"	,
	"qfm"	,"mw"	,"qif"	,"nno"	,"hm"	,"iw"	,"cco"	,"aH2O"	, "aO2"	,"aMgO"	,"aFeO"	,"aAl2O3"		,"aTiO2"								},
	{"liq"	,"fsp"	,"bi"	,"g"	,"ep"	,"ma"	,"mu"	,"opx"	,"sa"	,"cd"	,"st"	,"chl"	,"ctd"	,"sp"  ,"mt"  ,"ilm"  ,"ilmm" 	,"aq17" },
	
	{1		,1		,1		,1		,1		,1		,1		,1		,1 		,1 		,1 		,1 		,1 		,1 		,1 		,1 		,1 		,1		},  // allow solvus?
	{3150	,231 	,981	,757	,110 	,1875	,1877	,1277	,230	,343	,540	,2270	,216	,410 	,87 	,130 	,1430 	,1		},  // # of pseudocompound
	{0.15	,0.049	,0.19	,0.19	,0.049	,0.19	,0.19	,0.249	,0.19	,0.145	,0.19	,0.249	,0.19	,0.02 	,0.025 	,0.025 	,0.025 	,1.0	},  // discretization step

	6.0, 						/* max dG under which a phase is considered to be reintroduced  					*/
	473.15,						/* max temperature above which PGE solver is active 								*/
	773.15,						/** minimum temperature above which melt is considered 								*/

	4,							/** number of inner PGE iterations, this has to be made mass or dG dependent 		*/
	0.025,						/** maximum mol% phase change during one PGE iteration in wt% 						*/
	2.5,						/** maximum delta_G of reference change during PGE 									*/
	1.0,						/** maximum update factor during PGE under-relax < 0.0, over-relax > 0.0 	 		*/

	1e-1,						/** merge instances of solution phase if norm < val 								*/
	1e-4,						/** fraction of solution phase when re-introduced 									*/
	1e-6						/** objective function tolerance 				 									*/
};

metabasite_dataset metabasite_db = {
	62,							/* Endmember default dataset number */
	10,							/* number of oxides */			
	26,							/* number of pure phases */
	17,							/* number of solution phases */
	{"SiO2"	,"Al2O3","CaO"	,"MgO"	,"FeO"	,"K2O"	,"Na2O"	,"TiO2"	,"O"	,"H2O"													},
	{"q"	,"crst"	,"trd"	,"coe"	,"law"	,"ky"	,"sill"	,"and"	,"ru"	,"sph" 	,"ab"	,"H2O"	,"zo"	,
	"qfm"	,"mw" 	,"qif"	,"nno"	,"hm"	,"iw"	,"cco"	,"aH2O"	,"aO2"	,"aMgO"	,"aFeO"	,"aAl2O3"		,"aTiO2"								},
	{"sp"	,"opx"	,"fsp"	,"liq"	,"mu"	,"ilmm"	,"ilm"	,"ol"	,"amp"	,"ep"	,"g"	,"chl"	,"bi"	,"dio"	,"aug"	,"abc"  ,"spl"			},
	
	{1		,1		,1		,1		,1		,1		,1		,1		,1		,1 		,1 		,1 		,1 		,1 		,1 		,1 		,1				},  // allow solvus?
	{939	,1731 	,231	,3507	,4539 	,298	,422	,11		,7673	,110	,217	,3980	,1098	,1811	,2398	,21 	,196			},  // # of pseudocompound
	{0.04	,0.19	,0.02	,0.15	,0.049	,0.09	,0.049	,0.04	,0.15	,0.049	,0.19	,0.19	,0.125	,0.10	,0.10	,0.049 	,0.04			},  // discretization step

	6.0, 						/* max dG under which a phase is considered to be reintroduced  					*/
	473.15,						/* max temperature above which PGE solver is active 								*/
	873.15,						/** minimum temperature above which melt is considered 								*/

	4,							/** number of inner PGE iterations, this has to be made mass or dG dependent 		*/
	0.025,						/** maximum mol% phase change during one PGE iteration in wt% 						*/
	2.5,						/** maximum delta_G of reference change during PGE 									*/
	1.0,						/** maximum update factor during PGE under-relax < 0.0, over-relax > 0.0 	 		*/

	1e-1,						/** merge instances of solution phase if norm < val 								*/
	1e-4,						/** fraction of solution phase when re-introduced 									*/
	1e-6						/** objective function tolerance 				 									*/
};


metabasite_ext_dataset metabasite_ext_db = {
	62,							/* Endmember default dataset number */
	10,							/* number of oxides */			
	26,							/* number of pure phases */
	19,							/* number of solution phases */
	{"SiO2"	,"Al2O3","CaO"	,"MgO"	,"FeO"	,"K2O"	,"Na2O"	,"TiO2"	,"O"	,"H2O"													},
	{"q"	,"crst"	,"trd"	,"coe"	,"law"	,"ky"	,"sill"	,"and"	,"ru"	,"sph" 	,"ab"	,"H2O"	,"zo"	,
	"qfm"	,"mw" 	,"qif"	,"nno"	,"hm"	,"iw"	,"cco"	,"aH2O"	,"aO2"	,"aMgO"	,"aFeO"	,"aAl2O3"		,"aTiO2"								},
	{"sp"	,"opx"	,"fsp"	,"liq"	,"mu"	,"ilmm"	,"ilm"	,"ol"	,"amp"	,"ep"	,"g"	,"chl"	,"bi"	,"dio"	,"aug"	,"abc"  ,"spl"	,"ta"	,"oamp"		},
	
	{1		,1		,1		,1		,1		,1		,1		,1		,1		,1 		,1 		,1 		,1 		,1 		,1 		,1 		,1		,1 		,1			},  // allow solvus?
	{939	,1731 	,231	,3507	,4539 	,298	,422	,11		,7673	,110	,217	,3980	,1098	,1811	,2398	,21 	,196	,760	,3648		},  // # of pseudocompound
	{0.04	,0.19	,0.02	,0.15	,0.049	,0.09	,0.049	,0.04	,0.15	,0.049	,0.19	,0.19	,0.125	,0.10	,0.10	,0.049 	,0.04	,0.09	,0.2		},  // discretization step

	6.0, 						/* max dG under which a phase is considered to be reintroduced  					*/
	473.15,						/* max temperature above which PGE solver is active 								*/
	873.15,						/** minimum temperature above which melt is considered 								*/

	4,							/** number of inner PGE iterations, this has to be made mass or dG dependent 		*/
	0.025,						/** maximum mol% phase change during one PGE iteration in wt% 						*/
	2.5,						/** maximum delta_G of reference change during PGE 									*/
	1.0,						/** maximum update factor during PGE under-relax < 0.0, over-relax > 0.0 	 		*/

	1e-1,						/** merge instances of solution phase if norm < val 								*/
	1e-4,						/** fraction of solution phase when re-introduced 									*/
	1e-6						/** objective function tolerance 				 									*/
};

igneous_dataset igneous_db = {
	636,						/* Endmember default dataset number */
	11,							/* number of oxides */			
	26,							/* number of pure phases */
	16,							/* number of solution phases */
	{"SiO2"	,"Al2O3","CaO"	,"MgO"	,"FeO"	,"K2O"	,"Na2O"	,"TiO2"	,"O"	,"Cr2O3","H2O"														},
	{"ne"	,"q"	,"crst"	,"trd"	,"coe"	,"stv"	,"ky"	,"sill"	,"and"	,"ru"	,"sph"	,"O2"	,"H2O"	,
	"qfm"	,"mw"	,"qif"	,"nno"	,"hm"	,"iw"	,"cco"	,"aH2O"	, "aO2"	,"aMgO"	,"aFeO"	,"aAl2O3"		,"aTiO2"							},
	{1		, 1		, 1		, 1		, 1		, 1		, 1		, 1		, 1		, 1		, 1		, 1		, 0		,
	 1		, 1		, 1		, 1		, 1		, 1		, 1		, 1		, 1		, 1		, 1		, 1		, 1		,									},
	{"spl"	,"bi"	,"cd"	,"cpx"	,"ep"	,"g"	,"amp"	,"ilm"	,"liq"	,"ol"	,"opx"	,"fsp"	,"fl"	,"mu"	,"fper"	,"chl"				},
	
	{1		,1		,1		,1		,1		,1		,1		,1		,1 		,1 		,1 		,1 		,1		,1		,1		,1					}, // allow solvus?
	{1523	,3554	,121	,4127	,210	,2450	,5499	,1670	,3088	,381	,3413	,231	,2		,2376	,20		,3980				}, // # of pseudocompound
	{0.2	,0.124	,0.098	,0.20	,0.049	,0.145	,0.33	,0.10	,0.15	,0.04	,0.249	,0.049	,1.0 	,0.198	,0.05	,0.19				}, // discretization step

	6.0, 						/** max dG under which a phase is considered to be reintroduced  					*/
	673.15,						/** max temperature above which PGE solver is active 								*/
	773.15,						/** minimum temperature above which melt is considered 								*/

	4,							/** number of inner PGE iterations, this has to be made mass or dG dependent 		*/
	0.025,						/** maximum mol% phase change during one PGE iteration in wt% 						*/
	2.5,						/** maximum delta_G of reference change during PGE 									*/
	1.0,						/** maximum update factor during PGE under-relax < 0.0, over-relax > 0.0 	 		*/

	2e-1,						/** merge instances of solution phase if norm < val 								*/
	1e-4,						/** fraction of solution phase when re-introduced 									*/
	1e-6						/** objective function tolerance 				 									*/
};

igneous_igad_dataset igneous_igad_db = {
	636,						/* number of endmembers */
	10,							/* number of oxides */			
	24,							/* number of pure phases */
	12,							/* number of solution phases */
	{"SiO2"	,"Al2O3","CaO"	,"MgO"	,"FeO"	,"K2O"	,"Na2O"	,"TiO2"	,"O"	,"Cr2O3"														},
	{"q"	,"crst"	,"trd"	,"coe"	,"stv"	,"ky"	,"sill"	,"and"	,"ru"	,"sph"	,"O2" 	,
	"qfm"	,"mw"	,"qif"	,"nno"	,"hm"	,"iw"	,"cco"	,"aH2O"	, "aO2"	,"aMgO"	,"aFeO"	,"aAl2O3"		,"aTiO2"								},
	{"spl"	,"cpx"	,"g"	,"ilm"	,"liq"	,"ol"	,"opx"	,"fsp"	,"lct"	,"mel"	,"nph"	,"kals"		},
	
	{1		,1		,1		,1		,1 		,1 		,1 		,1 		,1		,1		,1		,1			}, // allow solvus?
	{3318	,4128	,2450	,1671	,4912	,381	,3413	,861	,21		,270	,1210	,21			}, // # of pseudocompound
	{0.195	,0.249	,0.145	,0.05	,0.1	,0.098	,0.249	,0.0249	,0.049	,0.19	,0.149 	,0.049		}, // discretization step

	6.0, 						/** max dG under which a phase is considered to be reintroduced  					*/
	673.15,						/** max temperature above which PGE solver is active 								*/
	873.15,						/** minimum temperature above which melt is considered 								*/

	4,							/** number of inner PGE iterations, this has to be made mass or dG dependent 		*/
	0.025,						/** maximum mol% phase change during one PGE iteration in wt% 						*/
	2.5,						/** maximum delta_G of reference change during PGE 									*/
	1.0,						/** maximum update factor during PGE under-relax < 0.0, over-relax > 0.0 	 		*/

	1e-1,						/** merge instances of solution phase if norm < val 								*/
	1e-4,						/** fraction of solution phase when re-introduced 									*/
	1e-6						/** objective function tolerance 				 									*/
};

ultramafic_dataset ultramafic_db = {
	633,						/* Endmember default dataset number */
	7,							/* number of oxides */			
	23,							/* number of pure phases */
	12,							/* number of solution phases */
	{"SiO2"	,"Al2O3","MgO"	,"FeO"	,"O"	,"H2O"	,"S"												},
	{"q"	,"crst"	,"trd"	,"coe"	,"stv"	,"ky"	,"sill"	,"and"	,"pyr"	,"O2"  	,
	"qfm"	,"qif"	,"nno"	,"hm"	,"mw"	,"iw"	,"cco"	,"aH2O"	, "aO2"	,"aMgO"	,"aFeO"	,"aAl2O3"		,"aTiO2"	},
	{"fl"	,"ol"  ,"br"	,"ch"	,"atg"	,"g"	,"ta"	,"chl"	,"spi"	,"opx"	,"po"	,"anth" 	},
	
	{1		,1		,1		,1		,1		,1		,1		,1		,1 		,1 		,1		,1			},  // allow solvus?
	{11  	,10  	,11 	,10 	,489 	,10  	,985 	,2691	,100	,196	,10		,274		},  // No. of pseudocompound
	{0.001	,0.1	,0.1	,0.1	,0.19	,0.1	,0.19	,0.19	,0.1	,0.19	,0.1	,0.249		},  // discretization step

	4.0, 						/* max dG under which a phase is considered to be reintroduced  					*/
	473.15,						/* max temperature above which PGE solver is active 								*/
	873.15,						/** minimum temperature above which melt is considered 								*/

	4,							/** number of inner PGE iterations, this has to be made mass or dG dependent 		*/
	0.025,						/** maximum mol% phase change during one PGE iteration in wt% 						*/
	2.5,						/** maximum delta_G of reference change during PGE 									*/
	1.0,						/** maximum update factor during PGE under-relax < 0.0, over-relax > 0.0 	 		*/

	1e-1,						/** merge instances of solution phase if norm < val 								*/
	1e-4,						/** fraction of solution phase when re-introduced 									*/
	1e-6						/** objective function tolerance 				 									*/
};

ultramafic_ext_dataset ultramafic_ext_db = {
	633,						/* Endmember default dataset number */
	9,							/* number of oxides */			
	23,							/* number of pure phases */
	15,							/* number of solution phases */
	{"SiO2"	,"Al2O3","MgO"	,"FeO"	,"O"	,"H2O"	,"S"	,"CaO"	,"Na2O"													},
	{"q"	,"crst"	,"trd"	,"coe"	,"stv"	,"ky"	,"sill"	,"and"	,"pyr"	,"O2"  	,
	"qfm"	,"qif"	,"nno"	,"hm"	,"mw"	,"iw"	,"cco"	,"aH2O"	, "aO2"	,"aMgO"	,"aFeO"	,"aAl2O3"		,"aTiO2"		},
	{"fl"	,"ol"  ,"br"	,"ch"	,"atg"	,"g"	,"ta"	,"chl"	,"spi"	,"opx"	,"po"	,"anth"	,"pl4tr","amp"	,"aug"	},
	
	{1		,1		,1		,1		,1		,1		,1		,1		,1 		,1 		,1		,1		,1		,1		,1		},  // allow solvus?
	{11  	,10  	,11 	,10 	,489 	,10  	,985 	,2691	,100	,196	,10		,274	,21		,5033	,3036	},  // No. of pseudocompound
	{0.001	,0.1	,0.1	,0.1	,0.19	,0.1	,0.19	,0.19	,0.1	,0.19	,0.1	,0.249	,0.049	,0.24	,0.24	},  // discretization step

	4.0, 						/* max dG under which a phase is considered to be reintroduced  					*/
	473.15,						/* max temperature above which PGE solver is active 								*/
	873.15,						/** minimum temperature above which melt is considered 								*/

	4,							/** number of inner PGE iterations, this has to be made mass or dG dependent 		*/
	0.025,						/** maximum mol% phase change during one PGE iteration in wt% 						*/
	2.5,						/** maximum delta_G of reference change during PGE 									*/
	1.0,						/** maximum update factor during PGE under-relax < 0.0, over-relax > 0.0 	 		*/

	1e-1,						/** merge instances of solution phase if norm < val 								*/
	1e-4,						/** fraction of solution phase when re-introduced 									*/
	1e-6						/** objective function tolerance 				 									*/
};


mantle_dataset mantle_db = {
	633,						/* Endmember default dataset number */
	6,							/* number of oxides */			
	8,							/* number of pure phases */
	14,							/* number of solution phases */
	{"SiO2"	,"Al2O3","CaO"	,"MgO"	,"FeO"	,"Na2O"																},
	{"q"	,"crst"	,"trd"	,"coe"	,"stv"	,"ky"	,"sill"	,"and"												},
	{"g"	,"fp"  ,"mpv"	,"cpv"	,"crn"	,"cf"	,"nal"	,"aki"	,"ol"	,"wad"	,"ring"	,"cpx"	,"opx"	,"hpx"		},
	
	{1		,1		,1		,1		,1		,1		,1		,1		,1		,1 		,1 		,1		,1		,1			},  // allow solvus?
	{800  	,21  	,336 	,336 	,144 	,1258 	,1672  	,144 	,21		,21		,21		,1005	,342	,343		},  // No. of pseudocompound
	{0.19	,0.049	,0.19	,0.10	,0.09	,0.19	,0.24	,0.9	,0.049	,0.049	,0.049	,0.10	,0.19	,0.10		},  // discretization step

	4.0, 						/* max dG under which a phase is considered to be reintroduced  					*/
	473.15,						/* max temperature above which PGE solver is active 								*/
	873.15,						/** minimum temperature above which melt is considered 								*/

	4,							/** number of inner PGE iterations, this has to be made mass or dG dependent 		*/
	0.025,						/** maximum mol% phase change during one PGE iteration in wt% 						*/
	2.5,						/** maximum delta_G of reference change during PGE 									*/
	1.0,						/** maximum update factor during PGE under-relax < 0.0, over-relax > 0.0 	 		*/

	1e-1,						/** merge instances of solution phase if norm < val 								*/
	1e-4,						/** fraction of solution phase when re-introduced 									*/
	1e-6						/** objective function tolerance 				 									*/
};


metapelite_dataset_ext metapelite_ext_db = {
	62,							/* Endmember default dataset number */
	13,							/* number of oxides */			
	31,							/* number of pure phases */
	24,							/* number of solution phases */
	{"SiO2"	,"Al2O3","CaO"	,"MgO"	,"FeO"	,"K2O"	,"Na2O"	,"TiO2"	,"O"	,"MnO"	,"H2O"	,"CO2"	,"S"										},
	{"q"	,"crst"	,"trd"	,"coe"	,"stv"	,"ky"	,"sill"	,"and"	,"ru"	,"sph"	,"O2" 	,"pyr"	,"gph"	,"law"	,"zo" ,"prl"  ,"mpm"   ,"pre"	,
	"qfm"	,"mw"	,"qif"	,"nno"	,"hm"	,"iw" 	,"cco"	,"aH2O"	, "aO2"	,"aMgO"	,"aFeO"	,"aAl2O3"		,"aTiO2"											},
	{"liq"	,"fsp"	,"bi"	,"g"	,"ep"	,"ma"	,"mu"	,"opx"	,"sa"	,"cd"	,"st"	,"chl"	,"ctd"	,"sp"  ,"mt"  ,"ilm"  ,"ilmm"  ,"occm"	,"fl"	,"po"	,"dio"	,"aug"	,"amp"	,"oamp"		},
	
	{1		,1		,1		,1		,1		,1		,1		,1		,1 		,1 		,1 		,1 		,1 		,1 		,1 		,1 		,1 		,1		,1		,1		,1		,1		,1		,1			},  // allow solvus?
	{3150	,231 	,981	,758	,110 	,1875	,1877	,1277	,230	,343	,540	,2270	,216	,407 	,87 	,131 	,1430 	,352	,12		,10		,1810	,2396,	7669	,3648		},  // # of pseudocompound
	{0.24	,0.049	,0.125	,0.049	,0.049	,0.19	,0.19	,0.249	,0.19	,0.145	,0.19	,0.249	,0.19	,0.124 	,0.099 	,0.09 	,0.099 	,0.10	,0.09	,0.1	,0.10	,0.10,	0.075	,0.2		},  // discretization step

	6.0, 						/* max dG under which a phase is considered to be reintroduced  					*/
	473.15,						/* max temperature above which PGE solver is active 								*/
	873.15,						/** minimum temperature above which melt is considered 								*/

	4,							/** number of inner PGE iterations, this has to be made mass or dG dependent 		*/
	0.025,						/** maximum mol% phase change during one PGE iteration in wt% 						*/
	2.5,						/** maximum delta_G of reference change during PGE 									*/
	1.0,						/** maximum update factor during PGE under-relax < 0.0, over-relax > 0.0 	 		*/

	1e-1,						/** merge instances of solution phase if norm < val 								*/
	1e-4,						/** fraction of solution phase when re-introduced 									*/
	1e-6						/** objective function tolerance 				 									*/
};


/* Function to allocate the memory of the data to be used/saved during PGE iterations */
global_variable global_variable_TC_init( 	global_variable  	 gv,
											bulk_info 			*z_b 	){
	int i, j;

	/* load database */
	if (gv.EM_database == 0){
		metapelite_dataset db 	= metapelite_db;
		if (gv.EM_dataset == -1){
			gv.EM_dataset = db.ds_version;	
		}
		gv.len_pp   		= db.n_pp;		
		gv.len_ss  			= db.n_ss;
		gv.len_ox  			= db.n_ox;

		gv.PC_df_add		= db.PC_df_add;					/** min value of df under which the PC is added 									*/
		gv.solver_switch_T  = db.solver_switch_T;
		gv.min_melt_T       = db.min_melt_T ;				/** minimum temperature above which melt is considered 								*/

		gv.inner_PGE_ite    = db.inner_PGE_ite;				/** number of inner PGE iterations, this has to be made mass or dG dependent 		*/
		gv.max_n_phase  	= db.max_n_phase;				/** maximum mol% phase change during one PGE iteration in wt% 						*/
		gv.max_g_phase  	= db.max_g_phase;				/** maximum delta_G of reference change during PGE 									*/
		gv.max_fac          = db.max_fac;					/** maximum update factor during PGE under-relax < 0.0, over-relax > 0.0 	 		*/

		gv.merge_value		= db.merge_value;				/** merge instances of solution phase if norm < val 								*/
		gv.re_in_n          = db.re_in_n;					/** fraction of phase when being reintroduce.  										*/
		gv.obj_tol 			= db.obj_tol;

		/* alloc */
		gv.ox 				= malloc (gv.len_ox * sizeof(char*)		);
		for (i = 0; i < gv.len_ox; i++){
			gv.ox[i] 		= malloc(20 * sizeof(char));	
			strcpy(gv.ox[i],db.ox[i]);
		}

		gv.PP_list 			= malloc (gv.len_pp * sizeof(char*)		);
		for (i = 0; i < (gv.len_pp); i++){	
			gv.PP_list[i] 	= malloc(20 * sizeof(char));
			strcpy(gv.PP_list[i],db.PP[i]);
		}
		gv.SS_list 			= malloc ((gv.len_ss) * sizeof (char*)	);
		gv.n_SS_PC     		= malloc ((gv.len_ss) * sizeof (int) 	);
		gv.verifyPC  		= malloc ((gv.len_ss) * sizeof (int) 	);
		gv.SS_PC_stp     	= malloc ((gv.len_ss) * sizeof (double) );
		for (i = 0; i < gv.len_ss; i++){ 
			gv.SS_list[i] 	= malloc(20 * sizeof(char)				);
			strcpy(gv.SS_list[i],db.SS[i]);
			gv.verifyPC[i]  = db.verifyPC[i]; 
			gv.n_SS_PC[i] 	= db.n_SS_PC[i]; 
			gv.SS_PC_stp[i] = db.SS_PC_stp[i]; 	
		}
		gv.act_PP     		= malloc ((gv.len_pp) * sizeof (int) 	);
		for (i = 0; i < gv.len_pp; i++){ gv.act_PP[i] = 1; 			};
	}
	else if (gv.EM_database == 1){
		metabasite_dataset db 	= metabasite_db;
		if (gv.EM_dataset == -1){
			gv.EM_dataset = db.ds_version;	
		}
		gv.len_pp   		= db.n_pp;		
		gv.len_ss  			= db.n_ss;
		gv.len_ox  			= db.n_ox;

		gv.PC_df_add		= db.PC_df_add;					/** min value of df under which the PC is added 									*/
		gv.solver_switch_T  = db.solver_switch_T;
		gv.min_melt_T       = db.min_melt_T ;				/** minimum temperature above which melt is considered 								*/

		gv.inner_PGE_ite    = db.inner_PGE_ite;				/** number of inner PGE iterations, this has to be made mass or dG dependent 		*/
		gv.max_n_phase  	= db.max_n_phase;				/** maximum mol% phase change during one PGE iteration in wt% 						*/
		gv.max_g_phase  	= db.max_g_phase;				/** maximum delta_G of reference change during PGE 									*/
		gv.max_fac          = db.max_fac;					/** maximum update factor during PGE under-relax < 0.0, over-relax > 0.0 	 		*/

		gv.merge_value		= db.merge_value;				/** merge instances of solution phase if norm < val 								*/
		gv.re_in_n          = db.re_in_n;					/** fraction of phase when being reintroduce.  										*/
		gv.obj_tol 			= db.obj_tol;

		/* alloc */
		gv.ox 				= malloc (gv.len_ox * sizeof(char*)		);
		for (i = 0; i < gv.len_ox; i++){
			gv.ox[i] 		= malloc(20 * sizeof(char));	
			strcpy(gv.ox[i],db.ox[i]);
		}

		gv.PP_list 			= malloc (gv.len_pp * sizeof(char*)		);
		for (i = 0; i < (gv.len_pp); i++){	
			gv.PP_list[i] 	= malloc(20 * sizeof(char));
			strcpy(gv.PP_list[i],db.PP[i]);
		}

		gv.SS_list 			= malloc ((gv.len_ss) * sizeof (char*)	);
		gv.n_SS_PC     		= malloc ((gv.len_ss) * sizeof (int) 	);

		gv.verifyPC  		= malloc ((gv.len_ss) * sizeof (int) 	);
		gv.SS_PC_stp     	= malloc ((gv.len_ss) * sizeof (double) );
		for (i = 0; i < gv.len_ss; i++){ 
			gv.SS_list[i] 	= malloc(20 * sizeof(char)				);
			strcpy(gv.SS_list[i],db.SS[i]);
			gv.verifyPC[i]  = db.verifyPC[i]; 
			gv.n_SS_PC[i] 	= db.n_SS_PC[i]; 
			gv.SS_PC_stp[i] = db.SS_PC_stp[i]; 	
		}
		gv.act_PP     		= malloc ((gv.len_pp) * sizeof (int) 	);
		for (i = 0; i < gv.len_pp; i++){ gv.act_PP[i] = 1; 			};
	}
	else if (gv.EM_database == 11){
		metabasite_ext_dataset db 	= metabasite_ext_db;
		if (gv.EM_dataset == -1){
			gv.EM_dataset = db.ds_version;	
		}
		gv.len_pp   		= db.n_pp;		
		gv.len_ss  			= db.n_ss;
		gv.len_ox  			= db.n_ox;

		gv.PC_df_add		= db.PC_df_add;					/** min value of df under which the PC is added 									*/
		gv.solver_switch_T  = db.solver_switch_T;
		gv.min_melt_T       = db.min_melt_T ;				/** minimum temperature above which melt is considered 								*/

		gv.inner_PGE_ite    = db.inner_PGE_ite;				/** number of inner PGE iterations, this has to be made mass or dG dependent 		*/
		gv.max_n_phase  	= db.max_n_phase;				/** maximum mol% phase change during one PGE iteration in wt% 						*/
		gv.max_g_phase  	= db.max_g_phase;				/** maximum delta_G of reference change during PGE 									*/
		gv.max_fac          = db.max_fac;					/** maximum update factor during PGE under-relax < 0.0, over-relax > 0.0 	 		*/

		gv.merge_value		= db.merge_value;				/** merge instances of solution phase if norm < val 								*/
		gv.re_in_n          = db.re_in_n;					/** fraction of phase when being reintroduce.  										*/
		gv.obj_tol 			= db.obj_tol;

		/* alloc */
		gv.ox 				= malloc (gv.len_ox * sizeof(char*)		);
		for (i = 0; i < gv.len_ox; i++){
			gv.ox[i] 		= malloc(20 * sizeof(char));	
			strcpy(gv.ox[i],db.ox[i]);
		}

		gv.PP_list 			= malloc (gv.len_pp * sizeof(char*)		);
		for (i = 0; i < (gv.len_pp); i++){	
			gv.PP_list[i] 	= malloc(20 * sizeof(char));
			strcpy(gv.PP_list[i],db.PP[i]);
		}

		gv.SS_list 			= malloc ((gv.len_ss) * sizeof (char*)	);
		gv.n_SS_PC     		= malloc ((gv.len_ss) * sizeof (int) 	);

		gv.verifyPC  		= malloc ((gv.len_ss) * sizeof (int) 	);
		gv.SS_PC_stp     	= malloc ((gv.len_ss) * sizeof (double) );
		for (i = 0; i < gv.len_ss; i++){ 
			gv.SS_list[i] 	= malloc(20 * sizeof(char)				);
			strcpy(gv.SS_list[i],db.SS[i]);
			gv.verifyPC[i]  = db.verifyPC[i]; 
			gv.n_SS_PC[i] 	= db.n_SS_PC[i]; 
			gv.SS_PC_stp[i] = db.SS_PC_stp[i]; 	
		}
		gv.act_PP     		= malloc ((gv.len_pp) * sizeof (int) 	);
		for (i = 0; i < gv.len_pp; i++){ gv.act_PP[i] = 1; 			};
	}
	else if (gv.EM_database == 2){
		igneous_dataset db 	= igneous_db;
		if (gv.EM_dataset == -1 || gv.EM_dataset == 62){
			if (gv.EM_dataset == 62 && gv.verbose == 1){
				printf(" dataset ds62 cannot be used with the igneous database (missing end-members), setting default dataset (ds634)\n");
			}
			gv.EM_dataset = db.ds_version;
		}
		gv.len_pp   		= db.n_pp;		
		gv.len_ss  			= db.n_ss;
		gv.len_ox  			= db.n_ox;

		gv.PC_df_add		= db.PC_df_add;					/** min value of df under which the PC is added 									*/
		gv.solver_switch_T  = db.solver_switch_T;
		gv.min_melt_T       = db.min_melt_T;				/** minimum temperature above which melt is considered 								*/

		gv.inner_PGE_ite    = db.inner_PGE_ite;				/** number of inner PGE iterations, this has to be made mass or dG dependent 		*/
		gv.max_n_phase  	= db.max_n_phase;				/** maximum mol% phase change during one PGE iteration in wt% 						*/
		gv.max_g_phase  	= db.max_g_phase;				/** maximum delta_G of reference change during PGE 									*/
		gv.max_fac          = db.max_fac;					/** maximum update factor during PGE under-relax < 0.0, over-relax > 0.0 	 		*/

		gv.merge_value		= db.merge_value;				/** merge instances of solution phase if norm < val 								*/
		gv.re_in_n          = db.re_in_n;					/** fraction of phase when being reintroduce.  										*/
		gv.obj_tol 			= db.obj_tol;

		gv.ox 				= malloc (gv.len_ox * sizeof(char*)		);
		for (i = 0; i < gv.len_ox; i++){
			gv.ox[i] 		= malloc(20 * sizeof(char));	
			strcpy(gv.ox[i],db.ox[i]);
		}

		gv.PP_list 			= malloc (gv.len_pp * sizeof(char*)		);
		for (i = 0; i < (gv.len_pp); i++){	
			gv.PP_list[i] 	= malloc(20 * sizeof(char));
			strcpy(gv.PP_list[i],db.PP[i]);
		}

		gv.SS_list 			= malloc ((gv.len_ss) * sizeof (char*)	);
		gv.n_SS_PC     		= malloc ((gv.len_ss) * sizeof (int) 	);
		gv.verifyPC  		= malloc ((gv.len_ss) * sizeof (int) 	);
		gv.SS_PC_stp     	= malloc ((gv.len_ss) * sizeof (double) );
		for (i = 0; i < gv.len_ss; i++){ 
			gv.SS_list[i] 	= malloc(20 * sizeof(char)				);
			strcpy(gv.SS_list[i],db.SS[i]);
			gv.verifyPC[i]  = db.verifyPC[i]; 
			gv.n_SS_PC[i] 	= db.n_SS_PC[i]; 
			gv.SS_PC_stp[i] = db.SS_PC_stp[i]; 	
		}
		gv.act_PP     		= malloc ((gv.len_pp) * sizeof (int) 	);
		for (i = 0; i < gv.len_pp; i++){ gv.act_PP[i]  = db.act_PP[i]; 			};
	}
	else if (gv.EM_database == 3){
		igneous_igad_dataset db 	= igneous_igad_db;
		if (gv.EM_dataset == -1 || gv.EM_dataset == 62){
			if (gv.EM_dataset == 62 && gv.verbose == 1){
				printf(" dataset ds62 cannot be used with the igneous database (missing end-members), setting default dataset (ds636)\n");
			}
			gv.EM_dataset = db.ds_version;
		}
		gv.len_pp   		= db.n_pp;		
		gv.len_ss  			= db.n_ss;
		gv.len_ox  			= db.n_ox;

		gv.PC_df_add		= db.PC_df_add;					/** min value of df under which the PC is added 									*/
		gv.solver_switch_T  = db.solver_switch_T;
		gv.min_melt_T       = db.min_melt_T;				/** minimum temperature above which melt is considered 								*/

		gv.inner_PGE_ite    = db.inner_PGE_ite;				/** number of inner PGE iterations, this has to be made mass or dG dependent 		*/
		gv.max_n_phase  	= db.max_n_phase;				/** maximum mol% phase change during one PGE iteration in wt% 						*/
		gv.max_g_phase  	= db.max_g_phase;				/** maximum delta_G of reference change during PGE 									*/
		gv.max_fac          = db.max_fac;					/** maximum update factor during PGE under-relax < 0.0, over-relax > 0.0 	 		*/

		gv.merge_value		= db.merge_value;				/** merge instances of solution phase if norm < val 								*/
		gv.re_in_n          = db.re_in_n;					/** fraction of phase when being reintroduce.  										*/
		gv.obj_tol 			= db.obj_tol;

		gv.ox 				= malloc (gv.len_ox * sizeof(char*)		);
		for (i = 0; i < gv.len_ox; i++){
			gv.ox[i] 		= malloc(20 * sizeof(char));	
			strcpy(gv.ox[i],db.ox[i]);
		}

		gv.PP_list 			= malloc (gv.len_pp * sizeof(char*)		);
		for (i = 0; i < (gv.len_pp); i++){	
			gv.PP_list[i] 	= malloc(20 * sizeof(char));
			strcpy(gv.PP_list[i],db.PP[i]);
		}

		gv.SS_list 			= malloc ((gv.len_ss) * sizeof (char*)	);
		gv.n_SS_PC     		= malloc ((gv.len_ss) * sizeof (int) 	);
		gv.verifyPC  		= malloc ((gv.len_ss) * sizeof (int) 	);
		gv.SS_PC_stp     	= malloc ((gv.len_ss) * sizeof (double) );
		for (i = 0; i < gv.len_ss; i++){ 
			gv.SS_list[i] 	= malloc(20 * sizeof(char)				);
			strcpy(gv.SS_list[i],db.SS[i]);
			gv.verifyPC[i]  = db.verifyPC[i]; 
			gv.n_SS_PC[i] 	= db.n_SS_PC[i]; 
			gv.SS_PC_stp[i] = db.SS_PC_stp[i]; 	
		}
		gv.act_PP     		= malloc ((gv.len_pp) * sizeof (int) 	);
		for (i = 0; i < gv.len_pp; i++){ gv.act_PP[i] = 1; 			};
	}
	else if (gv.EM_database == 4){
		ultramafic_dataset db = ultramafic_db;
		if (gv.EM_dataset == -1){
			gv.EM_dataset = db.ds_version;	
		}
		gv.len_pp   		= db.n_pp;		
		gv.len_ss  			= db.n_ss;
		gv.len_ox  			= db.n_ox;

		gv.PC_df_add		= db.PC_df_add;					/** min value of df under which the PC is added 									*/
		gv.solver_switch_T  = db.solver_switch_T;
		gv.min_melt_T       = db.min_melt_T;				/** minimum temperature above which melt is considered 								*/

		gv.inner_PGE_ite    = db.inner_PGE_ite;				/** number of inner PGE iterations, this has to be made mass or dG dependent 		*/
		gv.max_n_phase  	= db.max_n_phase;				/** maximum mol% phase change during one PGE iteration in wt% 						*/
		gv.max_g_phase  	= db.max_g_phase;				/** maximum delta_G of reference change during PGE 									*/
		gv.max_fac          = db.max_fac;					/** maximum update factor during PGE under-relax < 0.0, over-relax > 0.0 	 		*/

		gv.merge_value		= db.merge_value;				/** merge instances of solution phase if norm < val 								*/
		gv.re_in_n          = db.re_in_n;					/** fraction of phase when being reintroduce.  										*/
		gv.obj_tol 			= db.obj_tol;
		gv.ox 				= malloc (gv.len_ox * sizeof(char*)		);
		for (i = 0; i < gv.len_ox; i++){
			gv.ox[i] 		= malloc(20 * sizeof(char));	
			strcpy(gv.ox[i],db.ox[i]);
		}

		gv.PP_list 			= malloc (gv.len_pp * sizeof(char*)		);
		for (i = 0; i < (gv.len_pp); i++){	
			gv.PP_list[i] 	= malloc(20 * sizeof(char));
			strcpy(gv.PP_list[i],db.PP[i]);
		}

		gv.SS_list 			= malloc ((gv.len_ss) * sizeof (char*)	);
		gv.n_SS_PC     		= malloc ((gv.len_ss) * sizeof (int) 	);
		gv.verifyPC  		= malloc ((gv.len_ss) * sizeof (int) 	);
		gv.SS_PC_stp     	= malloc ((gv.len_ss) * sizeof (double) );
		for (i = 0; i < gv.len_ss; i++){ 
			gv.SS_list[i] 	= malloc(20 * sizeof(char)				);
			strcpy(gv.SS_list[i],db.SS[i]);
			gv.verifyPC[i]  = db.verifyPC[i]; 
			gv.n_SS_PC[i] 	= db.n_SS_PC[i]; 
			gv.SS_PC_stp[i] = db.SS_PC_stp[i]; 	
		}
		gv.act_PP     		= malloc ((gv.len_pp) * sizeof (int) 	);
		for (i = 0; i < gv.len_pp; i++){ gv.act_PP[i] = 1; 			};
	}
	else if (gv.EM_database == 5){
		ultramafic_ext_dataset db = ultramafic_ext_db;
		if (gv.EM_dataset == -1){
			gv.EM_dataset = db.ds_version;	
		}
		gv.len_pp   		= db.n_pp;		
		gv.len_ss  			= db.n_ss;
		gv.len_ox  			= db.n_ox;

		gv.PC_df_add		= db.PC_df_add;					/** min value of df under which the PC is added 									*/
		gv.solver_switch_T  = db.solver_switch_T;
		gv.min_melt_T       = db.min_melt_T;				/** minimum temperature above which melt is considered 								*/

		gv.inner_PGE_ite    = db.inner_PGE_ite;				/** number of inner PGE iterations, this has to be made mass or dG dependent 		*/
		gv.max_n_phase  	= db.max_n_phase;				/** maximum mol% phase change during one PGE iteration in wt% 						*/
		gv.max_g_phase  	= db.max_g_phase;				/** maximum delta_G of reference change during PGE 									*/
		gv.max_fac          = db.max_fac;					/** maximum update factor during PGE under-relax < 0.0, over-relax > 0.0 	 		*/

		gv.merge_value		= db.merge_value;				/** merge instances of solution phase if norm < val 								*/
		gv.re_in_n          = db.re_in_n;					/** fraction of phase when being reintroduce.  										*/
		gv.obj_tol 			= db.obj_tol;
		gv.ox 				= malloc (gv.len_ox * sizeof(char*)		);
		for (i = 0; i < gv.len_ox; i++){
			gv.ox[i] 		= malloc(20 * sizeof(char));	
			strcpy(gv.ox[i],db.ox[i]);
		}

		gv.PP_list 			= malloc (gv.len_pp * sizeof(char*)		);
		for (i = 0; i < (gv.len_pp); i++){	
			gv.PP_list[i] 	= malloc(20 * sizeof(char));
			strcpy(gv.PP_list[i],db.PP[i]);
		}

		gv.SS_list 			= malloc ((gv.len_ss) * sizeof (char*)	);
		gv.n_SS_PC     		= malloc ((gv.len_ss) * sizeof (int) 	);
		gv.verifyPC  		= malloc ((gv.len_ss) * sizeof (int) 	);
		gv.SS_PC_stp     	= malloc ((gv.len_ss) * sizeof (double) );
		for (i = 0; i < gv.len_ss; i++){ 
			gv.SS_list[i] 	= malloc(20 * sizeof(char)				);
			strcpy(gv.SS_list[i],db.SS[i]);
			gv.verifyPC[i]  = db.verifyPC[i]; 
			gv.n_SS_PC[i] 	= db.n_SS_PC[i]; 
			gv.SS_PC_stp[i] = db.SS_PC_stp[i]; 	
		}
		gv.act_PP     		= malloc ((gv.len_pp) * sizeof (int) 	);
		for (i = 0; i < gv.len_pp; i++){ gv.act_PP[i] = 1; 			};
	}
	else if (gv.EM_database == 6){
		mantle_dataset db = mantle_db;
		if (gv.EM_dataset == -1){
			gv.EM_dataset = db.ds_version;	
		}
		gv.len_pp   		= db.n_pp;		
		gv.len_ss  			= db.n_ss;
		gv.len_ox  			= db.n_ox;

		gv.PC_df_add		= db.PC_df_add;					/** min value of df under which the PC is added 									*/
		gv.solver_switch_T  = db.solver_switch_T;
		gv.min_melt_T       = db.min_melt_T;				/** minimum temperature above which melt is considered 								*/

		gv.inner_PGE_ite    = db.inner_PGE_ite;				/** number of inner PGE iterations, this has to be made mass or dG dependent 		*/
		gv.max_n_phase  	= db.max_n_phase;				/** maximum mol% phase change during one PGE iteration in wt% 						*/
		gv.max_g_phase  	= db.max_g_phase;				/** maximum delta_G of reference change during PGE 									*/
		gv.max_fac          = db.max_fac;					/** maximum update factor during PGE under-relax < 0.0, over-relax > 0.0 	 		*/

		gv.merge_value		= db.merge_value;				/** merge instances of solution phase if norm < val 								*/
		gv.re_in_n          = db.re_in_n;					/** fraction of phase when being reintroduce.  										*/
		gv.obj_tol 			= db.obj_tol;
		gv.ox 				= malloc (gv.len_ox * sizeof(char*)		);
		for (i = 0; i < gv.len_ox; i++){
			gv.ox[i] 		= malloc(20 * sizeof(char));	
			strcpy(gv.ox[i],db.ox[i]);
		}

		gv.PP_list 			= malloc (gv.len_pp * sizeof(char*)		);
		for (i = 0; i < (gv.len_pp); i++){	
			gv.PP_list[i] 	= malloc(20 * sizeof(char));
			strcpy(gv.PP_list[i],db.PP[i]);
		}

		gv.SS_list 			= malloc ((gv.len_ss) * sizeof (char*)	);
		gv.n_SS_PC     		= malloc ((gv.len_ss) * sizeof (int) 	);
		gv.verifyPC  		= malloc ((gv.len_ss) * sizeof (int) 	);
		gv.SS_PC_stp     	= malloc ((gv.len_ss) * sizeof (double) );
		for (i = 0; i < gv.len_ss; i++){ 
			gv.SS_list[i] 	= malloc(20 * sizeof(char)				);
			strcpy(gv.SS_list[i],db.SS[i]);
			gv.verifyPC[i]  = db.verifyPC[i]; 
			gv.n_SS_PC[i] 	= db.n_SS_PC[i]; 
			gv.SS_PC_stp[i] = db.SS_PC_stp[i]; 	
		}
		gv.act_PP     		= malloc ((gv.len_pp) * sizeof (int) 	);
		for (i = 0; i < gv.len_pp; i++){ gv.act_PP[i] = 1; 			};
	}
	else if (gv.EM_database == 7){
		metapelite_dataset_ext db 	= metapelite_ext_db;
		if (gv.EM_dataset == -1){
			gv.EM_dataset = db.ds_version;	
		}
		gv.len_pp   		= db.n_pp;		
		gv.len_ss  			= db.n_ss;
		gv.len_ox  			= db.n_ox;

		gv.PC_df_add		= db.PC_df_add;					/** min value of df under which the PC is added 									*/
		gv.solver_switch_T  = db.solver_switch_T;
		gv.min_melt_T       = db.min_melt_T ;				/** minimum temperature above which melt is considered 								*/

		gv.inner_PGE_ite    = db.inner_PGE_ite;				/** number of inner PGE iterations, this has to be made mass or dG dependent 		*/
		gv.max_n_phase  	= db.max_n_phase;				/** maximum mol% phase change during one PGE iteration in wt% 						*/
		gv.max_g_phase  	= db.max_g_phase;				/** maximum delta_G of reference change during PGE 									*/
		gv.max_fac          = db.max_fac;					/** maximum update factor during PGE under-relax < 0.0, over-relax > 0.0 	 		*/

		gv.merge_value		= db.merge_value;				/** merge instances of solution phase if norm < val 								*/
		gv.re_in_n          = db.re_in_n;					/** fraction of phase when being reintroduce.  										*/
		gv.obj_tol 			= db.obj_tol;

		/* alloc */
		gv.ox 				= malloc (gv.len_ox * sizeof(char*)		);
		for (i = 0; i < gv.len_ox; i++){
			gv.ox[i] 		= malloc(20 * sizeof(char));	
			strcpy(gv.ox[i],db.ox[i]);
		}

		gv.PP_list 			= malloc (gv.len_pp * sizeof(char*)		);
		for (i = 0; i < (gv.len_pp); i++){	
			gv.PP_list[i] 	= malloc(20 * sizeof(char));
			strcpy(gv.PP_list[i],db.PP[i]);
		}
		gv.SS_list 			= malloc ((gv.len_ss) * sizeof (char*)	);
		gv.n_SS_PC     		= malloc ((gv.len_ss) * sizeof (int) 	);
		gv.verifyPC  		= malloc ((gv.len_ss) * sizeof (int) 	);
		gv.SS_PC_stp     	= malloc ((gv.len_ss) * sizeof (double) );
		for (i = 0; i < gv.len_ss; i++){ 
			gv.SS_list[i] 	= malloc(20 * sizeof(char)				);
			strcpy(gv.SS_list[i],db.SS[i]);
			gv.verifyPC[i]  = db.verifyPC[i]; 
			gv.n_SS_PC[i] 	= db.n_SS_PC[i]; 
			gv.SS_PC_stp[i] = db.SS_PC_stp[i]; 	
		}
		gv.act_PP     		= malloc ((gv.len_pp) * sizeof (int) 	);
		for (i = 0; i < gv.len_pp; i++){ gv.act_PP[i] = 1; 			};
	}

	/**
	 ALLOCATE MEMORY OF OTHER GLOBAL VARIABLES
	*/
	gv.n_min     		= malloc ((gv.len_ss) * sizeof (int) 	);
	gv.n_ss_ph  		= malloc ((gv.len_ss) * sizeof (int) 	);
	gv.bulk_rock 		= malloc (gv.len_ox * sizeof(double)	);
	gv.PGE_mass_norm  	= malloc (gv.it_f*2 * sizeof (double) 	); 
	gv.Alg  			= malloc (gv.it_f*2 * sizeof (int) 		); 
	gv.gamma_norm  		= malloc (gv.it_f*2 * sizeof (double) 	); 
	gv.gibbs_ev  		= malloc (gv.it_f*2 * sizeof (double) 	); 
	gv.ite_time  		= malloc (gv.it_f*2 * sizeof (double) 	); 
	
	/* store values for numerical differentiation */
	/* The last entries MUST be [0-1][end] = 0.0  */
	gv.n_Diff = 8;
	gv.pdev = malloc (2 * sizeof(double*));			
	for (i = 0; i < 2; i++){
		gv.pdev[i] = malloc (gv.n_Diff * sizeof(double));
	}
	gv.pdev[0][0]  =  0.0;	gv.pdev[1][0]  =  1.0;	
	gv.pdev[0][1]  =  0.0;	gv.pdev[1][1]  = -1.0;	
	gv.pdev[0][2]  =  1.0;	gv.pdev[1][2]  =  1.0;	
	gv.pdev[0][3]  =  1.0;	gv.pdev[1][3]  = -1.0;	
	gv.pdev[0][4]  =  2.0;	gv.pdev[1][4]  =  0.0;
	gv.pdev[0][5]  =  1.0;	gv.pdev[1][5]  =  0.0;	
	gv.pdev[0][6]  =  3.0;	gv.pdev[1][6]  =  0.0;
	gv.pdev[0][7]  =  0.0;	gv.pdev[1][7]  =  0.0;

	gv.V_cor 	   		= malloc (2 * sizeof(double));

	/* declare size of chemical potential (gamma) vector */
	gv.dGamma 			= malloc (gv.len_ox * sizeof(double)	);
	gv.gam_tot  		= malloc (gv.len_ox * sizeof (double) 	); 
	gv.gam_tot_0		= malloc (gv.len_ox * sizeof (double) 	); 
	gv.delta_gam_tot  	= malloc (gv.len_ox * sizeof (double) 	); 
	gv.mass_residual 	= malloc (gv.len_ox * sizeof(double)	);	

	gv.lwork 			= 64;
	gv.ipiv     		= malloc ((gv.len_ox*3) * sizeof (int) 	);
	gv.work     		= malloc ((gv.len_ox*gv.lwork) * sizeof (double) 	);
	gv.n_solvi			= malloc ((gv.len_ss) * sizeof (int) 	);

	/* size of the flag array */
	gv.n_flags     = 5;

	/* allocate memory for pure and solution phase fractions */
	gv.pp_n    			= malloc (gv.len_pp * sizeof(double)	);									/** pure phase fraction vector */
	gv.pp_n_mol 		= malloc (gv.len_pp * sizeof(double)	);									/** pure phase fraction vector */
	gv.pp_n_wt  		= malloc (gv.len_pp * sizeof(double)	);									/** pure phase fraction vector */
	gv.pp_n_vol 		= malloc (gv.len_pp * sizeof(double)	);	
	gv.pp_xi    		= malloc (gv.len_pp * sizeof(double)	);									/** pure phase fraction vector */
	gv.delta_pp_n 		= malloc (gv.len_pp * sizeof(double)	);									/** pure phase fraction vector */
	gv.delta_pp_xi 		= malloc (gv.len_pp * sizeof(double)	);									/** pure phase fraction vector */
	gv.pp_flags 		= malloc (gv.len_pp * sizeof(int*)		);

	for (i = 0; i < (gv.len_pp); i++){	
		gv.pp_flags[i]   	= malloc (gv.n_flags  * sizeof(int));
	}
		
	/**
		PGE Matrix and RHS
	*/
	/* PGE method matrix and gradient arrays */
	gv.A_PGE  = malloc ((gv.len_ox*gv.len_ox*9) 	* sizeof(double));			
	gv.A0_PGE = malloc ((gv.len_ox*gv.len_ox*9) 	* sizeof(double));			
	gv.b_PGE  = malloc ((gv.len_ox*3) 				* sizeof(double));			

	gv.cp_id  = malloc ((gv.len_ox) 				* sizeof(int)	);			
	gv.pp_id  = malloc ((gv.len_ox) 				* sizeof(int)	);			

	gv.dn_cp  = malloc ((gv.len_ox) 				* sizeof(double));			
	gv.dn_pp  = malloc ((gv.len_ox) 				* sizeof(double));			

	/* stoechiometry matrix */
	gv.A  = malloc ((gv.len_ox) * sizeof(double*));		
	gv.A2 = malloc ((gv.len_ox) * sizeof(double*));			
	for (i = 0; i < (gv.len_ox); i++){
		gv.A[i]  = malloc ((gv.len_ox) * sizeof(double));
		gv.A2[i] = malloc ((gv.len_ox) * sizeof(double));
	
	gv.pc_id= malloc (gv.len_ox * sizeof(int));}
	gv.b 	= malloc (gv.len_ox * sizeof(double));	
	gv.b1 	= malloc (gv.len_ox * sizeof(double));	
	gv.tmp1 = malloc (gv.len_ox * sizeof(double));	
	gv.tmp2 = malloc (gv.len_ox * sizeof(double));	
	gv.tmp3 = malloc (gv.len_ox * sizeof(double));	
	gv.n_ss_array = malloc (gv.len_ss* sizeof(double));	
	/** 
		allocate oxides informations 						
	*/	
	z_b->apo     		= malloc (gv.len_ox * sizeof (double) ); 
	z_b->masspo     	= malloc (gv.len_ox * sizeof (double) );
	z_b->opo     		= malloc (gv.len_ox * sizeof (double) );
	z_b->cpo     		= malloc (gv.len_ox * sizeof (double) );
	z_b->ElEntropy     	= malloc (gv.len_ox * sizeof (double) );
	z_b->id     		= malloc (gv.len_ox * sizeof (int) 	  );
	z_b->elName     	= malloc (gv.len_ox * sizeof (char*) );
	for (i = 0; i < (gv.len_ox); i++){	
		z_b->elName[i] 	= malloc(20 * sizeof(char));
	}
	/**
		retrieve the right set of oxide and their informations 
	*/
	gv.H2O_id 	= -1;
	gv.CO2_id 	= -1;
	gv.CaO_id 	= -1;
	gv.Na2O_id 	= -1;
	gv.FeO_id 	= -1;
	gv.MgO_id 	= -1;
	gv.K2O_id 	= -1;
	gv.O_id 	= -1;
	gv.MnO_id 	= -1;
	
	oxide_data ox_in 	= oxide_info;
	for (i = 0; i < gv.len_ox; i++){
		for (j = 0; j < ox_in.n_ox; j++){
			if (strcmp( gv.ox[i], ox_in.oxName[j]) == 0){
				if (strcmp( gv.ox[i], "H2O") == 0){
					gv.H2O_id = i;
				}
				else if (strcmp( gv.ox[i], "Al2O3") == 0){
					gv.Al2O3_id = i;
				}
				else if (strcmp( gv.ox[i], "K2O") 	== 0){
					gv.K2O_id = i;
				}
				else if (strcmp( gv.ox[i], "TiO2") 	== 0){
					gv.TiO2_id = i;
				}
				else if (strcmp( gv.ox[i], "O") 	== 0){
					gv.O_id = i;
				}
				else if (strcmp( gv.ox[i], "S") 	== 0){
					gv.S_id = i;
				}
				else if (strcmp( gv.ox[i], "CO2") 	== 0){
					gv.CO2_id = i;
				}
				else if (strcmp( gv.ox[i], "Cr2O3") == 0){
					gv.Cr2O3_id = i;
				}
				else if (strcmp( gv.ox[i], "MnO") 	== 0){
					gv.MnO_id = i;
				}												
				z_b->apo[i]     	= ox_in.atPerOx[j];
				z_b->masspo[i]  	= ox_in.oxMass[j];
				z_b->opo[i]  		= ox_in.OPerOx[j];
				z_b->cpo[i]  		= ox_in.catPerOx[j];
				z_b->ElEntropy[i]   = ox_in.ElEntropy[j];
				strcpy(z_b->elName[i],ox_in.elName[j]);
				z_b->id[i]  		= j;
				break;
			}
		}
	}

	z_b->bulk_rock_cat  = malloc (gv.len_ox * sizeof (double) ); 
	z_b->bulk_rock  	= malloc (gv.len_ox * sizeof (double) ); 
	z_b->nzEl_array 	= malloc (gv.len_ox * sizeof (int) ); 
	z_b->zEl_array 		= malloc (gv.len_ox * sizeof (int) ); 

	/* sets end-member dataset information */
	if (gv.EM_dataset == 62){
			gv.n_em_db 			= 257;
	}
	else if (gv.EM_dataset == 633){
			gv.n_em_db 			= 289;
	}
	else if (gv.EM_dataset == 634){
			gv.n_em_db 			= 291;
	}		
	else if (gv.EM_dataset == 635){
			gv.n_em_db 			= 291;
	}	
	else if (gv.EM_dataset == 636){
			gv.n_em_db 			= 291;
	}	
	return gv;
}

/* Provide a list of test bulk-rock composition for the metapelite database (White et al., 2014)*/
global_variable get_bulk_metapelite( global_variable gv) {
 	if (gv.test != -1){
		if (gv.verbose == 1){
			printf("\n");
			printf("   - Minimization using in-built bulk-rock  : test %2d\n",gv.test);	
		}							
	}
	else{
		gv.test = 0;
		if (gv.verbose == 1){
			printf("\n");
			printf("   - No predefined bulk provided -> user custom bulk (if none provided, will run default FPWorldMedian pelite)\n");	
		}	
	}
	if (gv.test == 0){ 			//FPWorldMedian pelite
		/* SiO2 Al2O3 CaO MgO FeO K2O Na2O TiO2 O MnO H2O 	*/
		/* Forshaw, J. B., & Pattison, D. R. (2023) 		*/
		gv.bulk_rock[0]  = 70.999;		/** SiO2 	*/
		gv.bulk_rock[1]  = 12.8065;		/** Al2O2 	*/
		gv.bulk_rock[2]  = 0.771;		/** CaO  	*/
		gv.bulk_rock[3]  = 3.978;		/** MgO 	*/
		gv.bulk_rock[4]  = 6.342;		/** FeO 	*/
		gv.bulk_rock[5]  = 2.7895;		/** K2O	 	*/
		gv.bulk_rock[6]  = 1.481;		/** Na2O 	*/
		gv.bulk_rock[7]  = 0.758;		/** TiO2 	*/
		gv.bulk_rock[8]  = 0.72933;		/** O 		*/
		gv.bulk_rock[9]  = 0.075;		/** MnO 	*/
		gv.bulk_rock[10] = 30.000;		/** H2O 	*/
	}		
	else if (gv.test == 1){ 			//FPWorldMedian pelite !! WATER UNDER SATURATED!!
		/* SiO2 Al2O3 CaO MgO FeO K2O Na2O TiO2 O MnO H2O 	*/
		/* Forshaw, J. B., & Pattison, D. R. (2023) 		*/
		gv.bulk_rock[0]  = 70.999;		/** SiO2 	*/
		gv.bulk_rock[1]  = 12.8065;		/** Al2O2 	*/
		gv.bulk_rock[2]  = 0.771;		/** CaO  	*/
		gv.bulk_rock[3]  = 3.978;		/** MgO 	*/
		gv.bulk_rock[4]  = 6.342;		/** FeO 	*/
		gv.bulk_rock[5]  = 2.7895;		/** K2O	 	*/
		gv.bulk_rock[6]  = 1.481;		/** Na2O 	*/
		gv.bulk_rock[7]  = 0.758;		/** TiO2 	*/
		gv.bulk_rock[8]  = 0.72933;		/** O 		*/
		gv.bulk_rock[9]  = 0.075;		/** MnO 	*/
		gv.bulk_rock[10] = 5.000;		/** H2O 	*/
	}	
	else if (gv.test == 2){ 			//Pelite 
		/* SiO2 Al2O3 CaO MgO FeO K2O Na2O TiO2 O MnO H2O 	*/
		/* White et al., 2014, Fig 8. water oversaturated 	*/
		gv.bulk_rock[0]  = 64.578;		/** SiO2 	*/
		gv.bulk_rock[1]  = 13.651;		/** Al2O2 	*/
		gv.bulk_rock[2]  = 1.586;		/** CaO  	*/
		gv.bulk_rock[3]  = 5.529;		/** MgO 	*/
		gv.bulk_rock[4]  = 8.025;		/** FeO 	*/
		gv.bulk_rock[5]  = 2.943;		/** K2O	 	*/
		gv.bulk_rock[6]  = 2.000;		/** Na2O 	*/
		gv.bulk_rock[7]  = 0.907;		/** TiO2 	*/
		gv.bulk_rock[8]  = 0.65;		/** O 		*/
		gv.bulk_rock[9]  = 0.175;		/** MnO 	*/
		gv.bulk_rock[10] = 40.000;		/** H2O 	*/
	}	
	else if (gv.test == 3){ 			//Pelite 
		/* SiO2 Al2O3 CaO MgO FeO K2O Na2O TiO2 O MnO H2O 	*/
		/* White et al., 2014, Fig 8. water undersaturated 	*/
		gv.bulk_rock[0]  = 64.578;		/** SiO2 	*/
		gv.bulk_rock[1]  = 13.651;		/** Al2O2 	*/
		gv.bulk_rock[2]  = 1.586;		/** CaO  	*/
		gv.bulk_rock[3]  = 5.529;		/** MgO 	*/
		gv.bulk_rock[4]  = 8.025;		/** FeO 	*/
		gv.bulk_rock[5]  = 2.943;		/** K2O	 	*/
		gv.bulk_rock[6]  = 2.000;		/** Na2O 	*/
		gv.bulk_rock[7]  = 0.907;		/** TiO2 	*/
		gv.bulk_rock[8]  = 0.65;		/** O 		*/
		gv.bulk_rock[9]  = 0.175;		/** MnO 	*/
		gv.bulk_rock[10] = 6.244;		/** H2O 	*/
	}		
	else if (gv.test == 4){ 			//Pelite 
		/* SiO2 Al2O3 CaO MgO FeO K2O Na2O TiO2 O MnO H2O 	*/
		/* Garnet-Migmatite AV0832a (Riel et al., 2013) 	*/
		gv.bulk_rock[0]  = 73.9880;		/** SiO2 	*/
		gv.bulk_rock[1]  = 8.6143;		/** Al2O2 	*/
		gv.bulk_rock[2]  = 2.0146;		/** CaO  	*/
		gv.bulk_rock[3]  = 2.7401;		/** MgO 	*/
		gv.bulk_rock[4]  = 3.8451;		/** FeO 	*/
		gv.bulk_rock[5]  = 1.7686;		/** K2O	 	*/
		gv.bulk_rock[6]  = 2.4820;		/** Na2O 	*/
		gv.bulk_rock[7]  = 0.6393;		/** TiO2 	*/
		gv.bulk_rock[8]  = 0.1;			/** O 		*/
		gv.bulk_rock[9]  = 0.0630;		/** MnO 	*/
		gv.bulk_rock[10] = 10.0;		/** H2O 	*/
	}
	else{
		printf("Unknown test %i - please specify a different test! \n", gv.test);
	 	exit(EXIT_FAILURE);
	}
	return gv;
}

/* Provide a list of test bulk-rock composition for the metabasite database (Green et al., 2016)*/
global_variable get_bulk_metabasite( global_variable gv) {
 	if (gv.test != -1){
		if (gv.verbose == 1){
			printf("\n");
			printf("   - Minimization using in-built bulk-rock  : test %2d\n",gv.test);	
		}							
	}
	else{
		gv.test = 0;
		if (gv.verbose == 1){
			printf("\n");
			printf("   - No predefined bulk provided -> user custom bulk (if none provided, will run default SM89 oxidised average MORB)\n");	
		}	
	}
	if (gv.test == 0){
		/* SiO2 Al2O3 CaO MgO FeO K2O Na2O TiO2 O H2O 	*/
		/* SM89 oxidised average MORB composition of Sun & McDonough (1989)		*/
		gv.bulk_rock[0]  = 52.47;		/** SiO2 	*/
		gv.bulk_rock[1]  = 9.10;		/** Al2O2 	*/
		gv.bulk_rock[2]  = 12.21;		/** CaO  	*/
		gv.bulk_rock[3]  = 12.71;		/** MgO 	*/
		gv.bulk_rock[4]  = 8.15 ;		/** FeO 	*/
		gv.bulk_rock[5]  = 0.23;		/** K2O	 	*/
		gv.bulk_rock[6]  = 2.61;		/** Na2O 	*/
		gv.bulk_rock[7]  = 1.05;		/** TiO2 	*/
		gv.bulk_rock[8]  = 1.47;		/** O 		*/
		gv.bulk_rock[9]  = 20.0;		/** H2O 	*/
	}	
	else if (gv.test == 1){
		/* SiO2 Al2O3 CaO MgO FeO K2O Na2O TiO2 O H2O 	*/
		/* AG9: Natural amphibolites and low-temperature granulites (unpublished)		*/
		gv.bulk_rock[0]  = 51.08;		/** SiO2 	*/
		gv.bulk_rock[1]  = 9.68;		/** Al2O2 	*/
		gv.bulk_rock[2]  = 13.26;		/** CaO  	*/
		gv.bulk_rock[3]  = 11.21;		/** MgO 	*/
		gv.bulk_rock[4]  = 11.66 ;		/** FeO 	*/
		gv.bulk_rock[5]  = 0.16;		/** K2O	 	*/
		gv.bulk_rock[6]  = 0.79;		/** Na2O 	*/
		gv.bulk_rock[7]  = 1.37;		/** TiO2 	*/
		gv.bulk_rock[8]  = 0.80;		/** O 		*/
		gv.bulk_rock[9]  = 20.0;		/** H2O 	*/
	}
	else if (gv.test == 2){
		/* SiO2 Al2O3 CaO MgO FeO K2O Na2O TiO2 O H2O 	*/
		/* SQA: Synthetic amphibolite composition of Pati~no Douce & Beard(1995) (glass analysis)		*/
		gv.bulk_rock[0]  = 60.05;		/** SiO2 	*/
		gv.bulk_rock[1]  = 6.62;		/** Al2O2 	*/
		gv.bulk_rock[2]  = 8.31;		/** CaO  	*/
		gv.bulk_rock[3]  = 9.93;		/** MgO 	*/
		gv.bulk_rock[4]  = 6.57 ;		/** FeO 	*/
		gv.bulk_rock[5]  = 0.44;		/** K2O	 	*/
		gv.bulk_rock[6]  = 1.83;		/** Na2O 	*/
		gv.bulk_rock[7]  = 1.27;		/** TiO2 	*/
		gv.bulk_rock[8]  = 0.33;		/** O 		*/
		gv.bulk_rock[9]  = 4.64;		/** H2O 	*/
	}
	else if (gv.test == 3){
		/* SiO2 Al2O3 CaO MgO FeO K2O Na2O TiO2 O H2O 	*/
		/* BL478: Sample 478 of Beard & Lofgren (1991)		*/
		gv.bulk_rock[0]  = 52.73;		/** SiO2 	*/
		gv.bulk_rock[1]  = 9.57;		/** Al2O2 	*/
		gv.bulk_rock[2]  = 9.94;		/** CaO  	*/
		gv.bulk_rock[3]  = 6.76;		/** MgO 	*/
		gv.bulk_rock[4]  = 10.49;		/** FeO 	*/
		gv.bulk_rock[5]  = 0.24;		/** K2O	 	*/
		gv.bulk_rock[6]  = 3.28;		/** Na2O 	*/
		gv.bulk_rock[7]  = 1.2;			/** TiO2 	*/
		gv.bulk_rock[8]  = 1.3;			/** O 		*/
		gv.bulk_rock[9]  = 3.5;			/** H2O 	*/
	}
	else if (gv.test == 4){
		/* SiO2 Al2O3 CaO MgO FeO K2O Na2O TiO2 O H2O 	*/
		/* PDB95 */
		gv.bulk_rock[0]  = 60.0532;		/** SiO2 	*/
		gv.bulk_rock[1]  = 6.6231;		/** Al2O2 	*/
		gv.bulk_rock[2]  = 8.3095;		/** CaO  	*/
		gv.bulk_rock[3]  = 9.9281;		/** MgO 	*/
		gv.bulk_rock[4]  = 6.5693;		/** FeO 	*/
		gv.bulk_rock[5]  = 0.4435;		/** K2O	 	*/
		gv.bulk_rock[6]  = 1.8319;		/** Na2O 	*/
		gv.bulk_rock[7]  = 1.2708;		/** TiO2 	*/
		gv.bulk_rock[8]  = 0.3289;		/** O 		*/
		gv.bulk_rock[9]  = 4.6146;		/** H2O 	*/
	}
	else{
		printf("Unknown test %i - please specify a different test! \n", gv.test);
	 	exit(EXIT_FAILURE);
	}
	return gv;
}

/* Provide a list of test bulk-rock composition for the igneous database (Holland et al., 2018)*/
global_variable get_bulk_igneous( global_variable gv) {
 	if (gv.test != -1){
		if (gv.verbose == 1){
			printf("\n");
			printf("   - Minimization using in-built bulk-rock  : test %2d\n",gv.test);	
		}							
	}
	else{
		gv.test = 0;
		if (gv.verbose == 1){
			printf("\n");
			printf("   - No predefined bulk provided -> user custom bulk (if none provided, will run default KLB1)\n");	
		}	
	}
	if (gv.test == 0){ //KLB1
		/* SiO2 Al2O3 CaO MgO FeO K2O Na2O TiO2 O Cr2O3 H2O */
		/* Bulk rock composition of Peridotite from Holland et al., 2018, given by E. Green */
		gv.bulk_rock[0]  = 38.494 ;		/** SiO2 	*/
		gv.bulk_rock[1]  = 1.776;		/** Al2O2 	*/
		gv.bulk_rock[2]  = 2.824;		/** CaO  	*/
		gv.bulk_rock[3]  = 50.566;		/** MgO 	*/
		gv.bulk_rock[4]  = 5.886;		/** FeO 	*/
		gv.bulk_rock[5]  = 0.01;		/** K2O	 	*/
		gv.bulk_rock[6]  = 0.250;		/** Na2O 	*/
		gv.bulk_rock[7]  = 0.10;		/** TiO2 	*/
		gv.bulk_rock[8]  = 0.096;		/** O 		*/
		gv.bulk_rock[9]  = 0.109;		/** Cr2O3 	*/
		gv.bulk_rock[10] =	0.0;	
	}
	
	else if (gv.test == 1){ //RE46
		/* SiO2 Al2O3 CaO MgO FeO K2O Na2O TiO2 O Cr2O3 H2O */
		/* Bulk rock composition of RE46 - Icelandic basalt -Yang et al., 1996, given by E. Green */
		/*   50.72   9.16  15.21  16.25  7.06   0.01  1.47   0.39   0.35   0.01  */
		gv.bulk_rock[0] = 50.72;	
		gv.bulk_rock[1] = 9.16;	
		gv.bulk_rock[2] = 15.21;	
		gv.bulk_rock[3] = 16.25;	
		gv.bulk_rock[4] = 7.06;	
		gv.bulk_rock[5] = 0.01;	
		gv.bulk_rock[6]  = 1.47;
		gv.bulk_rock[7]  = 0.39;
		gv.bulk_rock[8]  = 0.35;
		gv.bulk_rock[9]  = 0.01;
		gv.bulk_rock[10] =	0.0;
	}
	else if (gv.test == 2){
		/* SiO2 Al2O3 CaO MgO FeO K2O Na2O TiO2 O Cr2O3 H2O */
		/* N_MORB, Gale et al., 2013, given by E. Green */
		gv.bulk_rock[0] = 53.21;	
		gv.bulk_rock[1] = 9.41;	
		gv.bulk_rock[2] = 12.21;	
		gv.bulk_rock[3] = 12.21;	
		gv.bulk_rock[4] = 8.65;	
		gv.bulk_rock[5] = 0.09;	
		gv.bulk_rock[6]  = 2.90;
		gv.bulk_rock[7]  = 1.21;
		gv.bulk_rock[8]  = 0.69;
		gv.bulk_rock[9]  = 0.02;
		gv.bulk_rock[10] =	0.0;
	}
	else if (gv.test == 3){ //MIX1G
		/* SiO2 Al2O3 CaO MgO FeO K2O Na2O TiO2 O Cr2O3 H2O */
		/* MIX1-G, Hirschmann et al., 2003, given by E. Green */
		gv.bulk_rock[0] = 45.25;	
		gv.bulk_rock[1] = 8.89;	
		gv.bulk_rock[2] = 12.22;	
		gv.bulk_rock[3] = 24.68;	
		gv.bulk_rock[4] = 6.45;	
		gv.bulk_rock[5] = 0.03;	
		gv.bulk_rock[6]  = 1.39;
		gv.bulk_rock[7]  = 0.67;
		gv.bulk_rock[8]  = 0.11;
		gv.bulk_rock[9]  = 0.02;
		gv.bulk_rock[10] =	0.0;
	}
	else if (gv.test == 4){
		/* SiO2 Al2O3 CaO MgO FeO K2O Na2O TiO2 O Cr2O3 H2O */
		/* High Al basalt Baker 1983 */
		gv.bulk_rock[0] = 54.40;	
		gv.bulk_rock[1] = 12.96;	
		gv.bulk_rock[2] = 11.31;	
		gv.bulk_rock[3] = 7.68;	
		gv.bulk_rock[4] = 8.63;	
		gv.bulk_rock[5] = 0.54;	
		gv.bulk_rock[6]  = 3.93;
		gv.bulk_rock[7]  = 0.79;
		gv.bulk_rock[8]  = 0.41;
		gv.bulk_rock[9]  = 0.01;
		gv.bulk_rock[10] =	0.0;
	}	
	else if (gv.test == 5){ //T101
		/* SiO2 Al2O3 CaO MgO FeO K2O Na2O TiO2 O Cr2O3 H2O */
		/* Tonalite 101 */
		gv.bulk_rock[0] = 66.01;	
		gv.bulk_rock[1] = 11.98;	
		gv.bulk_rock[2] = 7.06;	
		gv.bulk_rock[3] = 4.16;	
		gv.bulk_rock[4] = 5.30;	
		gv.bulk_rock[5] = 1.57;	
		gv.bulk_rock[6]  = 4.12;
		gv.bulk_rock[7]  = 0.66;
		gv.bulk_rock[8]  = 0.97;
		gv.bulk_rock[9]  = 0.01;
		gv.bulk_rock[10] =	50.0;
	}		
	else if (gv.test == 6){	//Wet Basalt
		/* SiO2 Al2O3 CaO MgO FeO K2O Na2O TiO2 O Cr2O3 H2O */
		/* Bulk rock composition of test 8 */
		gv.bulk_rock[0] = 50.0810;	
		gv.bulk_rock[1] = 8.6901;	
		gv.bulk_rock[2] = 11.6698;	
		gv.bulk_rock[3] = 12.1438;	
		gv.bulk_rock[4] = 7.7832;	
		gv.bulk_rock[5] = 0.2150;
		gv.bulk_rock[6]  = 2.4978;
		gv.bulk_rock[7]  = 1.0059;
		gv.bulk_rock[8]  = 0.4670;
		gv.bulk_rock[9]  = 0.0100;
		gv.bulk_rock[10] =	5.4364;
	}
	else if (gv.test == 7){	//BP002 harzburgite xenolith (Tomlinson & Holland, 2021)
		/* SiO2 Al2O3 CaO MgO FeO K2O Na2O TiO2 O Cr2O3 H2O */
		/* Bulk rock composition of test 8 */
		gv.bulk_rock[0] = 40.399;	
		gv.bulk_rock[1] = 0.923;	
		gv.bulk_rock[2] = 0.412;	
		gv.bulk_rock[3] = 54.091;	
		gv.bulk_rock[4] = 3.929;
		gv.bulk_rock[5] = 0.01;
		gv.bulk_rock[6] = 0.024;
		gv.bulk_rock[7]  = 0.01;
		gv.bulk_rock[8]  = 0.095;
		gv.bulk_rock[9]  = 0.122;
		gv.bulk_rock[10]  = 0.0;
	}              
	else{
		printf("Unknown test %i - please specify a different test! \n", gv.test);
	 	exit(EXIT_FAILURE);
	}
	return gv;
}

/* Get benchmark bulk rock composition given by Holland et al., 2018*/
global_variable get_bulk_igneous_igad( global_variable gv) {
 	if (gv.test != -1){
		if (gv.verbose == 1){
			printf("\n");
			printf("   - Minimization using in-built bulk-rock  : test %2d\n",gv.test);	
		}							
	}
	else{
		gv.test = 0;
		if (gv.verbose == 1){
			printf("\n");
			printf("   - No predefined bulk provided -> user custom bulk (if none provided, will run default KLB1)\n");	
		}	
	}
	if (gv.test == 0){ //Ne-syenite
		/* SiO2 Al2O3 CaO MgO FeO K2O Na2O TiO2 O Cr2O3 H2O */
		/* Weller et al., 2023: New thermodynamic models for alkaline systems */
		gv.bulk_rock[0]  = 63.84 ;		/** SiO2 	*/
		gv.bulk_rock[1]  = 13.72;		/** Al2O2 	*/
		gv.bulk_rock[2]  = 3.09;		/** CaO  	*/
		gv.bulk_rock[3]  = 1.55;		/** MgO 	*/
		gv.bulk_rock[4]  = 5.07;		/** FeOt 	*/
		gv.bulk_rock[5]  = 4.04;		/** K2O	 	*/
		gv.bulk_rock[6]  = 9.38;		/** Na2O 	*/
		gv.bulk_rock[7]  = 0.78;		/** TiO2 	*/
		gv.bulk_rock[8]  = 1.47;		/** O 		*/
		gv.bulk_rock[9]  = 0.01;		/** Cr2O3 	*/
	}
	else if (gv.test == 1){ // Syenite
		/* SiO2 Al2O3 CaO MgO FeO K2O Na2O TiO2 O Cr2O3 */
		/* Weller et al., 2023: New thermodynamic models for alkaline systems */
		gv.bulk_rock[0] = 70.06;	
		gv.bulk_rock[1] = 11.63;	
		gv.bulk_rock[2] = 2.76;	
		gv.bulk_rock[3] = 1.50;	
		gv.bulk_rock[4] = 4.30;	
		gv.bulk_rock[5] = 3.72;	
		gv.bulk_rock[6]  = 6.41;
		gv.bulk_rock[7]  = 0.51;
		gv.bulk_rock[8]  = 0.89;
		gv.bulk_rock[9]  = 0.01;
	}
	else if (gv.test == 2){ // Ijolite
		/* SiO2 Al2O3 CaO MgO FeO K2O Na2O TiO2 O Cr2O3 */
		/* Weller et al., 2023: New thermodynamic models for alkaline systems */
		gv.bulk_rock[0] = 48.97;	
		gv.bulk_rock[1] = 12.76;	
		gv.bulk_rock[2] = 12.87;	
		gv.bulk_rock[3] = 5.21;	
		gv.bulk_rock[4] = 7.97;	
		gv.bulk_rock[5] = 1.66;	
		gv.bulk_rock[6]  = 10.66;
		gv.bulk_rock[7]  = 1.36;
		gv.bulk_rock[8]  = 1.66;
		gv.bulk_rock[9]  = 0.01;
	}  
	else if (gv.test == 3){ // 9418-Fig 3c
		/* SiO2 Al2O3 CaO MgO FeO K2O Na2O TiO2 O Cr2O3 */
		/* Weller et al., 2023: New thermodynamic models for alkaline systems */
		gv.bulk_rock[0] = 53.221;	
		gv.bulk_rock[1] = 11.671;	
		gv.bulk_rock[2] = 10.009;	
		gv.bulk_rock[3] = 6.597;	
		gv.bulk_rock[4] = 7.053;	
		gv.bulk_rock[5] = 5.582;	
		gv.bulk_rock[6]  = 2.956;
		gv.bulk_rock[7]  = 0.825;
		gv.bulk_rock[8]  = 1.94;
		gv.bulk_rock[9]  = 0.146;
	}  
	else if (gv.test == 4){ //KLB1
		/* SiO2 Al2O3 CaO MgO FeO K2O Na2O TiO2 O Cr2O3 */
		/* Bulk rock composition of Peridotite from Holland et al., 2018, given by E. Green */
		gv.bulk_rock[0]  = 38.494 ;		/** SiO2 	*/
		gv.bulk_rock[1]  = 1.776;		/** Al2O2 	*/
		gv.bulk_rock[2]  = 2.824;		/** CaO  	*/
		gv.bulk_rock[3]  = 50.566;		/** MgO 	*/
		gv.bulk_rock[4]  = 5.886;		/** FeO 	*/
		gv.bulk_rock[5]  = 0.01;		/** K2O	 	*/
		gv.bulk_rock[6]  = 0.250;		/** Na2O 	*/
		gv.bulk_rock[7]  = 0.10;		/** TiO2 	*/
		gv.bulk_rock[8]  = 0.096;		/** O 		*/
		gv.bulk_rock[9]  = 0.109;		/** Cr2O3 	*/
	}   
	else{
		printf("Unknown test %i - please specify a different test! \n", gv.test);
	 	exit(EXIT_FAILURE);
	}
	return gv;
}

/* Provide a list of test bulk-rock composition for the ultramafic database (Evans & Frost, 2021)*/
global_variable get_bulk_ultramafic( global_variable gv) {
 	if (gv.test != -1){
		if (gv.verbose == 1){
			printf("\n");
			printf("   - Minimization using in-built bulk-rock  : test %2d\n",gv.test);	
		}							
	}
	else{
		gv.test = 0;
		if (gv.verbose == 1){
			printf("\n");
			printf("   - No predefined bulk provided -> user custom bulk (if none provided, will run default Serpentine oxidized)\n");	
		}	
	}
	if (gv.test == 0){ //Evans&Forst 2021, Serpentine oxidized
		/* SiO2 Al2O3 MgO FeO O H2O S */
		gv.bulk_rock[0]  = 20.044 ;		/** SiO2 	*/
		gv.bulk_rock[1]  = 0.6256;		/** Al2O2 	*/
		gv.bulk_rock[2]  = 29.24;		/** MgO 	*/
		gv.bulk_rock[3]  = 3.149;		/** FeO 	*/
		gv.bulk_rock[4]  = 0.7324;		/** O 		*/
		gv.bulk_rock[5]  = 46.755;		/** H2O 	*/
		gv.bulk_rock[6]  = 0.3;			/** S 		*/		
	}
	else if (gv.test == 1){ //Evans&Forst 2021, Serpentine reduced
		/* SiO2 Al2O3 MgO FeO O H2O S */
		gv.bulk_rock[0]  = 20.044 ;		/** SiO2 	*/
		gv.bulk_rock[1]  = 0.6256;		/** Al2O2 	*/
		gv.bulk_rock[2]  = 29.24;		/** MgO 	*/
		gv.bulk_rock[3]  = 3.149;		/** FeO 	*/
		gv.bulk_rock[4]  = 0.1324;		/** O 		*/
		gv.bulk_rock[5]  = 46.755;		/** H2O 	*/	
		gv.bulk_rock[6]  = 0.3;			/** S 		*/		
	}
                
	else{
		printf("Unknown test %i - please specify a different test! \n", gv.test);
	 	exit(EXIT_FAILURE);
	}
	return gv;
}

/* Provide a list of test bulk-rock composition for the ultramafic database (Evans & Frost, 2021)*/
global_variable get_bulk_ultramafic_ext( global_variable gv) {
 	if (gv.test != -1){
		if (gv.verbose == 1){
			printf("\n");
			printf("   - Minimization using in-built bulk-rock  : test %2d\n",gv.test);	
		}							
	}
	else{
		gv.test = 0;
		if (gv.verbose == 1){
			printf("\n");
			printf("   - No predefined bulk provided -> user custom bulk (if none provided, will run default Serpentine oxidized)\n");	
		}	
	}
	if (gv.test == 0){ //Barberton Komatiite, Tamblyn et al.,2022 
		/* SiO2 Al2O3 MgO FeO O H2O S */
		gv.bulk_rock[0]  = 38.51 ;		/** SiO2 	*/
		gv.bulk_rock[1]  = 2.25;		/** Al2O2 	*/
		gv.bulk_rock[2]  = 29.03;		/** MgO 	*/
		gv.bulk_rock[3]  = 4.65;		/** FeO 	*/
		gv.bulk_rock[4]  = 0.5;			/** O 		*/
		gv.bulk_rock[5]  = 16.0;		/** H2O 	*/	
		gv.bulk_rock[6]  = 0.1;			/** S 		*/		
		gv.bulk_rock[7]  = 6.92;		/** CaO		*/	
		gv.bulk_rock[8]  = 0.25;		/** Na2O	*/	
	}      
	else if (gv.test == 1){ //Evans&Forst 2021, Serpentine oxidized
		/* SiO2 Al2O3 MgO FeO O H2O S */
		gv.bulk_rock[0]  = 20.044 ;		/** SiO2 	*/
		gv.bulk_rock[1]  = 0.6256;		/** Al2O2 	*/
		gv.bulk_rock[2]  = 29.24;		/** MgO 	*/
		gv.bulk_rock[3]  = 3.149;		/** FeO 	*/
		gv.bulk_rock[4]  = 0.7324;		/** O 		*/
		gv.bulk_rock[5]  = 46.755;		/** H2O 	*/
		gv.bulk_rock[6]  = 0.3;			/** S 		*/		
		gv.bulk_rock[7]  = 2.0;			/** CaO		*/	
		gv.bulk_rock[8]  = 0.15;		/** Na2O	*/	
	}
	else if (gv.test == 2){ //Evans&Forst 2021, Serpentine reduced
		/* SiO2 Al2O3 MgO FeO O H2O S */
		gv.bulk_rock[0]  = 20.044 ;		/** SiO2 	*/
		gv.bulk_rock[1]  = 0.6256;		/** Al2O2 	*/
		gv.bulk_rock[2]  = 29.24;		/** MgO 	*/
		gv.bulk_rock[3]  = 3.149;		/** FeO 	*/
		gv.bulk_rock[4]  = 0.1324;		/** O 		*/
		gv.bulk_rock[5]  = 46.755;		/** H2O 	*/	
		gv.bulk_rock[6]  = 0.3;			/** S 		*/		
		gv.bulk_rock[7]  = 2.0;			/** CaO		*/	
		gv.bulk_rock[8]  = 0.15;			/** Na2O	*/	
	}
  
	else{
		printf("Unknown test %i - please specify a different test! \n", gv.test);
	 	exit(EXIT_FAILURE);
	}
	return gv;
}


/* Get benchmark bulk rock composition given by Holland et al., 2018*/
global_variable get_bulk_mantle( global_variable gv) {
 	if (gv.test != -1){
		if (gv.verbose == 1){
			printf("\n");
			printf("   - Minimization using in-built bulk-rock  : test %2d\n",gv.test);	
		}							
	}
	else{
		gv.test = 0;
		if (gv.verbose == 1){
			printf("\n");
			printf("   - No predefined bulk provided -> user custom bulk (if none provided, will run default KLB1)\n");	
		}	
	}
	if (gv.test == 0){ //KLB1
		/* SiO2 Al2O3 CaO MgO FeO Na2O */
		/* Bulk rock composition of Peridotite from Holland et al., 2018, given by E. Green */
		gv.bulk_rock[0]  = 38.494 ;		/** SiO2 	*/
		gv.bulk_rock[1]  = 1.776;		/** Al2O2 	*/
		gv.bulk_rock[2]  = 2.824;		/** CaO  	*/
		gv.bulk_rock[3]  = 50.566;		/** MgO 	*/
		gv.bulk_rock[4]  = 5.886;		/** FeO 	*/
		gv.bulk_rock[5]  = 0.250;		/** Na2O 	*/
	}  
	else if (gv.test == 1){ //Pyrolite
		/* SiO2 Al2O3 CaO MgO FeO Na2O */
		gv.bulk_rock[0]  = 38.89 ;		/** SiO2 	*/
		gv.bulk_rock[1]  = 2.2;			/** Al2O2 	*/
		gv.bulk_rock[2]  = 3.1;			/** CaO  	*/
		gv.bulk_rock[3]  = 50.0;		/** MgO 	*/
		gv.bulk_rock[4]  = 5.8;			/** FeO 	*/
		gv.bulk_rock[5]  = 0.01;		/** Na2O 	*/
	}  
	else if (gv.test == 1){ //Harzburgite
		/* SiO2 Al2O3 CaO MgO FeO Na2O */
		gv.bulk_rock[0]  = 36.39 ;		/** SiO2 	*/
		gv.bulk_rock[1]  = 0.7;			/** Al2O2 	*/
		gv.bulk_rock[2]  = 0.9;			/** CaO  	*/
		gv.bulk_rock[3]  = 56.6;		/** MgO 	*/
		gv.bulk_rock[4]  = 5.4;			/** FeO 	*/
		gv.bulk_rock[5]  = 0.01;		/** Na2O 	*/
	}  
	else{
		printf("Unknown test %i - please specify a different test! \n", gv.test);
	 	exit(EXIT_FAILURE);
	}
	return gv;
}


/* Provide a list of test bulk-rock composition for the metapelite database (White et al., 2014)*/
global_variable get_bulk_metapelite_ext( global_variable gv) {
 	if (gv.test != -1){
		if (gv.verbose == 1){
			printf("\n");
			printf("   - Minimization using in-built bulk-rock  : test %2d\n",gv.test);	
		}							
	}
	else{
		gv.test = 0;
		if (gv.verbose == 1){
			printf("\n");
			printf("   - No predefined bulk provided -> user custom bulk (if none provided, will run default FPWorldMedian pelite)\n");	
		}	
	}
	if (gv.test == 0){ 			//FPWorldMedian pelite
		/* SiO2 Al2O3 CaO MgO FeO K2O Na2O TiO2 O MnO H2O 	*/
		/* Forshaw, J. B., & Pattison, D. R. (2023) 		*/
		gv.bulk_rock[0]  = 70.999;		/** SiO2 	*/
		gv.bulk_rock[1]  = 12.8065;		/** Al2O2 	*/
		gv.bulk_rock[2]  = 0.771;		/** CaO  	*/
		gv.bulk_rock[3]  = 3.978;		/** MgO 	*/
		gv.bulk_rock[4]  = 6.342;		/** FeO 	*/
		gv.bulk_rock[5]  = 2.7895;		/** K2O	 	*/
		gv.bulk_rock[6]  = 1.481;		/** Na2O 	*/
		gv.bulk_rock[7]  = 0.758;		/** TiO2 	*/
		gv.bulk_rock[8]  = 0.72933;		/** O 		*/
		gv.bulk_rock[9]  = 0.075;		/** MnO 	*/
		gv.bulk_rock[10] = 30.000;		/** H2O 	*/
		gv.bulk_rock[11] = 10.0;			/** CO2  	*/
		gv.bulk_rock[12] = 0.2;			/** S 		*/
	}		
	else if (gv.test == 1){ 			//FPWorldMedian pelite !! WATER UNDER SATURATED!!
		/* SiO2 Al2O3 CaO MgO FeO K2O Na2O TiO2 O MnO H2O 	*/
		/* Forshaw, J. B., & Pattison, D. R. (2023) 		*/
		gv.bulk_rock[0]  = 70.999;		/** SiO2 	*/
		gv.bulk_rock[1]  = 12.8065;		/** Al2O2 	*/
		gv.bulk_rock[2]  = 0.771;		/** CaO  	*/
		gv.bulk_rock[3]  = 3.978;		/** MgO 	*/
		gv.bulk_rock[4]  = 6.342;		/** FeO 	*/
		gv.bulk_rock[5]  = 2.7895;		/** K2O	 	*/
		gv.bulk_rock[6]  = 1.481;		/** Na2O 	*/
		gv.bulk_rock[7]  = 0.758;		/** TiO2 	*/
		gv.bulk_rock[8]  = 0.72933;		/** O 		*/
		gv.bulk_rock[9]  = 0.075;		/** MnO 	*/
		gv.bulk_rock[10] = 5.000;		/** H2O 	*/
		gv.bulk_rock[11] = 10.0;		/** CO2  	*/
		gv.bulk_rock[12] = 0.2;			/** S 		*/
	}	
	else if (gv.test == 2){ 			//Pelite 
		/* SiO2 Al2O3 CaO MgO FeO K2O Na2O TiO2 O MnO H2O 	*/
		/* White et al., 2014, Fig 8. water oversaturated 	*/
		gv.bulk_rock[0]  = 64.578;		/** SiO2 	*/
		gv.bulk_rock[1]  = 13.651;		/** Al2O2 	*/
		gv.bulk_rock[2]  = 1.586;		/** CaO  	*/
		gv.bulk_rock[3]  = 5.529;		/** MgO 	*/
		gv.bulk_rock[4]  = 8.025;		/** FeO 	*/
		gv.bulk_rock[5]  = 2.943;		/** K2O	 	*/
		gv.bulk_rock[6]  = 2.000;		/** Na2O 	*/
		gv.bulk_rock[7]  = 0.907;		/** TiO2 	*/
		gv.bulk_rock[8]  = 0.65;		/** O 		*/
		gv.bulk_rock[9]  = 0.175;		/** MnO 	*/
		gv.bulk_rock[10] = 40.000;		/** H2O 	*/
		gv.bulk_rock[11] = 10.0;		/** CO2  	*/
		gv.bulk_rock[12] = 0.2;			/** S 		*/
	}	
	else if (gv.test == 3){ 			//Pelite 
		/* SiO2 Al2O3 CaO MgO FeO K2O Na2O TiO2 O MnO H2O 	*/
		/* White et al., 2014, Fig 8. water undersaturated 	*/
		gv.bulk_rock[0]  = 64.578;		/** SiO2 	*/
		gv.bulk_rock[1]  = 13.651;		/** Al2O2 	*/
		gv.bulk_rock[2]  = 1.586;		/** CaO  	*/
		gv.bulk_rock[3]  = 5.529;		/** MgO 	*/
		gv.bulk_rock[4]  = 8.025;		/** FeO 	*/
		gv.bulk_rock[5]  = 2.943;		/** K2O	 	*/
		gv.bulk_rock[6]  = 2.000;		/** Na2O 	*/
		gv.bulk_rock[7]  = 0.907;		/** TiO2 	*/
		gv.bulk_rock[8]  = 0.65;		/** O 		*/
		gv.bulk_rock[9]  = 0.175;		/** MnO 	*/
		gv.bulk_rock[10] = 6.244;		/** H2O 	*/
		gv.bulk_rock[11] = 10.0;		/** CO2  	*/
		gv.bulk_rock[12] = 0.2;			/** S 		*/
	}		
	else if (gv.test == 4){ 			//Pelite 
		/* SiO2 Al2O3 CaO MgO FeO K2O Na2O TiO2 O MnO H2O 	*/
		/* Garnet-Migmatite AV0832a (Riel et al., 2013) 	*/
		gv.bulk_rock[0]  = 73.9880;		/** SiO2 	*/
		gv.bulk_rock[1]  = 8.6143;		/** Al2O2 	*/
		gv.bulk_rock[2]  = 2.0146;		/** CaO  	*/
		gv.bulk_rock[3]  = 2.7401;		/** MgO 	*/
		gv.bulk_rock[4]  = 3.8451;		/** FeO 	*/
		gv.bulk_rock[5]  = 1.7686;		/** K2O	 	*/
		gv.bulk_rock[6]  = 2.4820;		/** Na2O 	*/
		gv.bulk_rock[7]  = 0.6393;		/** TiO2 	*/
		gv.bulk_rock[8]  = -0.5;		/** O 		*/
		gv.bulk_rock[9]  = 0.0630;		/** MnO 	*/
		gv.bulk_rock[10] = 10.0;		/** H2O 	*/
		gv.bulk_rock[11] = 10.0;		/** CO2  	*/
		gv.bulk_rock[12] = 0.2;			/** S 		*/
	}
	else{
		printf("Unknown test %i - please specify a different test! \n", gv.test);
	 	exit(EXIT_FAILURE);
	}
	return gv;
}