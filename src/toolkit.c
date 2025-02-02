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
#include <string.h>
#include <complex.h> 
#if __APPLE__
	/* dgetrf & dgetri routines */
	extern void dgetrf( int* M, int* N, double* A, int* lda, int* ipiv, int* info);
	extern void dgetri( int* N, double* A, int* lda, int* ipiv, double* work, int* lwork, int* info);
#else
	#include <lapacke.h> 
#endif 

#include "MAGEMin.h"

#include "gem_function.h"
#include "simplex_levelling.h"
#include "all_solution_phases.h"
#include "ss_min_function.h"
#include "nlopt.h"

#define FALSE  0
#define TRUE   1

/**
  function to print help and give command line parameters
*/
void print_help(	global_variable gv	){
	
	printf("\n");
	printf(" +---------------+-----------------------+\n");
	printf(" |    MAGEMin    |   %15s  |\n",gv.version);
	printf(" +---------------+-----------------------+\n");
	printf("\n");
	printf(" List of valid arguments:\n");
	printf(" ------------------------\n");	
	printf("\n");
	printf("  --version             : Sends back the called version of MAGEMin\n");	
	printf("  --Verb=       [int]   : Verbose option, 0. inactive, 1. active\n");	
	printf("  --File=       [str]   : File name containing multiple point calculation\n");
	printf("  --n_points=   [int]   : Number of points when using 'File' argument\n");
	printf("  --rg=         [str]   : ResearchGroup, can be 'tc' or 'sb (THERMOCALC, Stixrude-Lithgow-Bertelloni)\n");
	printf("  --db=         [str]   : Database, can be 'mp', 'ig', 'igad', 'um'or 'ume'* for TC and 2011/2024 for SB\n");
	printf("  --ds=         [int]   : TC End-member dataset, 62, 633 or 634 (stands for ds6xx)\n");
	printf("  --test=       [int]   : Number of points when using 'File' argument\n");
	printf("  --Pres=       [float] : Pressure in kilobar\n");
	printf("  --Temp=       [float] : Temperature in Celsius\n");
	printf("  --Bulk=       [float] : Bulk rock composition in [mol] or [wt] fraction**\n");
	printf("  --Gam=        [float] : Chemical potential of oxides (pure components)**\n");
	printf("  --sys_in=     [str]   : inputed system composition, [mol](default) or [wt]\n");
	printf("  --solver=     [int]   : solver: 0 for legacy and 1 for PGE (default)\n");
	printf("  --out_matlab= [int]   : Matlab text file output, 0. inactive, 1. active\n");
	printf("  --buffer= 	[str]   : choose among O2, qfm, mw, qif, nno, hm, cco, aH2O, aO2, aMgO, aFeO, aAl2O3, aTiO2\n");
	printf("  --buffer_n= 	[float] : multiplier with respect to qfm buffer\n");
	printf("  --mbCpx= 		[int]   : 0. omphacite, 1. augite (applies to metabasite database, see Green et al., 2016)\n");
	printf("  --mbIlm=      [int]   : 0. Ilm, 1. Ilmm see (Green et al., 2016)\n");
	printf("  --mpSp=       [int]   : 0. Ilm, 1. Ilmm see (White et al., 2014)\n");
	printf("  --mpIlm=      [int]   : 0. Sp, 1. Mt1 see (White et al., 2014)\n");
	printf("\n");
	printf(" *'mp': metapelite, 'mb': metabasite, 'ig': igneous H18->G23, 'igad': igneous alkaline dry, 'um': ultramafic, 'ume': ultramafic extended, 'mtl': mantle\n");
	printf("\n");
	printf(" **the list of oxides must be provided as follow:\n");
	printf("  'ig':               SiO2, Al2O3, CaO, MgO, FeOt, K2O, Na2O, TiO2, O, Cr2O3, H2O\n");
	printf("  'igad':             SiO2, Al2O3, CaO, MgO, FeOt, K2O, Na2O, TiO2, O, Cr2O3\n");
	printf("  'mp':               SiO2, Al2O3, CaO, MgO, FeOt, K2O, Na2O, TiO2, O, MnO, H2O\n");
	printf("  'mb':               SiO2, Al2O3, CaO, MgO, FeOt, K2O, Na2O, TiO2, O, H2O\n");
	printf("  'um':               SiO2, Al2O3, MgO, FeOt, O, H2O, S\n");
	printf("  'ume':              SiO2, Al2O3, MgO, FeOt, O, H2O, S, CaO, Na2O\n");
	printf("  'mtl':              SiO2, Al2O3, CaO, MgO, FeOt, Na2O\n");
	printf("\n");
	printf(" Note that FeOt (total iron) is used here!\n");	
	printf("\n\n");
	printf(" Examples of single point calculation:\n");
	printf(" ------------------------------------\n");	
	printf("\n");
    printf("  ./MAGEMin --Verb=1 --db=ig --Temp=718.750 --Pres=30.5000 --test=0 >&log.txt\n");
	printf("  ./MAGEMin --Verb=1 --db=mp --ds=633 --Temp=518.750 --Pres=3.5000 --test=2 >&log.txt\n");
    printf("\n");
	printf(" Here, the verbose is active and the bulk rock composition of 'test 0' is selected. The output of the verbose is saved as a log file 'log.txt'\n");
    printf(" Note that you don't have to use a test bulk composition, you can provide you own using arg '--Bulk='\n");
	printf("\n\n");
	printf(" Example multiple points calculation:\n");
	printf(" ------------------------------------\n");	
    printf("\n");
	printf(" To run multiple points at once you can pass an input file containing the list of points such as\n");
    printf("\n");
    printf("  ./MAGEMin --Verb=1 --db=ig  --File='path_to_file' --n_points=x\n");
    printf("\n");
	printf(" where 'path_to_file' is the location of the file and 'x' is an integer corresponding to the total number of points contained in the file. The file must have one point per line using the following structure\n");
    printf("\n");
	printf("  Mode(0-1), Pressure(kbar), Temperature(C), Gam1, Gam2, ..., Gamn\n");
    printf("\n");
	printf(" *Mode = 0 for global minimization\n");
    printf("\n");
	printf(" A valid list of points is for instance (MAGEMin_input.dat):\n");
    printf("\n");	
    printf("  0 0.0 800.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0\n");
    printf("  0 4.0 800.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0\n");
    printf("  0 7.0 800.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0\n");
    printf("  0 7.0 700.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0\n");	
    printf("\n");	
	printf(" Example of parallel points calculation:\n");
	printf(" ---------------------------------------\n");	
    printf("\n");	
    printf(" Simply call 'MPI' before MAGEMin and give the number of cores you want to use. Valid calls using previously defined input file are for instance\n");	
    printf("\n");	
    printf("  mpirun -np 8 ./MAGEMin --db=ig --File=/path_to_file/MAGEMin_input.dat --n_points=4\n");
    printf("  mpiexec -n 8 ./MAGEMin --db=ig --File=/path_to_file/MAGEMin_input.dat --n_points=4\n");	
    printf("\n");	
	printf(" Other useful information:\n");
	printf(" -------------------------\n");	
    printf("\n");	 
	printf(" Bulk-rock composition (wt or mol) expressed as FeO and Fe2O3 can be converted using the Matlab GUI to MAGEMin format ( [wt,mol] -> [mol], [FeO and Fe2O3] -> [FeOt and O])\n");
    printf("\n");	 
    printf("\n");	 

}


void print_2D_double_array(double nx, double ny, double **array, char *title){
	int i,j;
	printf(" %s:\n",title);
	for (i = 0; i < nx; i++){
		for (j = 0; j < ny; j++){
			printf(" %+10f",array[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}

void print_1D_double_array(double nx, double *array, char *title){
	int j;
	printf(" %s:\n",title);
	for (j = 0; j < nx; j++){
		printf(" %+10f",array[j]);
	}
	printf("\n");
}


void print_1D_int_array(double nx, int *array, char *title){
	int j;
	printf(" %s:\n",title);
	for (j = 0; j < nx; j++){
		printf(" %d",array[j]);
	}
	printf("\n");
}

/* generate a random double between 0 and a*/
double rnd(double a){
 	double b = ((double)rand()/(double)(RAND_MAX)) * a;
	return a;
}


double SUPCRT_to_HSC(double *ElH, double *comp, int size){
	double cor = 0.0;
	for (int i = 0; i < size; i++){
		cor -= comp[i]*ElH[i];
	}
	return cor;
}
double HSC_to_SUPCRT(double *ElH, double *comp, int size){
	double cor = 0.0;
	for (int i = 0; i < size; i++){
		cor += comp[i]*ElH[i];
	}
	return cor;
}


/* Normalize array to sum to 1 */
double sum_norm_xipi(double *xi, double *pi, int size) {
	int i;
	double norm = 0.0;
	for (i = 0; i < size; i++) {
		norm += fabs(xi[i]*pi[i] - pi[i]);
	}

	return norm;
}


/* Normalize array to sum to 1 */
double* norm_array(double *array, int size) {
	int i;
	double sum = 0.0;
	for (i = 0; i < size; i++) {
		sum += array[i];
	}
	for (i = 0; i < size; i++) {
		array[i] /= sum;
	}	
	return array;
}

/**
  retrieve bulk rock composition and PT compositions
*/
bulk_info retrieve_bulk_PT(				global_variable      gv,
										io_data 		    *input_data,
										int					 sgleP,
										bulk_info 			 z_b			){

	/* bulk from command line arguments */
	if (gv.arg_bulk[0] > 0.0) {
		if (gv.verbose == 1){
			printf("\n");
			printf("  - Minimization using bulk-rock composition from arg\n");	
		}	
		for (int i = 0; i < gv.len_ox; i++){ 
			gv.bulk_rock[i] = gv.arg_bulk[i];
		}
	}
	/* bulk from file */
	if (strcmp( gv.File, "none") != 0){

		z_b.P = input_data[sgleP].P;
		z_b.T = input_data[sgleP].T + 273.15;					/** K to C 		*/

		if (input_data[sgleP].in_bulk[0] > 0.0){
			if (gv.verbose == 1){
				printf("\n");
				printf("   - Minimization using bulk-rock composition from input file\n");	
			}	
			for (int i = 0; i < gv.len_ox; i++){
				gv.bulk_rock[i] = input_data[sgleP].in_bulk[i];					
			}
		}
	}

	/* transform bulk from wt% to mol% for minimiation */
	if (strcmp( gv.sys_in, "wt") == 0){	
		for (int i = 0; i < gv.len_ox; i++){ gv.bulk_rock[i] /= z_b.masspo[i];}
	}

	if (gv.verbose == 1){	

		if (gv.EM_database 		== 0){
			printf("  - Database                  : Metapelite (White et al., 2014)\n"	);
		}
		if (gv.EM_database 		== 1){
			printf("  - Database                  : Metabasite (Green et al., 2016)\n"	);
		}
		else if (gv.EM_database == 2){
			printf("  - Database                  : Igneous (Holland et al., 2018 -> Green et al., 2024)\n"	);
		}
		else if (gv.EM_database == 3){
			printf("  - Database                  : Igneous alkaline dry (Weller et al., 2024)\n"	);
		}
		else if (gv.EM_database == 4 ){
			printf("  - Database                  : Ultramafic (Evans & Frost, 2021)\n"	);
		}
		else if (gv.EM_database == 5 ){
			printf("  - Database                  : Ultramafic extended (Evans & Frost, 2021 + pl, amp and aug from Green et al., 2016)\n"	);
		}
		else if (gv.EM_database == 6 ){
			printf("  - Database                  : Uppermost lower mantle to upper mantle database (Holland et al., 2013)\n"	);
		}
		else if (gv.EM_database == 7 ){
			printf("  - Database                  : Metapelite extended (White et al., 2014; po from Evans & Frost, 2021;  amp, dio and aug from Green et al., 2016)\n"	);
		}

		if (strcmp( gv.sys_in, "mol") == 0){	
			printf("  - input system composition  : mol fraction\n"	);
		}
		else if (strcmp( gv.sys_in, "wt") == 0){	
			printf("  - input system composition  : wt fraction\n"	);
		}
		else{
			printf("  - input system composition  : unknown! [has to be mol or wt]\n");
		}
		printf("  - Buffer                    : %s\n",gv.buffer	);

	}	

	/** Normalize composition to sum to 1. 										*/
	norm_array(							gv.bulk_rock,
										gv.len_ox					);		

	/** here we check if the normalized mol fraction is < 1e-4 for oxides != H2O */
	/** if it is, then the fraction is set to 1e-4 -> this is a current limitation of system component reduction */
	int renorm = 0;
	for (int i = 0; i < gv.len_ox; i++){ 
		if (gv.bulk_rock[i] < 1.0e-4){

			if (gv.EM_database == 0){ 				// metapelite database
				if(strcmp( gv.ox[i], "H2O") != 0 && strcmp( gv.ox[i], "MnO") != 0  && strcmp( gv.ox[i], "O") != 0  && strcmp( gv.ox[i], "TiO2") != 0){
					gv.bulk_rock[i] = 1.0e-4;
					renorm = 1;
					if (gv.verbose == 1){
						printf("  - mol of %4s = %+.5f < 1e-4        : set back to 1e-4 to avoid minimization issues\n",gv.ox[i],gv.bulk_rock[i]);
					}	
				}
			}
			else if (gv.EM_database == 1){ 			// metabasite database
				if(strcmp( gv.ox[i], "TiO2") != 0  && strcmp( gv.ox[i], "O") != 0){
					gv.bulk_rock[i] = 1.0e-4;
					renorm = 1;
					if (gv.verbose == 1){
						printf("  - mol of %4s = %+.5f < 1e-4        : set back to 1e-4 to avoid minimization issues\n",gv.ox[i],gv.bulk_rock[i]);
					}	
				}
			}
			else if (gv.EM_database == 2){ 			// igneous database
				if(strcmp( gv.ox[i], "H2O") != 0  && strcmp( gv.ox[i], "TiO2") != 0 && strcmp( gv.ox[i], "Cr2O3") != 0 && strcmp( gv.ox[i], "O")  != 0 && strcmp( gv.ox[i], "K2O") != 0){
					gv.bulk_rock[i] = 1.0e-4;
					renorm = 1;
					if (gv.verbose == 1){
						printf("  - mol of %4s = %+.5f < 1e-4        : set back to 1e-4 to avoid minimization issues\n",gv.ox[i],gv.bulk_rock[i]);
					}	
				}
			}
			else if (gv.EM_database == 3){ 			// igneous database
				if(strcmp( gv.ox[i], "TiO2") != 0 && strcmp( gv.ox[i], "Cr2O3") != 0 && strcmp( gv.ox[i], "O")  != 0 ){
					gv.bulk_rock[i] = 1.0e-4;
					renorm = 1;
					if (gv.verbose == 1){
						printf("  - mol of %4s = %+.5f < 1e-4        : set back to 1e-4 to avoid minimization issues\n",gv.ox[i],gv.bulk_rock[i]);
					}	
				}
			}
			else if (gv.EM_database == 4){ 			// ultramafic database
				if(strcmp( gv.ox[i], "H2O") != 0 && strcmp( gv.ox[i], "S") != 0  && strcmp( gv.ox[i], "O") != 0){
					gv.bulk_rock[i] = 1.0e-4;
					renorm = 1;
					if (gv.verbose == 1){
						printf("  - mol of %4s = %+.5f < 1e-4        : set back to 1e-4 to avoid minimization issues\n",gv.ox[i],gv.bulk_rock[i]);
					}	
				}
			}
			else if (gv.EM_database == 5){ 			// ultramafic database
				if(strcmp( gv.ox[i], "H2O") != 0 && strcmp( gv.ox[i], "S") != 0  && strcmp( gv.ox[i], "O") != 0){
					gv.bulk_rock[i] = 1.0e-4;
					renorm = 1;
					if (gv.verbose == 1){
						printf("  - mol of %4s = %+.5f < 1e-4        : set back to 1e-4 to avoid minimization issues\n",gv.ox[i],gv.bulk_rock[i]);
					}	
				}
			}
			else if (gv.EM_database == 7){ 			// ultramafic database
				if(strcmp( gv.ox[i], "H2O") != 0 && strcmp( gv.ox[i], "S") != 0  && strcmp( gv.ox[i], "O") != 0 && strcmp( gv.ox[i], "MnO") != 0 && strcmp( gv.ox[i], "TiO2") != 0){
					gv.bulk_rock[i] = 1.0e-4;
					renorm = 1;
					if (gv.verbose == 1){
						printf("  - mol of %4s = %+.5f < 1e-4        : set back to 1e-4 to avoid minimization issues\n",gv.ox[i],gv.bulk_rock[i]);
					}	
				}
			}


		}
	}
	if (gv.verbose == 1){
		printf("\n");
	}
	if (renorm == 1){
		norm_array(							gv.bulk_rock,
											gv.len_ox					);						
	}

	return z_b;
};


/**
  retrieve bulk rock composition and PT compositions
  This function is not used in the C version of MAGEMin, but can be called via the Julia wrapper MAGEMin_C to normalize the composition
 */
void convert_system_comp(				global_variable      gv,
										char 				*sys_in,
										bulk_info 			 z_b,		
										double 				*bulk_rock			){

	/* transform bulk from wt% to mol% for minimization */
	if (strcmp( sys_in, "wt") == 0){	
		for (int i = 0; i < gv.len_ox; i++){ bulk_rock[i] /= z_b.masspo[i];}
	}

};


/**
  test function for Brent method
*/
double AFunction(int mode, double v, double *data) {
	double val = 0.0;
	if (mode == 0){
		double r 	 = 1.0/v;
		double R1    = data[0];
		double T     = data[1];
		double c1    = data[2];
		double c2    = data[3];
		double c3    = data[4];
		double c4    = data[5];
		double c5    = data[6];
		double c6    = data[7];
		double c7    = data[8];
		double c8    = data[9];
		double c9    = data[10];
		double c10   = data[11];
		double p_bar = data[12];

		val =  R1*T*( (r) + c1*pow(r,2.0) - pow(r,2.0)*( (c3 + 2.0*c4*(r) + 3.0*c5*pow(r,2.0) + 4.0*c6*pow(r,3.0) ) / pow(c2 + c3*(r) + c4*pow(r,2.0) + c5*pow(r,3.0) + c6*pow(r,4.0),2.0)) + c7*pow(r,2.0)*exp(-c8*(r)) + c9*pow(r,2.0)*exp(-c10*(r))) - p_bar;
	}
	else if (mode == 1){
		
		double sfdh    = data[0];
		double P       = data[1];
		double sfdhv   = data[2];
		double sfw     = data[3];
		double T       = data[4];
		double sfwv    = data[5];
		double sfn     = data[6];
		double R       = data[7];
		double sffac   = data[8];

		val = ( sfdh + P*sfdhv + (sfw + P*sfwv)*(2.*v - 1.) + sfn/(sfn + 1.)*R*T * (log(sfn*(1. - v)/(1. + sfn*v)) - sffac*log((1. - v)/(sfn + v))) );
	}
	else if (mode == 2) {

		double sfdh    = data[0];
		double P       = data[1];
		double sfdhv   = data[2];
		double sfw     = data[3];
		double sfwv    = data[4];
		double sffac   = data[5];
		double sfn     = data[6];
		double R       = data[7];
		double T       = data[8];
		
		val = ( sfdh + P*sfdhv + (sfw + P*sfwv)*(2.*v - 1.) + sffac*sfn/(sfn + 1.)*R*T * log(sfn*pow(1. - v,2.0)/((1. + sfn*v)*(sfn + v))) );
	}
	else{
		printf("Mode is not implemented!");
	}
	
  return val;
}

/** 
  TRUE if x1*x2 negative
*/
int RootBracketed(double x1,double x2) {
  int result;
  if ((x1 > 0. && x2 > 0.) || (x1 < 0. && x2 < 0.)) 
    result = FALSE;
  else
    result = TRUE;
  return result;
}

/** 
  returns the minimum of two real numbers
*/
double Minimum(double a,double b) {
  return ((a) < (b) ? (a) : (b));
}

/** 
  returns the maimum of two real numbers
*/
double Maximum(double a,double b) {
  return ((a) > (b) ? (a) : (b));
}

/**
  main Brent root finding routine 
*/
double BrentRoots(  double  x1, 
					double  x2, 
					double *data,
					double  Tolerance,
					
					int 	mode,
					int 	maxIterations,
					double *valueAtRoot,
					
					int 	*niter, 
					int 	*error 			
){

  double FPP = 1e-11, nearzero = 1e-40;

  double result, AA, BB, CC, DD, EE, FA, FB, FC, Tol1, PP, QQ, RR, SS, xm;
  int i, done;
  result 	= 0.0;
  EE 		= 0.0;
  CC 		= 0.0;

  i = 0; done = FALSE;   *error = 0;
  AA = x1;  BB = x2;  FA = AFunction(mode,AA,data); FB = AFunction(mode,BB,data);

  if (!(RootBracketed(FA,FB))) 
    *error = 1;
  else {
    FC = FB;
    do {
      if (!(RootBracketed(FC,FB))) {
        CC = AA; FC = FA; DD = BB - AA; EE = DD;
      }
      if (fabs(FC) < fabs(FB)) {
        AA = BB; BB = CC; CC = AA;
        FA = FB; FB = FC; FC = FA;
      }
      Tol1 = 2.0 * FPP * fabs(BB) + 0.5 * Tolerance;
      xm = 0.5 * (CC-BB);
      if ((fabs(xm) <= Tol1) || (fabs(FA) < nearzero)) {
        result = BB;
        done = TRUE;
        *valueAtRoot = AFunction(mode,result,data);
      } // A root has been found
      else {
        if ((fabs(EE) >= Tol1) && (fabs(FA) > fabs(FB))) {
          SS = FB/ FA;
          if (fabs(AA - CC) < nearzero) {
            PP = 2.0 * xm * SS;
            QQ = 1.0 - SS;
          }
          else {
            QQ = FA/FC;
            RR = FB /FC;
            PP = SS * (2.0 * xm * QQ * (QQ - RR) - (BB-AA) * (RR - 1.0));
            QQ = (QQ - 1.0) * (RR - 1.0) * (SS - 1.0);
          }
          if (PP > nearzero) QQ = -QQ;
          PP = fabs(PP);
          if ((2.0 * PP) < Minimum(3.0*xm *QQ-fabs(Tol1 * QQ), fabs(EE * QQ))) {
            EE = DD;  DD = PP/QQ;
          }
          else {
            DD = xm;   EE = DD;
          }
        }
        else {
          DD = xm;
          EE = DD;
        }
        AA = BB;
        FA = FB;
        if (fabs(DD) > Tol1) 
          BB = BB + DD;
        else {
          if (xm > 0.0) BB = BB + fabs(Tol1);
          else BB = BB - fabs(Tol1);
        }
        FB = AFunction(mode,BB,data);
        i++;
      }
	}  while ((!done) && (i < maxIterations));
    if (i >= maxIterations) *error = 2;
  }
  *niter = i;
  return result;
}


/**
  function to calculate norm of a vector
*/
double norm_vector(double *array ,int n){
	double norm = 0.0;
	for (int i = 0; i < n; i++){
		norm += array[i]*array[i];
	}
	norm = pow(norm,0.5);
	return norm;
}

/**
  function to calculate the Euclidean distance between two normalized vectors
  - used to decypher which phase to add during phase update
*/
double euclidean_distance(double *array1 ,double *array2 ,int n){
	double norm = 0.0;
	
	for (int i = 0; i < n; i++){
		norm += (array1[i]-array2[i])*(array1[i]-array2[i]);
	}
	norm = pow(norm,0.5);
	return norm;
}

double partial_euclidean_distance(double *array1 ,double *array2 ,int n){
	double norm = 0.0;
	
	for (int i = 0; i < n; i++){
		norm += (array1[i]-array2[i])*(array1[i]-array2[i]);
	}

	return norm;
}
/**
  inverse a matrix using LAPACKE dgetrf and dgetri
*/	
void inverseMatrix(int *ipiv, double *A1, int n, double *work, int lwork){	
	int    info;

	/* call lapacke to inverse Matrix */
	#if __APPLE__	
			dgetrf(&n, &n, A1, &n, ipiv, &info); 
			dgetri(&n, A1, &n, ipiv, work, &lwork, &info);

	#else
		info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, A1, n, ipiv); 
		info = LAPACKE_dgetri_work(LAPACK_ROW_MAJOR, n, A1, n, ipiv, work, lwork);
	#endif 

};

// Function to perform vector-matrix multiplication
void vector_matrix_multiplication(double* v, double** M, double* result, int vector_size, int matrix_cols) {
    for (int j = 0; j < matrix_cols; j++) {
        result[j] = 0.0;
        for (int i = 0; i < vector_size; i++) {
            result[j] += v[i] * M[i][j];
        }
    }
}
// Function to perform matrix-vector multiplication
void matrix_vector_multiplication(double** M, double* v, double* result, int matrix_rows, int matrix_cols) {
    for (int i = 0; i < matrix_rows; i++) {
        result[i] = 0.0;
        for (int j = 0; j < matrix_cols; j++) {
            result[i] += M[i][j] * v[j];
        }
    }
}

/**
  vector vectorT multiplication
*/	
double VecVecMul(double *B0, double *B1, int n){
	double result = 0.0;
	for (int i = 0; i < n; i++){
		result += (B0[i]*B1[i]);
	}
	return result;
};

/**
  vector matrix multiplication
*/	
void VecMatMul(double *B1, double *A1, double *B, int n){
	int i,j,k;
	for (i = 0; i < n; i++){
		B1[i] = 0.0;
		for (j = 0; j < n; j++){
			k = j + i*n;
			B1[i] += B[j]*A1[k];
		}
	}
};

/**
  matrix vector multiplication
*/	
void MatVecMul(double *A1, double *br, double *n_vec, int n){
	
	int i,j,k;
	for (i = 0; i < n; i++){
		n_vec[i] = 0.0;
		for (j = 0; j < n; j++){
			k = j + i*n;
			n_vec[i] += br[j]*A1[k];
		}
		if ( n_vec[i] < 1e-15){
			n_vec[i] = 1e-15; 
		}
	}
};


/**
  get active endmember
*/
int get_active_em(double *array, int n){
	int sum = 0.;
	for (int i = 0; i < n; i++){
		sum += array[i];
	}
	return sum;
};

/** 
  compare last character of a string with a reference character
  - function to find out the Liquid endmembers, and remove them from considerations after levelling
*/
int EndsWithTail(char *name, char* tail) {
    if (strlen(tail) > strlen(name))
        return 0;

    int len = strlen(name);

    if (strcmp(&name[len-strlen(tail)],tail) == 0)
        return 1;
    return 0;
}

/**
  function to calculate the sum of an array
*/
double sum_array(double *array, int size) {
	int i;
	double sum = 0.0;
	for (i = 0; i < size; i++) {
		sum += array[i];
	}
	return sum;
}

/**
  function to compare the sign of two variables
*/
int check_sign(double v1, double v2) {
	return ( (v1/fabs(v1) + v2/fabs(v2) == 0.0) ? (1) : (0));
}

/**
  function to compare the sign of two variables
*/
double sign(double x) {
	//return ( fabs(x)/x );
	return ( 1.0 );
}

/* Auxiliary routine: printing norms of matrix columns */
void print_vector_norm( char* desc, int m, int n, double* a, int lda ) {
        int i, j;
        double norm;
        printf( "\n %s\n", desc );
        for( j = 0; j < n; j++ ) {
                norm = 0.0;
                for( i = 0; i < m; i++ ) norm += a[i*lda+j] * a[i*lda+j];
                printf( " %g", norm );
        }
        printf( "\n" );
}

/*****************************************************************************************/


/**
	function to print out considered phases structure 
*/
void print_cp(		global_variable	 gv,
					csd_phase_set  	*cp
){
	printf("PRINT CONSIDERED PHASES\n");
	printf("------------------------\n\n");
	
	printf(" N_solvi %d: \n",gv.len_cp);
	for (int i = 0; i < gv.len_ss; i++){	
		printf(" %4s %d | ",gv.SS_list[i],gv.n_solvi[i]);
		printf("\n");
	}
	printf("\n");
	
	for (int i = 0; i < gv.len_cp; i++){		
		printf("[ #%d ]\n",i);
		printf(" SS name:  %4s\n",   cp[i].name);
		printf(" SS id:     %d\n",   cp[i].id);
		printf(" SS_nxeos:  %d\n",   cp[i].n_xeos);
		printf(" SS_nem:    %d\n",   cp[i].n_em);
		printf(" SS_df:    %+10f\n", cp[i].df*cp[i].factor);
		printf(" SS_factor:%+10f\n", cp[i].factor);
		printf(" SS_min_time:%+10f\n", cp[i].min_time);
	 	printf(" SS_flags: ");
		for (int ii = 0; ii < gv.n_flags; ii++){
			printf(" %d",cp[i].ss_flags[ii]);
		}
		printf("\n");
		printf(" SS_mode:  %+10f\n", cp[i].ss_n);
		printf("\n");
	 	printf(" SS_p_em:  ");
		for (int ii = 0; ii < cp[i].n_em; ii++){
			printf("%+10f ",cp[i].p_em[ii]);
		}
		printf("\n");
	 	printf(" SS_mu:  ");
		for (int ii = 0; ii < cp[i].n_em; ii++){
			printf("%+10f ",cp[i].mu[ii]);
		}
		printf("\n");		
	 	printf(" SS_xi_em:  ");
		for (int ii = 0; ii < cp[i].n_em; ii++){
			printf("%+10f ",cp[i].xi_em[ii]*cp[i].p_em[ii]);
		}
		printf("\n");
	 	printf(" SS_dgss:  ");
		for (int ii = 0; ii < cp[i].n_xeos; ii++){
			printf("%+10f ",cp[i].dguess[ii]);
		}
		printf("\n");
	 	printf(" SS_xgss:  ");
		for (int ii = 0; ii < cp[i].n_xeos; ii++){
			printf("%+10f ",cp[i].xeos[ii]);
		}
		printf("\n");
		printf("\n");
	}
};


global_variable get_tests_bulks(	global_variable  	 gv		
){

	if ( strcmp(gv.research_group, "tc") 	== 0 ){
		if 		(gv.EM_database == 0){
			gv = get_bulk_metapelite( 		gv );
		}
		else if (gv.EM_database == 1){
			gv = get_bulk_metabasite( 		gv );
		}
		else if (gv.EM_database == 2){
			gv = get_bulk_igneous( 			gv );
		}
		else if (gv.EM_database == 3){
			gv = get_bulk_igneous_igad( 	gv );
		}
		else if (gv.EM_database == 4){
			gv = get_bulk_ultramafic( 		gv );
		}
		else if (gv.EM_database == 5){
			gv = get_bulk_ultramafic_ext( 	gv );
		}
		else if (gv.EM_database == 6){
			gv = get_bulk_mantle( 			gv );
		}
		else if (gv.EM_database == 7){
			gv = get_bulk_metapelite_ext( 		gv );
		}
		else{
			printf(" Wrong database...\n");
		}
	}
	else if ( strcmp(gv.research_group, "sb") 	== 0 ){
		if 		(gv.EM_database == 0){
			gv = get_bulk_stx11( 		gv );
		}
		else{
			printf(" Wrong database...\n");
		}
	}
	else{
		printf(" Wrong research group...\n");
	}


	return gv;
}

/**
   rotate G-hyperplane using Gamma
*/
void print_SS_informations(		global_variable gv,
								SS_ref 			SS_ref_db,
								int				iss		
){
	printf(" %4s  | %+10f | %2d | %+10f | %+10f | ",gv.SS_list[iss],SS_ref_db.df,SS_ref_db.sf_ok,SS_ref_db.sum_xi,SS_ref_db.LM_time);
	for (int k = 0; k < SS_ref_db.n_xeos; k++) {
		printf(" %+6f",SS_ref_db.xeos[k]);
	}
	printf("\n");
	// for (int k = 0; k < SS_ref_db.n_xeos; k++) {
	// 	printf(" %+10f",SS_ref_db.dfx[k]);
	// }
	// for (int k = SS_ref_db.n_xeos; k < 12; k++){
	// 	printf(" %10s","-");
	// }
	// printf("\n");
}

/**
   rotate G-hyperplane using Gamma
*/
SS_ref rotate_hyperplane(	global_variable gv,
							SS_ref 			SS_ref_db		){
	
	/** rotate gbase with respect to the G-hyperplane (change of base) */
	for (int k = 0; k < SS_ref_db.n_em; k++) {
		SS_ref_db.gb_lvl[k] = SS_ref_db.gbase[k];
		for (int j = 0; j < gv.len_ox; j++) {
			SS_ref_db.gb_lvl[k] -= SS_ref_db.Comp[k][j]*gv.gam_tot[j];
		}
	}	

	return SS_ref_db;
}

/**
   non rotated G-hyperplane
*/
SS_ref non_rot_hyperplane(	global_variable gv,
							SS_ref 			SS_ref_db		){
	
	/** rotate gbase with respect to the G-hyperplane (change of base) */
	for (int k = 0; k < SS_ref_db.n_em; k++) {
		SS_ref_db.gb_lvl[k] = SS_ref_db.gbase[k];
	}	
	
	return SS_ref_db;
}

/**
   raw G-hyperplane using Gamma
*/
SS_ref raw_hyperplane(		global_variable  gv,
							SS_ref 			 SS_ref_db,
							double 			*gb				){
	
	/** rotate gbase with respect to the G-hyperplane (change of base) */
	for (int k = 0; k < SS_ref_db.n_em; k++) {
		SS_ref_db.gb_lvl[k] = gb[k];
	}	
	
	return SS_ref_db;
}


/**
   restrict solution phase hyper volume for local minimization
*/
SS_ref restrict_SS_HyperVolume(		global_variable gv, 
									SS_ref 			SS_ref_db,
									double 			box_size		){
									
	for (int j = 0; j < SS_ref_db.n_xeos; j++){
		SS_ref_db.bounds[j][0] = SS_ref_db.iguess[j] - box_size;
		SS_ref_db.bounds[j][1] = SS_ref_db.iguess[j] + box_size;
		
		if (SS_ref_db.bounds[j][0] < SS_ref_db.bounds_ref[j][0]){
			SS_ref_db.bounds[j][0] = SS_ref_db.bounds_ref[j][0];
		}
		if (SS_ref_db.bounds[j][1] > SS_ref_db.bounds_ref[j][1]){
			SS_ref_db.bounds[j][1] = SS_ref_db.bounds_ref[j][1];
		}
	}
						
	return SS_ref_db;										
}

/**
   check bounds
*/
SS_ref check_SS_bounds(		global_variable gv, 
							SS_ref SS_ref_db					){
									
	for (int j = 0; j < SS_ref_db.n_xeos; j++){
		if (SS_ref_db.iguess[j] < SS_ref_db.bounds_ref[j][0]){
			SS_ref_db.iguess[j] = SS_ref_db.bounds_ref[j][0];
		}
		if (SS_ref_db.iguess[j] > SS_ref_db.bounds_ref[j][1]){
			SS_ref_db.iguess[j] = SS_ref_db.bounds_ref[j][1];
		}
	}
						
	return SS_ref_db;										
}

/**
   retrieve the number of solution phase that are active 
*/
int getActiveSPhaseN(	global_variable gv,
						PP_ref *PP_ref_db,
						SS_ref *SS_ref_db		){
	
	int n = 0;
	for (int i = 0; i < gv.len_ss; i++){
		if (SS_ref_db[i].ss_flags[1] == 1){	n += 1;	}
	}
	
	return n;
}

/**
   retrieve the number of phases that are active 
*/
int getActivePhaseN(	global_variable gv,
						PP_ref *PP_ref_db,
						SS_ref *SS_ref_db		){
	
	int n = 0;
	for (int i = 0; i < gv.len_ss; i++){
		if (SS_ref_db[i].ss_flags[1] == 1){	n += 1;	}
	}
	for (int i = 0; i < gv.len_pp; i++){
		if (gv.pp_flags[i][1] == 1){	n += 1;	}
	}	
	return n;
}

int get_act_sf(double *A, int n){
	int n_act_sf = 0;						
	for (int j = 0; j < n; j++){
		if (A[j] < 0.0){
			n_act_sf += 1;
		}
	}
	return n_act_sf;
}

void get_act_sf_id(int *result, double *A, int n){
	int ix = 0;
	for (int i = 0; i < n; i++){
		if (A[i] < 0.0){
			result[ix] = i;
			ix += 1;
		}
	}	
}

/**
	get id of active pure phases
*/
global_variable get_pp_id(	global_variable  	 gv		
){
	/* extact pure phases to take into account */
	int k = 0;
	for (int i = 0; i < gv.len_pp; i++){
		if (gv.pp_flags[i][1] == 1){
			//gv.nEntry_sys += 1;
			gv.pp_id[k] = i;
			k       += 1;
		}
	}
	
	if (gv.n_pp_phase != k){
		printf("\n   !WARNING! inconsistent number of active phases (n_pp_phase vs sum(pp_flag[1])\n");
		printf("   !WARNING! n_pp_phase %i; sum(pp_flag[1]) %i;\n\n",gv.n_pp_phase,k);
	}
	
	return gv;
}

/**
	get id of active solution phases
*/
global_variable get_ss_id(	global_variable  	 gv,
							csd_phase_set  		*cp		
){
	/* extact solution phases to take into account */
	int k = 0;
	for (int i = 0; i < gv.len_cp; i++){
		if (cp[i].ss_flags[1] == 1){
			//gv.nEntry_sys += cp[i].n_xeos;
			gv.cp_id[k] = i;
			k += 1;
		}
	}
	
	if (gv.n_cp_phase != k){
		printf("\n   !WARNING! inconsistent number of active phases (n_ss_phase vs sum(ss_flag[1])\n");
		printf("   !WARNING! n_ss_phase %i; sum(ss_flag[1]) %i;\n\n",gv.n_cp_phase,k);
	}
	
	return gv;
}

/**
  Melt-fraction correction for P-wave and S-wave velocities
  The routine uses the reduction formulation of Clark et al., (2017) and is based on the equilibrium geometry model for the solid skeleton of Takei et al., 1997.
  * aspectRatio:  Coefficient defining the geometry of the solid framework (contiguity): 0.0 (layered melt distributed) < 0.1 (grain boundary melt) < 1.0 (melt in separated bubble pockets)
  
	printf(" Vp_Sol      (harm) : %+12.5f\t [km/s]\n",gv.solid_Vp);
	printf(" Vs_Sol      (harm) : %+12.5f\t [km/s]\n",gv.solid_Vs);
	gv.solid_Vs = anelastic_correction( 0,
										gv.solid_Vs,
										z_b.P,
										z_b.T 		);
	printf(" Vs_Sol      (anel) : %+12.5f\t [km/s]\n",gv.solid_Vs);
	if (gv.melt_fraction > 0.0){
		wave_melt_correction(  	gv.melt_bulkModulus,
								gv.solid_bulkModulus,
								gv.solid_shearModulus,
								gv.melt_density,
								gv.solid_density,
								gv.solid_Vp,	
								gv.solid_Vs,
								gv.melt_fraction,
								0.1,
								gv.V_cor				);
		printf("\n Vp_Melt_cor        : %+12.5f\t [km/s]\n",  gv.V_cor[0]);
		printf(" Vs_Melt_cor        : %+12.5f\t [km/s]\n",	gv.V_cor[1]);
		printf(" Vp/Vs_Melt_cor     : %+12.5f\t [km/s]\n",gv.V_cor[0]/gv.V_cor[1]);
	}
	printf("\n");

*/
global_variable wave_melt_correction( 	global_variable     gv,
										bulk_info 			z_b,	
										double  			aspectRatio		)
{

	if (gv.melt_fraction > 0.0 && gv.V_cor[1] > 0.0){
		double sum 	 =  gv.melt_fraction + gv.solid_fraction;
		gv.solid_fraction		/= sum;
		gv.melt_fraction	/= sum;

		double poisson = 0.25;

		double aij[3][4] ={ {0.318, 6.780, 57.560, 0.182},
							{0.164, 4.290, 26.658, 0.464},
							{1.549, 4.814, 8.777, -0.290}   };

		double bij[2][2] ={  {-0.3238, 0.2341},
							{-0.1819, 0.5103}  			};

		double a[3];
		double b[2];

		for (int i = 0; i < 3; i++){
			a[i] = aij[i][0]*exp( aij[i][1]*(poisson-0.25) + aij[i][2]*pow(poisson - 0.25,3.0) ) + aij[i][3];
		}
		for (int i = 0; i < 2; i++){
			b[i] = bij[i][0]*poisson + bij[i][1];
		}
		
		double nk      = a[0]*aspectRatio + a[1]*(1.0 - aspectRatio) + a[2]*aspectRatio*(1.0 - aspectRatio)*(0.5 - aspectRatio);
		double nmu     = b[0]*aspectRatio + b[1]*(1.0 - aspectRatio);

		double ksk_k   = pow(aspectRatio,nk);
		double musk_mu = pow(aspectRatio,nmu);

		double ksk     = ksk_k*gv.solid_bulkModulus;
		double musk    = musk_mu*gv.solid_shearModulus;

		double kb      = (1.0-gv.melt_fraction)*ksk;
		double mu      = (1.0-gv.melt_fraction)*musk;

		double LambdaK = gv.solid_bulkModulus/kb;
		double LambdaG = gv.solid_shearModulus/mu;

		double beta    = gv.solid_bulkModulus/gv.melt_bulkModulus;
		double gamma   = gv.solid_shearModulus/gv.solid_bulkModulus;

		double deltaVp = (((((beta -1.0)*LambdaK) / ((beta-1.0) + LambdaK) + 4.0/3.0*gamma*LambdaG ) / ( 1.0 + 4.0/3.0*gamma)) - (1.0 - gv.melt_density/gv.solid_density) )*(gv.melt_fraction/2.0);
		double deltaVs = ( LambdaG - (1.0 - gv.melt_density/gv.solid_density) )*(gv.melt_fraction/2.0);

		gv.V_cor[0] = gv.solid_Vp - deltaVp*gv.solid_Vp;
		gv.V_cor[1] = gv.solid_Vs - deltaVs*gv.solid_Vs;

		if (gv.V_cor[0] < 0.0){ gv.V_cor[0] = 0.0;}
		if (gv.V_cor[1] < 0.0){ gv.V_cor[1] = 0.0;}

	}


	if (gv.melt_fraction == 0.0){
		//aspect ratio has been choosen by hand to fit Vs = 2.0km/s at surface.
		aspectRatio 				= 0.25;
		double poisson 			 	= 0.25;

		double depth_m 				= z_b.P*1e5/(2600.0*9.81);
		double porosity_fraction 	= 0.474/pow(1.0+0.071*depth_m,5.989);

		double solid_fraction    	= 1.0 - porosity_fraction;
		double porosity_density  	= 1000.0;
		double porosity_bulkModulus = 2.9;
		double aij[3][4] ={ {0.318, 6.780, 57.560, 0.182},
							{0.164, 4.290, 26.658, 0.464},
							{1.549, 4.814, 8.777, -0.290}   };

		double bij[2][2] ={  {-0.3238, 0.2341},
							{-0.1819, 0.5103}  			};

		double a[3];
		double b[2];

		for (int i = 0; i < 3; i++){
			a[i] = aij[i][0]*exp( aij[i][1]*(poisson-0.25) + aij[i][2]*pow(poisson - 0.25,3.0) ) + aij[i][3];
		}
		for (int i = 0; i < 2; i++){
			b[i] = bij[i][0]*poisson + bij[i][1];
		}
	
		double nk      = a[0]*aspectRatio + a[1]*(1.0 - aspectRatio) + a[2]*aspectRatio*(1.0 - aspectRatio)*(0.5 - aspectRatio);
		double nmu     = b[0]*aspectRatio + b[1]*(1.0 - aspectRatio);

		double ksk_k   = pow(aspectRatio,nk);
		double musk_mu = pow(aspectRatio,nmu);

		double ksk     = ksk_k*gv.solid_bulkModulus;
		double musk    = musk_mu*gv.solid_shearModulus;

		double kb      = (1.0-porosity_fraction)*ksk;
		double mu      = (1.0-porosity_fraction)*musk;

		double LambdaK = gv.solid_bulkModulus/kb;
		double LambdaG = gv.solid_shearModulus/mu;

		double beta    = gv.solid_bulkModulus/porosity_bulkModulus;
		double gamma   = gv.solid_shearModulus/gv.solid_bulkModulus;

		// double deltaVp = (((((beta -1.0)*LambdaK) / ((beta-1.0) + LambdaK) + 4.0/3.0*gamma*LambdaG ) / ( 1.0 + 4.0/3.0*gamma)) - (1.0 - porosity_density/gv.solid_density) )*(porosity_fraction/2.0);
		double deltaVs = ( LambdaG - (1.0 - porosity_density/gv.solid_density) )*(porosity_fraction/2.0);

		// gv.V_cor[0] = gv.solid_Vp - deltaVp*gv.solid_Vp;
		gv.V_cor[1] = gv.solid_Vs - deltaVs*gv.solid_Vs;

		// if (gv.V_cor[0] < 0.0){ gv.V_cor[0] = 0.0;}
		if (gv.V_cor[1] < 0.0){ gv.V_cor[1] = 0.0;}	

	}
	return gv;

}	

double anelastic_correction(  	int 	water,
								double 	Vs0,
								double 	P,
								double 	T 		)
{
	double rH = 0.0, COH = 0.0;
    double kbar2pa = 100.0e3;

    double Pref    = P*kbar2pa;            	//pa
    double R       = 8.31446261815324;

    // values based on fitting experimental constraints (Behn et al., 2009)
    double alpha   = 0.27;
    double B0      = 1.28e8;               	//m/s
    double dref    = 1.24e-5;             	//m
    double COHref  = 50.0/1e6;             	//50H/1e6Si

    double Gref    = 1.09;
    double Eref    = 505.0e3;              	//J/mol
    double Vref    = 1.2e-5;               	//m3*mol

    double G       = 1.00;
    double E       = 420.0e3;              	//J/mol (activation energy)
    double V       = 1.2e-5;               	//m3*mol (activation volume)

    // using remaining values from Cobden et al., 2018
    // double omega   = 0.01;                 	//Hz (frequency to match for studied seismic system)
    double omega   = 3.0;                 	//Hz (frequency for Toba)
    double d       = 1e-2;                 	//m (grain size)
    
    if (water == 0){
        COH     = 50.0/1e6;             	//for dry mantle
        rH      = 0.0;                  	//for dry mantle
	}
	else if(water == 1){
        COH     = 1000.0/1e6;          	 	//for damp mantle    
        rH      = 1.0;                	  	//for damp mantle
	}
	else if(water == 2){
        COH     = 3000.0/1e6;           	//for wet mantle (saturated water)
        rH      = 2.0;                  	//for wet mantle
	}
	else{
        printf("WARN: water mode is not implemented. Valid values are 0 (dry),1 (dampened) and 2 (wet)\n");	
	}

    double B       = B0*pow(dref,G-Gref)*pow(COH/COHref,rH) * exp(((E+Pref*V)-(Eref + Pref*Vref))/(R*T));
    double Qinv    = pow( B*pow(d,(-G))*(1.0/omega) * exp(- (E+Pref*V)/(R*T)) ,alpha);
    double Vs_anel = Vs0*(1.0 - (Qinv)/(2.0*tan(3.141592*alpha/2.0) ) );

    return Vs_anel;

}


// char* retrieve_phase_name(		global_variable 	 gv,
// 								SS_ref  			*SS_ref_db,
// 								int 				 iss,		){


// 		if (gv.EM_database 		== 0){
// 			printf("  - Database                  : Metapelite (White et al., 2014)\n"	);
// 		}
// 		if (gv.EM_database 		== 1){
// 			printf("  - Database                  : Metabasite (Green et al., 2016)\n"	);
// 		}
// 		else if (gv.EM_database == 2){
// 			printf("  - Database                  : Igneous (Holland et al., 2018 -> Green et al., 2024)\n"	);
// 		}
// 		else if (gv.EM_database == 3){
// 			printf("  - Database                  : Igneous alkaline dry (Weller et al., 2024)\n"	);
// 		}
// 		else if (gv.EM_database == 4 ){
// 			printf("  - Database                  : Ultramafic (Evans & Frost, 2021)\n"	);
// 		}
// 		else if (gv.EM_database == 5 ){
// 			printf("  - Database                  : Ultramafic extended (Evans & Frost, 2021 + pl, amp and aug from Green et al., 2016)\n"	);
// 		}
// 		else if (gv.EM_database == 6 ){
// 			printf("  - Database                  : Uppermost lower mantle to upper mantle database (Holland et al., 2013)\n"	);
// 		}
// 		else if (gv.EM_database == 7 ){
// 			printf("  - Database                  : Metapelite extended (White et al., 2014; po from Evans & Frost, 2021;  amp, dio and aug from Green et al., 2016)\n"	);
// 		}

// 	// char   *ss_name[];
// 	// ss_name = 


// 	return ss_name;
// }



/** 
   This routine convert the molar fraction on 1 atom basis to mol fraction 
*/
global_variable compute_phase_mol_fraction(			global_variable 	 gv,
													bulk_info 	 		 z_b,
													PP_ref  			*PP_ref_db,
													SS_ref  			*SS_ref_db,
													csd_phase_set  		*cp					){

	int nox  = gv.len_ox;
	double n_at_bulk = 0.0;
	/* number of atom of bulk-rock composition */
	for (int i = 0; i < nox; i++){
		n_at_bulk += z_b.bulk_rock[i] * z_b.apo[i];
	}

	/* solution phases */
	double sum, n_at_ph, sum_mol_tot, sum_wt, sum_wt_tot;
	sum_mol_tot = 0.0;

	sum_wt_tot	= 0.0;
	for (int i = 0; i < gv.len_cp; i++){
		if (cp[i].ss_flags[1] == 1){

			sum 		= 0.0;
			for (int j = 0; j < gv.len_ox; j++){
				sum += cp[i].ss_comp[j];
			}

			n_at_ph 	= 0.0;
			for (int j = 0; j < gv.len_ox; j++){
				n_at_ph += cp[i].ss_comp[j]/sum * z_b.apo[j];
			}

			for (int j = 0; j < gv.len_ox; j++){
				cp[i].ss_comp_mol[j] = cp[i].ss_comp[j]/sum;
			}
			sum_wt  	= 0.0;
			for (int j = 0; j < gv.len_ox; j++){
				cp[i].ss_comp_wt[j] = cp[i].ss_comp_mol[j] * z_b.masspo[j];
				sum_wt += cp[i].ss_comp_wt[j];
			}	
		
			for (int j = 0; j < gv.len_ox; j++){
				cp[i].ss_comp_wt[j] /= sum_wt;
			}	

			cp[i].factor_norm = n_at_bulk/n_at_ph;

			cp[i].ss_n_mol   = cp[i].ss_n * cp[i].factor_norm;
			cp[i].ss_n_wt 	 = cp[i].ss_n_mol * sum_wt;
			sum_mol_tot 	+= cp[i].ss_n_mol;
			sum_wt_tot 		+= cp[i].ss_n_wt;
		}
	}

	/* pure phases */
	for (int i = 0; i < gv.len_pp; i++){
		if (gv.pp_flags[i][1] == 1 && gv.pp_flags[i][4] == 0){

			sum 		= 0.0;
			for (int j = 0; j < gv.len_ox; j++){
				if (PP_ref_db[i].Comp[j] > 0.0 ){
					sum += PP_ref_db[i].Comp[j];
				}
			}

			n_at_ph 	= 0.0;
			for (int j = 0; j < gv.len_ox; j++){
				n_at_ph += PP_ref_db[i].Comp[j]/sum * z_b.apo[j];
			}

			for (int j = 0; j < gv.len_ox; j++){
				PP_ref_db[i].Comp_mol[j] = PP_ref_db[i].Comp[j]/sum;
			}
			sum_wt  	= 0.0;
			for (int j = 0; j < gv.len_ox; j++){
				PP_ref_db[i].Comp_wt[j] = PP_ref_db[i].Comp_mol[j] * z_b.masspo[j];
				sum_wt += PP_ref_db[i].Comp_wt[j];
			}	
			for (int j = 0; j < gv.len_ox; j++){
				PP_ref_db[i].Comp_wt[j] /= sum_wt;
			}	

			PP_ref_db[i].factor_norm = n_at_bulk/n_at_ph;

			gv.pp_n_mol[i]   = gv.pp_n[i] * PP_ref_db[i].factor_norm;
			gv.pp_n_wt[i] 	 = gv.pp_n_mol[i] * sum_wt;
			sum_mol_tot 	+= gv.pp_n_mol[i];
			sum_wt_tot 		+= gv.pp_n_wt[i];
		}
	}

	/* normalize mol fractions */
	for (int i = 0; i < gv.len_cp; i++){
		if (cp[i].ss_flags[1] == 1){
			cp[i].ss_n_mol   /= sum_mol_tot;
			cp[i].ss_n_wt    /= sum_wt_tot;
		}
	}
	for (int i = 0; i < gv.len_pp; i++){
		if (gv.pp_flags[i][1] == 1){
			gv.pp_n_mol[i]   /= sum_mol_tot;
			gv.pp_n_wt[i]    /= sum_wt_tot;
		}
	}

	return gv;										
}


global_variable compute_density_volume_modulus(				int 				 EM_database,
															bulk_info 	 		 z_b,
															global_variable 	 gv,
															PP_ref  			*PP_ref_db,
															SS_ref  			*SS_ref_db,
															csd_phase_set  		*cp					){

							
	PP_ref PP_db;					
								
	double muC, muN, muE, muW, muNN, muNNN, muNE, muNW, G, muN0, muC0;
	
	double P 			  = z_b.P;					/** PC function uses the z_b structure this is why the Pressure is saved here */
	double T 			  = z_b.T;					/** PC function uses the z_b structure this is why the Pressure is saved here */
	double sum_volume     = 0.0;
	double sum_volume_sol = 0.0;
	double dGdTPP, dGdTMP, dG2dT2, dGdP, dG2dP2, dG2dP2_N, dGdP_N;// dGdP_P0,
	double mut, mut_N;
	double phase_isoTbulkModulus_P1;

	double density[gv.len_ox];
	int not_only_liq = 0;
	int ss;


	/** calculate mass, volume and densities */
	for (int i = 0; i < gv.len_cp; i++){
		if (cp[i].ss_flags[1] == 1){

			if (strcmp( cp[i].name, "liq") != 0){
				not_only_liq = 1;
			}

			ss = cp[i].id;	
			
			for (int k = 0; k < cp[i].n_xeos; k++) {
				SS_ref_db[ss].iguess[k] = cp[i].xeos[k];
			}
										
			/** calculate Molar Mass of solution phase */
			cp[i].mass = 0.0;
			for (int k = 0; k < gv.len_ox; k++){
				cp[i].mass	+= cp[i].ss_comp[k]*z_b.masspo[k];
			}

			/** calculate cp of solution phase */
			cp[i].phase_cp 	 		 = 0.0;
			cp[i].volume 	 		 = 0.0;
			cp[i].phase_expansivity	 = 0.0;
			cp[i].phase_bulkModulus  = 0.0;
			cp[i].phase_shearModulus = 0.0;
			cp[i].phase_entropy  	 = 0.0;
			cp[i].phase_enthalpy 	 = 0.0;
			cp[i].phase_isoTbulkModulus = 0.0;
			// cp[i].volume_P0 		 = 0.0;
			cp[i].thetaExp 			 = 0.0;
			phase_isoTbulkModulus_P1 = 0.0;

			for (int j = 0; j < cp[i].n_em; j++){ 
				if (SS_ref_db[ss].z_em[j] == 1.0){
					dG2dT2 					 = (SS_ref_db[ss].mu_array[0][j]-2.0*SS_ref_db[ss].mu_array[7][j]+SS_ref_db[ss].mu_array[1][j])/(gv.gb_T_eps*gv.gb_T_eps);
					dG2dP2 					 = (SS_ref_db[ss].mu_array[4][j]-2.0*SS_ref_db[ss].mu_array[5][j]+SS_ref_db[ss].mu_array[7][j])/(gv.gb_P_eps*gv.gb_P_eps);
					dG2dP2_N				 = (SS_ref_db[ss].mu_array[6][j]-2.0*SS_ref_db[ss].mu_array[4][j]+SS_ref_db[ss].mu_array[5][j])/(gv.gb_P_eps*gv.gb_P_eps);
					dGdTPP 					 = (SS_ref_db[ss].mu_array[2][j]-SS_ref_db[ss].mu_array[3][j])/(2.0*gv.gb_T_eps);
					dGdTMP 					 = (SS_ref_db[ss].mu_array[0][j]-SS_ref_db[ss].mu_array[1][j])/(2.0*gv.gb_T_eps);
					dGdP					 = (SS_ref_db[ss].mu_array[5][j]-SS_ref_db[ss].mu_array[7][j])/(gv.gb_P_eps);
					dGdP_N 					 = (SS_ref_db[ss].mu_array[4][j]-SS_ref_db[ss].mu_array[5][j])/(gv.gb_P_eps);
					// dGdP_P0 				 = (SS_ref_db[ss].mu_array[8][j]-SS_ref_db[ss].mu_array[9][j])/(gv.gb_P_eps);
					/* heat capacity 	*/
					cp[i].phase_cp    		+= -T*(dG2dT2)*cp[i].p_em[j];
					
					/* volume 			*/
					cp[i].volume    		+= (dGdP)*cp[i].p_em[j];

					/* volume 			*/
					// cp[i].volume_P0    		+= (dGdP_P0)*cp[i].p_em[j];

					/* entropy   		*/
					cp[i].phase_entropy 	+= -(dGdTMP)*cp[i].p_em[j]*cp[i].factor;

					/* expansivity 		*/
					cp[i].phase_expansivity += (1.0/(dGdP)*((dGdTPP-dGdTMP)/(gv.gb_P_eps)))*cp[i].p_em[j];
					
					/* bulk modulus	*/
					if ( strcmp(gv.research_group, "sb") 	== 0 ){
						cp[i].phase_bulkModulus += SS_ref_db[ss].ElBulkMod[j] * cp[i].p_em[j];
					}
					else{
						cp[i].phase_bulkModulus += -dGdP/( dG2dP2 + pow(((dGdTPP-dGdTMP)/(gv.gb_P_eps)),2.0)/dG2dT2 ) * cp[i].p_em[j];
					}

					/* iso bulk modulus	*/
					cp[i].phase_isoTbulkModulus += -dGdP/( dG2dP2 ) 	* cp[i].p_em[j];
					phase_isoTbulkModulus_P1	+= -dGdP_N/( dG2dP2_N ) * cp[i].p_em[j];
							
					/* shear modulus	*/
					cp[i].phase_shearModulus += SS_ref_db[ss].ElShearMod[j] * cp[i].p_em[j];
				}
			}	

			G = 0.0;
			for (int j = 0; j < gv.len_ox; j++){
				G += cp[i].ss_comp[j]*gv.gam_tot[j];
			}

			/* enthalpy   		*/
			cp[i].phase_enthalpy = cp[i].phase_entropy*T + G*cp[i].factor;
	
			/** calculate density from volume */
			cp[i].phase_density  = (cp[i].mass*1000.0)/(cp[i].volume*10.0);

			mut 				 = (3.0*cp[i].phase_isoTbulkModulus - 6.0*cp[i].phase_isoTbulkModulus*gv.poisson_ratio) / (2. + 2.0*gv.poisson_ratio)/10.0;
			mut_N 				 = (3.0*phase_isoTbulkModulus_P1 	- 6.0*phase_isoTbulkModulus_P1*gv.poisson_ratio) 	/ (2. + 2.0*gv.poisson_ratio)/10.0;

			cp[i].thetaExp 		 = (mut_N - mut)/gv.gb_P_eps - (cp[i].phase_bulkModulus*cp[i].phase_expansivity)/(cp[i].phase_cp*cp[i].phase_density);


			if (strcmp( cp[i].name, "liq") == 0){
				gv.melt_density   	= cp[i].phase_density;
				gv.melt_fraction  	= cp[i].ss_n_mol;
				gv.melt_bulkModulus = cp[i].phase_bulkModulus/10.0;
			}

			/** get sum of volume*fraction*factor to calculate vol% from mol% */
			sum_volume += cp[i].volume*cp[i].ss_n_mol*cp[i].factor;

			if (strcmp( cp[i].name, "liq") != 0 && strcmp( cp[i].name, "fl") != 0){
				sum_volume_sol 		+= cp[i].volume*cp[i].ss_n_mol*cp[i].factor;
				gv.solid_fraction 	+= cp[i].ss_n_mol;
			}

		}
	}

	for (int i = 0; i < gv.len_pp; i++){
		/* if pure phase is active or on hold (PP cannot be removed from consideration */
		if (gv.pp_flags[i][1] == 1 && gv.pp_flags[i][4] == 0){

			/* calculate phase volume as V = dG/dP */
			PP_db    	 = G_EM_function(gv.research_group, gv.EM_dataset, gv.len_ox,z_b.id,z_b.bulk_rock, z_b.apo, 1., z_b.T, gv.PP_list[i], "equilibrium");
			muC0	 	 = PP_db.gbase;

			PP_db    	 = G_EM_function(gv.research_group, gv.EM_dataset, gv.len_ox,z_b.id,z_b.bulk_rock, z_b.apo, 1. + gv.gb_P_eps, z_b.T, gv.PP_list[i], "equilibrium");
			muN0 	 	 = PP_db.gbase;

			PP_db    	 = G_EM_function(gv.research_group, gv.EM_dataset, gv.len_ox,z_b.id,z_b.bulk_rock, z_b.apo, z_b.P, z_b.T, gv.PP_list[i], "equilibrium");
			muC	 	 	 = PP_db.gbase;

			PP_db    	 = G_EM_function(gv.research_group, gv.EM_dataset, gv.len_ox,z_b.id,z_b.bulk_rock, z_b.apo, z_b.P + gv.gb_P_eps, z_b.T, gv.PP_list[i], "equilibrium");
			muN 	 	 = PP_db.gbase;

			PP_db    	 = G_EM_function(gv.research_group, gv.EM_dataset, gv.len_ox,z_b.id,z_b.bulk_rock, z_b.apo, z_b.P + gv.gb_P_eps*2.0, z_b.T, gv.PP_list[i], "equilibrium");
			muNN 	 	 = PP_db.gbase;

			PP_db    	 = G_EM_function(gv.research_group, gv.EM_dataset, gv.len_ox,z_b.id,z_b.bulk_rock, z_b.apo, z_b.P + gv.gb_P_eps*3.0, z_b.T, gv.PP_list[i], "equilibrium");
			muNNN 	 	 = PP_db.gbase;

			PP_db    	 = G_EM_function(gv.research_group, gv.EM_dataset, gv.len_ox,z_b.id,z_b.bulk_rock, z_b.apo, z_b.P, z_b.T + gv.gb_T_eps, gv.PP_list[i], "equilibrium");
			muE	 		 = PP_db.gbase;

			PP_db    	 = G_EM_function(gv.research_group, gv.EM_dataset, gv.len_ox,z_b.id,z_b.bulk_rock, z_b.apo, z_b.P, z_b.T - gv.gb_T_eps, gv.PP_list[i], "equilibrium");			
			muW	 		 = PP_db.gbase;

			PP_db    	 = G_EM_function(gv.research_group, gv.EM_dataset, gv.len_ox,z_b.id,z_b.bulk_rock, z_b.apo, z_b.P + gv.gb_P_eps, z_b.T + gv.gb_T_eps, gv.PP_list[i], "equilibrium");
			muNE	 	 = PP_db.gbase;

			PP_db    	 = G_EM_function(gv.research_group, gv.EM_dataset, gv.len_ox,z_b.id,z_b.bulk_rock, z_b.apo, z_b.P + gv.gb_P_eps, z_b.T - gv.gb_T_eps, gv.PP_list[i], "equilibrium");			
			muNW	 	 = PP_db.gbase;

			/* Calculate mass per pure phase */
			PP_ref_db[i].mass = 0.0;
			for (int j = 0; j< gv.len_ox; j++){
				PP_ref_db[i].mass += PP_ref_db[i].Comp[j]*z_b.masspo[j];
			}

			dG2dT2 		= (muE-2.0*muC+muW)		/(gv.gb_T_eps*gv.gb_T_eps);
			dG2dP2 		= (muNN-2.0*muN+muC)	/(gv.gb_P_eps*gv.gb_P_eps);
			dG2dP2_N	= (muNNN-2.0*muNN+muN)	/(gv.gb_P_eps*gv.gb_P_eps);
			dGdTPP 		= (muNE-muNW)			/(2.0*gv.gb_T_eps);
			dGdTMP 		= (muE-muW)				/(2.0*gv.gb_T_eps);
			dGdP		= (muN-muC)				/(gv.gb_P_eps);
			dGdP_N		= (muNN-muN)			/(gv.gb_P_eps);
			// dGdP_P0 	= (muN0-muC0)			/(gv.gb_P_eps);
			
			/* Calculate volume  per pure phase */
			PP_ref_db[i].volume  	   		= dGdP; 

			/* Calculate volume  per pure phase */
			// PP_ref_db[i].volume_P0  	    = dGdP_P0; 
			
			/* Calculate density per pure phase */
			PP_ref_db[i].phase_density 		= (1000.0*PP_ref_db[i].mass)/(PP_ref_db[i].volume*10.0);
			
			/* calculate cp of pure phase */
			PP_ref_db[i].phase_cp 			= -T*(dG2dT2);
			
			/* expansivity 		*/
			PP_ref_db[i].phase_expansivity 	= 1.0/(dGdP)*((dGdTPP-dGdTMP)/(gv.gb_P_eps));
			
			/* entropy 		*/
			PP_ref_db[i].phase_entropy 		= -dGdTMP*PP_ref_db[i].factor;
			
			/* enthalpy   		*/
			PP_ref_db[i].phase_enthalpy 	= PP_ref_db[i].phase_entropy*T + PP_ref_db[i].gbase*PP_ref_db[i].factor;;
	
			/* bulk modulus	*/
			if ( strcmp(gv.research_group, "sb") 	== 0 ){
				
			}
			else{
				PP_ref_db[i].phase_bulkModulus	= -dGdP/( dG2dP2 + pow(((dGdTPP-dGdTMP)/(gv.gb_P_eps)),2.0)/dG2dT2 );
			}

			/* shear modulus	*/
			PP_ref_db[i].phase_isoTbulkModulus	= -dGdP/( dG2dP2  );
			phase_isoTbulkModulus_P1			= -dGdP_N/( dG2dP2_N  );
	
			mut 				 = (3.0*PP_ref_db[i].phase_isoTbulkModulus - 6.0*PP_ref_db[i].phase_isoTbulkModulus*gv.poisson_ratio) / (2. + 2.0*gv.poisson_ratio)/10.0;
			mut_N 				 = (3.0*phase_isoTbulkModulus_P1 	- 6.0*phase_isoTbulkModulus_P1*gv.poisson_ratio) 	/ (2. + 2.0*gv.poisson_ratio)/10.0;

			PP_ref_db[i].thetaExp= (mut_N - mut)/gv.gb_P_eps - (PP_ref_db[i].phase_bulkModulus*PP_ref_db[i].phase_expansivity)/(PP_ref_db[i].phase_cp*PP_ref_db[i].phase_density);

			/** get sum of volume*fraction*factor to calculate vol% from mol% */
			sum_volume 			+= PP_ref_db[i].volume*gv.pp_n_mol[i]*PP_ref_db[i].factor;

			if (strcmp( gv.PP_list[i], "H2O") != 0){
				sum_volume_sol 		+= PP_ref_db[i].volume*gv.pp_n_mol[i]*PP_ref_db[i].factor;
				gv.solid_fraction 	+= gv.pp_n_mol[i];
			}
		}
	}




	/** 
	 calculate the bulk and shear modulus of the aggregate using the Voigt-Reuss-Hill averaging scheme with a weighting factor of 0.5 
	*/
	double s1 = 0.0; double b1 = 0.0;
	double s2 = 0.0; double b2 = 0.0;
	double s1S = 0.0; double b1S = 0.0;
	double s2S = 0.0; double b2S = 0.0;

	for (int i = 0; i < gv.len_cp; i++){
		if (cp[i].ss_flags[1] == 1){
			s1 +=  cp[i].volume*cp[i].ss_n_mol*cp[i].factor/sum_volume *  (cp[i].phase_shearModulus/10.0);
			s2 += (cp[i].volume*cp[i].ss_n_mol*cp[i].factor/sum_volume) / (cp[i].phase_shearModulus/10.0);
			b1 +=  cp[i].volume*cp[i].ss_n_mol*cp[i].factor/sum_volume *  (cp[i].phase_bulkModulus /10.0);
			b2 += (cp[i].volume*cp[i].ss_n_mol*cp[i].factor/sum_volume) / (cp[i].phase_bulkModulus /10.0);
			if (strcmp( cp[i].name, "liq") != 0 && strcmp( cp[i].name, "fl") != 0){
				s1S +=  cp[i].volume*cp[i].ss_n_mol*cp[i].factor/sum_volume_sol *  (cp[i].phase_shearModulus/10.0);
				s2S += (cp[i].volume*cp[i].ss_n_mol*cp[i].factor/sum_volume_sol) / (cp[i].phase_shearModulus/10.0);
				b1S +=  cp[i].volume*cp[i].ss_n_mol*cp[i].factor/sum_volume_sol *  (cp[i].phase_bulkModulus /10.0);
				b2S += (cp[i].volume*cp[i].ss_n_mol*cp[i].factor/sum_volume_sol) / (cp[i].phase_bulkModulus /10.0);
			}
		}
	}
	for (int i = 0; i < gv.len_pp; i++){
		if (gv.pp_flags[i][1] == 1  && gv.pp_flags[i][4] == 0){
			s1 +=  PP_ref_db[i].volume*gv.pp_n_mol[i]*PP_ref_db[i].factor/sum_volume *  (PP_ref_db[i].phase_shearModulus/10.0);
			s2 += (PP_ref_db[i].volume*gv.pp_n_mol[i]*PP_ref_db[i].factor/sum_volume) / (PP_ref_db[i].phase_shearModulus/10.0);
			b1 +=  PP_ref_db[i].volume*gv.pp_n_mol[i]*PP_ref_db[i].factor/sum_volume *  (PP_ref_db[i].phase_bulkModulus /10.0);
			b2 += (PP_ref_db[i].volume*gv.pp_n_mol[i]*PP_ref_db[i].factor/sum_volume) / (PP_ref_db[i].phase_bulkModulus /10.0);

			s1S +=  PP_ref_db[i].volume*gv.pp_n_mol[i]*PP_ref_db[i].factor/sum_volume_sol *  (PP_ref_db[i].phase_shearModulus/10.0);
			s2S += (PP_ref_db[i].volume*gv.pp_n_mol[i]*PP_ref_db[i].factor/sum_volume_sol) / (PP_ref_db[i].phase_shearModulus/10.0);
			b1S +=  PP_ref_db[i].volume*gv.pp_n_mol[i]*PP_ref_db[i].factor/sum_volume_sol *  (PP_ref_db[i].phase_bulkModulus /10.0);
			b2S += (PP_ref_db[i].volume*gv.pp_n_mol[i]*PP_ref_db[i].factor/sum_volume_sol) / (PP_ref_db[i].phase_bulkModulus /10.0);
		}
	}

	// Voight-Reuss-Hill averaging
	gv.system_shearModulus 	= 0.50 * s1 + 0.50 * (1.0/(s2));
	gv.system_bulkModulus  	= 0.50 * b1 + 0.50 * (1.0/(b2));

	gv.solid_shearModulus 	= 0.50 * s1S + 0.50 * (1.0/(s2S));
	gv.solid_bulkModulus  	= 0.50 * b1S + 0.50 * (1.0/(b2S));

	/* calculate density of the system */
	for (int i = 0; i < gv.len_cp; i++){
		if (cp[i].ss_flags[1] == 1){
			gv.system_density += cp[i].phase_density*cp[i].ss_n_wt;
			gv.system_entropy += cp[i].phase_entropy*cp[i].ss_n_mol;//*cp[i].factor;
			gv.system_cp 	  += cp[i].phase_cp*cp[i].ss_n_mol;
			gv.system_expansivity 	  += cp[i].phase_expansivity*cp[i].ss_n_mol;
			if (strcmp( cp[i].name, "liq") != 0 && strcmp( cp[i].name, "fl") != 0){
				gv.solid_density += cp[i].phase_density*cp[i].ss_n_wt;
			}
		}
	}
	for (int i = 0; i < gv.len_pp; i++){
		if (gv.pp_flags[i][1] == 1 && gv.pp_flags[i][4] == 0){
			gv.system_density += PP_ref_db[i].phase_density*gv.pp_n_wt[i];
			gv.system_entropy += PP_ref_db[i].phase_entropy*gv.pp_n_mol[i];//*PP_ref_db[i].factor;		
			gv.system_cp 	  += PP_ref_db[i].phase_cp*gv.pp_n_mol[i];	
			gv.system_expansivity 	  += PP_ref_db[i].phase_expansivity*gv.pp_n_mol[i];	
			if (strcmp( gv.PP_list[i], "H2O") != 0){		
				gv.solid_density  += PP_ref_db[i].phase_density*gv.pp_n_wt[i];
			}
		}
	}

	gv.system_volume = sum_volume;
	G = 0.0;
	for (int j = 0; j < gv.len_ox; j++){
		G += z_b.bulk_rock[j]*gv.gam_tot[j];
	}
	gv.system_enthalpy = gv.system_entropy*T + G;

	gv.system_Vp 	= sqrt((gv.system_bulkModulus +4.0/3.0*gv.system_shearModulus)/(gv.system_density/1e3));
	gv.system_Vs 	= sqrt(gv.system_shearModulus/(gv.system_density/1e3));
	gv.solid_Vp 	= sqrt((gv.solid_bulkModulus +4.0/3.0*gv.solid_shearModulus)/(gv.solid_density/1e3));
	gv.solid_Vs 	= sqrt(gv.solid_shearModulus/(gv.solid_density/1e3));

	// gv.V_cor[0] 	= gv.solid_Vp;
	// gv.V_cor[1] 	= gv.solid_Vs;
	// if (gv.calc_seismic_cor == 1){
	// 	gv.solid_Vs 	= anelastic_correction( 0,
	// 											gv.solid_Vs,
	// 											z_b.P,
	// 											z_b.T 		);
	// 	gv.V_cor[0] 	= gv.solid_Vp;
	// 	gv.V_cor[1] 	= gv.solid_Vs;
	// 	gv = wave_melt_correction(  	gv,
	// 									z_b,
	// 									0.1				);
	// }

	return gv;
}


global_variable compute_activities(			int					 EM_database,	
											global_variable 	 gv,
											PP_ref  			*PP_ref_db,
											bulk_info 			 z_b			){

	PP_ref PP_db;	

	/** calculate oxygen fugacity: mu_O2 = G0_O2 + RTlog(fO2) */
	/* get O2 Gibbs energy of reference */
	double G0_O = 0.0;
	for (int i = 0; i < gv.len_pp; i++){
		if	(strcmp( gv.PP_list[i], "O2") == 0){
			G0_O = PP_ref_db[i].gbase;//*PP_ref_db[i].factor;
			break;
		}
	}
	/* get chemical potential of Oxygen (index)*/
	int O_ix = -1;
	for (int i = 0; i < gv.len_ox; i++){
		if	(strcmp( gv.ox[i], "O") == 0){
			O_ix = i;
			break;
		}
	}
	if (O_ix != -1){
		gv.system_fO2 = exp( (gv.gam_tot[O_ix]*2.0 - G0_O) / (z_b.R*z_b.T)) ;

		PP_ref q 	= G_EM_function(	gv.research_group, gv.EM_dataset, 
										gv.len_ox,
										z_b.id,
										z_b.bulk_rock, 
										z_b.apo, 
										z_b.P, 
										z_b.T, 
										"q", 
										"equilibrium"				);
		PP_ref fa 	= G_EM_function(	gv.research_group, gv.EM_dataset, 
										gv.len_ox,
										z_b.id,
										z_b.bulk_rock, 
										z_b.apo, 
										z_b.P, 
										z_b.T, 
										"fa", 
										"equilibrium"				);

		PP_ref mt 	= G_EM_function(	gv.research_group, gv.EM_dataset, 
										gv.len_ox,
										z_b.id,
										z_b.bulk_rock, 
										z_b.apo, 
										z_b.P, 
										z_b.T, 
										"mt", 
										"equilibrium"				);

		double gO_qfm   	=  -3.0 * fa.gbase + 3.0*q.gbase + 2.0*mt.gbase;
		gv.system_deltaQFM 	= exp( (gv.gam_tot[O_ix]*2.0 - gO_qfm) / (z_b.R*z_b.T));

		// double  QFM  =  42.14743 - 27792.8/z_b.T + 104.59*z_b.P/z_b.T - 4.7113*log(z_b.T) + 0.00180*z_b.T;
		// printf("fO2_Dqfm_an %g\n", log10(gv.system_fO2)-QFM);

	}
	else {
		if (gv.verbose == 1){
			printf(" INFO: no fO2 -> O oxide not part of the database\n");
		}
	}

	/* compute activities for pure component phases */
	/* get chemical potential of pure components (index)*/
	int H2O_ix 	= -1;
	int TiO2_ix = -1;
	int SiO2_ix = -1;
	int Al2O3_ix = -1;
	int MgO_ix = -1;
	int FeO_ix = -1;
	for (int i = 0; i < gv.len_ox; i++){
		if	(strcmp( gv.ox[i], "H2O") 	   == 0 && z_b.bulk_rock[i] > 0.0){
			H2O_ix = i;
		}
		else if(strcmp( gv.ox[i], "TiO2")  == 0 && z_b.bulk_rock[i] > 0.0){
			TiO2_ix = i;
		}
		else if(strcmp( gv.ox[i], "SiO2")  == 0 && z_b.bulk_rock[i] > 0.0){
			SiO2_ix = i;
		}
		else if(strcmp( gv.ox[i], "Al2O3") == 0 && z_b.bulk_rock[i] > 0.0){
			Al2O3_ix = i;
		}
		else if(strcmp( gv.ox[i], "FeO")  == 0 && z_b.bulk_rock[i] > 0.0){
			FeO_ix = i;
		}
		else if(strcmp( gv.ox[i], "MgO") == 0 && z_b.bulk_rock[i] > 0.0){
			MgO_ix = i;
		}
	}

	/* if we can compute the activity of MgO (if Gamma MgO exists i.e., if the MgO is taken into account) */
	if (MgO_ix != -1){
		double G0_per = 0.0;
		if (strcmp(gv.research_group, "tc") 	== 0 ){
			PP_db  			= G_EM_function(gv.research_group, gv.EM_dataset, gv.len_ox,z_b.id,z_b.bulk_rock, z_b.apo, z_b.P, z_b.T , "per", "equilibrium");
		}
		else if (strcmp(gv.research_group, "sb") 	== 0 ){
			PP_db  			= G_EM_function(gv.research_group, gv.EM_dataset, gv.len_ox,z_b.id,z_b.bulk_rock, z_b.apo, z_b.P, z_b.T , "pe", "equilibrium");
		}
		G0_per  		= PP_db.gbase;//*PP_db.factor;
		gv.system_aMgO = exp( (gv.gam_tot[MgO_ix] - G0_per) / (z_b.R*z_b.T));
	}

	/* if we can compute the activity of FeO (if Gamma FeO exists i.e., if the FeO is taken into account) */
	if (FeO_ix != -1){
		double G0_fper  = 0.0;
		if (strcmp(gv.research_group, "tc") 	== 0 ){
			PP_db  			= G_EM_function(gv.research_group, gv.EM_dataset, gv.len_ox,z_b.id,z_b.bulk_rock, z_b.apo, z_b.P, z_b.T , "fper", "equilibrium");
		}
		else if (strcmp(gv.research_group, "sb") 	== 0 ){
			PP_db  			= G_EM_function(gv.research_group, gv.EM_dataset, gv.len_ox,z_b.id,z_b.bulk_rock, z_b.apo, z_b.P, z_b.T , "wu", "equilibrium");
		}
		G0_fper  		= PP_db.gbase;//*PP_db.factor;
		gv.system_aFeO  = exp( (gv.gam_tot[FeO_ix] - G0_fper) / (z_b.R*z_b.T));
	}


	/* if we can compute the activity of Al2O3 (if Gamma Al2O3 exists i.e., if the Al2O3 is taken into account) */
	if (Al2O3_ix != -1){
		double G0_cor = 0.0;
		if (strcmp(gv.research_group, "tc") 	== 0 ){
		PP_db  			= G_EM_function(gv.research_group, gv.EM_dataset, gv.len_ox,z_b.id,z_b.bulk_rock, z_b.apo, z_b.P, z_b.T , "cor", "equilibrium");
		}
		else if (strcmp(gv.research_group, "sb") 	== 0 ){
			PP_db  			= G_EM_function(gv.research_group, gv.EM_dataset, gv.len_ox,z_b.id,z_b.bulk_rock, z_b.apo, z_b.P, z_b.T , "co", "equilibrium");
		}
		G0_cor  		= PP_db.gbase;//*PP_db.factor;
		gv.system_aAl2O3 = exp( (gv.gam_tot[Al2O3_ix] - G0_cor) / (z_b.R*z_b.T));
	}

	/* if we can compute the activity of TiO2 (if Gamma TiO2 exists i.e., if the TiO2 is taken into account) */
	if (TiO2_ix != -1){
		double G0_ru = 0.0;
		PP_db  			= G_EM_function(gv.research_group, gv.EM_dataset, gv.len_ox,z_b.id,z_b.bulk_rock, z_b.apo, z_b.P, z_b.T , "ru", "equilibrium");
		G0_ru  			= PP_db.gbase;//*PP_db.factor;
		gv.system_aTiO2 = exp( (gv.gam_tot[TiO2_ix] - G0_ru) / (z_b.R*z_b.T));
	}

	/* if we can compute the activity of H2O (if Gamma H2O exists i.e., if the H2O is taken into account) */
	if (H2O_ix != -1){
		double G0_H2O = 0.0;
		PP_db  			= G_EM_function(gv.research_group, gv.EM_dataset, gv.len_ox,z_b.id,z_b.bulk_rock, z_b.apo, z_b.P, z_b.T , "H2O", "equilibrium");
		G0_H2O  		= PP_db.gbase;//*PP_db.factor;
		gv.system_aH2O  = exp( (gv.gam_tot[H2O_ix] - G0_H2O) / (z_b.R*z_b.T));
	}

	/* if we can compute the activity of H2O (if Gamma H2O exists i.e., if the H2O is taken into account) */
	if (SiO2_ix != -1){
		double G0_q 	= 0.0;
		double G0_coe 	= 0.0;
		double G0_st 	= 0.0;

		if (strcmp(gv.research_group, "tc") 	== 0 ){
			PP_db  			= G_EM_function(gv.research_group, gv.EM_dataset, gv.len_ox,z_b.id,z_b.bulk_rock, z_b.apo, z_b.P, z_b.T , "q", "equilibrium");
			G0_q  			= PP_db.gbase;//*PP_db.factor;
			PP_db  			= G_EM_function(gv.research_group, gv.EM_dataset, gv.len_ox,z_b.id,z_b.bulk_rock, z_b.apo, z_b.P, z_b.T , "coe", "equilibrium");
			G0_coe  		= PP_db.gbase;//*PP_db.factor;
		}
		else if (strcmp(gv.research_group, "sb") 	== 0 ){
			PP_db  			= G_EM_function(gv.research_group, gv.EM_dataset, gv.len_ox,z_b.id,z_b.bulk_rock, z_b.apo, z_b.P, z_b.T , "qtz", "equilibrium");
			G0_q  			= PP_db.gbase;//*PP_db.factor;
			PP_db  			= G_EM_function(gv.research_group, gv.EM_dataset, gv.len_ox,z_b.id,z_b.bulk_rock, z_b.apo, z_b.P, z_b.T , "coe", "equilibrium");
			G0_coe  		= PP_db.gbase;//*PP_db.factor;
			PP_db  			= G_EM_function(gv.research_group, gv.EM_dataset, gv.len_ox,z_b.id,z_b.bulk_rock, z_b.apo, z_b.P, z_b.T , "st", "equilibrium");
			G0_st  		= PP_db.gbase;//*PP_db.factor;
		}
		double G0_SiO2 	= G0_q;
		if (G0_coe < G0_SiO2){
			G0_SiO2 = G0_coe;
			if (G0_st < G0_coe){
				G0_SiO2 = G0_st;
			}
		}
		gv.system_aSiO2 = exp( (gv.gam_tot[SiO2_ix] - G0_SiO2) / (z_b.R*z_b.T));
	}

	return gv;
}

