#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <complex.h> 
#include <lapacke.h> 

#include "MAGEMin.h"
#include "gem_function.h"
#include "gss_function.h"
#include "NLopt_opt_function.h"
#include "simplex_levelling.h"
#include "objective_functions.h"
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
	printf("  --Verb=     [int]   : Verbose option, 0. inactive, 1. active\n");	
	printf("  --File=     [str]   : Given file for multiple point calculation\n");
	printf("  --n_points= [int]   : Number of points when using 'File' argument\n");
	printf("  --test=     [int]   : Number of points when using 'File' argument\n");
	printf("  --Pres=     [float] : Pressure in kilobar\n");
	printf("  --Temp=     [float] : Temperature in Celsius\n");
	printf("  --Bulk=     [float] : Bulk rock composition in molar amount*\n");
	printf("  --Gam=      [float] : Chemical potential of oxides (pure components)*\n");
	printf("\n");
	printf(" *the list of oxides must be given as follow:\n");
	printf("  SiO2, Al2O3, CaO, MgO, FeO, K2O, Na2O, TiO2, O, Cr2O3, H2O\n");
	printf("\n\n");
	printf(" Example of single point calculation:\n");
	printf(" ------------------------------------\n");	
	printf("\n");
    printf("  ./MAGEMin --Verb=1 --Temp=718.750 --Pres=30.5000 --test=0 >&log.txt\n");
    printf("\n");
	printf(" Here, the verbose is active and the bulk rock composition of 'test 0' is selected. The output of the verbose is saved as a log file 'log.txt'\n");
	printf("\n\n");
	printf(" Example multiple points calculation:\n");
	printf(" ------------------------------------\n");	
    printf("\n");
	printf(" To run multiple points at once you can pass an input file containing the list of points such as\n");
    printf("\n");
    printf("  ./MAGEMin --Verb=1 --File='path_to_file' --n_points=x\n");
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
    printf("  mpirun -np 8 ./MAGEMin --Verb=1 --File=/path_to_file/MAGEMin_input.dat --n_points=4\n");
    printf("  mpiexec -n 8 ./MAGEMin --Verb=1 --File=/path_to_file/MAGEMin_input.dat --n_points=4\n");	   
    printf("\n");	 
    printf("\n");	 

}

/**
  Get the number of max_pc
*/
int get_max_n_pc(int tot_pc, int n_pc){
	return ((tot_pc >= n_pc) ? (n_pc) : (tot_pc));
};

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

  i = 0; done = FALSE;   error = 0;
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
  create matrix typedef
*/
typedef struct TMatrix {
	double **m;
	int nRows; int nCols;
} TMATRIX;


/** 
  create matrix
*/
TMATRIX createMatrix (int nRows, int nCols) {
	TMATRIX oMatrix;
 
	oMatrix.nRows = nRows;
	oMatrix.nCols = nCols;
 
	oMatrix.m = malloc (nRows * sizeof (double*) ); 
	for (int i = 0; i < (nRows); i++){
		oMatrix.m[i] = malloc (nCols * sizeof (double) );
	}
	
	for (int i = 0; i < (nRows); i++){
		for (int j = 0; j < (nCols); j++){	
			oMatrix.m[i][j] = 0.0;
		}
	}
 
	return oMatrix;
}

/**
  Remove values less than tolerance
*/
void freeMatrix (TMATRIX oMatrix) {
	// Removes all numbers close to zero, i.e between -tol and +tol 
	for (int i = 0; i < oMatrix.nRows; i++){
		free(oMatrix.m[i]);
	}
	free(oMatrix.m);
}
 
/**
  Remove values less than tolerance
*/
void cleanUpMatrix (TMATRIX oMatrix,  double tolerance) {
	
	// Removes all numbers close to zero, i.e between -tol and +tol 
	for (int i = 0; i < oMatrix.nRows; i++){
		for (int j = 0; j < oMatrix.nCols; j++){
			if (fabs (oMatrix.m[i][j]) < tolerance){
				oMatrix.m[i][j] = 0;
			}
		}
	}
}
 
/**
  exchange rows in the matrix
*/
void exchangeRows (TMATRIX oMatrix, int r1, int r2) {
 
	double t = 0;
	for (int i = 0; i < oMatrix.nCols; i++) {
		t = oMatrix.m[r1][i];
		oMatrix.m[r1][i] = oMatrix.m[r2][i];
		oMatrix.m[r2][i] = t;
	}
}
 
/**
  Mainly for debugging, helps spot out of range errors when indexing elements
*/
double getValue (TMATRIX oMatrix, int i, int j) {
	
	if ((i < 0) || (j < 0)) {
		printf ("Error in indexing\n");
		getchar ();
		exit (0);
	}
 
	if ((i >= oMatrix.nRows) || (j >= oMatrix.nCols)) {
		printf ("Error in indexing: %d, %d\n", i, j);
		getchar ();
		exit (0);
	}
 
	return oMatrix.m[i][j];
}


/**
  main rref routine (reduced row echelon form)
  - used to calculate the potential reactions between endmembers satisfying the mass constraint 
*/
TMATRIX rref(TMATRIX oMatrix, int *pivot, double tolerance) {
	int currentRow; double factor;
 
	TMATRIX oEchelon = createMatrix (oMatrix.nRows, oMatrix.nCols);
 
	// Make a copy and work on that.
	for (int i = 0; i < oMatrix.nRows; i++){
		for (int j = 0; j < oMatrix.nCols; j++){
			oEchelon.m[i][j] = oMatrix.m[i][j];
		}
	}
 
	int Arow = 0; int Acol = 0; int pvt = 0;
	while ((Arow < oEchelon.nRows) && (Acol < oEchelon.nCols)) {
		// locate a nonzero column
		if (abs (getValue (oEchelon, Arow, Acol) < tolerance)) {
			// If the entry is zero work our way down the matrix
			// looking for a nonzero entry, when found, swap it for Arow 
			currentRow = Arow;
			do {
				// next row
				currentRow++;
				// Have we reached the end of the rows but we've still got columns left to scan?
				if ((currentRow >= oEchelon.nRows) && (Acol <= oEchelon.nCols)) {
					// reset row counter back to where it was and try next column 
					currentRow = Arow; Acol++;
				}
 
				// If we've scanned the whole matrix, then lets get out... 
				if (currentRow >= oEchelon.nRows) {
					cleanUpMatrix (oEchelon, tolerance);
					return oEchelon;
				}
			} while (fabs (getValue (oEchelon, currentRow, Acol)) < tolerance);
 
			// We've found a nonzero row entry so swap it with 'Arow' which did have a zero as its entry 
			exchangeRows (oEchelon, Arow, currentRow);
		}
		// Arow now holds the row of interest }
		factor = 1.0 / getValue (oEchelon, Arow, Acol);
		pivot[pvt] = Acol;
		pvt 	  += 1;
		//printf(" pivot %d\n",Acol);
		// reduce all the entries along the column by the factor 
		for (int i = Acol; i < oEchelon.nCols; i++)
			oEchelon.m[Arow][i] = getValue (oEchelon, Arow, i) * factor;
 
		// now eliminate all entries above and below Arow, this generates the reduced form 
		for (int i = 0; i < oEchelon.nRows; i++) {
			// miss out Arow itself 
			if ((i != Arow) && (fabs (getValue (oEchelon, i, Acol)) > tolerance)) {
				factor = getValue (oEchelon, i, Acol);
				// work your way along the column doing the same operation 
				for (int j = Acol; j < oEchelon.nCols; j++) {
					oEchelon.m[i][j] = getValue (oEchelon, i, j) - factor * getValue (oEchelon, Arow, j);
				}
			}
		}
		Arow++; Acol++;
	}
	cleanUpMatrix (oEchelon, tolerance);
	
	// Make a copy and work on that.
	for (int i = 0; i < oMatrix.nRows; i++){
		for (int j = 0; j < oMatrix.nCols; j++){
			oMatrix.m[i][j] = oEchelon.m[i][j];
		}
	}
	freeMatrix(oEchelon);
	
	return oMatrix;
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
void inverseMatrix(double *A1, int n){
	int    ipiv[n];					
	int    info;
	
	/* call lapacke to inverse Matrix */
	info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, A1, n, ipiv); 
	info = LAPACKE_dgetri(LAPACK_ROW_MAJOR, n, A1, n, ipiv);
};

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
	int k;
	for (int i = 0; i < n; i++){
		B1[i] = 0.0;
		for (int j = 0; j < n; j++){
			k = j + i*n;
			B1[i] += B[j]*A1[k];
		}
	}
};

/**
  matrix vector multiplication
*/	
void MatVecMul(double *A1, double *br, double *n_vec, int n){
	
	int k;
	for (int i = 0; i < n; i++){
		n_vec[i] = 0.0;
		for (int j = 0; j < n; j++){
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
		for (int j = 0; j < gv.n_solvi[i]; j++){	
			printf(" %4s %d",cp[gv.id_solvi[i][j]].name,gv.id_solvi[i][j]);
		}
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


/**
   rotate G-hyperplane using Gamma
*/
void print_SS_informations(		global_variable gv,
								SS_ref SS_ref_db,
								int		iss					){
	printf(" %4s  | %+10f | %2d | %+10f | %+10f | ",gv.SS_list[iss],SS_ref_db.df,SS_ref_db.sf_ok,SS_ref_db.sum_xi,SS_ref_db.LM_time);
	for (int k = 0; k < SS_ref_db.n_xeos; k++) {
		printf(" %+10f",SS_ref_db.xeos[k]);
	}
	for (int k = SS_ref_db.n_xeos; k < 11; k++){
		printf(" %10s","-");
	}

	printf(" | ");
	for (int k = 0; k < SS_ref_db.n_xeos; k++) {
		printf(" %+10f",SS_ref_db.dfx[k]);
	}
	printf("\n");
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
									SS_ref SS_ref_db,
									double box_size		){
									
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
void wave_melt_correction( 	 double  Kb_L,
							 double  Kb_S,
							 double  Ks_S,
							 double  rhoL,
							 double  rhoS,
							 double  Vp0,
							 double  Vs0,
							 double  meltFrac,
							 double  aspectRatio,
							 double *V_cor		)
{
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

    double ksk     = ksk_k*Kb_S;
    double musk    = musk_mu*Ks_S;

    double kb      = (1.0-meltFrac)*ksk;
    double mu      = (1.0-meltFrac)*musk;

    double LambdaK = Kb_S/kb;
    double LambdaG = Ks_S/mu;

    double beta    = Kb_S/Kb_L;
    double gamma   = Ks_S/Kb_S;

	double deltaVp = (((((beta -1.0)*LambdaK) / ((beta-1.0) + LambdaK) + 4.0/3.0*gamma*LambdaG ) / ( 1.0 + 4.0/3.0*gamma)) - (1.0 - rhoL/rhoS) )*(meltFrac/2.0);
	double deltaVs = ( LambdaG - (1.0 - rhoL/rhoS) )*(meltFrac/2.0);

    V_cor[0] = Vp0 - deltaVp*Vp0;
    V_cor[1] = Vs0 - deltaVs*Vs0;

	if (V_cor[0] < 0.0){ V_cor[0] = 0.0;}
	if (V_cor[1] < 0.0){ V_cor[1] = 0.0;}

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
    double omega   = 0.01;                 	//Hz (frequency to match for studied seismic system)
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


/**
  function to calculate delta G (position of potential phases with current G hyperplane*)
*/
void update_dG(	simplex_data *splx_data
){

	simplex_data *d  = (simplex_data *) splx_data;
	double F;
	double minF = 0.0;

	VecMatMul(	d->B1,
				d->A1,
				d->B,
				d->n_Ox	);
	
	d->dG_B = d->g0_B;
	for (int j = 0; j < d->n_Ox; j++){
		d->dG_B -= d->B1[j]*d->g0_A[j];
	}
	
	d->ph2swp = -1;	
	if (d->dG_B < -1e-6){	
		d->min_F  =  1e6;												/** max value for F, tentative here because F can tend to +inf */
		for (int i = 0; i < d->n_Ox; i++){
			F = (d->n_vec[i])/d->B1[i];
			if (F < d->min_F && F > minF){
				d->min_F  = F;
				d->ph2swp = i;
			}
		}
	}
}



/**
  update Gamma using LAPACKE dgesv
*/	
void update_local_gamma(	double *A1, 
							double *g0_A, 
							double *gam, 
							int 	n			){
	int k;
	for (int i = 0; i < n; i++){
		gam[i] = 0.0;
		for (int j = 0; j < n; j++){
			k = i + j*n;
			gam[i] += g0_A[j]*A1[k];
		}
	}
};

/**
  update global Gamma 
*/	
void update_global_gamma( 				struct bulk_info 	z_b,
										simplex_data 	   *splx_data
){

	simplex_data *d  = (simplex_data *) splx_data;

	/** update gam_tot using solution phase levelling gamma */
	for (int i = 0; i < d->n_Ox; i++){
		d->gamma_delta[z_b.nzEl_array[i]] 	= d->gamma_ss[i] - d->gamma_tot[z_b.nzEl_array[i]];
		d->gamma_tot[z_b.nzEl_array[i]]     = d->gamma_ss[i];
	}	
	
};


/**
  function to swp pure phases
*/	
void swap_pure_phases(				struct bulk_info 	 z_b,
									simplex_data 		*splx_data,
									global_variable 	 gv,
									
									PP_ref 				*PP_ref_db,
									SS_ref 				*SS_ref_db
){
	simplex_data *d  = (simplex_data *) splx_data;
	int     k;
	
	for (int i = 0; i < gv.len_pp; i++){
		if (gv.pp_flags[i][0] == 1){
			
			d->g0_B 		= PP_ref_db[i].gbase*PP_ref_db[i].factor;
			d->ph_id_B[0]  	= 1;															/** added phase is a pure species */
			d->ph_id_B[1]	= i;															/** save pure species index */
		
			/* retrieve the composition in the right (reduced) chemical space */
			for (int j = 0; j < z_b.nzEl_val; j++){
				d->B[j] = PP_ref_db[i].Comp[z_b.nzEl_array[j]]*PP_ref_db[i].factor; 
			}
			
			/** update deltaG with respect to G hyperplane */
			update_dG(splx_data);
			
			/** swap phase */
			if (d->ph2swp != -1){															/** if the phase can be added */
				d->swp 						   = 1;
				d->n_swp 					  += 1;
				d->ph_id_A[d->ph2swp][0] = d->ph_id_B[0];
				d->ph_id_A[d->ph2swp][1] = d->ph_id_B[1];
				d->ph_id_A[d->ph2swp][2] = d->ph_id_B[2];
				d->g0_A[d->ph2swp] 	   = d->g0_B;
				
				for (int j = 0; j < d->n_Ox; j++){				
					k = d->ph2swp + j*d->n_Ox;
					d->A[k] = d->B[j];
				}
				
				for (int k = 0; k < d->n_Ox*d->n_Ox; k++){ d->A1[k] = d->A[k];}

				/** inverse guessed assemblage stoechiometry matrix */
				inverseMatrix(	d->A1,
								d->n_Ox		);
				
				/** update phase fractions */
				MatVecMul(		d->A1,
								// br,
								z_b.bulk_rock_cat,
								d->n_vec,
								d->n_Ox		);	
			}
		}
	}
	
}

/**
  function to swp pure endmembers
*/	
void swap_pure_endmembers(				struct bulk_info 	 z_b,
										simplex_data 		*splx_data,
										global_variable 	 gv,
										
										PP_ref 				*PP_ref_db,
										SS_ref 				*SS_ref_db
){

	simplex_data *d  = (simplex_data *) splx_data;

	int     k;
	double factor;

	for (int i = 0; i < gv.len_ss; i++){												/**loop to pass informations from active endmembers */
		if (SS_ref_db[i].ss_flags[0] == 1){														/** if SS is not filtered out then continue */

			for (int l = 0; l < SS_ref_db[i].n_em; l++){	
				/** if bulk-rock satisfy the compositions of endmembers, retrieve their informations */
				if (SS_ref_db[i].z_em[l] == 1.0){ 
					
					/* update normalizing factor for solution models than need it */
					factor 	= z_b.fbc/SS_ref_db[i].ape[l];	


					d->g0_B		 = SS_ref_db[i].gbase[l]*factor;
					d->ph_id_B[0] = 2;														/** added phase is a pure species */
					d->ph_id_B[1] = i;														/** save pure species index */
					d->ph_id_B[2] = 0;														/** save used initial guess */	
					
					/** retrieve the composition in the right (reduced) chemical space */
					for (int j = 0; j < z_b.nzEl_val; j++){
						d->B[j] 	 = SS_ref_db[i].Comp[l][z_b.nzEl_array[j]]*factor; 
					}

					/** update deltaG with respect to G hyperplane */
					update_dG(splx_data);
									
					/** swap phase */
					if (d->ph2swp != -1){													/** if the phase can be added */
						d->swp 						   = 1;
						d->n_swp 					  += 1;
						d->ph_id_A[d->ph2swp][0] = d->ph_id_B[0];
						d->ph_id_A[d->ph2swp][1] = d->ph_id_B[1];
						d->ph_id_A[d->ph2swp][2] = d->ph_id_B[2];
						d->ph_id_A[d->ph2swp][3] = l;	
						d->g0_A[d->ph2swp] 	   = d->g0_B;
						
						for (int j = 0; j < d->n_Ox; j++){				
							int k = d->ph2swp + j*d->n_Ox;
							d->A[k] = d->B[j];
						}
						for (int k = 0; k < d->n_Ox*d->n_Ox; k++){ d->A1[k] = d->A[k];}

						/** inverse guessed assemblage stoechiometry matrix */
						inverseMatrix(	d->A1,
										d->n_Ox		);
						
						/** update phase fractions */
						MatVecMul(		d->A1,
										// br,
										z_b.bulk_rock_cat,
										d->n_vec,
										d->n_Ox		);	
					}
				}
			}
		}
	}
	
}

/**
  function to swp solution phase pseudocompounds
*/	
void swap_pseudocompounds(				struct bulk_info 	 z_b,
										simplex_data 		*splx_data,
										global_variable 	 gv,
										
										PP_ref 				*PP_ref_db,
										SS_ref 				*SS_ref_db
){
	simplex_data *d  = (simplex_data *) splx_data;

	int     k, max_n_pc;
	
	for (int i = 0; i < gv.len_ss; i++){										/**loop to pass informations from active endmembers */
		if (SS_ref_db[i].ss_flags[0] == 1){

			max_n_pc = get_max_n_pc(SS_ref_db[i].tot_pc, SS_ref_db[i].n_pc);

			for (int l = 0; l < max_n_pc; l++){

				d->g0_B		  = SS_ref_db[i].G_pc[l];
				d->ph_id_B[0] = 3;														/** added phase is a pure species */
				d->ph_id_B[1] = i;														/** save solution phase index */
				d->ph_id_B[2] = 0;														/** save pseudocompound index */
			
				/* retrieve the composition in the right (reduced) chemical space */
				for (int j = 0; j < z_b.nzEl_val; j++){
					d->B[j] = SS_ref_db[i].comp_pc[l][z_b.nzEl_array[j]]; 
				}

				/** update deltaG with respect to G hyperplane */
				update_dG(splx_data);
				
				/** save updated driving force */
				SS_ref_db[i].DF_pc[l] = d->dG_B;
				
				/** swap phase */
				if (d->ph2swp != -1){													/** if the phase can be added */
					SS_ref_db[i].n_swap[l]   = d->n_swp;
					d->swp 				 	 = 1;
					d->n_swp 				+= 1;
					d->ph_id_A[d->ph2swp][0] = d->ph_id_B[0];
					d->ph_id_A[d->ph2swp][1] = d->ph_id_B[1];
					d->ph_id_A[d->ph2swp][2] = d->ph_id_B[2];
					d->ph_id_A[d->ph2swp][3] = l;									/** save pseudocompound number */
					d->g0_A[d->ph2swp] 	     = d->g0_B;
					
					for (int j = 0; j < d->n_Ox; j++){				
						int k = d->ph2swp + j*d->n_Ox;
						d->A[k] = d->B[j];
					}
					
					for (int k = 0; k < d->n_Ox*d->n_Ox; k++){ d->A1[k] = d->A[k];}

					/** inverse guessed assemblage stoechiometry matrix */
					inverseMatrix(	d->A1,
									d->n_Ox		);
					
					/** update phase fractions */
					MatVecMul(		d->A1,
									// br,
									z_b.bulk_rock_cat,
									d->n_vec,
									d->n_Ox		);
				}
			}
		}
	}
}

/**
  function to swp solution phase pseudocompounds during PGE iterations
*/	
void swap_PGE_pseudocompounds(			struct bulk_info 	 z_b,
										simplex_data 		*splx_data,
										global_variable 	 gv,
										
										PP_ref 				*PP_ref_db,
										SS_ref 				*SS_ref_db
){
	simplex_data *d  = (simplex_data *) splx_data;

	int     k, max_n_Ppc;
	
	for (int i = 0; i < gv.len_ss; i++){										/**loop to pass informations from active endmembers */
		if (SS_ref_db[i].ss_flags[0] == 1){

			max_n_Ppc = get_max_n_pc(SS_ref_db[i].tot_Ppc, SS_ref_db[i].n_Ppc);

			for (int l = 0; l < max_n_Ppc; l++){

				d->g0_B		  = SS_ref_db[i].G_Ppc[l];
				d->ph_id_B[0] = 3;														/** added phase is a pure species */
				d->ph_id_B[1] = i;														/** save solution phase index */
				d->ph_id_B[2] = 0;														/** save pseudocompound index */
			
				/* retrieve the composition in the right (reduced) chemical space */
				for (int j = 0; j < z_b.nzEl_val; j++){
					d->B[j] = SS_ref_db[i].comp_Ppc[l][z_b.nzEl_array[j]]; 
				}

				/** update deltaG with respect to G hyperplane */
				update_dG(splx_data);
				
				/** save updated driving force */
				SS_ref_db[i].DF_Ppc[l] = d->dG_B;
				
				/** swap phase */
				if (d->ph2swp != -1){													/** if the phase can be added */
					SS_ref_db[i].n_swap[l]   = d->n_swp;
					d->swp 				 	 = 1;
					d->n_swp 				+= 1;
					d->ph_id_A[d->ph2swp][0] = d->ph_id_B[0];
					d->ph_id_A[d->ph2swp][1] = d->ph_id_B[1];
					d->ph_id_A[d->ph2swp][2] = d->ph_id_B[2];
					d->ph_id_A[d->ph2swp][3] = l;									/** save pseudocompound number */
					d->g0_A[d->ph2swp] 	     = d->g0_B;
					
					for (int j = 0; j < d->n_Ox; j++){				
						int k = d->ph2swp + j*d->n_Ox;
						d->A[k] = d->B[j];
					}
					
					for (int k = 0; k < d->n_Ox*d->n_Ox; k++){ d->A1[k] = d->A[k];}

					/** inverse guessed assemblage stoechiometry matrix */
					inverseMatrix(	d->A1,
									d->n_Ox		);
					
					/** update phase fractions */
					MatVecMul(		d->A1,
									// br,
									z_b.bulk_rock_cat,
									d->n_vec,
									d->n_Ox		);
				}
			}
		}
	}
}


/**
  function to run simplex linear programming during PGE with pseudocompounds 
*/	
void run_simplex_PGE_pseudocompounds(	struct bulk_info 	 z_b,
										simplex_data 		*splx_data,
										global_variable 	 gv,
										
										PP_ref 				*PP_ref_db,
										SS_ref 				*SS_ref_db
){
	simplex_data *d  = (simplex_data *) splx_data;

	int     k = 0;

	d->swp = 1;
	while (d->swp == 1){					/** as long as a phase can be added to the guessed assemblage, go on */
		k 		  += 1;
		d->swp     = 0;
		
		swap_pure_phases(					z_b,
											splx_data,
											gv,
											PP_ref_db,
											SS_ref_db		);	

		swap_PGE_pseudocompounds(			z_b,
											splx_data,
											gv,
											PP_ref_db,
											SS_ref_db		);				
	}

	/* update gamma of SS */
	update_local_gamma(						d->A1,
											d->g0_A,
											d->gamma_ss,
											d->n_Ox			);

	/* update global variable gamma */
	update_global_gamma(					z_b,
											splx_data		);	

	if (gv.verbose == 1){
		printf("    (# iterations %d)",k);	
	}

}
