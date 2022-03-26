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
							SS_ref SS_ref_db		){
	
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
