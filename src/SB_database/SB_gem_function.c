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
/**
  Function to calculate chemical potential of endmembers/pure phases for Stixrude database
*/

#include <math.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <complex.h> 

#include "../MAGEMin.h"
#include "../toolkit.h"
#include "SB_gem_function.h"

#define eps 1e-10


double plg(double t) {
	double p4, dinc;

    double P0 	= exp(-t);
    double p1 	= 1.0;
    double p2 	= t * t;
    double p3 	= 2.0 * t;

    double plg 	= -2.1646464674222763831;

    int i 	= 1;
    while (i < 100000) {

        p4 		= (double)i;
        p1 		= p1 * P0;
        dinc 	= (p2 + (p3 + 2.0 / p4) / p4) * p1 / p4 / p4;
        plg 	= plg + dinc;

        if (fabs(dinc / (1.0 + fabs(plg))) < eps) {
            return plg;
        }

        i += 1;
    }

    return plg;
}


PP_ref SB_G_EM_function(	int 		 EM_dataset, 
							int 		 len_ox,
							int         *id,
							double 		*bulk_rock, 
							double 		*apo, 
							double 		 Pkbar, 
							double 		 T, 
							char 		*name, 
							char		*state			
){

	/* Get thermodynamic data */
	EM_db_sb EM_return;
	int i, p_id = find_EM_id(name);
	EM_return   = Access_SB_EM_DB(p_id, EM_dataset);
	
	/* Get composition (in molar amount) */
	double composition[len_ox];
	for (i = 0; i < len_ox; i ++){
		composition[i] = EM_return.Comp[id[i]];
	}

	double kbar2bar = 1e3;
	double T0 		= 298.15;
	double P0 		= 0.001;
	double R  		= 8.31446261815324;
 	double P       = Pkbar * kbar2bar;
	/* declare the variables */
	double nr9, nr9T0, c1, c2, c3, aii, aiikk2, aii2;
	double r23, r59, t1, t2, nr9t, tht, thT0;
	double dfth, dfth0, root, V, V23;
	double f,df,d2f,dfc,d2fc,z,a2f,da,dtht,d2tht,dthT0,d2thT0,fpoly,fpoly0,etht,letht,d2fth,ethT0,lethT0,d2fth0,f1,df1,dv;
	double a,gbase;

	int    max_ite, itic, ibad, bad;

	max_ite 		= 100;

	double F0		= EM_return.input_1[0];
	double n		= EM_return.input_1[1];
	double V0 		= -EM_return.input_1[2];
	double K0 		= EM_return.input_1[3];
	double Kp 		= EM_return.input_1[4];
	double z00		= EM_return.input_1[5];
	double gamma0 	= EM_return.input_1[6];
	double q0 		= EM_return.input_1[7];
	double etaS0	= EM_return.input_1[8];
	double cme		= EM_return.input_1[9];
    
    nr9 	= -9.0 * n * R;
    nr9T0 	= nr9 * T0;
    c1 		= -9.0 * (-V0) * K0;
    c2 		= Kp / 2.0 - 2.0;
    c3 		= 3.0 * c1 * c2;
    aii 	= 6.0 * gamma0;
    aiikk2 	= 0.5 * aii * (-2.0 + 6.0 * gamma0 - 3.0 * q0);
    aii2 	= 3.0 * gamma0;

    r23 = 2.0 / 3.0;
    r59 = 5.0 / 9.0;

    t1 		= z00 / T;
    t2 		= T / T0;
    nr9t 	= nr9 * T;

	tht 	= t1;
	thT0 	= tht * t2;

	dfth 	= nr9t * gamma0 / 	V0 * (3.0 * plg(tht) 	/ pow(tht,3.0) 	- log(1.0 - exp(-tht)));
	dfth0 	= nr9T0 * gamma0 / 	V0 * (3.0 * plg(thT0) 	/ pow(thT0,3.0) - log(1.0 - exp(-thT0)));

	root 	= K0 * ((2.0 + 2.0 * Kp) * (P + dfth - dfth0) + K0);

	V 		= 0.0;
	if (root > 0.0){
		V = ((2.0 + Kp) - sqrt(root) / K0) * V0 / (1.0 + Kp);
		if (V < V0 / 10.0 || V > V0 * 10.0){
			V = V0;
		}
	}
	else{
		V = V0;
	}

	ibad 	= 4;
	bad 	= 1;

	itic 	= 0;
	while (itic < max_ite){

		itic   += 1;

		V23 	= pow((V0 / V),r23);
		f 		= 0.5 * V23 - 0.5;
		df 		= -V23 / V / 3.0;
		d2f 	= r59 * V23 / pow(V,2.0);
		dfc 	= (c3 * f + c1) * f * df;
		d2fc 	= (2.0 * c3 * f + c1) * pow(df,2.0) + (c3 * f + c1) * f * d2f;

		z 		= 1.0 + (aii + aiikk2 * f) * f;

		if (z < 0.0 || V / V0 > 100.0 || V / V0 < 1e-2){
			printf(" ERROR z or V/V0\n");
		}

		root 	= sqrt(z);
		tht 	= t1 * root;
		thT0 	= tht * T / T0;
		a2f 	= aii2 + aiikk2 * f;
		da 		= a2f / root;
		dtht 	= t1 * da * df;
		d2tht 	= t1 * ((aiikk2 / root - pow(a2f,2.0) / pow(z,1.5)) * pow(df,2.0) + da * d2f);

		dthT0 	= dtht * t2;
		d2thT0 	= d2tht * t2;

		fpoly 	= 3.0 * plg(tht) / pow(tht,3.0);
		fpoly0 	= 3.0 * plg(thT0) / pow(thT0,3.0);
		etht 	= exp(-tht);

		if (1.0 - etht < 0.0){
			printf("ERROR 1-etht\n");
		}

		letht 	= log(1.0 - etht);

		dfth 	= (letht - fpoly) * nr9t * dtht / tht;
		d2fth 	= ((4.0 * pow(dtht,2.0) / tht - d2tht) * (fpoly - letht) + pow(dtht,2.0) * etht / (1.0 - etht)) * nr9t / tht;
		ethT0 	= exp(-thT0);

		if (1.0 - ethT0 < 0.0){
			printf("ERROR 1-thT0\n");
		}

		lethT0 	= log(1.0 - ethT0);

		dfth0 	= (lethT0 - fpoly0) * nr9T0 * dthT0 / thT0;
		d2fth0 	= ((4.0 * pow(dthT0,2.0) / thT0 - d2thT0) * (fpoly0 - lethT0) + pow(dthT0,2.0) * ethT0 / (1.0 - ethT0)) * nr9T0 / thT0;

		f1 		= -dfc - dfth + dfth0 - P;
		df1 	= -d2fc - d2fth + d2fth0;
		dv 		= f1 / df1;

		if (V - dv < 0.0){
			dv = V / 2.0;
		}

		V -= dv;
		
		if (itic >= max_ite || fabs(f1) > 1e40){

			if (fabs(f1 / P) < 0.0){
				ibad = 5;
				printf("ERROR abs(f1/p)\n");
			}
		}
		else if( fabs(dv / (1.0 + V)) < eps ){
				bad = 0;
				break;
		}
		
	}

	if (bad == 1){
		printf("ERROR bad\n");
	}

	/* get helmoltz energy:*/
	f 		= 0.5 * pow((V0 / V),r23) - 0.5;
	z 		= 1.0 + (aii + aiikk2 * f) * f;
	root 	= sqrt(z);

	tht 	= t1 * root;
	thT0 	= tht * t2;

	/* helmholtz energy */
	a 		= F0 + c1 * pow(f,2.0) * (0.5 + c2 * f) + nr9 * (T / pow(tht,3.0) * plg(tht) - T0 / pow(thT0,3.0) * plg(thT0));
	gbase 	= a + P * V - T * cme;



	/* fill structure to send back to main */
	PP_ref PP_ref_db;
	
	/* Calculate normalizing factor using bulk-rock composition */
	double factor  = 0.0;
	
	/* Calculate the number of atoms in the bulk-rock composition */
	double fbc     = 0.0;
	for (i = 0; i < len_ox; i++){
		fbc += bulk_rock[i]*apo[i];
	}
	
	/* Calculate the number of atom in the solution */
	double ape = 0.0;
	for (i = 0; i < len_ox; i++){
		ape += composition[i]*apo[i];
	}
	
	/* Calculate normalizing factor */
	factor = fbc/ape;

	strcpy(PP_ref_db.Name, name);
	for (i = 0; i < len_ox; i++){
		PP_ref_db.Comp[i] = composition[i];
	}
	PP_ref_db.gbase   =  gbase/kbar2bar;
	PP_ref_db.factor  =  factor;
	PP_ref_db.phase_shearModulus  =  (EM_return.input_2[0] + (P - P0)*(EM_return.input_2[1]) + (T - T0)*(EM_return.input_2[2]))/kbar2bar;

    //   SHM(PHMAX+1) = (1D0 + 2D0*F)**(2.5D0)*                         &
    //        (SHM0*(1D0 - 5D0*F) + F*SHMP*3D0*K0)  &
    //      -  D4*VOLUM/V0R*((DFT0-DFT)/D2/VQ)
	
	// printf("gbase %4s %+10f\n",name,PP_ref_db.gbase);
	// printf("phase_shearModulus %4s %+10f\n",name,PP_ref_db.phase_shearModulus);
	// for (i = 0; i < len_ox; i++){
	// 	printf("%+10f",PP_ref_db.Comp[i]*PP_ref_db.factor); 
	// }
	// printf("\n");

	return (PP_ref_db);
}
