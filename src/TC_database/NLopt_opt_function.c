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
	Function to perform local minimization of the solution phases
	-> using NLopt external library
	-> penalty on SF inequality constraints seems good with a small value of 1e-12 (1e-10 was giving problems)                                         
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <complex.h> 

#include "nlopt.h"                  // requires specifying this in the makefile
#include "../MAGEMin.h"
#include "NLopt_opt_function.h"
#include "../all_solution_phases.h"
#include "../toolkit.h"


/**************************************************************************************/
/**************************************************************************************/
/**********************METABASITE DATABASE (Green et al., 2016)************************/
/**************************************************************************************/
/**************************************************************************************/

/**
    Inequality constraints for L
*/
void liq_mb_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( x[7] -(x[7] + 1.0)*(x[0] + x[1] + x[3] + x[4] + x[5]));
    result[1] = ( -x[0]*(x[7] + 1.0));
    result[2] = ( -x[1]*x[2]*(x[7] + 1.0));
    result[3] = ( -x[1]*(1.0 -x[2])*(x[7] + 1.0));
    result[4] = ( -x[3]*(x[7] + 1.0) + x[7]);
    result[5] = ( -x[4]*(x[7] + 1.0) + x[7]);
    result[6] = ( -x[7] -(x[7] + 1.0)*(- x[0] -x[1] -x[3] -x[4] -x[5]) - 1.0);
    result[7] = ( -x[7]);
    result[8] = ( -x[5]*(x[7] + 1.0));
    result[9] = ( -x[6]);
    result[10] = ( x[6] - 1.0);

    if (grad) {
        grad[0] = -x[7] - 1.0;
        grad[1] = -x[7] - 1.0;
        grad[2] = 0.0;
        grad[3] = -x[7] - 1.0;
        grad[4] = -x[7] - 1.0;
        grad[5] = -x[7] - 1.0;
        grad[6] = 0.0;
        grad[7] = -x[0] -x[1] -x[3] -x[4] -x[5] + 1.0;
        grad[8] = -x[7] - 1.0;
        grad[9] = 0.0;
        grad[10] = 0.0;
        grad[11] = 0.0;
        grad[12] = 0.0;
        grad[13] = 0.0;
        grad[14] = 0.0;
        grad[15] = -x[0];
        grad[16] = 0.0;
        grad[17] = -x[2]*(x[7] + 1.0);
        grad[18] = -x[1]*(x[7] + 1.0);
        grad[19] = 0.0;
        grad[20] = 0.0;
        grad[21] = 0.0;
        grad[22] = 0.0;
        grad[23] = -x[1]*x[2];
        grad[24] = 0.0;
        grad[25] = -(1.0 -x[2])*(x[7] + 1.0);
        grad[26] = x[1]*(x[7] + 1.0);
        grad[27] = 0.0;
        grad[28] = 0.0;
        grad[29] = 0.0;
        grad[30] = 0.0;
        grad[31] = -x[1]*(1.0 -x[2]);
        grad[32] = 0.0;
        grad[33] = 0.0;
        grad[34] = 0.0;
        grad[35] = -x[7] - 1.0;
        grad[36] = 0.0;
        grad[37] = 0.0;
        grad[38] = 0.0;
        grad[39] = 1.0 -x[3];
        grad[40] = 0.0;
        grad[41] = 0.0;
        grad[42] = 0.0;
        grad[43] = 0.0;
        grad[44] = -x[7] - 1.0;
        grad[45] = 0.0;
        grad[46] = 0.0;
        grad[47] = 1.0 -x[4];
        grad[48] = x[7] + 1.0;
        grad[49] = x[7] + 1.0;
        grad[50] = 0.0;
        grad[51] = x[7] + 1.0;
        grad[52] = x[7] + 1.0;
        grad[53] = x[7] + 1.0;
        grad[54] = 0.0;
        grad[55] = x[0] + x[1] + x[3] + x[4] + x[5] - 1.0;
        grad[56] = 0.0;
        grad[57] = 0.0;
        grad[58] = 0.0;
        grad[59] = 0.0;
        grad[60] = 0.0;
        grad[61] = 0.0;
        grad[62] = 0.0;
        grad[63] = -1.00;
        grad[64] = 0.0;
        grad[65] = 0.0;
        grad[66] = 0.0;
        grad[67] = 0.0;
        grad[68] = 0.0;
        grad[69] = -x[7] - 1.0;
        grad[70] = 0.0;
        grad[71] = -x[5];
        grad[72] = 0.0;
        grad[73] = 0.0;
        grad[74] = 0.0;
        grad[75] = 0.0;
        grad[76] = 0.0;
        grad[77] = 0.0;
        grad[78] = -1.00;
        grad[79] = 0.0;
        grad[80] = 0.0;
        grad[81] = 0.0;
        grad[82] = 0.0;
        grad[83] = 0.0;
        grad[84] = 0.0;
        grad[85] = 0.0;
        grad[86] = 1.00;
        grad[87] = 0.0;
    }

    return;
};

/**
    Inequality constraints for hb
*/
void hb_mb_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( x[3] - 1.0);
    result[1] = ( x[3]*x[4] -x[3]);
    result[2] = ( -x[3]*x[4]);
    result[3] = ( x[0] -x[8] - 1.0);
    result[4] = ( -x[0] + x[8]);
    result[5] = ( -x[0]*x[1] -x[0]*x[6] -x[0]*x[7] + x[0] + x[1]*x[9] + x[1] + x[6]*x[9] + x[6] + x[7]*x[9] + x[7] -x[9] - 1.0);
    result[6] = ( x[0]*x[1] + x[0]*x[6] + x[0]*x[7] -x[0] -x[1]*x[9] -x[6]*x[9] -x[7]*x[9] + x[9]);
    result[7] = ( -x[1]);
    result[8] = ( -x[6]);
    result[9] = ( -x[7]);
    result[10] = ( -x[5]);
    result[11] = ( -x[0]*x[2] -x[0]*x[5] + x[0] -x[1]*x[9] + x[2] + x[5] -x[6]*x[9] -x[7]*x[9] + 1.5*x[8] + x[9] - 1.0);
    result[12] = ( x[0]*x[2] + x[0]*x[5] -x[0] + x[1]*x[9] + x[6]*x[9] + x[7]*x[9] - 1.5*x[8] -x[9]);
    result[13] = ( -x[2]);
    result[14] = ( 0.5*x[1] - 0.5*x[2] + 0.25*x[3] + 0.5*x[6] + 0.5*x[7] - 1.0);
    result[15] = ( -0.5*x[1] + 0.5*x[2] - 0.25*x[3] - 0.5*x[6] - 0.5*x[7]);
    result[16] = ( x[7] - 1.0);
    result[17] = ( -x[7]);

    if (grad) {
        grad[0] = 0.0;
        grad[1] = 0.0;
        grad[2] = 0.0;
        grad[3] = 1.00;
        grad[4] = 0.0;
        grad[5] = 0.0;
        grad[6] = 0.0;
        grad[7] = 0.0;
        grad[8] = 0.0;
        grad[9] = 0.0;
        grad[10] = 0.0;
        grad[11] = 0.0;
        grad[12] = 0.0;
        grad[13] = x[4] - 1.0;
        grad[14] = x[3];
        grad[15] = 0.0;
        grad[16] = 0.0;
        grad[17] = 0.0;
        grad[18] = 0.0;
        grad[19] = 0.0;
        grad[20] = 0.0;
        grad[21] = 0.0;
        grad[22] = 0.0;
        grad[23] = -x[4];
        grad[24] = -x[3];
        grad[25] = 0.0;
        grad[26] = 0.0;
        grad[27] = 0.0;
        grad[28] = 0.0;
        grad[29] = 0.0;
        grad[30] = 1.00;
        grad[31] = 0.0;
        grad[32] = 0.0;
        grad[33] = 0.0;
        grad[34] = 0.0;
        grad[35] = 0.0;
        grad[36] = 0.0;
        grad[37] = 0.0;
        grad[38] = -1.00;
        grad[39] = 0.0;
        grad[40] = -1.00;
        grad[41] = 0.0;
        grad[42] = 0.0;
        grad[43] = 0.0;
        grad[44] = 0.0;
        grad[45] = 0.0;
        grad[46] = 0.0;
        grad[47] = 0.0;
        grad[48] = 1.00;
        grad[49] = 0.0;
        grad[50] = -x[1] -x[6] -x[7] + 1.0;
        grad[51] = -x[0] + x[9] + 1.0;
        grad[52] = 0.0;
        grad[53] = 0.0;
        grad[54] = 0.0;
        grad[55] = 0.0;
        grad[56] = -x[0] + x[9] + 1.0;
        grad[57] = -x[0] + x[9] + 1.0;
        grad[58] = 0.0;
        grad[59] = x[1] + x[6] + x[7] - 1.0;
        grad[60] = x[1] + x[6] + x[7] - 1.0;
        grad[61] = x[0] -x[9];
        grad[62] = 0.0;
        grad[63] = 0.0;
        grad[64] = 0.0;
        grad[65] = 0.0;
        grad[66] = x[0] -x[9];
        grad[67] = x[0] -x[9];
        grad[68] = 0.0;
        grad[69] = -x[1] -x[6] -x[7] + 1.0;
        grad[70] = 0.0;
        grad[71] = -1.00;
        grad[72] = 0.0;
        grad[73] = 0.0;
        grad[74] = 0.0;
        grad[75] = 0.0;
        grad[76] = 0.0;
        grad[77] = 0.0;
        grad[78] = 0.0;
        grad[79] = 0.0;
        grad[80] = 0.0;
        grad[81] = 0.0;
        grad[82] = 0.0;
        grad[83] = 0.0;
        grad[84] = 0.0;
        grad[85] = 0.0;
        grad[86] = -1.00;
        grad[87] = 0.0;
        grad[88] = 0.0;
        grad[89] = 0.0;
        grad[90] = 0.0;
        grad[91] = 0.0;
        grad[92] = 0.0;
        grad[93] = 0.0;
        grad[94] = 0.0;
        grad[95] = 0.0;
        grad[96] = 0.0;
        grad[97] = -1.00;
        grad[98] = 0.0;
        grad[99] = 0.0;
        grad[100] = 0.0;
        grad[101] = 0.0;
        grad[102] = 0.0;
        grad[103] = 0.0;
        grad[104] = 0.0;
        grad[105] = -1.00;
        grad[106] = 0.0;
        grad[107] = 0.0;
        grad[108] = 0.0;
        grad[109] = 0.0;
        grad[110] = -x[2] -x[5] + 1.0;
        grad[111] = -x[9];
        grad[112] = 1.0 -x[0];
        grad[113] = 0.0;
        grad[114] = 0.0;
        grad[115] = 1.0 -x[0];
        grad[116] = -x[9];
        grad[117] = -x[9];
        grad[118] = 1.50;
        grad[119] = -x[1] -x[6] -x[7] + 1.0;
        grad[120] = x[2] + x[5] - 1.0;
        grad[121] = x[9];
        grad[122] = x[0];
        grad[123] = 0.0;
        grad[124] = 0.0;
        grad[125] = x[0];
        grad[126] = x[9];
        grad[127] = x[9];
        grad[128] = -1.50;
        grad[129] = x[1] + x[6] + x[7] - 1.0;
        grad[130] = 0.0;
        grad[131] = 0.0;
        grad[132] = -1.00;
        grad[133] = 0.0;
        grad[134] = 0.0;
        grad[135] = 0.0;
        grad[136] = 0.0;
        grad[137] = 0.0;
        grad[138] = 0.0;
        grad[139] = 0.0;
        grad[140] = 0.0;
        grad[141] = 0.500;
        grad[142] = -0.500;
        grad[143] = 0.250;
        grad[144] = 0.0;
        grad[145] = 0.0;
        grad[146] = 0.500;
        grad[147] = 0.500;
        grad[148] = 0.0;
        grad[149] = 0.0;
        grad[150] = 0.0;
        grad[151] = -0.500;
        grad[152] = 0.500;
        grad[153] = -0.250;
        grad[154] = 0.0;
        grad[155] = 0.0;
        grad[156] = -0.500;
        grad[157] = -0.500;
        grad[158] = 0.0;
        grad[159] = 0.0;
        grad[160] = 0.0;
        grad[161] = 0.0;
        grad[162] = 0.0;
        grad[163] = 0.0;
        grad[164] = 0.0;
        grad[165] = 0.0;
        grad[166] = 0.0;
        grad[167] = 1.00;
        grad[168] = 0.0;
        grad[169] = 0.0;
        grad[170] = 0.0;
        grad[171] = 0.0;
        grad[172] = 0.0;
        grad[173] = 0.0;
        grad[174] = 0.0;
        grad[175] = 0.0;
        grad[176] = 0.0;
        grad[177] = -1.00;
        grad[178] = 0.0;
        grad[179] = 0.0;
    }

    return;
};

/**
    Inequality constraints for aug
*/
void aug_mb_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( -x[0]*x[1] -x[0]*x[4] + x[0] + x[1] + 0.5*x[3]*x[5] + 0.5*x[4]*x[5] + x[4] - 0.5*x[5] - 1.0);
    result[1] = ( x[0]*x[1] + x[0]*x[4] -x[0] - 0.5*x[3]*x[5] - 0.5*x[4]*x[5] + 0.5*x[5]);
    result[2] = ( -x[1] + x[2] -x[4]);
    result[3] = ( -x[2]);
    result[4] = ( -x[0]*x[3] -x[0]*x[4] + x[0] - 0.5*x[3]*x[5] + x[3] - 0.5*x[4]*x[5] + x[4] + 0.5*x[5] - 1.0);
    result[5] = ( x[0]*x[3] + x[0]*x[4] -x[0] + 0.5*x[3]*x[5] + 0.5*x[4]*x[5] - 0.5*x[5]);
    result[6] = ( -x[3]);
    result[7] = ( -x[4]);
    result[8] = ( 0.5*x[1] - 0.5*x[6] - 1.0);
    result[9] = ( -0.5*x[1] + 0.5*x[6]);
    result[10] = ( 0.5*x[1] + 0.5*x[6] - 1.0);
    result[11] = ( -0.5*x[1] - 0.5*x[6]);

    if (grad) {
        grad[0] = -x[1] -x[4] + 1.0;
        grad[1] = 1.0 -x[0];
        grad[2] = 0.0;
        grad[3] = 0.5*x[5];
        grad[4] = -x[0] + 0.5*x[5] + 1.0;
        grad[5] = 0.5*x[3] + 0.5*x[4] - 0.5;
        grad[6] = 0.0;
        grad[7] = x[1] + x[4] - 1.0;
        grad[8] = x[0];
        grad[9] = 0.0;
        grad[10] = -0.5*x[5];
        grad[11] = x[0] - 0.5*x[5];
        grad[12] = -0.5*x[3] - 0.5*x[4] + 0.5;
        grad[13] = 0.0;
        grad[14] = 0.0;
        grad[15] = -1.00;
        grad[16] = 1.00;
        grad[17] = 0.0;
        grad[18] = -1.00;
        grad[19] = 0.0;
        grad[20] = 0.0;
        grad[21] = 0.0;
        grad[22] = 0.0;
        grad[23] = -1.00;
        grad[24] = 0.0;
        grad[25] = 0.0;
        grad[26] = 0.0;
        grad[27] = 0.0;
        grad[28] = -x[3] -x[4] + 1.0;
        grad[29] = 0.0;
        grad[30] = 0.0;
        grad[31] = -x[0] - 0.5*x[5] + 1.0;
        grad[32] = -x[0] - 0.5*x[5] + 1.0;
        grad[33] = -0.5*x[3] - 0.5*x[4] + 0.5;
        grad[34] = 0.0;
        grad[35] = x[3] + x[4] - 1.0;
        grad[36] = 0.0;
        grad[37] = 0.0;
        grad[38] = x[0] + 0.5*x[5];
        grad[39] = x[0] + 0.5*x[5];
        grad[40] = 0.5*x[3] + 0.5*x[4] - 0.5;
        grad[41] = 0.0;
        grad[42] = 0.0;
        grad[43] = 0.0;
        grad[44] = 0.0;
        grad[45] = -1.00;
        grad[46] = 0.0;
        grad[47] = 0.0;
        grad[48] = 0.0;
        grad[49] = 0.0;
        grad[50] = 0.0;
        grad[51] = 0.0;
        grad[52] = 0.0;
        grad[53] = -1.00;
        grad[54] = 0.0;
        grad[55] = 0.0;
        grad[56] = 0.0;
        grad[57] = 0.500;
        grad[58] = 0.0;
        grad[59] = 0.0;
        grad[60] = 0.0;
        grad[61] = 0.0;
        grad[62] = -0.500;
        grad[63] = 0.0;
        grad[64] = -0.500;
        grad[65] = 0.0;
        grad[66] = 0.0;
        grad[67] = 0.0;
        grad[68] = 0.0;
        grad[69] = 0.500;
        grad[70] = 0.0;
        grad[71] = 0.500;
        grad[72] = 0.0;
        grad[73] = 0.0;
        grad[74] = 0.0;
        grad[75] = 0.0;
        grad[76] = 0.500;
        grad[77] = 0.0;
        grad[78] = -0.500;
        grad[79] = 0.0;
        grad[80] = 0.0;
        grad[81] = 0.0;
        grad[82] = 0.0;
        grad[83] = -0.500;
    }

    return;
};

/**
    Inequality constraints for dio
*/
void dio_mb_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( -x[0]*x[1] + x[0]*x[3] + x[0] + x[1]*x[5] + x[1] + x[3]*x[5] -x[3] -x[5] - 1.0);
    result[1] = ( x[0]*x[1] -x[0]*x[3] -x[0] -x[1]*x[5] -x[3]*x[5] + x[5]);
    result[2] = ( -x[1]*x[2] + x[4]);
    result[3] = ( x[1]*x[2] -x[1] + x[3] -x[4]);
    result[4] = ( -x[0]*x[1] -x[0]*x[3] + x[0] -x[1]*x[5] + x[1] -x[3]*x[5] + x[3] + x[5] - 1.0);
    result[5] = ( x[0]*x[1] + x[0]*x[3] -x[0] + x[1]*x[5] + x[3]*x[5] -x[5]);
    result[6] = ( -x[1]*x[2] -x[4]);
    result[7] = ( x[1]*x[2] -x[1] -x[3] + x[4]);
    result[8] = ( -x[1] + x[3]);
    result[9] = ( x[1] -x[3] - 1.0);
    result[10] = ( -x[1] -x[3]);
    result[11] = ( x[1] + x[3] - 1.0);

    if (grad) {
        grad[0] = -x[1] + x[3] + 1.0;
        grad[1] = -x[0] + x[5] + 1.0;
        grad[2] = 0.0;
        grad[3] = x[0] + x[5] - 1.0;
        grad[4] = 0.0;
        grad[5] = x[1] + x[3] - 1.0;
        grad[6] = x[1] -x[3] - 1.0;
        grad[7] = x[0] -x[5];
        grad[8] = 0.0;
        grad[9] = -x[0] -x[5];
        grad[10] = 0.0;
        grad[11] = -x[1] -x[3] + 1.0;
        grad[12] = 0.0;
        grad[13] = -x[2];
        grad[14] = -x[1];
        grad[15] = 0.0;
        grad[16] = 1.00;
        grad[17] = 0.0;
        grad[18] = 0.0;
        grad[19] = x[2] - 1.0;
        grad[20] = x[1];
        grad[21] = 1.00;
        grad[22] = -1.00;
        grad[23] = 0.0;
        grad[24] = -x[1] -x[3] + 1.0;
        grad[25] = -x[0] -x[5] + 1.0;
        grad[26] = 0.0;
        grad[27] = -x[0] -x[5] + 1.0;
        grad[28] = 0.0;
        grad[29] = -x[1] -x[3] + 1.0;
        grad[30] = x[1] + x[3] - 1.0;
        grad[31] = x[0] + x[5];
        grad[32] = 0.0;
        grad[33] = x[0] + x[5];
        grad[34] = 0.0;
        grad[35] = x[1] + x[3] - 1.0;
        grad[36] = 0.0;
        grad[37] = -x[2];
        grad[38] = -x[1];
        grad[39] = 0.0;
        grad[40] = -1.00;
        grad[41] = 0.0;
        grad[42] = 0.0;
        grad[43] = x[2] - 1.0;
        grad[44] = x[1];
        grad[45] = -1.00;
        grad[46] = 1.00;
        grad[47] = 0.0;
        grad[48] = 0.0;
        grad[49] = -1.00;
        grad[50] = 0.0;
        grad[51] = 1.00;
        grad[52] = 0.0;
        grad[53] = 0.0;
        grad[54] = 0.0;
        grad[55] = 1.00;
        grad[56] = 0.0;
        grad[57] = -1.00;
        grad[58] = 0.0;
        grad[59] = 0.0;
        grad[60] = 0.0;
        grad[61] = -1.00;
        grad[62] = 0.0;
        grad[63] = -1.00;
        grad[64] = 0.0;
        grad[65] = 0.0;
        grad[66] = 0.0;
        grad[67] = 1.00;
        grad[68] = 0.0;
        grad[69] = 1.00;
        grad[70] = 0.0;
        grad[71] = 0.0;
    }

    return;
};

/**
    Inequality constraints for opx
*/
void opx_mb_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( -x[0]*x[1] -x[0]*x[2] + x[0] + x[1] + x[2] + 0.5*x[3]*x[4] - 0.5*x[4] - 1.0);
    result[1] = ( x[0]*x[1] + x[0]*x[2] -x[0] - 0.5*x[3]*x[4] + 0.5*x[4]);
    result[2] = ( -x[2]);
    result[3] = ( -x[1]);
    result[4] = ( -x[0]*x[3] + x[0] - 0.5*x[3]*x[4] + x[3] + 0.5*x[4] - 1.0);
    result[5] = ( x[0]*x[3] -x[0] + 0.5*x[3]*x[4] - 0.5*x[4]);
    result[6] = ( -x[3]);
    result[7] = ( -0.5*x[1] - 0.5*x[2]);
    result[8] = ( 0.5*x[1] + 0.5*x[2] - 1.0);

    if (grad) {
        grad[0] = -x[1] -x[2] + 1.0;
        grad[1] = 1.0 -x[0];
        grad[2] = 1.0 -x[0];
        grad[3] = 0.5*x[4];
        grad[4] = 0.5*x[3] - 0.5;
        grad[5] = x[1] + x[2] - 1.0;
        grad[6] = x[0];
        grad[7] = x[0];
        grad[8] = -0.5*x[4];
        grad[9] = 0.5 - 0.5*x[3];
        grad[10] = 0.0;
        grad[11] = 0.0;
        grad[12] = -1.00;
        grad[13] = 0.0;
        grad[14] = 0.0;
        grad[15] = 0.0;
        grad[16] = -1.00;
        grad[17] = 0.0;
        grad[18] = 0.0;
        grad[19] = 0.0;
        grad[20] = 1.0 -x[3];
        grad[21] = 0.0;
        grad[22] = 0.0;
        grad[23] = -x[0] - 0.5*x[4] + 1.0;
        grad[24] = 0.5 - 0.5*x[3];
        grad[25] = x[3] - 1.0;
        grad[26] = 0.0;
        grad[27] = 0.0;
        grad[28] = x[0] + 0.5*x[4];
        grad[29] = 0.5*x[3] - 0.5;
        grad[30] = 0.0;
        grad[31] = 0.0;
        grad[32] = 0.0;
        grad[33] = -1.00;
        grad[34] = 0.0;
        grad[35] = 0.0;
        grad[36] = -0.500;
        grad[37] = -0.500;
        grad[38] = 0.0;
        grad[39] = 0.0;
        grad[40] = 0.0;
        grad[41] = 0.500;
        grad[42] = 0.500;
        grad[43] = 0.0;
        grad[44] = 0.0;
    }

    return;
};

/**
    Inequality constraints for g
*/
void g_mb_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( -x[0]*x[1] + x[0] + x[1] - 1.0);
    result[1] = ( x[0]*x[1] -x[0]);
    result[2] = ( -x[1]);
    result[3] = ( x[2] - 1.0);
    result[4] = ( -x[2]);

    if (grad) {
        grad[0] = 1.0 -x[1];
        grad[1] = 1.0 -x[0];
        grad[2] = 0.0;
        grad[3] = x[1] - 1.0;
        grad[4] = x[0];
        grad[5] = 0.0;
        grad[6] = 0.0;
        grad[7] = -1.00;
        grad[8] = 0.0;
        grad[9] = 0.0;
        grad[10] = 0.0;
        grad[11] = 1.00;
        grad[12] = 0.0;
        grad[13] = 0.0;
        grad[14] = -1.00;
    }

    return;
};

/**
    Inequality constraints for ol
*/
void ol_mb_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( x[0] - 1.0);
    result[1] = ( -x[0]);

    if (grad) {
        grad[0] = 1.00;
        grad[1] = -1.00;
    }

    return;
};

/**
    Inequality constraints for fsp
*/
void fsp_mb_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( x[0] + x[1] - 1.0);
    result[1] = ( -x[0]);
    result[2] = ( -x[1]);
    result[3] = ( -0.25*x[0] - 0.25);
    result[4] = ( 0.25*x[0] - 0.75);

    if (grad) {
        grad[0] = 1.00;
        grad[1] = 1.00;
        grad[2] = -1.00;
        grad[3] = 0.0;
        grad[4] = 0.0;
        grad[5] = -1.00;
        grad[6] = -0.250;
        grad[7] = 0.0;
        grad[8] = 0.250;
        grad[9] = 0.0;
    }

    return;
};

/**
    Inequality constraints for abc
*/
void abc_mb_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( x[0] - 1.0);
    result[1] = ( -x[0]);

    if (grad) {
        grad[0] = 1.00;
        grad[1] = -1.00;
    }

    return;
};

/**
    Inequality constraints for k4tr
*/
void k4tr_mb_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( -x[0]);
    result[1] = ( -x[1]);
    result[2] = ( x[0] + x[1] - 1.0);
    result[3] = ( -0.25*x[1] - 0.25);
    result[4] = ( 0.25*x[1] - 0.75);

    if (grad) {
        grad[0] = -1.00;
        grad[1] = 0.0;
        grad[2] = 0.0;
        grad[3] = -1.00;
        grad[4] = 1.00;
        grad[5] = 1.00;
        grad[6] = 0.0;
        grad[7] = -0.250;
        grad[8] = 0.0;
        grad[9] = 0.250;
    }

    return;
};

/**
    Inequality constraints for spl
*/
void spl_mb_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = (-x[1]);
    result[1] = (x[1] - 1.0);
    result[2] = (x[0] - 1.0);
    result[3] = (-x[0]);

    if (grad) {
        grad[0] = 0.0;
        grad[1] = -1.0;
        grad[2] = 0.0;
        grad[3] = 1.0;
        grad[4] = 1.0;
        grad[5] = 0.0;
        grad[6] = -1.0;
        grad[7] = 0.0;
    }

    return;
};
/**
    Inequality constraints for sp
*/
void sp_mb_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( -x[1]);
    result[1] = ( x[1] + x[2] - 1.0);
    result[2] = ( -x[2]);
    result[3] = ( x[0] - 1.0);
    result[4] = ( -x[0]);

    if (grad) {
        grad[0] = 0.0;
        grad[1] = -1.00;
        grad[2] = 0.0;
        grad[3] = 0.0;
        grad[4] = 1.00;
        grad[5] = 1.00;
        grad[6] = 0.0;
        grad[7] = 0.0;
        grad[8] = -1.00;
        grad[9] = 1.00;
        grad[10] = 0.0;
        grad[11] = 0.0;
        grad[12] = -1.00;
        grad[13] = 0.0;
        grad[14] = 0.0;
    }

    return;
};

/**
    Inequality constraints for ilm
*/
void ilm_mb_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( -0.5*x[0] - 0.5*x[1]);
    result[1] = ( -0.5*x[0] + 0.5*x[1]);
    result[2] = ( x[0] - 1.0);
    result[3] = ( -0.5*x[0] + 0.5*x[1]);
    result[4] = ( -0.5*x[0] - 0.5*x[1]);
    result[5] = ( x[0] - 1.0);

    if (grad) {
        grad[0] = -0.500;
        grad[1] = -0.500;
        grad[2] = -0.500;
        grad[3] = 0.500;
        grad[4] = 1.00;
        grad[5] = 0.0;
        grad[6] = -0.500;
        grad[7] = 0.500;
        grad[8] = -0.500;
        grad[9] = -0.500;
        grad[10] = 1.00;
        grad[11] = 0.0;
    }

    return;
};

/**
    Inequality constraints for ilmm
*/
void ilmm_mb_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( -0.5*x[0] + 0.5*x[1] - 0.5*x[2]);
    result[1] = ( -0.5*x[0] + 0.5*x[1] + 0.5*x[2]);
    result[2] = ( -x[1]);
    result[3] = ( x[0] - 1.0);
    result[4] = ( -0.5*x[0] + 0.5*x[1] + 0.5*x[2]);
    result[5] = ( -0.5*x[0] - 0.5*x[1] - 0.5*x[2]);
    result[6] = ( x[0] - 1.0);

    if (grad) {
        grad[0] = -0.500;
        grad[1] = 0.500;
        grad[2] = -0.500;
        grad[3] = -0.500;
        grad[4] = 0.500;
        grad[5] = 0.500;
        grad[6] = 0.0;
        grad[7] = -1.00;
        grad[8] = 0.0;
        grad[9] = 1.00;
        grad[10] = 0.0;
        grad[11] = 0.0;
        grad[12] = -0.500;
        grad[13] = 0.500;
        grad[14] = 0.500;
        grad[15] = -0.500;
        grad[16] = -0.500;
        grad[17] = -0.500;
        grad[18] = 1.00;
        grad[19] = 0.0;
        grad[20] = 0.0;
    }

    return;
};

/**
    Inequality constraints for ep
*/
void ep_mb_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( -x[0] + x[1]);
    result[1] = ( x[0] -x[1] - 1.0);
    result[2] = ( -x[0] -x[1]);
    result[3] = ( x[0] + x[1] - 1.0);

    if (grad) {
        grad[0] = -1.00;
        grad[1] = 1.00;
        grad[2] = 1.00;
        grad[3] = -1.00;
        grad[4] = -1.00;
        grad[5] = -1.00;
        grad[6] = 1.00;
        grad[7] = 1.00;
    }

    return;
};

/**
    Inequality constraints for bi
*/
void bi_mb_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( -x[0]*x[1] -x[0]*x[2] -x[0]*x[3] + x[0] + x[1] + x[2] + x[3] + 2.0/3.0*x[4] - 1.0);
    result[1] = ( x[0]*x[1] + x[0]*x[2] + x[0]*x[3] -x[0] - 2.0/3.0*x[4]);
    result[2] = ( -x[2]);
    result[3] = ( -x[3]);
    result[4] = ( -x[1]);
    result[5] = ( x[0] - 1.0/3.0*x[4] - 1.0);
    result[6] = ( -x[0] + 1.0/3.0*x[4]);
    result[7] = ( 0.5*x[1] + 0.5*x[2] - 0.5);
    result[8] = ( -0.5*x[1] - 0.5*x[2] - 0.5);
    result[9] = ( x[3] - 1.0);
    result[10] = ( -x[3]);

    if (grad) {
        grad[0] = -x[1] -x[2] -x[3] + 1.0;
        grad[1] = 1.0 -x[0];
        grad[2] = 1.0 -x[0];
        grad[3] = 1.0 -x[0];
        grad[4] = 2.0/3.0;
        grad[5] = x[1] + x[2] + x[3] - 1.0;
        grad[6] = x[0];
        grad[7] = x[0];
        grad[8] = x[0];
        grad[9] = -2.0/3.0;
        grad[10] = 0.0;
        grad[11] = 0.0;
        grad[12] = -1.00;
        grad[13] = 0.0;
        grad[14] = 0.0;
        grad[15] = 0.0;
        grad[16] = 0.0;
        grad[17] = 0.0;
        grad[18] = -1.00;
        grad[19] = 0.0;
        grad[20] = 0.0;
        grad[21] = -1.00;
        grad[22] = 0.0;
        grad[23] = 0.0;
        grad[24] = 0.0;
        grad[25] = 1.00;
        grad[26] = 0.0;
        grad[27] = 0.0;
        grad[28] = 0.0;
        grad[29] = -1.0/3.0;
        grad[30] = -1.00;
        grad[31] = 0.0;
        grad[32] = 0.0;
        grad[33] = 0.0;
        grad[34] = 1.0/3.0;
        grad[35] = 0.0;
        grad[36] = 0.500;
        grad[37] = 0.500;
        grad[38] = 0.0;
        grad[39] = 0.0;
        grad[40] = 0.0;
        grad[41] = -0.500;
        grad[42] = -0.500;
        grad[43] = 0.0;
        grad[44] = 0.0;
        grad[45] = 0.0;
        grad[46] = 0.0;
        grad[47] = 0.0;
        grad[48] = 1.00;
        grad[49] = 0.0;
        grad[50] = 0.0;
        grad[51] = 0.0;
        grad[52] = 0.0;
        grad[53] = -1.00;
        grad[54] = 0.0;
    }

    return;
};

/**
    Inequality constraints for mu
*/
void mu_mb_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( x[3] + x[4] - 1.0);
    result[1] = ( -x[3]);
    result[2] = ( -x[4]);
    result[3] = ( -x[0]*x[1] + x[0] + x[1] - 1.0);
    result[4] = ( x[0]*x[1] -x[0]);
    result[5] = ( -x[1]);
    result[6] = ( x[2] - 1.0);
    result[7] = ( -x[2]);
    result[8] = ( 0.5*x[1] + 0.5*x[4] - 1.0);
    result[9] = ( -0.5*x[1] - 0.5*x[4]);

    if (grad) {
        grad[0] = 0.0;
        grad[1] = 0.0;
        grad[2] = 0.0;
        grad[3] = 1.00;
        grad[4] = 1.00;
        grad[5] = 0.0;
        grad[6] = 0.0;
        grad[7] = 0.0;
        grad[8] = -1.00;
        grad[9] = 0.0;
        grad[10] = 0.0;
        grad[11] = 0.0;
        grad[12] = 0.0;
        grad[13] = 0.0;
        grad[14] = -1.00;
        grad[15] = 1.0 -x[1];
        grad[16] = 1.0 -x[0];
        grad[17] = 0.0;
        grad[18] = 0.0;
        grad[19] = 0.0;
        grad[20] = x[1] - 1.0;
        grad[21] = x[0];
        grad[22] = 0.0;
        grad[23] = 0.0;
        grad[24] = 0.0;
        grad[25] = 0.0;
        grad[26] = -1.00;
        grad[27] = 0.0;
        grad[28] = 0.0;
        grad[29] = 0.0;
        grad[30] = 0.0;
        grad[31] = 0.0;
        grad[32] = 1.00;
        grad[33] = 0.0;
        grad[34] = 0.0;
        grad[35] = 0.0;
        grad[36] = 0.0;
        grad[37] = -1.00;
        grad[38] = 0.0;
        grad[39] = 0.0;
        grad[40] = 0.0;
        grad[41] = 0.500;
        grad[42] = 0.0;
        grad[43] = 0.0;
        grad[44] = 0.500;
        grad[45] = 0.0;
        grad[46] = -0.500;
        grad[47] = 0.0;
        grad[48] = 0.0;
        grad[49] = -0.500;
    }

    return;
};

/**
    Inequality constraints for chl
*/
void chl_mb_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( -x[0]*x[1] + x[0]*x[3] + x[0] + x[1]*x[4] + x[1] -x[3]*x[4] -x[3] -x[4] - 1.0);
    result[1] = ( x[0]*x[1] -x[0]*x[3] -x[0] -x[1]*x[4] + x[3]*x[4] + x[4]);
    result[2] = ( -x[1] + x[3]);
    result[3] = ( x[0] - 0.25*x[1]*x[4] - 0.25*x[1]*x[5] - 0.25*x[2]*x[5] + 0.25*x[3]*x[4] - 0.25*x[3]*x[5] + 0.25*x[4] + 0.25*x[5] - 1.0);
    result[4] = ( -x[0] + 0.25*x[1]*x[4] + 0.25*x[1]*x[5] + 0.25*x[2]*x[5] - 0.25*x[3]*x[4] + 0.25*x[3]*x[5] - 0.25*x[4] - 0.25*x[5]);
    result[5] = ( -x[0]*x[1] -x[0]*x[2] -x[0]*x[3] + x[0] + x[1]*x[5] + x[1] + x[2]*x[5] + x[2] + x[3]*x[5] + x[3] -x[5] - 1.0);
    result[6] = ( x[0]*x[1] + x[0]*x[2] + x[0]*x[3] -x[0] -x[1]*x[5] -x[2]*x[5] -x[3]*x[5] + x[5]);
    result[7] = ( -x[2]);
    result[8] = ( -x[1] -x[3]);
    result[9] = ( x[1] + 0.5*x[2] - 1.0);
    result[10] = ( -x[1] - 0.5*x[2]);

    if (grad) {
        grad[0] = -x[1] + x[3] + 1.0;
        grad[1] = -x[0] + x[4] + 1.0;
        grad[2] = 0.0;
        grad[3] = x[0] -x[4] - 1.0;
        grad[4] = x[1] -x[3] - 1.0;
        grad[5] = 0.0;
        grad[6] = x[1] -x[3] - 1.0;
        grad[7] = x[0] -x[4];
        grad[8] = 0.0;
        grad[9] = -x[0] + x[4];
        grad[10] = -x[1] + x[3] + 1.0;
        grad[11] = 0.0;
        grad[12] = 0.0;
        grad[13] = -1.00;
        grad[14] = 0.0;
        grad[15] = 1.00;
        grad[16] = 0.0;
        grad[17] = 0.0;
        grad[18] = 1.00;
        grad[19] = -0.25*x[4] - 0.25*x[5];
        grad[20] = -0.25*x[5];
        grad[21] = 0.25*x[4] - 0.25*x[5];
        grad[22] = -0.25*x[1] + 0.25*x[3] + 0.25;
        grad[23] = -0.25*x[1] - 0.25*x[2] - 0.25*x[3] + 0.25;
        grad[24] = -1.00;
        grad[25] = 0.25*x[4] + 0.25*x[5];
        grad[26] = 0.25*x[5];
        grad[27] = -0.25*x[4] + 0.25*x[5];
        grad[28] = 0.25*x[1] - 0.25*x[3] - 0.25;
        grad[29] = 0.25*x[1] + 0.25*x[2] + 0.25*x[3] - 0.25;
        grad[30] = -x[1] -x[2] -x[3] + 1.0;
        grad[31] = -x[0] + x[5] + 1.0;
        grad[32] = -x[0] + x[5] + 1.0;
        grad[33] = -x[0] + x[5] + 1.0;
        grad[34] = 0.0;
        grad[35] = x[1] + x[2] + x[3] - 1.0;
        grad[36] = x[1] + x[2] + x[3] - 1.0;
        grad[37] = x[0] -x[5];
        grad[38] = x[0] -x[5];
        grad[39] = x[0] -x[5];
        grad[40] = 0.0;
        grad[41] = -x[1] -x[2] -x[3] + 1.0;
        grad[42] = 0.0;
        grad[43] = 0.0;
        grad[44] = -1.00;
        grad[45] = 0.0;
        grad[46] = 0.0;
        grad[47] = 0.0;
        grad[48] = 0.0;
        grad[49] = -1.00;
        grad[50] = 0.0;
        grad[51] = -1.00;
        grad[52] = 0.0;
        grad[53] = 0.0;
        grad[54] = 0.0;
        grad[55] = 1.00;
        grad[56] = 0.500;
        grad[57] = 0.0;
        grad[58] = 0.0;
        grad[59] = 0.0;
        grad[60] = 0.0;
        grad[61] = -1.00;
        grad[62] = -0.500;
        grad[63] = 0.0;
        grad[64] = 0.0;
        grad[65] = 0.0;
    }

    return;
};

SS_ref NLopt_opt_mb_liq_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_mb_liq, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, liq_mb_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    
    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};
SS_ref NLopt_opt_mb_hb_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_mb_hb, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, hb_mb_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    
    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};
SS_ref NLopt_opt_mb_aug_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_mb_aug, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, aug_mb_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    
    double minf;
   SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);
    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};
SS_ref NLopt_opt_mb_dio_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_mb_dio, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, dio_mb_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    
    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};
SS_ref NLopt_opt_mb_opx_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_mb_opx, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, opx_mb_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    
    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};
SS_ref NLopt_opt_mb_g_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_mb_g, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, g_mb_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    
    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};
SS_ref NLopt_opt_mb_ol_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_mb_ol, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, ol_mb_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    
    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};
SS_ref NLopt_opt_mb_fsp_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_mb_fsp, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, fsp_mb_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    
    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};
SS_ref NLopt_opt_mb_abc_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_mb_abc, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, abc_mb_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    
    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};
SS_ref NLopt_opt_mb_k4tr_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_mb_k4tr, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, k4tr_mb_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    
    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};
SS_ref NLopt_opt_mb_spl_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_mb_spl, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, spl_mb_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    
    double minf;
    if (gv.maxeval==1){  
       // we are only interested in evaluating the objective function  
       minf = obj_mb_spl(n, x, NULL, &SS_ref_db);
    }
    else{
      // do optimization
      SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);
    }
    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};
SS_ref NLopt_opt_mb_sp_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_mb_sp, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, sp_mb_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    
    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};
SS_ref NLopt_opt_mb_ilm_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_mb_ilm, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, ilm_mb_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    
    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};
SS_ref NLopt_opt_mb_ilmm_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_mb_ilmm, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, ilmm_mb_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    
    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};
SS_ref NLopt_opt_mb_ep_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_mb_ep, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, ep_mb_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    
    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};
SS_ref NLopt_opt_mb_bi_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_mb_bi, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, bi_mb_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    
    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};
SS_ref NLopt_opt_mb_mu_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_mb_mu, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, mu_mb_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    
    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};
SS_ref NLopt_opt_mb_chl_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_mb_chl, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, chl_mb_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    
    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};



/**
    Inequality constraints for liq_mp
*/
void liq_mp_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( x[6] - 1.0);
    result[1] = ( -x[0]);
    result[2] = ( -x[1]*x[2]);
    result[3] = ( -x[1]*(1.0 - x[2]));
    result[4] = ( -x[3]);
    result[5] = ( x[3] + x[1] + x[6] + x[4] + x[0] - 1.0);
    result[6] = ( -x[4]);
    result[7] = ( -x[5]);
    result[8] = ( x[5] - 1.0);
    result[9] = ( -x[6]);

    if (grad) {
        grad[0] = 0.0;
        grad[1] = 0.0;
        grad[2] = 0.0;
        grad[3] = 0.0;
        grad[4] = 0.0;
        grad[5] = 0.0;
        grad[6] = 1.0;
        grad[7] = -1.0;
        grad[8] = 0.0;
        grad[9] = 0.0;
        grad[10] = 0.0;
        grad[11] = 0.0;
        grad[12] = 0.0;
        grad[13] = 0.0;
        grad[14] = 0.0;
        grad[15] = -x[2];
        grad[16] = -x[1];
        grad[17] = 0.0;
        grad[18] = 0.0;
        grad[19] = 0.0;
        grad[20] = 0.0;
        grad[21] = 0.0;
        grad[22] = x[2] - 1.0;
        grad[23] = x[1];
        grad[24] = 0.0;
        grad[25] = 0.0;
        grad[26] = 0.0;
        grad[27] = 0.0;
        grad[28] = 0.0;
        grad[29] = 0.0;
        grad[30] = 0.0;
        grad[31] = -1.0;
        grad[32] = 0.0;
        grad[33] = 0.0;
        grad[34] = 0.0;
        grad[35] = 1.0;
        grad[36] = 1.0;
        grad[37] = 0.0;
        grad[38] = 1.0;
        grad[39] = 1.0;
        grad[40] = 0.0;
        grad[41] = 1.0;
        grad[42] = 0.0;
        grad[43] = 0.0;
        grad[44] = 0.0;
        grad[45] = 0.0;
        grad[46] = -1.0;
        grad[47] = 0.0;
        grad[48] = 0.0;
        grad[49] = 0.0;
        grad[50] = 0.0;
        grad[51] = 0.0;
        grad[52] = 0.0;
        grad[53] = 0.0;
        grad[54] = -1.0;
        grad[55] = 0.0;
        grad[56] = 0.0;
        grad[57] = 0.0;
        grad[58] = 0.0;
        grad[59] = 0.0;
        grad[60] = 0.0;
        grad[61] = 1.0;
        grad[62] = 0.0;
        grad[63] = 0.0;
        grad[64] = 0.0;
        grad[65] = 0.0;
        grad[66] = 0.0;
        grad[67] = 0.0;
        grad[68] = 0.0;
        grad[69] = -1.0;
    }

    return;
};

/**
    Inequality constraints for st_mp
*/
void st_mp_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( -x[1]*x[0] + x[1] + x[0] - 1.0);
    result[1] = ( x[1]*x[0] - x[0]);
    result[2] = ( -x[1]);
    result[3] = ( x[2] + 4.0/3.0*x[3] - 1.0);
    result[4] = ( -x[2]);
    result[5] = ( -x[3]);
    result[6] = ( -1./3.*x[3]);

    if (grad) {
        grad[0] = 1.0 - x[1];
        grad[1] = 1.0 - x[0];
        grad[2] = 0.0;
        grad[3] = 0.0;
        grad[4] = x[1] - 1.0;
        grad[5] = x[0];
        grad[6] = 0.0;
        grad[7] = 0.0;
        grad[8] = 0.0;
        grad[9] = -1.0;
        grad[10] = 0.0;
        grad[11] = 0.0;
        grad[12] = 0.0;
        grad[13] = 0.0;
        grad[14] = 1.0;
        grad[15] = 4.0/3.0;
        grad[16] = 0.0;
        grad[17] = 0.0;
        grad[18] = -1.0;
        grad[19] = 0.0;
        grad[20] = 0.0;
        grad[21] = 0.0;
        grad[22] = 0.0;
        grad[23] = -1.0;
        grad[24] = 0.0;
        grad[25] = 0.0;
        grad[26] = 0.0;
        grad[27] = -1./3.;
    }

    return;
};

/**
    Inequality constraints for sp_mp
*/
void sp_mp_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( -x[1]);
    result[1] = ( x[1] + x[2] - 1.0);
    result[2] = ( -x[2]);
    result[3] = ( x[0] - 1.0);
    result[4] = ( -x[0]);

    if (grad) {
        grad[0] = 0.0;
        grad[1] = -1.0;
        grad[2] = 0.0;
        grad[3] = 0.0;
        grad[4] = 1.0;
        grad[5] = 1.0;
        grad[6] = 0.0;
        grad[7] = 0.0;
        grad[8] = -1.0;
        grad[9] = 1.0;
        grad[10] = 0.0;
        grad[11] = 0.0;
        grad[12] = -1.0;
        grad[13] = 0.0;
        grad[14] = 0.0;
    }

    return;
};

/**
    Inequality constraints for sa_mp
*/
void sa_mp_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( -x[2]*x[0] + x[2] - 0.75*x[3] - x[0]*x[1] + x[0] + x[1] - 1.0);
    result[1] = ( x[2]*x[0] + 0.75*x[3] + x[0]*x[1] - x[0]);
    result[2] = ( -x[2]);
    result[3] = ( -x[1]);
    result[4] = ( 0.25*x[3] + x[0] - 1.0);
    result[5] = ( -0.25*x[3] - x[0]);
    result[6] = ( x[2] + x[1] - 1.0);
    result[7] = ( -x[2] - x[1]);

    if (grad) {
        grad[0] = -x[2] - x[1] + 1.0;
        grad[1] = 1.0 - x[0];
        grad[2] = 1.0 - x[0];
        grad[3] = -0.75;
        grad[4] = x[2] + x[1] - 1.0;
        grad[5] = x[0];
        grad[6] = x[0];
        grad[7] = 0.75;
        grad[8] = 0.0;
        grad[9] = 0.0;
        grad[10] = -1.0;
        grad[11] = 0.0;
        grad[12] = 0.0;
        grad[13] = -1.0;
        grad[14] = 0.0;
        grad[15] = 0.0;
        grad[16] = 1.0;
        grad[17] = 0.0;
        grad[18] = 0.0;
        grad[19] = 0.25;
        grad[20] = -1.0;
        grad[21] = 0.0;
        grad[22] = 0.0;
        grad[23] = -0.25;
        grad[24] = 0.0;
        grad[25] = 1.0;
        grad[26] = 1.0;
        grad[27] = 0.0;
        grad[28] = 0.0;
        grad[29] = -1.0;
        grad[30] = -1.0;
        grad[31] = 0.0;
    }

    return;
};

/**
    Inequality constraints for fsp_mp
*/
void fsp_mp_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( x[0] + x[1] - 1.0);
    result[1] = ( -x[0]);
    result[2] = ( -x[1]);
    result[3] = ( -0.25*x[0] - 0.25);
    result[4] = ( 0.25*x[0] - 0.75);

    if (grad) {
        grad[0] = 1.0;
        grad[1] = 1.0;
        grad[2] = -1.0;
        grad[3] = 0.0;
        grad[4] = 0.0;
        grad[5] = -1.0;
        grad[6] = -0.25;
        grad[7] = 0.0;
        grad[8] = 0.25;
        grad[9] = 0.0;
    }

    return;
};

/**
    Inequality constraints for opx_mp
*/
void opx_mp_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( 0.5*x[4]*x[5] - x[3]*x[0] + x[3] + 0.5*x[1]*x[5] - x[1]*x[0] + x[1] - 0.5*x[5] - x[0]*x[2] + x[0] + x[2] - 1.0);
    result[1] = ( -0.5*x[4]*x[5] + x[3]*x[0] - 0.5*x[1]*x[5] + x[1]*x[0] + 0.5*x[5] + x[0]*x[2] - x[0]);
    result[2] = ( -x[1]);
    result[3] = ( -x[3]);
    result[4] = ( -x[2]);
    result[5] = ( -0.5*x[4]*x[5] - x[4]*x[0] + x[4] - 0.5*x[1]*x[5] - x[1]*x[0] + x[1] + 0.5*x[5] + x[0] - 1.0);
    result[6] = ( 0.5*x[4]*x[5] + x[4]*x[0] + 0.5*x[1]*x[5] + x[1]*x[0] - 0.5*x[5] - x[0]);
    result[7] = ( -x[1]);
    result[8] = ( -x[4]);
    result[9] = ( -0.5*x[3] - 0.5*x[2]);
    result[10] = ( 0.5*x[3] + 0.5*x[2] - 1.0);

    if (grad) {
        grad[0] = -x[3] - x[1] - x[2] + 1.0;
        grad[1] = 0.5*x[5] - x[0] + 1.0;
        grad[2] = 1.0 - x[0];
        grad[3] = 1.0 - x[0];
        grad[4] = 0.5*x[5];
        grad[5] = 0.5*x[4] + 0.5*x[1] - 0.5;
        grad[6] = x[3] + x[1] + x[2] - 1.0;
        grad[7] = -0.5*x[5] + x[0];
        grad[8] = x[0];
        grad[9] = x[0];
        grad[10] = -0.5*x[5];
        grad[11] = -0.5*x[4] - 0.5*x[1] + 0.5;
        grad[12] = 0.0;
        grad[13] = -1.0;
        grad[14] = 0.0;
        grad[15] = 0.0;
        grad[16] = 0.0;
        grad[17] = 0.0;
        grad[18] = 0.0;
        grad[19] = 0.0;
        grad[20] = 0.0;
        grad[21] = -1.0;
        grad[22] = 0.0;
        grad[23] = 0.0;
        grad[24] = 0.0;
        grad[25] = 0.0;
        grad[26] = -1.0;
        grad[27] = 0.0;
        grad[28] = 0.0;
        grad[29] = 0.0;
        grad[30] = -x[4] - x[1] + 1.0;
        grad[31] = -0.5*x[5] - x[0] + 1.0;
        grad[32] = 0.0;
        grad[33] = 0.0;
        grad[34] = -0.5*x[5] - x[0] + 1.0;
        grad[35] = -0.5*x[4] - 0.5*x[1] + 0.5;
        grad[36] = x[4] + x[1] - 1.0;
        grad[37] = 0.5*x[5] + x[0];
        grad[38] = 0.0;
        grad[39] = 0.0;
        grad[40] = 0.5*x[5] + x[0];
        grad[41] = 0.5*x[4] + 0.5*x[1] - 0.5;
        grad[42] = 0.0;
        grad[43] = -1.0;
        grad[44] = 0.0;
        grad[45] = 0.0;
        grad[46] = 0.0;
        grad[47] = 0.0;
        grad[48] = 0.0;
        grad[49] = 0.0;
        grad[50] = 0.0;
        grad[51] = 0.0;
        grad[52] = -1.0;
        grad[53] = 0.0;
        grad[54] = 0.0;
        grad[55] = 0.0;
        grad[56] = -0.50;
        grad[57] = -0.50;
        grad[58] = 0.0;
        grad[59] = 0.0;
        grad[60] = 0.0;
        grad[61] = 0.0;
        grad[62] = 0.50;
        grad[63] = 0.50;
        grad[64] = 0.0;
        grad[65] = 0.0;
    }

    return;
};

/**
    Inequality constraints for mu_mp
*/
void mu_mp_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( x[4] + x[3] - 1.0);
    result[1] = ( -x[3]);
    result[2] = ( -x[4]);
    result[3] = ( -x[0]*x[1] + x[0] + x[1] - 1.0);
    result[4] = ( x[0]*x[1] - x[0]);
    result[5] = ( -x[1]);
    result[6] = ( x[2] - 1.0);
    result[7] = ( -x[2]);
    result[8] = ( 0.5*x[4] + 0.5*x[1] - 1.0);
    result[9] = ( -0.5*x[4] - 0.5*x[1]);

    if (grad) {
        grad[0] = 0.0;
        grad[1] = 0.0;
        grad[2] = 0.0;
        grad[3] = 1.0;
        grad[4] = 1.0;
        grad[5] = 0.0;
        grad[6] = 0.0;
        grad[7] = 0.0;
        grad[8] = -1.0;
        grad[9] = 0.0;
        grad[10] = 0.0;
        grad[11] = 0.0;
        grad[12] = 0.0;
        grad[13] = 0.0;
        grad[14] = -1.0;
        grad[15] = 1.0 - x[1];
        grad[16] = 1.0 - x[0];
        grad[17] = 0.0;
        grad[18] = 0.0;
        grad[19] = 0.0;
        grad[20] = x[1] - 1.0;
        grad[21] = x[0];
        grad[22] = 0.0;
        grad[23] = 0.0;
        grad[24] = 0.0;
        grad[25] = 0.0;
        grad[26] = -1.0;
        grad[27] = 0.0;
        grad[28] = 0.0;
        grad[29] = 0.0;
        grad[30] = 0.0;
        grad[31] = 0.0;
        grad[32] = 1.0;
        grad[33] = 0.0;
        grad[34] = 0.0;
        grad[35] = 0.0;
        grad[36] = 0.0;
        grad[37] = -1.0;
        grad[38] = 0.0;
        grad[39] = 0.0;
        grad[40] = 0.0;
        grad[41] = 0.50;
        grad[42] = 0.0;
        grad[43] = 0.0;
        grad[44] = 0.50;
        grad[45] = 0.0;
        grad[46] = -0.50;
        grad[47] = 0.0;
        grad[48] = 0.0;
        grad[49] = -0.50;
    }

    return;
};

/**
    Inequality constraints for mt_mp
*/
void mt_mp_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( 0.5*x[0] - 0.5);
    result[1] = ( 0.5*x[1] - x[0]);
    result[2] = ( -0.5*x[1] + 0.5*x[0] - 0.5);
    result[3] = ( -x[1]);
    result[4] = ( x[1] - 1.0);

    if (grad) {
        grad[0] = 0.50;
        grad[1] = 0.0;
        grad[2] = -1.0;
        grad[3] = 0.50;
        grad[4] = 0.50;
        grad[5] = -0.50;
        grad[6] = 0.0;
        grad[7] = -1.0;
        grad[8] = 0.0;
        grad[9] = 1.0;
    }

    return;
};

/**
    Inequality constraints for ma_mp
*/
void ma_mp_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( x[4] + x[3] - 1.0);
    result[1] = ( -x[3]);
    result[2] = ( -x[4]);
    result[3] = ( -x[0]*x[1] + x[0] + x[1] - 1.0);
    result[4] = ( x[0]*x[1] - x[0]);
    result[5] = ( -x[1]);
    result[6] = ( x[2] - 1.0);
    result[7] = ( -x[2]);
    result[8] = ( 0.5*x[4] + 0.5*x[1] - 1.0);
    result[9] = ( -0.5*x[4] - 0.5*x[1]);

    if (grad) {
        grad[0] = 0.0;
        grad[1] = 0.0;
        grad[2] = 0.0;
        grad[3] = 1.0;
        grad[4] = 1.0;
        grad[5] = 0.0;
        grad[6] = 0.0;
        grad[7] = 0.0;
        grad[8] = -1.0;
        grad[9] = 0.0;
        grad[10] = 0.0;
        grad[11] = 0.0;
        grad[12] = 0.0;
        grad[13] = 0.0;
        grad[14] = -1.0;
        grad[15] = 1.0 - x[1];
        grad[16] = 1.0 - x[0];
        grad[17] = 0.0;
        grad[18] = 0.0;
        grad[19] = 0.0;
        grad[20] = x[1] - 1.0;
        grad[21] = x[0];
        grad[22] = 0.0;
        grad[23] = 0.0;
        grad[24] = 0.0;
        grad[25] = 0.0;
        grad[26] = -1.0;
        grad[27] = 0.0;
        grad[28] = 0.0;
        grad[29] = 0.0;
        grad[30] = 0.0;
        grad[31] = 0.0;
        grad[32] = 1.0;
        grad[33] = 0.0;
        grad[34] = 0.0;
        grad[35] = 0.0;
        grad[36] = 0.0;
        grad[37] = -1.0;
        grad[38] = 0.0;
        grad[39] = 0.0;
        grad[40] = 0.0;
        grad[41] = 0.50;
        grad[42] = 0.0;
        grad[43] = 0.0;
        grad[44] = 0.50;
        grad[45] = 0.0;
        grad[46] = -0.50;
        grad[47] = 0.0;
        grad[48] = 0.0;
        grad[49] = -0.50;
    }

    return;
};

/**
    Inequality constraints for ilm
*/
void ilm_mp_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( -0.5*x[0] - 0.5*x[1]);
    result[1] = ( -0.5*x[0] + 0.5*x[1]);
    result[2] = ( x[0] - 1.0);
    result[3] = ( -0.5*x[0] + 0.5*x[1]);
    result[4] = ( -0.5*x[0] - 0.5*x[1]);
    result[5] = ( x[0] - 1.0);

    if (grad) {
        grad[0] = -0.500;
        grad[1] = -0.500;
        grad[2] = -0.500;
        grad[3] = 0.500;
        grad[4] = 1.00;
        grad[5] = 0.0;
        grad[6] = -0.500;
        grad[7] = 0.500;
        grad[8] = -0.500;
        grad[9] = -0.500;
        grad[10] = 1.00;
        grad[11] = 0.0;
    }

    return;
};
/**
    Inequality constraints for ilmm_mp
*/
void ilmm_mp_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( -0.5*x[0] + 0.5*x[1] + 0.5*x[2] - 0.5*x[3]);
    result[1] = ( -0.5*x[0] + 0.5*x[1] + 0.5*x[2] + 0.5*x[3]);
    result[2] = ( - x[1]);
    result[3] = ( - x[2]);
    result[4] = (  x[0] - 1.0);
    result[5] = ( -0.5*x[0] + 0.5*x[1] + 0.5*x[2] + 0.5*x[3]);
    result[6] = ( -0.5*x[0] - 0.5*x[1] - 0.5*x[2] - 0.5*x[3]);

    if (grad) {
        grad[0] = -0.50;
        grad[1] = 0.50;
        grad[2] = 0.50;
        grad[3] = -0.50;
        grad[4] = -0.50;
        grad[5] = 0.50;
        grad[6] = 0.50;
        grad[7] = 0.50;
        grad[8] = 0.0;
        grad[9] = -1.0;
        grad[10] = 0.0;
        grad[11] = 0.0;
        grad[12] = 0.0;
        grad[13] = 0.0;
        grad[14] = -1.0;
        grad[15] = 0.0;
        grad[16] = 1.0;
        grad[17] = 0.0;
        grad[18] = 0.0;
        grad[19] = 0.0;
        grad[20] = -0.50;
        grad[21] = 0.50;
        grad[22] = 0.50;
        grad[23] = 0.50;
        grad[24] = -0.50;
        grad[25] = -0.50;
        grad[26] = -0.50;
        grad[27] = -0.50;
    }

    return;
};

/**
    Inequality constraints for g_mp
*/
void g_mp_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( - x[2]*x[0] + x[2] - x[0]*x[1] + x[0] + x[1] - 1.0);
    result[1] = ( x[2]*x[0] + x[0]*x[1] - x[0]);
    result[2] = ( - x[2]);
    result[3] = ( - x[1]);
    result[4] = ( - 1.0 + x[3]);
    result[5] = ( - x[3]);

    if (grad) {
        grad[0] = -x[2] - x[1] + 1.0;
        grad[1] = 1.0 - x[0];
        grad[2] = 1.0 - x[0];
        grad[3] = 0.0;
        grad[4] = x[2] + x[1] - 1.0;
        grad[5] = x[0];
        grad[6] = x[0];
        grad[7] = 0.0;
        grad[8] = 0.0;
        grad[9] = 0.0;
        grad[10] = -1.0;
        grad[11] = 0.0;
        grad[12] = 0.0;
        grad[13] = -1.0;
        grad[14] = 0.0;
        grad[15] = 0.0;
        grad[16] = 0.0;
        grad[17] = 0.0;
        grad[18] = 0.0;
        grad[19] = 1.0;
        grad[20] = 0.0;
        grad[21] = 0.0;
        grad[22] = 0.0;
        grad[23] = -1.0;
    }

    return;
};

/**
    Inequality constraints for ep_mp
*/
void ep_mp_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( -x[0] + x[1]);
    result[1] = ( x[0] - x[1] - 1.0);
    result[2] = ( -x[0] - x[1]);
    result[3] = ( x[0] + x[1] - 1.0);

    if (grad) {
        grad[0] = -1.0;
        grad[1] = 1.0;
        grad[2] = 1.0;
        grad[3] = -1.0;
        grad[4] = -1.0;
        grad[5] = -1.0;
        grad[6] = 1.0;
        grad[7] = 1.0;
    }

    return;
};

/**
    Inequality constraints for ctd_mp
*/
void ctd_mp_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( x[2] - 1.0);
    result[1] = ( -x[2]);
    result[2] = ( x[1]*x[0] - x[0]);
    result[3] = ( -x[1]*x[0] + x[1] + x[0] - 1.0);
    result[4] = ( -x[1]);

    if (grad) {
        grad[0] = 0.0;
        grad[1] = 0.0;
        grad[2] = 1.0;
        grad[3] = 0.0;
        grad[4] = 0.0;
        grad[5] = -1.0;
        grad[6] = x[1] - 1.0;
        grad[7] = x[0];
        grad[8] = 0.0;
        grad[9] = 1.0 - x[1];
        grad[10] = 1.0 - x[0];
        grad[11] = 0.0;
        grad[12] = 0.0;
        grad[13] = -1.0;
        grad[14] = 0.0;
    }

    return;
}

/**
    Inequality constraints for chl_mp
*/
void chl_mp_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( x[3]*x[5] - x[3]*x[0] + x[3] - x[5]*x[4] + x[5]*x[1] - x[5] + x[4]*x[0] - x[4] - x[0]*x[1] + x[0] + x[1] - 1.0);
    result[1] = ( -x[3]*x[5] + x[3]*x[0] + x[5]*x[4] - x[5]*x[1] + x[5] - x[4]*x[0] + x[0]*x[1] - x[0]);
    result[2] = ( x[4] - x[1]);
    result[3] = ( -0.25*x[2]*x[6] - 0.25*x[3]*x[5] - x[3]*x[0] + x[3] + 0.25*x[5]*x[4] - 0.25*x[5]*x[1] + 0.25*x[5] - 0.25*x[6]*x[4] - 0.25*x[6]*x[1] + 0.25*x[6] + x[0] - 1.0);
    result[4] = ( -x[3]);
    result[5] = ( 0.25*x[2]*x[6] + 0.25*x[3]*x[5] + x[3]*x[0] - 0.25*x[5]*x[4] + 0.25*x[5]*x[1] - 0.25*x[5] + 0.25*x[6]*x[4] + 0.25*x[6]*x[1] - 0.25*x[6] - x[0]);
    result[6] = ( x[2]*x[6] - x[2]*x[0] + x[2] + x[6]*x[4] + x[6]*x[1] - x[6] - x[4]*x[0] + x[4] - x[0]*x[1] + x[0] + x[1] - 1.0);
    result[7] = ( -x[2]*x[6] + x[2]*x[0] - x[6]*x[4] - x[6]*x[1] + x[6] + x[4]*x[0] + x[0]*x[1] - x[0]);
    result[8] = ( -x[2]);
    result[9] = ( -x[4] - x[1]);
    result[10] = ( 0.5*x[2] + x[1] - 1.0);
    result[11] = ( -0.5*x[2] - x[1]);

    if (grad) {
        grad[0] = -x[3] + x[4] - x[1] + 1.0;
        grad[1] = x[5] - x[0] + 1.0;
        grad[2] = 0.0;
        grad[3] = x[5] - x[0] + 1.0;
        grad[4] = -x[5] + x[0] - 1.0;
        grad[5] = x[3] - x[4] + x[1] - 1.0;
        grad[6] = 0.0;
        grad[7] = x[3] - x[4] + x[1] - 1.0;
        grad[8] = -x[5] + x[0];
        grad[9] = 0.0;
        grad[10] = -x[5] + x[0];
        grad[11] = x[5] - x[0];
        grad[12] = -x[3] + x[4] - x[1] + 1.0;
        grad[13] = 0.0;
        grad[14] = 0.0;
        grad[15] = -1.0;
        grad[16] = 0.0;
        grad[17] = 0.0;
        grad[18] = 1.0;
        grad[19] = 0.0;
        grad[20] = 0.0;
        grad[21] = 1.0 - x[3];
        grad[22] = -0.25*x[5] - 0.25*x[6];
        grad[23] = -0.25*x[6];
        grad[24] = -0.25*x[5] - x[0] + 1.0;
        grad[25] = 0.25*x[5] - 0.25*x[6];
        grad[26] = -0.25*x[3] + 0.25*x[4] - 0.25*x[1] + 0.25;
        grad[27] = -0.25*x[2] - 0.25*x[4] - 0.25*x[1] + 0.25;
        grad[28] = 0.0;
        grad[29] = 0.0;
        grad[30] = 0.0;
        grad[31] = -1.0;
        grad[32] = 0.0;
        grad[33] = 0.0;
        grad[34] = 0.0;
        grad[35] = x[3] - 1.0;
        grad[36] = 0.25*x[5] + 0.25*x[6];
        grad[37] = 0.25*x[6];
        grad[38] = 0.25*x[5] + x[0];
        grad[39] = -0.25*x[5] + 0.25*x[6];
        grad[40] = 0.25*x[3] - 0.25*x[4] + 0.25*x[1] - 0.25;
        grad[41] = 0.25*x[2] + 0.25*x[4] + 0.25*x[1] - 0.25;
        grad[42] = -x[2] - x[4] - x[1] + 1.0;
        grad[43] = x[6] - x[0] + 1.0;
        grad[44] = x[6] - x[0] + 1.0;
        grad[45] = 0.0;
        grad[46] = x[6] - x[0] + 1.0;
        grad[47] = 0.0;
        grad[48] = x[2] + x[4] + x[1] - 1.0;
        grad[49] = x[2] + x[4] + x[1] - 1.0;
        grad[50] = -x[6] + x[0];
        grad[51] = -x[6] + x[0];
        grad[52] = 0.0;
        grad[53] = -x[6] + x[0];
        grad[54] = 0.0;
        grad[55] = -x[2] - x[4] - x[1] + 1.0;
        grad[56] = 0.0;
        grad[57] = 0.0;
        grad[58] = -1.0;
        grad[59] = 0.0;
        grad[60] = 0.0;
        grad[61] = 0.0;
        grad[62] = 0.0;
        grad[63] = 0.0;
        grad[64] = -1.0;
        grad[65] = 0.0;
        grad[66] = 0.0;
        grad[67] = -1.0;
        grad[68] = 0.0;
        grad[69] = 0.0;
        grad[70] = 0.0;
        grad[71] = 1.0;
        grad[72] = 0.50;
        grad[73] = 0.0;
        grad[74] = 0.0;
        grad[75] = 0.0;
        grad[76] = 0.0;
        grad[77] = 0.0;
        grad[78] = -1.0;
        grad[79] = -0.50;
        grad[80] = 0.0;
        grad[81] = 0.0;
        grad[82] = 0.0;
        grad[83] = 0.0;
    }

    return;
};

/**
    Inequality constraints for cd_mp
*/
void cd_mp_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( x[0]*x[1] - x[0]);
    result[1] = ( -x[0]*x[1] + x[0] + x[1] - 1.0);
    result[2] = ( -x[1]);
    result[3] = ( -x[2]);
    result[4] = ( x[2] - 1.0);

    if (grad) {
        grad[0] = x[1] - 1.0;
        grad[1] = x[0];
        grad[2] = 0.0;
        grad[3] = 1.0 - x[1];
        grad[4] = 1.0 - x[0];
        grad[5] = 0.0;
        grad[6] = 0.0;
        grad[7] = -1.0;
        grad[8] = 0.0;
        grad[9] = 0.0;
        grad[10] = 0.0;
        grad[11] = -1.0;
        grad[12] = 0.0;
        grad[13] = 0.0;
        grad[14] = 1.0;
    }

    return;
};

/**
    Inequality constraints for bi_mp
*/
void bi_mp_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = (  - x[3]*x[0] + x[3] - 3.0*x[1]*x[0] + x[1] + 2.0/3.0*x[5] - x[4]*x[0] + x[4] - x[0]*x[2] + x[0] + x[2] - 1.0);
    result[1] = (  - x[1]);
    result[2] = (  + x[3]*x[0] + 3.0*x[1]*x[0] - 2.0/3.0*x[5] + x[4]*x[0] + x[0]*x[2] - x[0]);
    result[3] = (  - x[3]);
    result[4] = (  - x[4]);
    result[5] = (  - x[2]);
    result[6] = (  + x[1] - 1.0/3.0*x[5] + x[0] - 1.0);
    result[7] = (  - x[1]);
    result[8] = (  + 1.0/3.0*x[5] - x[0]);
    result[9] = (  + 0.5*x[3] + 0.5*x[2] - 0.5);
    result[10] = ( - 0.5*x[3] - 0.5*x[2] - 0.5);
    result[11] = ( x[4] - 1.0);
    result[12] = ( - x[4]);

    if (grad) {
        grad[0] = -x[3] - 3.0*x[1] - x[4] - x[2] + 1.0;
        grad[1] = 1.0 - 3.0*x[0];
        grad[2] = 1.0 - x[0];
        grad[3] = 1.0 - x[0];
        grad[4] = 1.0 - x[0];
        grad[5] = 2.0/3.0;
        grad[6] = 0.0;
        grad[7] = -1.0;
        grad[8] = 0.0;
        grad[9] = 0.0;
        grad[10] = 0.0;
        grad[11] = 0.0;
        grad[12] = x[3] + 3.0*x[1] + x[4] + x[2] - 1.0;
        grad[13] = 3.0*x[0];
        grad[14] = x[0];
        grad[15] = x[0];
        grad[16] = x[0];
        grad[17] = -2.0/3.0;
        grad[18] = 0.0;
        grad[19] = 0.0;
        grad[20] = 0.0;
        grad[21] = -1.0;
        grad[22] = 0.0;
        grad[23] = 0.0;
        grad[24] = 0.0;
        grad[25] = 0.0;
        grad[26] = 0.0;
        grad[27] = 0.0;
        grad[28] = -1.0;
        grad[29] = 0.0;
        grad[30] = 0.0;
        grad[31] = 0.0;
        grad[32] = -1.0;
        grad[33] = 0.0;
        grad[34] = 0.0;
        grad[35] = 0.0;
        grad[36] = 1.0;
        grad[37] = 1.0;
        grad[38] = 0.0;
        grad[39] = 0.0;
        grad[40] = 0.0;
        grad[41] = -1.0/3.0;
        grad[42] = 0.0;
        grad[43] = -1.0;
        grad[44] = 0.0;
        grad[45] = 0.0;
        grad[46] = 0.0;
        grad[47] = 0.0;
        grad[48] = -1.0;
        grad[49] = 0.0;
        grad[50] = 0.0;
        grad[51] = 0.0;
        grad[52] = 0.0;
        grad[53] = 1.0/3.0;
        grad[54] = 0.0;
        grad[55] = 0.0;
        grad[56] = 0.50;
        grad[57] = 0.50;
        grad[58] = 0.0;
        grad[59] = 0.0;
        grad[60] = 0.0;
        grad[61] = 0.0;
        grad[62] = -0.50;
        grad[63] = -0.50;
        grad[64] = 0.0;
        grad[65] = 0.0;
        grad[66] = 0.0;
        grad[67] = 0.0;
        grad[68] = 0.0;
        grad[69] = 0.0;
        grad[70] = 1.0;
        grad[71] = 0.0;
        grad[72] = 0.0;
        grad[73] = 0.0;
        grad[74] = 0.0;
        grad[75] = 0.0;
        grad[76] = -1.0;
        grad[77] = 0.0;
    }

    return;
};

/**
    Inequality constraints for fper_S11
*/
void fper_ig_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( - x[0]);
    result[1] = ( x[0] - 1.0);

    if (grad) {
        grad[0] = -1.0;
        grad[1] = 1.0;
    }

    return;
};

/** 
  local minimization for clinopyroxene
*/
void cpx_ig_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( -x[8]*x[4] - x[8]*x[0] + x[8] - x[3]*x[4] - x[3]*x[0] + x[3] + x[4]*x[7] - x[4]*x[1] + x[4] + x[7]*x[0] - x[7] - x[0]*x[1] + x[0] + x[1] - 1.0);
    result[1] = ( x[8]*x[4] + x[8]*x[0] + x[3]*x[4] + x[3]*x[0] - x[4]*x[7] + x[4]*x[1] - x[4] - x[7]*x[0] + x[0]*x[1] - x[0]);
    result[2] = ( x[6] + x[5] - x[8] - x[3] + 2.0*x[7] - x[1]);
    result[3] = ( -x[5]);
    result[4] = ( -x[6]);
    result[5] = ( -x[7]);
    result[6] = ( x[8]*x[4] + x[3]*x[4] + x[2]*x[0] - x[2] - x[4]*x[7] + x[4]*x[1] - x[4]);
    result[7] = ( -x[8]*x[4] - x[3]*x[4] - x[2]*x[0] + x[4]*x[7] - x[4]*x[1] + x[4]);
    result[8] = ( x[8] + x[3] + x[2] - 1.0);
    result[9] = ( -x[3]);
    result[10] = ( -x[8]);
    result[11] = ( 0.5*x[1] - 1.0);
    result[12] = ( -0.5*x[1]);

    if (grad) {
        grad[0] = -x[8] - x[3] + x[7] - x[1] + 1.0;
        grad[1] = -x[4] - x[0] + 1.0;
        grad[2] = 0.0;
        grad[3] = -x[4] - x[0] + 1.0;
        grad[4] = -x[8] - x[3] + x[7] - x[1] + 1.0;
        grad[5] = 0.0;
        grad[6] = 0.0;
        grad[7] = x[4] + x[0] - 1.0;
        grad[8] = -x[4] - x[0] + 1.0;
        grad[9] = x[8] + x[3] - x[7] + x[1] - 1.0;
        grad[10] = x[4] + x[0];
        grad[11] = 0.0;
        grad[12] = x[4] + x[0];
        grad[13] = x[8] + x[3] - x[7] + x[1] - 1.0;
        grad[14] = 0.0;
        grad[15] = 0.0;
        grad[16] = -x[4] - x[0];
        grad[17] = x[4] + x[0];
        grad[18] = 0.0;
        grad[19] = -1.0;
        grad[20] = 0.0;
        grad[21] = -1.0;
        grad[22] = 0.0;
        grad[23] = 1.0;
        grad[24] = 1.0;
        grad[25] = 2.0;
        grad[26] = -1.0;
        grad[27] = 0.0;
        grad[28] = 0.0;
        grad[29] = 0.0;
        grad[30] = 0.0;
        grad[31] = 0.0;
        grad[32] = -1.0;
        grad[33] = 0.0;
        grad[34] = 0.0;
        grad[35] = 0.0;
        grad[36] = 0.0;
        grad[37] = 0.0;
        grad[38] = 0.0;
        grad[39] = 0.0;
        grad[40] = 0.0;
        grad[41] = 0.0;
        grad[42] = -1.0;
        grad[43] = 0.0;
        grad[44] = 0.0;
        grad[45] = 0.0;
        grad[46] = 0.0;
        grad[47] = 0.0;
        grad[48] = 0.0;
        grad[49] = 0.0;
        grad[50] = 0.0;
        grad[51] = 0.0;
        grad[52] = -1.0;
        grad[53] = 0.0;
        grad[54] = x[2];
        grad[55] = x[4];
        grad[56] = x[0] - 1.0;
        grad[57] = x[4];
        grad[58] = x[8] + x[3] - x[7] + x[1] - 1.0;
        grad[59] = 0.0;
        grad[60] = 0.0;
        grad[61] = -x[4];
        grad[62] = x[4];
        grad[63] = -x[2];
        grad[64] = -x[4];
        grad[65] = -x[0];
        grad[66] = -x[4];
        grad[67] = -x[8] - x[3] + x[7] - x[1] + 1.0;
        grad[68] = 0.0;
        grad[69] = 0.0;
        grad[70] = x[4];
        grad[71] = -x[4];
        grad[72] = 0.0;
        grad[73] = 0.0;
        grad[74] = 1.0;
        grad[75] = 1.0;
        grad[76] = 0.0;
        grad[77] = 0.0;
        grad[78] = 0.0;
        grad[79] = 0.0;
        grad[80] = 1.0;
        grad[81] = 0.0;
        grad[82] = 0.0;
        grad[83] = 0.0;
        grad[84] = -1.0;
        grad[85] = 0.0;
        grad[86] = 0.0;
        grad[87] = 0.0;
        grad[88] = 0.0;
        grad[89] = 0.0;
        grad[90] = 0.0;
        grad[91] = 0.0;
        grad[92] = 0.0;
        grad[93] = 0.0;
        grad[94] = 0.0;
        grad[95] = 0.0;
        grad[96] = 0.0;
        grad[97] = 0.0;
        grad[98] = -1.0;
        grad[99] = 0.0;
        grad[100] = 0.50;
        grad[101] = 0.0;
        grad[102] = 0.0;
        grad[103] = 0.0;
        grad[104] = 0.0;
        grad[105] = 0.0;
        grad[106] = 0.0;
        grad[107] = 0.0;
        grad[108] = 0.0;
        grad[109] = -0.50;
        grad[110] = 0.0;
        grad[111] = 0.0;
        grad[112] = 0.0;
        grad[113] = 0.0;
        grad[114] = 0.0;
        grad[115] = 0.0;
        grad[116] = 0.0;
    }

    return;
};  


/** 
  local minimization for epidote
*/
void ep_ig_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( -x[0] + x[1]);
    result[1] = ( x[0] - x[1] - 1.0);
    result[2] = ( -x[0] - x[1]);
    result[3] = ( x[0] + x[1] - 1.0);

    if (grad) {
        grad[0] = -1.0;
        grad[1] = 1.0;
        grad[2] = 1.0;
        grad[3] = -1.0;
        grad[4] = -1.0;
        grad[5] = -1.0;
        grad[6] = 1.0;
        grad[7] = 1.0;
    }

    return;
};


/** 
  local minimization for fluid
*/
void fl_ig_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( x[0] + x[1] + x[2] + x[3] + x[4] + x[5] + x[6] + x[7] + x[8] + x[9] - 1.0);
    result[1] = ( -x[1]);
    result[2] = ( -x[0]);
    result[3] = ( -x[2]);
    result[4] = ( -x[3]);
    result[5] = ( -x[4]);
    result[6] = ( -x[5]);
    result[7] = ( -x[6]);
    result[8] = ( -x[7]);
    result[9] = ( -x[8]);
    result[10] = ( -x[9]);
    result[11] = ( x[9] - 1.0);

    if (grad) {
        grad[0] = 1.00;
        grad[1] = 1.00;
        grad[2] = 1.00;
        grad[3] = 1.00;
        grad[4] = 1.00;
        grad[5] = 1.00;
        grad[6] = 1.00;
        grad[7] = 1.00;
        grad[8] = 1.00;
        grad[9] = 1.00;
        grad[10] = 0.0;
        grad[11] = -1.00;
        grad[12] = 0.0;
        grad[13] = 0.0;
        grad[14] = 0.0;
        grad[15] = 0.0;
        grad[16] = 0.0;
        grad[17] = 0.0;
        grad[18] = 0.0;
        grad[19] = 0.0;
        grad[20] = -1.00;
        grad[21] = 0.0;
        grad[22] = 0.0;
        grad[23] = 0.0;
        grad[24] = 0.0;
        grad[25] = 0.0;
        grad[26] = 0.0;
        grad[27] = 0.0;
        grad[28] = 0.0;
        grad[29] = 0.0;
        grad[30] = 0.0;
        grad[31] = 0.0;
        grad[32] = -1.00;
        grad[33] = 0.0;
        grad[34] = 0.0;
        grad[35] = 0.0;
        grad[36] = 0.0;
        grad[37] = 0.0;
        grad[38] = 0.0;
        grad[39] = 0.0;
        grad[40] = 0.0;
        grad[41] = 0.0;
        grad[42] = 0.0;
        grad[43] = -1.00;
        grad[44] = 0.0;
        grad[45] = 0.0;
        grad[46] = 0.0;
        grad[47] = 0.0;
        grad[48] = 0.0;
        grad[49] = 0.0;
        grad[50] = 0.0;
        grad[51] = 0.0;
        grad[52] = 0.0;
        grad[53] = 0.0;
        grad[54] = -1.00;
        grad[55] = 0.0;
        grad[56] = 0.0;
        grad[57] = 0.0;
        grad[58] = 0.0;
        grad[59] = 0.0;
        grad[60] = 0.0;
        grad[61] = 0.0;
        grad[62] = 0.0;
        grad[63] = 0.0;
        grad[64] = 0.0;
        grad[65] = -1.00;
        grad[66] = 0.0;
        grad[67] = 0.0;
        grad[68] = 0.0;
        grad[69] = 0.0;
        grad[70] = 0.0;
        grad[71] = 0.0;
        grad[72] = 0.0;
        grad[73] = 0.0;
        grad[74] = 0.0;
        grad[75] = 0.0;
        grad[76] = -1.00;
        grad[77] = 0.0;
        grad[78] = 0.0;
        grad[79] = 0.0;
        grad[80] = 0.0;
        grad[81] = 0.0;
        grad[82] = 0.0;
        grad[83] = 0.0;
        grad[84] = 0.0;
        grad[85] = 0.0;
        grad[86] = 0.0;
        grad[87] = -1.00;
        grad[88] = 0.0;
        grad[89] = 0.0;
        grad[90] = 0.0;
        grad[91] = 0.0;
        grad[92] = 0.0;
        grad[93] = 0.0;
        grad[94] = 0.0;
        grad[95] = 0.0;
        grad[96] = 0.0;
        grad[97] = 0.0;
        grad[98] = -1.00;
        grad[99] = 0.0;
        grad[100] = 0.0;
        grad[101] = 0.0;
        grad[102] = 0.0;
        grad[103] = 0.0;
        grad[104] = 0.0;
        grad[105] = 0.0;
        grad[106] = 0.0;
        grad[107] = 0.0;
        grad[108] = 0.0;
        grad[109] = -1.00;
        grad[110] = 0.0;
        grad[111] = 0.0;
        grad[112] = 0.0;
        grad[113] = 0.0;
        grad[114] = 0.0;
        grad[115] = 0.0;
        grad[116] = 0.0;
        grad[117] = 0.0;
        grad[118] = 0.0;
        grad[119] = 1.00;
    }

    return;
};



/** 
  local minimization for garnet
*/
void g_ig_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = (-x[0]*x[1] + x[0] + x[1] - 1.0);
    result[1] = (x[0]*x[1] -x[0]);
    result[2] = (-x[1]);
    result[3] = (x[2] + x[3] + 2.0*x[4] - 1.0);
    result[4] = (-x[3]);
    result[5] = (-x[2]);
    result[6] = (-x[4]);
    result[7] = (-x[4]);

    if (grad) {
        grad[0] = 1.0 -x[1];
        grad[1] = 1.0 -x[0];
        grad[2] = 0.0;
        grad[3] = 0.0;
        grad[4] = 0.0;
        grad[5] = x[1] - 1.0;
        grad[6] = x[0];
        grad[7] = 0.0;
        grad[8] = 0.0;
        grad[9] = 0.0;
        grad[10] = 0.0;
        grad[11] = -1.0;
        grad[12] = 0.0;
        grad[13] = 0.0;
        grad[14] = 0.0;
        grad[15] = 0.0;
        grad[16] = 0.0;
        grad[17] = 1.0;
        grad[18] = 1.0;
        grad[19] = 2.0;
        grad[20] = 0.0;
        grad[21] = 0.0;
        grad[22] = 0.0;
        grad[23] = -1.0;
        grad[24] = 0.0;
        grad[25] = 0.0;
        grad[26] = 0.0;
        grad[27] = -1.0;
        grad[28] = 0.0;
        grad[29] = 0.0;
        grad[30] = 0.0;
        grad[31] = 0.0;
        grad[32] = 0.0;
        grad[33] = 0.0;
        grad[34] = -1.0;
        grad[35] = 0.0;
        grad[36] = 0.0;
        grad[37] = 0.0;
        grad[38] = 0.0;
        grad[39] = -1.0;
    }

    return;
};


/** 
  local minimization for hornblende
*/
void hb_ig_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( x[3] - 1.0);
    result[1] = ( x[3]*x[4] - x[3]);
    result[2] = ( -x[3]*x[4]);
    result[3] = ( x[0] - x[8] - 1.0);
    result[4] = ( -x[0] + x[8]);
    result[5] = ( -x[0]*x[1] - x[0]*x[6] - x[0]*x[7] + x[0] + x[1]*x[9] + x[1] + x[6]*x[9] + x[6] + x[7]*x[9] + x[7] - x[9] - 1.0);
    result[6] = ( x[0]*x[1] + x[0]*x[6] + x[0]*x[7] - x[0] - x[1]*x[9] - x[6]*x[9] - x[7]*x[9] + x[9]);
    result[7] = ( -x[1]);
    result[8] = ( -x[6]);
    result[9] = ( -x[7]);
    result[10] = ( -x[5]);
    result[11] = ( -x[0]*x[2] - x[0]*x[5] + x[0] - x[1]*x[9] + x[2] + x[5] - x[6]*x[9] - x[7]*x[9] + 1.5*x[8] + x[9] - 1.0);
    result[12] = ( x[0]*x[2] + x[0]*x[5] - x[0] + x[1]*x[9] + x[6]*x[9] + x[7]*x[9] - 1.5*x[8] - x[9]);
    result[13] = ( -x[2]);
    result[14] = ( 0.5*x[1] - 0.5*x[2] + 0.25*x[3] + 0.5*x[6] + 0.5*x[7] - 1.0);
    result[15] = ( -0.5*x[1] + 0.5*x[2] - 0.25*x[3] - 0.5*x[6] - 0.5*x[7]);
    result[16] = ( x[7] - 1.0);
    result[17] = ( -x[7]);

    if (grad) {
        grad[0] = 0.0;
        grad[1] = 0.0;
        grad[2] = 0.0;
        grad[3] = 1.0;
        grad[4] = 0.0;
        grad[5] = 0.0;
        grad[6] = 0.0;
        grad[7] = 0.0;
        grad[8] = 0.0;
        grad[9] = 0.0;
        grad[10] = 0.0;
        grad[11] = 0.0;
        grad[12] = 0.0;
        grad[13] = x[4] - 1.0;
        grad[14] = x[3];
        grad[15] = 0.0;
        grad[16] = 0.0;
        grad[17] = 0.0;
        grad[18] = 0.0;
        grad[19] = 0.0;
        grad[20] = 0.0;
        grad[21] = 0.0;
        grad[22] = 0.0;
        grad[23] = -x[4];
        grad[24] = -x[3];
        grad[25] = 0.0;
        grad[26] = 0.0;
        grad[27] = 0.0;
        grad[28] = 0.0;
        grad[29] = 0.0;
        grad[30] = 1.0;
        grad[31] = 0.0;
        grad[32] = 0.0;
        grad[33] = 0.0;
        grad[34] = 0.0;
        grad[35] = 0.0;
        grad[36] = 0.0;
        grad[37] = 0.0;
        grad[38] = -1.0;
        grad[39] = 0.0;
        grad[40] = -1.0;
        grad[41] = 0.0;
        grad[42] = 0.0;
        grad[43] = 0.0;
        grad[44] = 0.0;
        grad[45] = 0.0;
        grad[46] = 0.0;
        grad[47] = 0.0;
        grad[48] = 1.0;
        grad[49] = 0.0;
        grad[50] = -x[1] - x[6] - x[7] + 1.0;
        grad[51] = -x[0] + x[9] + 1.0;
        grad[52] = 0.0;
        grad[53] = 0.0;
        grad[54] = 0.0;
        grad[55] = 0.0;
        grad[56] = -x[0] + x[9] + 1.0;
        grad[57] = -x[0] + x[9] + 1.0;
        grad[58] = 0.0;
        grad[59] = x[1] + x[6] + x[7] - 1.0;
        grad[60] = x[1] + x[6] + x[7] - 1.0;
        grad[61] = x[0] - x[9];
        grad[62] = 0.0;
        grad[63] = 0.0;
        grad[64] = 0.0;
        grad[65] = 0.0;
        grad[66] = x[0] - x[9];
        grad[67] = x[0] - x[9];
        grad[68] = 0.0;
        grad[69] = -x[1] - x[6] - x[7] + 1.0;
        grad[70] = 0.0;
        grad[71] = -1.0;
        grad[72] = 0.0;
        grad[73] = 0.0;
        grad[74] = 0.0;
        grad[75] = 0.0;
        grad[76] = 0.0;
        grad[77] = 0.0;
        grad[78] = 0.0;
        grad[79] = 0.0;
        grad[80] = 0.0;
        grad[81] = 0.0;
        grad[82] = 0.0;
        grad[83] = 0.0;
        grad[84] = 0.0;
        grad[85] = 0.0;
        grad[86] = -1.0;
        grad[87] = 0.0;
        grad[88] = 0.0;
        grad[89] = 0.0;
        grad[90] = 0.0;
        grad[91] = 0.0;
        grad[92] = 0.0;
        grad[93] = 0.0;
        grad[94] = 0.0;
        grad[95] = 0.0;
        grad[96] = 0.0;
        grad[97] = -1.0;
        grad[98] = 0.0;
        grad[99] = 0.0;
        grad[100] = 0.0;
        grad[101] = 0.0;
        grad[102] = 0.0;
        grad[103] = 0.0;
        grad[104] = 0.0;
        grad[105] = -1.0;
        grad[106] = 0.0;
        grad[107] = 0.0;
        grad[108] = 0.0;
        grad[109] = 0.0;
        grad[110] = -x[2] - x[5] + 1.0;
        grad[111] = -x[9];
        grad[112] = 1.0 - x[0];
        grad[113] = 0.0;
        grad[114] = 0.0;
        grad[115] = 1.0 - x[0];
        grad[116] = -x[9];
        grad[117] = -x[9];
        grad[118] = 1.500;
        grad[119] = -x[1] - x[6] - x[7] + 1.0;
        grad[120] = x[2] + x[5] - 1.0;
        grad[121] = x[9];
        grad[122] = x[0];
        grad[123] = 0.0;
        grad[124] = 0.0;
        grad[125] = x[0];
        grad[126] = x[9];
        grad[127] = x[9];
        grad[128] = -1.500;
        grad[129] = x[1] + x[6] + x[7] - 1.0;
        grad[130] = 0.0;
        grad[131] = 0.0;
        grad[132] = -1.0;
        grad[133] = 0.0;
        grad[134] = 0.0;
        grad[135] = 0.0;
        grad[136] = 0.0;
        grad[137] = 0.0;
        grad[138] = 0.0;
        grad[139] = 0.0;
        grad[140] = 0.0;
        grad[141] = 0.50;
        grad[142] = -0.50;
        grad[143] = 0.2500;
        grad[144] = 0.0;
        grad[145] = 0.0;
        grad[146] = 0.50;
        grad[147] = 0.50;
        grad[148] = 0.0;
        grad[149] = 0.0;
        grad[150] = 0.0;
        grad[151] = -0.50;
        grad[152] = 0.50;
        grad[153] = -0.2500;
        grad[154] = 0.0;
        grad[155] = 0.0;
        grad[156] = -0.50;
        grad[157] = -0.50;
        grad[158] = 0.0;
        grad[159] = 0.0;
        grad[160] = 0.0;
        grad[161] = 0.0;
        grad[162] = 0.0;
        grad[163] = 0.0;
        grad[164] = 0.0;
        grad[165] = 0.0;
        grad[166] = 0.0;
        grad[167] = 1.0;
        grad[168] = 0.0;
        grad[169] = 0.0;
        grad[170] = 0.0;
        grad[171] = 0.0;
        grad[172] = 0.0;
        grad[173] = 0.0;
        grad[174] = 0.0;
        grad[175] = 0.0;
        grad[176] = 0.0;
        grad[177] = -1.0;
        grad[178] = 0.0;
        grad[179] = 0.0;
    }

    return;
};



/**
    Inequality constraints for ilm
*/
void ilm_ig_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = (0.5*x[0]*x[1] - 0.5*x[0] - 0.5*x[2]);
    result[1] = (-0.5*x[0] + 0.5*x[3]);
    result[2] = (x[0] - 1.0);
    result[3] = (-0.5*x[0]*x[1] + 0.5*x[2] - 0.5*x[3]);
    result[4] = (0.5*x[0]*x[1] - 0.5*x[0] + 0.5*x[2]);
    result[5] = (-0.5*x[0] - 0.5*x[3]);
    result[6] = (x[0] - 1.0);
    result[7] = (-0.5*x[0]*x[1] - 0.5*x[2] + 0.5*x[3]);

    if (grad) {
        grad[0] = 0.5*x[1] - 0.5;
        grad[1] = 0.5*x[0];
        grad[2] = -0.50;
        grad[3] = 0.0;
        grad[4] = -0.50;
        grad[5] = 0.0;
        grad[6] = 0.0;
        grad[7] = 0.50;
        grad[8] = 1.0;
        grad[9] = 0.0;
        grad[10] = 0.0;
        grad[11] = 0.0;
        grad[12] = -0.5*x[1];
        grad[13] = -0.5*x[0];
        grad[14] = 0.50;
        grad[15] = -0.50;
        grad[16] = 0.5*x[1] - 0.5;
        grad[17] = 0.5*x[0];
        grad[18] = 0.50;
        grad[19] = 0.0;
        grad[20] = -0.50;
        grad[21] = 0.0;
        grad[22] = 0.0;
        grad[23] = -0.50;
        grad[24] = 1.0;
        grad[25] = 0.0;
        grad[26] = 0.0;
        grad[27] = 0.0;
        grad[28] = -0.5*x[1];
        grad[29] = -0.5*x[0];
        grad[30] = -0.50;
        grad[31] = 0.50;
    }

    return;
};



/**
    Inequality constraints for liqHw
*/
void liq_ig_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( x[6] +x[3] +x[2] +x[10] +x[5] +x[4] +x[8] +x[1] +x[7] +x[0] - 0.25*x[9]*(-3.0*x[6] - 3.0*x[3] - 3.0*x[2] - 3.0*x[10] - 3.0*x[5] - 3.0*x[4] - 3.0*x[8] - 3.0*x[1] - 3.0*x[7] - 3.0*x[0] + 4.0) - 1.0);
    result[1] = ( -0.75*x[1]*x[9] - x[1] +x[9]);
    result[2] = ( -0.75*x[0]*x[9] - x[0] +x[9]);
    result[3] = ( -0.75*x[4]*x[9] - x[4]);
    result[4] = ( -0.75*x[5]*x[9] - x[5]);
    result[5] = ( -0.75*x[6]*x[9] - x[6]);
    result[6] = ( -0.75*x[7]*x[9] - x[7]);
    result[7] = ( -0.75*x[8]*x[9] - x[8]);
    result[8] = ( -x[9]);
    result[9] = ( -x[3] - x[2] - 0.75*x[9]*(x[3] + x[2]));
    result[10] = ( 0.75*x[10]*x[9] +x[10] - 1.0);
    result[11] = ( -4.0*x[2]*(0.75*x[9] + 1.0));
    result[12] = ( -4.0*x[3]*(0.75*x[9] + 1.0));
    result[13] = ( -x[0]*(0.75*x[9] + 1.0) +x[9]);
    result[14] = ( -x[1]*(0.75*x[9] + 1.0) +x[9]);
    result[15] = ( 2.0*x[9] - (0.75*x[9] + 1.0)*(4.0*x[3] + 4.0*x[2] + x[1] + x[0]));
    result[16] = ( -x[10]*(0.75*x[9] + 1.0));
    result[17] = ( 0.75*x[10]*x[9] +x[10] - 1.0);

    if (grad) {
        grad[0] = 0.75*x[9] + 1.0;
        grad[1] = 0.75*x[9] + 1.0;
        grad[2] = 0.75*x[9] + 1.0;
        grad[3] = 0.75*x[9] + 1.0;
        grad[4] = 0.75*x[9] + 1.0;
        grad[5] = 0.75*x[9] + 1.0;
        grad[6] = 0.75*x[9] + 1.0;
        grad[7] = 0.75*x[9] + 1.0;
        grad[8] = 0.75*x[9] + 1.0;
        grad[9] = 0.75*x[6] + 0.75*x[3] + 0.75*x[2] + 0.75*x[10] + 0.75*x[5] + 0.75*x[4] + 0.75*x[8] + 0.75*x[1] + 0.75*x[7] + 0.75*x[0] - 1.0;
        grad[10] = 0.75*x[9] + 1.0;
        grad[11] = 0.0;
        grad[12] = -0.75*x[9] - 1.0;
        grad[13] = 0.0;
        grad[14] = 0.0;
        grad[15] = 0.0;
        grad[16] = 0.0;
        grad[17] = 0.0;
        grad[18] = 0.0;
        grad[19] = 0.0;
        grad[20] = 1.0 - 0.75*x[1];
        grad[21] = 0.0;
        grad[22] = -0.75*x[9] - 1.0;
        grad[23] = 0.0;
        grad[24] = 0.0;
        grad[25] = 0.0;
        grad[26] = 0.0;
        grad[27] = 0.0;
        grad[28] = 0.0;
        grad[29] = 0.0;
        grad[30] = 0.0;
        grad[31] = 1.0 - 0.75*x[0];
        grad[32] = 0.0;
        grad[33] = 0.0;
        grad[34] = 0.0;
        grad[35] = 0.0;
        grad[36] = 0.0;
        grad[37] = -0.75*x[9] - 1.0;
        grad[38] = 0.0;
        grad[39] = 0.0;
        grad[40] = 0.0;
        grad[41] = 0.0;
        grad[42] = -0.75*x[4];
        grad[43] = 0.0;
        grad[44] = 0.0;
        grad[45] = 0.0;
        grad[46] = 0.0;
        grad[47] = 0.0;
        grad[48] = 0.0;
        grad[49] = -0.75*x[9] - 1.0;
        grad[50] = 0.0;
        grad[51] = 0.0;
        grad[52] = 0.0;
        grad[53] = -0.75*x[5];
        grad[54] = 0.0;
        grad[55] = 0.0;
        grad[56] = 0.0;
        grad[57] = 0.0;
        grad[58] = 0.0;
        grad[59] = 0.0;
        grad[60] = 0.0;
        grad[61] = -0.75*x[9] - 1.0;
        grad[62] = 0.0;
        grad[63] = 0.0;
        grad[64] = -0.75*x[6];
        grad[65] = 0.0;
        grad[66] = 0.0;
        grad[67] = 0.0;
        grad[68] = 0.0;
        grad[69] = 0.0;
        grad[70] = 0.0;
        grad[71] = 0.0;
        grad[72] = 0.0;
        grad[73] = -0.75*x[9] - 1.0;
        grad[74] = 0.0;
        grad[75] = -0.75*x[7];
        grad[76] = 0.0;
        grad[77] = 0.0;
        grad[78] = 0.0;
        grad[79] = 0.0;
        grad[80] = 0.0;
        grad[81] = 0.0;
        grad[82] = 0.0;
        grad[83] = 0.0;
        grad[84] = 0.0;
        grad[85] = -0.75*x[9] - 1.0;
        grad[86] = -0.75*x[8];
        grad[87] = 0.0;
        grad[88] = 0.0;
        grad[89] = 0.0;
        grad[90] = 0.0;
        grad[91] = 0.0;
        grad[92] = 0.0;
        grad[93] = 0.0;
        grad[94] = 0.0;
        grad[95] = 0.0;
        grad[96] = 0.0;
        grad[97] = -1.00;
        grad[98] = 0.0;
        grad[99] = 0.0;
        grad[100] = 0.0;
        grad[101] = -0.75*x[9] - 1.0;
        grad[102] = -0.75*x[9] - 1.0;
        grad[103] = 0.0;
        grad[104] = 0.0;
        grad[105] = 0.0;
        grad[106] = 0.0;
        grad[107] = 0.0;
        grad[108] = -0.75*x[3] - 0.75*x[2];
        grad[109] = 0.0;
        grad[110] = 0.0;
        grad[111] = 0.0;
        grad[112] = 0.0;
        grad[113] = 0.0;
        grad[114] = 0.0;
        grad[115] = 0.0;
        grad[116] = 0.0;
        grad[117] = 0.0;
        grad[118] = 0.0;
        grad[119] = 0.75*x[10];
        grad[120] = 0.75*x[9] + 1.0;
        grad[121] = 0.0;
        grad[122] = 0.0;
        grad[123] = -3.0*x[9] - 4.0;
        grad[124] = 0.0;
        grad[125] = 0.0;
        grad[126] = 0.0;
        grad[127] = 0.0;
        grad[128] = 0.0;
        grad[129] = 0.0;
        grad[130] = -3.0*x[2];
        grad[131] = 0.0;
        grad[132] = 0.0;
        grad[133] = 0.0;
        grad[134] = 0.0;
        grad[135] = -3.0*x[9] - 4.0;
        grad[136] = 0.0;
        grad[137] = 0.0;
        grad[138] = 0.0;
        grad[139] = 0.0;
        grad[140] = 0.0;
        grad[141] = -3.0*x[3];
        grad[142] = 0.0;
        grad[143] = -0.75*x[9] - 1.0;
        grad[144] = 0.0;
        grad[145] = 0.0;
        grad[146] = 0.0;
        grad[147] = 0.0;
        grad[148] = 0.0;
        grad[149] = 0.0;
        grad[150] = 0.0;
        grad[151] = 0.0;
        grad[152] = 1.0 - 0.75*x[0];
        grad[153] = 0.0;
        grad[154] = 0.0;
        grad[155] = -0.75*x[9] - 1.0;
        grad[156] = 0.0;
        grad[157] = 0.0;
        grad[158] = 0.0;
        grad[159] = 0.0;
        grad[160] = 0.0;
        grad[161] = 0.0;
        grad[162] = 0.0;
        grad[163] = 1.0 - 0.75*x[1];
        grad[164] = 0.0;
        grad[165] = -0.75*x[9] - 1.0;
        grad[166] = -0.75*x[9] - 1.0;
        grad[167] = -3.0*x[9] - 4.0;
        grad[168] = -3.0*x[9] - 4.0;
        grad[169] = 0.0;
        grad[170] = 0.0;
        grad[171] = 0.0;
        grad[172] = 0.0;
        grad[173] = 0.0;
        grad[174] = -3.0*x[3] - 3.0*x[2] - 0.75*x[1] - 0.75*x[0] + 2.0;
        grad[175] = 0.0;
        grad[176] = 0.0;
        grad[177] = 0.0;
        grad[178] = 0.0;
        grad[179] = 0.0;
        grad[180] = 0.0;
        grad[181] = 0.0;
        grad[182] = 0.0;
        grad[183] = 0.0;
        grad[184] = 0.0;
        grad[185] = -0.75*x[10];
        grad[186] = -0.75*x[9] - 1.0;
        grad[187] = 0.0;
        grad[188] = 0.0;
        grad[189] = 0.0;
        grad[190] = 0.0;
        grad[191] = 0.0;
        grad[192] = 0.0;
        grad[193] = 0.0;
        grad[194] = 0.0;
        grad[195] = 0.0;
        grad[196] = 0.75*x[10];
        grad[197] = 0.75*x[9] + 1.0;
    }

    return;
};

/** 
  local minimization for muscovite
*/
void mu_ig_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( x[4] + x[3] - 1.0);
    result[1] = ( -x[3]);
    result[2] = ( -x[4]);
    result[3] = ( -x[0]*x[1] + x[0] + x[1] - 1.0);
    result[4] = ( x[0]*x[1] - x[0]);
    result[5] = ( -x[1]);
    result[6] = ( x[2] - 1.0);
    result[7] = ( -x[2]);
    result[8] = ( 0.5*x[4] + 0.5*x[1] - 1.0);
    result[9] = ( -0.5*x[4] - 0.5*x[1]);

    if (grad) {
        grad[0] = 0.0;
        grad[1] = 0.0;
        grad[2] = 0.0;
        grad[3] = 1.0;
        grad[4] = 1.0;
        grad[5] = 0.0;
        grad[6] = 0.0;
        grad[7] = 0.0;
        grad[8] = -1.0;
        grad[9] = 0.0;
        grad[10] = 0.0;
        grad[11] = 0.0;
        grad[12] = 0.0;
        grad[13] = 0.0;
        grad[14] = -1.0;
        grad[15] = 1.0 - x[1];
        grad[16] = 1.0 - x[0];
        grad[17] = 0.0;
        grad[18] = 0.0;
        grad[19] = 0.0;
        grad[20] = x[1] - 1.0;
        grad[21] = x[0];
        grad[22] = 0.0;
        grad[23] = 0.0;
        grad[24] = 0.0;
        grad[25] = 0.0;
        grad[26] = -1.0;
        grad[27] = 0.0;
        grad[28] = 0.0;
        grad[29] = 0.0;
        grad[30] = 0.0;
        grad[31] = 0.0;
        grad[32] = 1.0;
        grad[33] = 0.0;
        grad[34] = 0.0;
        grad[35] = 0.0;
        grad[36] = 0.0;
        grad[37] = -1.0;
        grad[38] = 0.0;
        grad[39] = 0.0;
        grad[40] = 0.0;
        grad[41] = 0.50;
        grad[42] = 0.0;
        grad[43] = 0.0;
        grad[44] = 0.50;
        grad[45] = 0.0;
        grad[46] = -0.50;
        grad[47] = 0.0;
        grad[48] = 0.0;
        grad[49] = -0.50;
    }

    return;
};


/** 
  local minimization for olivine
*/
void ol_ig_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( -x[2] + x[0] - 1.0);
    result[1] = ( x[2] - x[0]);
    result[2] = ( -x[1]*x[0] + x[1] + x[2] + x[0] - 1.0);
    result[3] = ( x[1]*x[0] - x[2] - x[0]);
    result[4] = ( -x[1]);

    if (grad) {
        grad[0] = 1.0;
        grad[1] = 0.0;
        grad[2] = -1.0;
        grad[3] = -1.0;
        grad[4] = 0.0;
        grad[5] = 1.0;
        grad[6] = 1.0 - x[1];
        grad[7] = 1.0 - x[0];
        grad[8] = 1.0;
        grad[9] = x[1] - 1.0;
        grad[10] = x[0];
        grad[11] = -1.0;
        grad[12] = 0.0;
        grad[13] = -1.0;
        grad[14] = 0.0;
    }

    return;
};


/** 
  local minimization for orthopyroxene
*/
void opx_ig_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( -x[7]*x[3] - x[7]*x[0] + x[7] + x[3]*x[5] - x[3]*x[1] + x[3] + x[5]*x[0] - x[5] - x[0]*x[1] + x[0] + x[1] - 1.0);
    result[1] = ( x[7]*x[3] + x[7]*x[0] - x[3]*x[5] + x[3]*x[1] - x[3] - x[5]*x[0] + x[0]*x[1] - x[0]);
    result[2] = ( x[6] + x[4] - x[7] + 2.0*x[5] - x[1]);
    result[3] = ( -x[4]);
    result[4] = ( -x[6]);
    result[5] = ( -x[5]);
    result[6] = ( -x[2]*x[0] + x[2] + x[7]*x[3] - x[7]*x[0] + x[7] - x[3]*x[5] + x[3]*x[1] - x[3] + x[0] - 1.0);
    result[7] = ( x[2]*x[0] - x[7]*x[3] + x[7]*x[0] + x[3]*x[5] - x[3]*x[1] + x[3] - x[0]);
    result[8] = ( -x[2]);
    result[9] = ( -x[7]);
    result[10] = ( 0.5*x[1] - 1.0);
    result[11] = ( -0.5*x[1]);

    if (grad) {
        grad[0] = -x[7] + x[5] - x[1] + 1.0;
        grad[1] = -x[3] - x[0] + 1.0;
        grad[2] = 0.0;
        grad[3] = -x[7] + x[5] - x[1] + 1.0;
        grad[4] = 0.0;
        grad[5] = x[3] + x[0] - 1.0;
        grad[6] = 0.0;
        grad[7] = -x[3] - x[0] + 1.0;
        grad[8] = x[7] - x[5] + x[1] - 1.0;
        grad[9] = x[3] + x[0];
        grad[10] = 0.0;
        grad[11] = x[7] - x[5] + x[1] - 1.0;
        grad[12] = 0.0;
        grad[13] = -x[3] - x[0];
        grad[14] = 0.0;
        grad[15] = x[3] + x[0];
        grad[16] = 0.0;
        grad[17] = -1.0;
        grad[18] = 0.0;
        grad[19] = 0.0;
        grad[20] = 1.0;
        grad[21] = 2.0;
        grad[22] = 1.0;
        grad[23] = -1.0;
        grad[24] = 0.0;
        grad[25] = 0.0;
        grad[26] = 0.0;
        grad[27] = 0.0;
        grad[28] = -1.0;
        grad[29] = 0.0;
        grad[30] = 0.0;
        grad[31] = 0.0;
        grad[32] = 0.0;
        grad[33] = 0.0;
        grad[34] = 0.0;
        grad[35] = 0.0;
        grad[36] = 0.0;
        grad[37] = 0.0;
        grad[38] = -1.0;
        grad[39] = 0.0;
        grad[40] = 0.0;
        grad[41] = 0.0;
        grad[42] = 0.0;
        grad[43] = 0.0;
        grad[44] = 0.0;
        grad[45] = -1.0;
        grad[46] = 0.0;
        grad[47] = 0.0;
        grad[48] = -x[2] - x[7] + 1.0;
        grad[49] = x[3];
        grad[50] = 1.0 - x[0];
        grad[51] = x[7] - x[5] + x[1] - 1.0;
        grad[52] = 0.0;
        grad[53] = -x[3];
        grad[54] = 0.0;
        grad[55] = x[3] - x[0] + 1.0;
        grad[56] = x[2] + x[7] - 1.0;
        grad[57] = -x[3];
        grad[58] = x[0];
        grad[59] = -x[7] + x[5] - x[1] + 1.0;
        grad[60] = 0.0;
        grad[61] = x[3];
        grad[62] = 0.0;
        grad[63] = -x[3] + x[0];
        grad[64] = 0.0;
        grad[65] = 0.0;
        grad[66] = -1.0;
        grad[67] = 0.0;
        grad[68] = 0.0;
        grad[69] = 0.0;
        grad[70] = 0.0;
        grad[71] = 0.0;
        grad[72] = 0.0;
        grad[73] = 0.0;
        grad[74] = 0.0;
        grad[75] = 0.0;
        grad[76] = 0.0;
        grad[77] = 0.0;
        grad[78] = 0.0;
        grad[79] = -1.0;
        grad[80] = 0.0;
        grad[81] = 0.50;
        grad[82] = 0.0;
        grad[83] = 0.0;
        grad[84] = 0.0;
        grad[85] = 0.0;
        grad[86] = 0.0;
        grad[87] = 0.0;
        grad[88] = 0.0;
        grad[89] = -0.50;
        grad[90] = 0.0;
        grad[91] = 0.0;
        grad[92] = 0.0;
        grad[93] = 0.0;
        grad[94] = 0.0;
        grad[95] = 0.0;
    }

    return;
};


/** 
  local minimization for felsdpar
*/
void fsp_ig_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( x[0] + x[1] - 1.0);
    result[1] = ( -x[0]);
    result[2] = ( -x[1]);
    result[3] = ( -0.25*x[0] - 0.25);
    result[4] = ( 0.25*x[0] - 0.75);

    if (grad) {
        grad[0] = 1.0;
        grad[1] = 1.0;
        grad[2] = -1.0;
        grad[3] = 0.0;
        grad[4] = 0.0;
        grad[5] = -1.0;
        grad[6] = -0.250;
        grad[7] = 0.0;
        grad[8] = 0.250;
        grad[9] = 0.0;
    }

    return;
};


/** 
  local minimization for spinel
*/
void spl_ig_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *SS_ref_db){
    result[0] = ( 1./3.*x[0]*x[3] + 1./3.*x[0] - 1./3.*x[3] - 2./3.*x[4] - 1./3.);
    result[1] = ( -1./3.*x[0]*x[3] - 1./3.*x[0] - 2./3.*x[5]);
    result[2] = ( -2./3.*x[1]*x[2] - 2./3.*x[1]*x[3] + 2./3.*x[1] + 1./3.*x[3] + 2./3.*x[4] + 2./3.*x[5] + 2./3.*x[6] - 2./3.);
    result[3] = ( 2./3.*x[1]*x[2] + 2./3.*x[1]*x[3] - 2./3.*x[1] - 2./3.*x[6]);
    result[4] = ( 1./3.*x[0]*x[3] + 1./3.*x[0] - 1./3.*x[3] + 1./3.*x[4] - 1./3.);
    result[5] = ( -1./3.*x[0]*x[3] - 1./3.*x[0] + 1./3.*x[5]);
    result[6] = ( -2./3.*x[1]*x[2] - 2./3.*x[1]*x[3] + 2./3.*x[1] + x[2] + 5.0/6.0*x[3] - 1./3.*x[4] - 1./3.*x[5] - 1./3.*x[6] - 2./3.);
    result[7] = ( 2./3.*x[1]*x[2] + 2./3.*x[1]*x[3] - 2./3.*x[1] + 1./3.*x[6]);
    result[8] = ( -x[2]);
    result[9] = ( -0.5*x[3]);

   if (grad) {
        grad[0] = 1./3.*x[3] + 1./3.;
        grad[1] = 0.0;
        grad[2] = 0.0;
        grad[3] = 1./3.*x[0] - 1./3.;
        grad[4] = -2./3.;
        grad[5] = 0.0;
        grad[6] = 0.0;
        grad[7] = -1./3.*x[3] - 1./3.;
        grad[8] = 0.0;
        grad[9] = 0.0;
        grad[10] = -1./3.*x[0];
        grad[11] = 0.0;
        grad[12] = -2./3.;
        grad[13] = 0.0;
        grad[14] = 0.0;
        grad[15] = -2./3.*x[2] - 2./3.*x[3] + 2./3.;
        grad[16] = -2./3.*x[1];
        grad[17] = 1./3. - 2./3.*x[1];
        grad[18] = 2./3.;
        grad[19] = 2./3.;
        grad[20] = 2./3.;
        grad[21] = 0.0;
        grad[22] = 2./3.*x[2] + 2./3.*x[3] - 2./3.;
        grad[23] = 2./3.*x[1];
        grad[24] = 2./3.*x[1];
        grad[25] = 0.0;
        grad[26] = 0.0;
        grad[27] = -2./3.;
        grad[28] = 1./3.*x[3] + 1./3.;
        grad[29] = 0.0;
        grad[30] = 0.0;
        grad[31] = 1./3.*x[0] - 1./3.;
        grad[32] = 1./3.;
        grad[33] = 0.0;
        grad[34] = 0.0;
        grad[35] = -1./3.*x[3] - 1./3.;
        grad[36] = 0.0;
        grad[37] = 0.0;
        grad[38] = -1./3.*x[0];
        grad[39] = 0.0;
        grad[40] = 1./3.;
        grad[41] = 0.0;
        grad[42] = 0.0;
        grad[43] = -2./3.*x[2] - 2./3.*x[3] + 2./3.;
        grad[44] = 1.0 - 2./3.*x[1];
        grad[45] = 5.0/6.0 - 2./3.*x[1];
        grad[46] = -1./3.;
        grad[47] = -1./3.;
        grad[48] = -1./3.;
        grad[49] = 0.0;
        grad[50] = 2./3.*x[2] + 2./3.*x[3] - 2./3.;
        grad[51] = 2./3.*x[1];
        grad[52] = 2./3.*x[1];
        grad[53] = 0.0;
        grad[54] = 0.0;
        grad[55] = 1./3.;
        grad[56] = 0.0;
        grad[57] = 0.0;
        grad[58] = -1.0;
        grad[59] = 0.0;
        grad[60] = 0.0;
        grad[61] = 0.0;
        grad[62] = 0.0;
        grad[63] = 0.0;
        grad[64] = 0.0;
        grad[65] = 0.0;
        grad[66] = -0.50;
        grad[67] = 0.0;
        grad[68] = 0.0;
        grad[69] = 0.0;
    }

    return;
};

/** 
  local minimization for biotite
*/
void bi_ig_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( -x[2]*x[0] + x[2] + 2.0/3.0*x[4] - x[3]*x[0] + x[3] - x[0]*x[1] + x[0] + x[1] - 1.0);
    result[1] = ( x[2]*x[0] - 2.0/3.0*x[4] + x[3]*x[0] + x[0]*x[1] - x[0]);
    result[2] = ( -x[2]);
    result[3] = ( -x[3]);
    result[4] = ( -x[1]);
    result[5] = ( -1.0/3.0*x[4] + x[0] - 1.0);
    result[6] = ( 1.0/3.0*x[4] - x[0]);
    result[7] = ( 0.5*x[2] + 0.5*x[1] - 0.5);
    result[8] = ( -0.5*x[2] - 0.5*x[1] - 0.5);
    result[9] = ( x[3] - 1.0);
    result[10] = ( -x[3]);

    if (grad) {
        grad[0] = -x[2] - x[3] - x[1] + 1.0;
        grad[1] = 1.0 - x[0];
        grad[2] = 1.0 - x[0];
        grad[3] = 1.0 - x[0];
        grad[4] = 2.0/3.0;
        grad[5] = x[2] + x[3] + x[1] - 1.0;
        grad[6] = x[0];
        grad[7] = x[0];
        grad[8] = x[0];
        grad[9] = -2.0/3.0;
        grad[10] = 0.0;
        grad[11] = 0.0;
        grad[12] = -1.0;
        grad[13] = 0.0;
        grad[14] = 0.0;
        grad[15] = 0.0;
        grad[16] = 0.0;
        grad[17] = 0.0;
        grad[18] = -1.0;
        grad[19] = 0.0;
        grad[20] = 0.0;
        grad[21] = -1.0;
        grad[22] = 0.0;
        grad[23] = 0.0;
        grad[24] = 0.0;
        grad[25] = 1.0;
        grad[26] = 0.0;
        grad[27] = 0.0;
        grad[28] = 0.0;
        grad[29] = -1.0/3.0;
        grad[30] = -1.0;
        grad[31] = 0.0;
        grad[32] = 0.0;
        grad[33] = 0.0;
        grad[34] = 1.0/3.0;
        grad[35] = 0.0;
        grad[36] = 0.5;
        grad[37] = 0.5;
        grad[38] = 0.0;
        grad[39] = 0.0;
        grad[40] = 0.0;
        grad[41] = -0.5;
        grad[42] = -0.5;
        grad[43] = 0.0;
        grad[44] = 0.0;
        grad[45] = 0.0;
        grad[46] = 0.0;
        grad[47] = 0.0;
        grad[48] = 1.0;
        grad[49] = 0.0;
        grad[50] = 0.0;
        grad[51] = 0.0;
        grad[52] = 0.0;
        grad[53] = -1.0;
        grad[54] = 0.0;
    }

    return;
};

/** 
  local minimization for cordierite
*/
void cd_ig_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( -x[0]);
    result[1] = ( x[0] - 1.0);
    result[2] = ( -x[1]);
    result[3] = ( x[1] - 1.0);

    if (grad) {
        grad[0] = -1.0;
        grad[1] = 0.0;
        grad[2] = 1.0;
        grad[3] = 0.0;
        grad[4] = 0.0;
        grad[5] = -1.0;
        grad[6] = 0.0;
        grad[7] = 1.0;
    }

    return;
};


/**************************************************************************************/
/**************************************************************************************/
/********************IGNEOUS ALKALINE DATABASE (Weller et al., 2023)*******************/
/**************************************************************************************/
/**************************************************************************************/

/**
    Inequality constraints for liq
*/
void liq_igad_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = (x[0] - 1.0/6.0*x[10]*(-7.0*x[0] - 7.0*x[1] - 7.0*x[2] - 7.0*x[3] - 7.0*x[4] - 7.0*x[5] - 7.0*x[6] - 7.0*x[7] - 7.0*x[8] + 3.0) - 1.0/6.0*x[11]*(x[0] + x[1] + x[2] + x[3] + x[4] + x[5] + x[6] + x[7] + x[8] - 3.0) - 1.0/6.0*x[12]*(-7.0*x[0] - 7.0*x[1] - 7.0*x[2] - 7.0*x[3] - 7.0*x[4] - 7.0*x[5] - 7.0*x[6] - 7.0*x[7] - 7.0*x[8] + 3.0) + x[1] + x[2] + x[3] + x[4] + x[5] + x[6] + x[7] + x[8] -x[9]*(-x[0] -x[1] -x[2] -x[3] -x[4] -x[5] -x[6] -x[7] -x[8] + 1.0) - 1.0);
    result[1] = (-7.0/6.0*x[10]*x[1] + 0.5*x[10] + 1.0/6.0*x[11]*x[1] - 7.0/6.0*x[12]*x[1] + 0.5*x[12] -x[1]*x[9] -x[1] + x[9]);
    result[2] = (-7.0/6.0*x[0]*x[10] + 1.0/6.0*x[0]*x[11] - 7.0/6.0*x[0]*x[12] -x[0]*x[9] -x[0] + x[9]);
    result[3] = (-7.0/6.0*x[10]*x[4] + x[10] + 1.0/6.0*x[11]*x[4] - 7.0/6.0*x[12]*x[4] -x[4]*x[9] -x[4]);
    result[4] = (-7.0/6.0*x[10]*x[5] + 1.0/6.0*x[11]*x[5] - 7.0/6.0*x[12]*x[5] -x[5]*x[9] -x[5]);
    result[5] = (-7.0/6.0*x[10]*x[6] + 1.0/6.0*x[11]*x[6] - 7.0/6.0*x[12]*x[6] -x[6]*x[9] -x[6]);
    result[6] = (-7.0/6.0*x[10]*x[7] + 1.0/6.0*x[11]*x[7] - 7.0/6.0*x[12]*x[7] -x[7]*x[9] -x[7]);
    result[7] = (-7.0/6.0*x[10]*x[8] + 1.0/6.0*x[11]*x[8] - 7.0/6.0*x[12]*x[8] + x[12] -x[8]*x[9] -x[8]);
    result[8] = (-x[10]);
    result[9] = (-x[9]);
    result[10] = (-x[11]);
    result[11] = (-x[12]);
    result[12] = (-7.0/6.0*x[10]*(x[2] + x[3]) + 1.0/6.0*x[11]*(x[2] + x[3]) + 0.5*x[11] - 7.0/6.0*x[12]*(x[2] + x[3]) -x[2] -x[3] -x[9]*(x[2] + x[3]));
    result[13] = (2.0*x[11] - 4.0*x[2]*(7.0/6.0*x[10] - 1.0/6.0*x[11] + 7.0/6.0*x[12] + x[9] + 1.0));
    result[14] = (-4.0*x[3]*(7.0/6.0*x[10] - 1.0/6.0*x[11] + 7.0/6.0*x[12] + x[9] + 1.0));
    result[15] = (-x[0]*(7.0/6.0*x[10] - 1.0/6.0*x[11] + 7.0/6.0*x[12] + x[9] + 1.0) + x[9]);
    result[16] = (0.5*x[10] + 0.5*x[12] -x[1]*(7.0/6.0*x[10] - 1.0/6.0*x[11] + 7.0/6.0*x[12] + x[9] + 1.0) + x[9]);
    result[17] = (0.5*x[10] + 2.0*x[11] + 0.5*x[12] + 2.0*x[9] -(x[0] + x[1] + 4.0*x[2] + 4.0*x[3])*(7.0/6.0*x[10] - 1.0/6.0*x[11] + 7.0/6.0*x[12] + x[9] + 1.0));

    if (grad) {
        grad[0] = 7.0/6.0*x[10] - 1.0/6.0*x[11] + 7.0/6.0*x[12] + x[9] + 1.0;
        grad[1] = 7.0/6.0*x[10] - 1.0/6.0*x[11] + 7.0/6.0*x[12] + x[9] + 1.0;
        grad[2] = 7.0/6.0*x[10] - 1.0/6.0*x[11] + 7.0/6.0*x[12] + x[9] + 1.0;
        grad[3] = 7.0/6.0*x[10] - 1.0/6.0*x[11] + 7.0/6.0*x[12] + x[9] + 1.0;
        grad[4] = 7.0/6.0*x[10] - 1.0/6.0*x[11] + 7.0/6.0*x[12] + x[9] + 1.0;
        grad[5] = 7.0/6.0*x[10] - 1.0/6.0*x[11] + 7.0/6.0*x[12] + x[9] + 1.0;
        grad[6] = 7.0/6.0*x[10] - 1.0/6.0*x[11] + 7.0/6.0*x[12] + x[9] + 1.0;
        grad[7] = 7.0/6.0*x[10] - 1.0/6.0*x[11] + 7.0/6.0*x[12] + x[9] + 1.0;
        grad[8] = 7.0/6.0*x[10] - 1.0/6.0*x[11] + 7.0/6.0*x[12] + x[9] + 1.0;
        grad[9] = x[0] + x[1] + x[2] + x[3] + x[4] + x[5] + x[6] + x[7] + x[8] - 1.0;
        grad[10] = 7.0/6.0*x[0] + 7.0/6.0*x[1] + 7.0/6.0*x[2] + 7.0/6.0*x[3] + 7.0/6.0*x[4] + 7.0/6.0*x[5] + 7.0/6.0*x[6] + 7.0/6.0*x[7] + 7.0/6.0*x[8] - 0.5;
        grad[11] = -1.0/6.0*x[0] - 1.0/6.0*x[1] - 1.0/6.0*x[2] - 1.0/6.0*x[3] - 1.0/6.0*x[4] - 1.0/6.0*x[5] - 1.0/6.0*x[6] - 1.0/6.0*x[7] - 1.0/6.0*x[8] + 0.5;
        grad[12] = 7.0/6.0*x[0] + 7.0/6.0*x[1] + 7.0/6.0*x[2] + 7.0/6.0*x[3] + 7.0/6.0*x[4] + 7.0/6.0*x[5] + 7.0/6.0*x[6] + 7.0/6.0*x[7] + 7.0/6.0*x[8] - 0.5;
        grad[13] = 0.0;
        grad[14] = -7.0/6.0*x[10] + 1.0/6.0*x[11] - 7.0/6.0*x[12] -x[9] - 1.0;
        grad[15] = 0.0;
        grad[16] = 0.0;
        grad[17] = 0.0;
        grad[18] = 0.0;
        grad[19] = 0.0;
        grad[20] = 0.0;
        grad[21] = 0.0;
        grad[22] = 1.0 -x[1];
        grad[23] = 0.5 - 7.0/6.0*x[1];
        grad[24] = 1.0/6.0*x[1];
        grad[25] = 0.5 - 7.0/6.0*x[1];
        grad[26] = -7.0/6.0*x[10] + 1.0/6.0*x[11] - 7.0/6.0*x[12] -x[9] - 1.0;
        grad[27] = 0.0;
        grad[28] = 0.0;
        grad[29] = 0.0;
        grad[30] = 0.0;
        grad[31] = 0.0;
        grad[32] = 0.0;
        grad[33] = 0.0;
        grad[34] = 0.0;
        grad[35] = 1.0 -x[0];
        grad[36] = -7.0/6.0*x[0];
        grad[37] = 1.0/6.0*x[0];
        grad[38] = -7.0/6.0*x[0];
        grad[39] = 0.0;
        grad[40] = 0.0;
        grad[41] = 0.0;
        grad[42] = 0.0;
        grad[43] = -7.0/6.0*x[10] + 1.0/6.0*x[11] - 7.0/6.0*x[12] -x[9] - 1.0;
        grad[44] = 0.0;
        grad[45] = 0.0;
        grad[46] = 0.0;
        grad[47] = 0.0;
        grad[48] = -x[4];
        grad[49] = 1.0 - 7.0/6.0*x[4];
        grad[50] = 1.0/6.0*x[4];
        grad[51] = -7.0/6.0*x[4];
        grad[52] = 0.0;
        grad[53] = 0.0;
        grad[54] = 0.0;
        grad[55] = 0.0;
        grad[56] = 0.0;
        grad[57] = -7.0/6.0*x[10] + 1.0/6.0*x[11] - 7.0/6.0*x[12] -x[9] - 1.0;
        grad[58] = 0.0;
        grad[59] = 0.0;
        grad[60] = 0.0;
        grad[61] = -x[5];
        grad[62] = -7.0/6.0*x[5];
        grad[63] = 1.0/6.0*x[5];
        grad[64] = -7.0/6.0*x[5];
        grad[65] = 0.0;
        grad[66] = 0.0;
        grad[67] = 0.0;
        grad[68] = 0.0;
        grad[69] = 0.0;
        grad[70] = 0.0;
        grad[71] = -7.0/6.0*x[10] + 1.0/6.0*x[11] - 7.0/6.0*x[12] -x[9] - 1.0;
        grad[72] = 0.0;
        grad[73] = 0.0;
        grad[74] = -x[6];
        grad[75] = -7.0/6.0*x[6];
        grad[76] = 1.0/6.0*x[6];
        grad[77] = -7.0/6.0*x[6];
        grad[78] = 0.0;
        grad[79] = 0.0;
        grad[80] = 0.0;
        grad[81] = 0.0;
        grad[82] = 0.0;
        grad[83] = 0.0;
        grad[84] = 0.0;
        grad[85] = -7.0/6.0*x[10] + 1.0/6.0*x[11] - 7.0/6.0*x[12] -x[9] - 1.0;
        grad[86] = 0.0;
        grad[87] = -x[7];
        grad[88] = -7.0/6.0*x[7];
        grad[89] = 1.0/6.0*x[7];
        grad[90] = -7.0/6.0*x[7];
        grad[91] = 0.0;
        grad[92] = 0.0;
        grad[93] = 0.0;
        grad[94] = 0.0;
        grad[95] = 0.0;
        grad[96] = 0.0;
        grad[97] = 0.0;
        grad[98] = 0.0;
        grad[99] = -7.0/6.0*x[10] + 1.0/6.0*x[11] - 7.0/6.0*x[12] -x[9] - 1.0;
        grad[100] = -x[8];
        grad[101] = -7.0/6.0*x[8];
        grad[102] = 1.0/6.0*x[8];
        grad[103] = 1.0 - 7.0/6.0*x[8];
        grad[104] = 0.0;
        grad[105] = 0.0;
        grad[106] = 0.0;
        grad[107] = 0.0;
        grad[108] = 0.0;
        grad[109] = 0.0;
        grad[110] = 0.0;
        grad[111] = 0.0;
        grad[112] = 0.0;
        grad[113] = 0.0;
        grad[114] = -1.00;
        grad[115] = 0.0;
        grad[116] = 0.0;
        grad[117] = 0.0;
        grad[118] = 0.0;
        grad[119] = 0.0;
        grad[120] = 0.0;
        grad[121] = 0.0;
        grad[122] = 0.0;
        grad[123] = 0.0;
        grad[124] = 0.0;
        grad[125] = 0.0;
        grad[126] = -1.00;
        grad[127] = 0.0;
        grad[128] = 0.0;
        grad[129] = 0.0;
        grad[130] = 0.0;
        grad[131] = 0.0;
        grad[132] = 0.0;
        grad[133] = 0.0;
        grad[134] = 0.0;
        grad[135] = 0.0;
        grad[136] = 0.0;
        grad[137] = 0.0;
        grad[138] = 0.0;
        grad[139] = 0.0;
        grad[140] = 0.0;
        grad[141] = -1.00;
        grad[142] = 0.0;
        grad[143] = 0.0;
        grad[144] = 0.0;
        grad[145] = 0.0;
        grad[146] = 0.0;
        grad[147] = 0.0;
        grad[148] = 0.0;
        grad[149] = 0.0;
        grad[150] = 0.0;
        grad[151] = 0.0;
        grad[152] = 0.0;
        grad[153] = 0.0;
        grad[154] = 0.0;
        grad[155] = -1.00;
        grad[156] = 0.0;
        grad[157] = 0.0;
        grad[158] = -7.0/6.0*x[10] + 1.0/6.0*x[11] - 7.0/6.0*x[12] -x[9] - 1.0;
        grad[159] = -7.0/6.0*x[10] + 1.0/6.0*x[11] - 7.0/6.0*x[12] -x[9] - 1.0;
        grad[160] = 0.0;
        grad[161] = 0.0;
        grad[162] = 0.0;
        grad[163] = 0.0;
        grad[164] = 0.0;
        grad[165] = -x[2] -x[3];
        grad[166] = -7.0/6.0*x[2] - 7.0/6.0*x[3];
        grad[167] = 1.0/6.0*x[2] + 1.0/6.0*x[3] + 0.5;
        grad[168] = -7.0/6.0*x[2] - 7.0/6.0*x[3];
        grad[169] = 0.0;
        grad[170] = 0.0;
        grad[171] = -28.0/6.0*x[10] + 2.0/3.0*x[11] - 28.0/6.0*x[12] - 4.0*x[9] - 4.0;
        grad[172] = 0.0;
        grad[173] = 0.0;
        grad[174] = 0.0;
        grad[175] = 0.0;
        grad[176] = 0.0;
        grad[177] = 0.0;
        grad[178] = -4.0*x[2];
        grad[179] = -28.0/6.0*x[2];
        grad[180] = 2.0/3.0*x[2] + 2.0;
        grad[181] = -28.0/6.0*x[2];
        grad[182] = 0.0;
        grad[183] = 0.0;
        grad[184] = 0.0;
        grad[185] = -28.0/6.0*x[10] + 2.0/3.0*x[11] - 28.0/6.0*x[12] - 4.0*x[9] - 4.0;
        grad[186] = 0.0;
        grad[187] = 0.0;
        grad[188] = 0.0;
        grad[189] = 0.0;
        grad[190] = 0.0;
        grad[191] = -4.0*x[3];
        grad[192] = -28.0/6.0*x[3];
        grad[193] = 2.0/3.0*x[3];
        grad[194] = -28.0/6.0*x[3];
        grad[195] = -7.0/6.0*x[10] + 1.0/6.0*x[11] - 7.0/6.0*x[12] -x[9] - 1.0;
        grad[196] = 0.0;
        grad[197] = 0.0;
        grad[198] = 0.0;
        grad[199] = 0.0;
        grad[200] = 0.0;
        grad[201] = 0.0;
        grad[202] = 0.0;
        grad[203] = 0.0;
        grad[204] = 1.0 -x[0];
        grad[205] = -7.0/6.0*x[0];
        grad[206] = 1.0/6.0*x[0];
        grad[207] = -7.0/6.0*x[0];
        grad[208] = 0.0;
        grad[209] = -7.0/6.0*x[10] + 1.0/6.0*x[11] - 7.0/6.0*x[12] -x[9] - 1.0;
        grad[210] = 0.0;
        grad[211] = 0.0;
        grad[212] = 0.0;
        grad[213] = 0.0;
        grad[214] = 0.0;
        grad[215] = 0.0;
        grad[216] = 0.0;
        grad[217] = 1.0 -x[1];
        grad[218] = 0.5 - 7.0/6.0*x[1];
        grad[219] = 1.0/6.0*x[1];
        grad[220] = 0.5 - 7.0/6.0*x[1];
        grad[221] = -7.0/6.0*x[10] + 1.0/6.0*x[11] - 7.0/6.0*x[12] -x[9] - 1.0;
        grad[222] = -7.0/6.0*x[10] + 1.0/6.0*x[11] - 7.0/6.0*x[12] -x[9] - 1.0;
        grad[223] = -28.0/6.0*x[10] + 2.0/3.0*x[11] - 28.0/6.0*x[12] - 4.0*x[9] - 4.0;
        grad[224] = -28.0/6.0*x[10] + 2.0/3.0*x[11] - 28.0/6.0*x[12] - 4.0*x[9] - 4.0;
        grad[225] = 0.0;
        grad[226] = 0.0;
        grad[227] = 0.0;
        grad[228] = 0.0;
        grad[229] = 0.0;
        grad[230] = -x[0] -x[1] - 4.0*x[2] - 4.0*x[3] + 2.0;
        grad[231] = -7.0/6.0*x[0] - 7.0/6.0*x[1] - 28.0/6.0*x[2] - 28.0/6.0*x[3] + 0.5;
        grad[232] = 1.0/6.0*x[0] + 1.0/6.0*x[1] + 2.0/3.0*x[2] + 2.0/3.0*x[3] + 2.0;
        grad[233] = -7.0/6.0*x[0] - 7.0/6.0*x[1] - 28.0/6.0*x[2] - 28.0/6.0*x[3] + 0.5;
    }

    return;
};


/**
    Inequality constraints for fsp
*/
void fsp_igad_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( x[0] + x[1] - 1.0);
    result[1] = ( -x[0]);
    result[2] = ( -x[1]);
    result[3] = ( -0.25*x[0] - 0.25);
    result[4] = ( 0.25*x[0] - 0.75);

    if (grad) {
        grad[0] = 1.00;
        grad[1] = 1.00;
        grad[2] = -1.00;
        grad[3] = 0.0;
        grad[4] = 0.0;
        grad[5] = -1.00;
        grad[6] = -0.250;
        grad[7] = 0.0;
        grad[8] = 0.250;
        grad[9] = 0.0;
    }

    return;
};

/**
    Inequality constraints for spl
*/
void spl_igad_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( 1.0/3.0*x[0]*x[3] + 1.0/3.0*x[0] - 1.0/3.0*x[3] - 2.0/3.0*x[4] - 1.0/3.0);
    result[1] = ( -1.0/3.0*x[0]*x[3] - 1.0/3.0*x[0] - 2.0/3.0*x[5]);
    result[2] = ( -2.0/3.0*x[1]*x[2] - 2.0/3.0*x[1]*x[3] + 2.0/3.0*x[1] + 1.0/3.0*x[3] + 2.0/3.0*x[4] + 2.0/3.0*x[5] + 2.0/3.0*x[6] - 2.0/3.0);
    result[3] = ( 2.0/3.0*x[1]*x[2] + 2.0/3.0*x[1]*x[3] - 2.0/3.0*x[1] - 2.0/3.0*x[6]);
    result[4] = ( 1.0/3.0*x[0]*x[3] + 1.0/3.0*x[0] - 1.0/3.0*x[3] + 1.0/3.0*x[4] - 1.0/3.0);
    result[5] = ( -1.0/3.0*x[0]*x[3] - 1.0/3.0*x[0] + 1.0/3.0*x[5]);
    result[6] = ( -2.0/3.0*x[1]*x[2] - 2.0/3.0*x[1]*x[3] + 2.0/3.0*x[1] + x[2] + 5.0/6.0*x[3] - 1.0/3.0*x[4] - 1.0/3.0*x[5] - 1.0/3.0*x[6] - 2.0/3.0);
    result[7] = ( 2.0/3.0*x[1]*x[2] + 2.0/3.0*x[1]*x[3] - 2.0/3.0*x[1] + 1.0/3.0*x[6]);
    result[8] = ( -x[2]);
    result[9] = ( -0.5*x[3]);

    if (grad) {
        grad[0] = 1.0/3.0*x[3] + 1.0/3.0;
        grad[1] = 0.0;
        grad[2] = 0.0;
        grad[3] = 1.0/3.0*x[0] - 1.0/3.0;
        grad[4] = -2.0/3.0;
        grad[5] = 0.0;
        grad[6] = 0.0;
        grad[7] = -1.0/3.0*x[3] - 1.0/3.0;
        grad[8] = 0.0;
        grad[9] = 0.0;
        grad[10] = -1.0/3.0*x[0];
        grad[11] = 0.0;
        grad[12] = -2.0/3.0;
        grad[13] = 0.0;
        grad[14] = 0.0;
        grad[15] = -2.0/3.0*x[2] - 2.0/3.0*x[3] + 2.0/3.0;
        grad[16] = -2.0/3.0*x[1];
        grad[17] = 1.0/3.0 - 2.0/3.0*x[1];
        grad[18] = 2.0/3.0;
        grad[19] = 2.0/3.0;
        grad[20] = 2.0/3.0;
        grad[21] = 0.0;
        grad[22] = 2.0/3.0*x[2] + 2.0/3.0*x[3] - 2.0/3.0;
        grad[23] = 2.0/3.0*x[1];
        grad[24] = 2.0/3.0*x[1];
        grad[25] = 0.0;
        grad[26] = 0.0;
        grad[27] = -2.0/3.0;
        grad[28] = 1.0/3.0*x[3] + 1.0/3.0;
        grad[29] = 0.0;
        grad[30] = 0.0;
        grad[31] = 1.0/3.0*x[0] - 1.0/3.0;
        grad[32] = 1.0/3.0;
        grad[33] = 0.0;
        grad[34] = 0.0;
        grad[35] = -1.0/3.0*x[3] - 1.0/3.0;
        grad[36] = 0.0;
        grad[37] = 0.0;
        grad[38] = -1.0/3.0*x[0];
        grad[39] = 0.0;
        grad[40] = 1.0/3.0;
        grad[41] = 0.0;
        grad[42] = 0.0;
        grad[43] = -2.0/3.0*x[2] - 2.0/3.0*x[3] + 2.0/3.0;
        grad[44] = 1.0 - 2.0/3.0*x[1];
        grad[45] = 5.0/6.0 - 2.0/3.0*x[1];
        grad[46] = -1.0/3.0;
        grad[47] = -1.0/3.0;
        grad[48] = -1.0/3.0;
        grad[49] = 0.0;
        grad[50] = 2.0/3.0*x[2] + 2.0/3.0*x[3] - 2.0/3.0;
        grad[51] = 2.0/3.0*x[1];
        grad[52] = 2.0/3.0*x[1];
        grad[53] = 0.0;
        grad[54] = 0.0;
        grad[55] = 1.0/3.0;
        grad[56] = 0.0;
        grad[57] = 0.0;
        grad[58] = -1.00;
        grad[59] = 0.0;
        grad[60] = 0.0;
        grad[61] = 0.0;
        grad[62] = 0.0;
        grad[63] = 0.0;
        grad[64] = 0.0;
        grad[65] = 0.0;
        grad[66] = -0.500;
        grad[67] = 0.0;
        grad[68] = 0.0;
        grad[69] = 0.0;
    }

    return;
};

/**
    Inequality constraints for g
*/
void g_igad_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( -x[0]*x[1] + x[0] + x[1] - 1.0);
    result[1] = ( x[0]*x[1] - x[0]);
    result[2] = ( -x[1]);
    result[3] = ( x[2] + x[3] + 2.0*x[4] - 1.0);
    result[4] = ( -x[3]);
    result[5] = ( -x[2]);
    result[6] = ( -x[4]);
    result[7] = ( -x[4]);

    if (grad) {
        grad[0] = 1.0 - x[1];
        grad[1] = 1.0 - x[0];
        grad[2] = 0.0;
        grad[3] = 0.0;
        grad[4] = 0.0;
        grad[5] = x[1] - 1.0;
        grad[6] = x[0];
        grad[7] = 0.0;
        grad[8] = 0.0;
        grad[9] = 0.0;
        grad[10] = 0.0;
        grad[11] = -1.00;
        grad[12] = 0.0;
        grad[13] = 0.0;
        grad[14] = 0.0;
        grad[15] = 0.0;
        grad[16] = 0.0;
        grad[17] = 1.00;
        grad[18] = 1.00;
        grad[19] = 2.00;
        grad[20] = 0.0;
        grad[21] = 0.0;
        grad[22] = 0.0;
        grad[23] = -1.00;
        grad[24] = 0.0;
        grad[25] = 0.0;
        grad[26] = 0.0;
        grad[27] = -1.00;
        grad[28] = 0.0;
        grad[29] = 0.0;
        grad[30] = 0.0;
        grad[31] = 0.0;
        grad[32] = 0.0;
        grad[33] = 0.0;
        grad[34] = -1.00;
        grad[35] = 0.0;
        grad[36] = 0.0;
        grad[37] = 0.0;
        grad[38] = 0.0;
        grad[39] = -1.00;
    }

    return;
};

/**
    Inequality constraints for ol
*/
void ol_igad_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( x[0] - x[2] - 1.0);
    result[1] = ( -x[0] + x[2]);
    result[2] = ( -x[0]*x[1] + x[0] + x[1] + x[2] - 1.0);
    result[3] = ( x[0]*x[1] - x[0] - x[2]);
    result[4] = ( -x[1]);

    if (grad) {
        grad[0] = 1.00;
        grad[1] = 0.0;
        grad[2] = -1.00;
        grad[3] = -1.00;
        grad[4] = 0.0;
        grad[5] = 1.00;
        grad[6] = 1.0 - x[1];
        grad[7] = 1.0 - x[0];
        grad[8] = 1.00;
        grad[9] = x[1] - 1.0;
        grad[10] = x[0];
        grad[11] = -1.00;
        grad[12] = 0.0;
        grad[13] = -1.00;
        grad[14] = 0.0;
    }

    return;
};

/**
    Inequality constraints for opx
*/
void opx_igad_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = (-x[0]*x[1] + x[0]*x[5] -x[0]*x[7] + x[0] -x[1]*x[3] + x[1] + x[3]*x[5] -x[3]*x[7] + x[3] -x[5] + x[7] - 1.0);
    result[1] = (x[0]*x[1] -x[0]*x[5] + x[0]*x[7] -x[0] + x[1]*x[3] -x[3]*x[5] + x[3]*x[7] -x[3]);
    result[2] = (-x[1] + x[4] + 2.0*x[5] + x[6] -x[7]);
    result[3] = (-x[4]);
    result[4] = (-x[6]);
    result[5] = (-x[5]);
    result[6] = (-x[0]*x[2] -x[0]*x[7] + x[0] + x[1]*x[3] + x[2] -x[3]*x[5] + x[3]*x[7] -x[3] + x[7] - 1.0);
    result[7] = (x[0]*x[2] + x[0]*x[7] -x[0] -x[1]*x[3] + x[3]*x[5] -x[3]*x[7] + x[3]);
    result[8] = (-x[2]);
    result[9] = (-x[7]);
    result[10] = (0.5*x[1] - 1.0);
    result[11] = (-0.5*x[1]);

    if (grad) {
        grad[0] = -x[1] + x[5] -x[7] + 1.0;
        grad[1] = -x[0] -x[3] + 1.0;
        grad[2] = 0.0;
        grad[3] = -x[1] + x[5] -x[7] + 1.0;
        grad[4] = 0.0;
        grad[5] = x[0] + x[3] - 1.0;
        grad[6] = 0.0;
        grad[7] = -x[0] -x[3] + 1.0;
        grad[8] = x[1] -x[5] + x[7] - 1.0;
        grad[9] = x[0] + x[3];
        grad[10] = 0.0;
        grad[11] = x[1] -x[5] + x[7] - 1.0;
        grad[12] = 0.0;
        grad[13] = -x[0] -x[3];
        grad[14] = 0.0;
        grad[15] = x[0] + x[3];
        grad[16] = 0.0;
        grad[17] = -1.0;
        grad[18] = 0.0;
        grad[19] = 0.0;
        grad[20] = 1.0;
        grad[21] = 2.0;
        grad[22] = 1.0;
        grad[23] = -1.0;
        grad[24] = 0.0;
        grad[25] = 0.0;
        grad[26] = 0.0;
        grad[27] = 0.0;
        grad[28] = -1.0;
        grad[29] = 0.0;
        grad[30] = 0.0;
        grad[31] = 0.0;
        grad[32] = 0.0;
        grad[33] = 0.0;
        grad[34] = 0.0;
        grad[35] = 0.0;
        grad[36] = 0.0;
        grad[37] = 0.0;
        grad[38] = -1.0;
        grad[39] = 0.0;
        grad[40] = 0.0;
        grad[41] = 0.0;
        grad[42] = 0.0;
        grad[43] = 0.0;
        grad[44] = 0.0;
        grad[45] = -1.0;
        grad[46] = 0.0;
        grad[47] = 0.0;
        grad[48] = -x[2] -x[7] + 1.0;
        grad[49] = x[3];
        grad[50] = 1.0 -x[0];
        grad[51] = x[1] -x[5] + x[7] - 1.0;
        grad[52] = 0.0;
        grad[53] = -x[3];
        grad[54] = 0.0;
        grad[55] = -x[0] + x[3] + 1.0;
        grad[56] = x[2] + x[7] - 1.0;
        grad[57] = -x[3];
        grad[58] = x[0];
        grad[59] = -x[1] + x[5] -x[7] + 1.0;
        grad[60] = 0.0;
        grad[61] = x[3];
        grad[62] = 0.0;
        grad[63] = x[0] -x[3];
        grad[64] = 0.0;
        grad[65] = 0.0;
        grad[66] = -1.0;
        grad[67] = 0.0;
        grad[68] = 0.0;
        grad[69] = 0.0;
        grad[70] = 0.0;
        grad[71] = 0.0;
        grad[72] = 0.0;
        grad[73] = 0.0;
        grad[74] = 0.0;
        grad[75] = 0.0;
        grad[76] = 0.0;
        grad[77] = 0.0;
        grad[78] = 0.0;
        grad[79] = -1.0;
        grad[80] = 0.0;
        grad[81] = 0.50;
        grad[82] = 0.0;
        grad[83] = 0.0;
        grad[84] = 0.0;
        grad[85] = 0.0;
        grad[86] = 0.0;
        grad[87] = 0.0;
        grad[88] = 0.0;
        grad[89] = -0.50;
        grad[90] = 0.0;
        grad[91] = 0.0;
        grad[92] = 0.0;
        grad[93] = 0.0;
        grad[94] = 0.0;
        grad[95] = 0.0;
    }

    return;
};


/**
    Inequality constraints for cpx
*/
void cpx_igad_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = (-x[0]*x[1] -x[0]*x[3] + x[0]*x[7] -x[0]*x[8] + x[0] -x[1]*x[4] + x[1] -x[3]*x[4] + x[3] + x[4]*x[7] -x[4]*x[8] + x[4] -x[7] + x[8] - 1.0);
    result[1] = (x[0]*x[1] + x[0]*x[3] -x[0]*x[7] + x[0]*x[8] -x[0] + x[1]*x[4] + x[3]*x[4] -x[4]*x[7] + x[4]*x[8] -x[4]);
    result[2] = (-x[1] -x[3] + x[5] + x[6] + 2.0*x[7] -x[8]);
    result[3] = (-x[5]);
    result[4] = (-x[6]);
    result[5] = (-x[7]);
    result[6] = (x[0]*x[2] + x[1]*x[4] -x[2] + x[3]*x[4] -x[4]*x[7] + x[4]*x[8] -x[4]);
    result[7] = (-x[0]*x[2] -x[1]*x[4] -x[3]*x[4] + x[4]*x[7] -x[4]*x[8] + x[4]);
    result[8] = (x[2] + x[3] + x[8] - 1.0);
    result[9] = (-x[3]);
    result[10] = (-x[8]);
    result[11] = (0.5*x[1] - 1.0);
    result[12] = (-0.5*x[1]);

    if (grad) {
        grad[0] = -x[1] -x[3] + x[7] -x[8] + 1.0;
        grad[1] = -x[0] -x[4] + 1.0;
        grad[2] = 0.0;
        grad[3] = -x[0] -x[4] + 1.0;
        grad[4] = -x[1] -x[3] + x[7] -x[8] + 1.0;
        grad[5] = 0.0;
        grad[6] = 0.0;
        grad[7] = x[0] + x[4] - 1.0;
        grad[8] = -x[0] -x[4] + 1.0;
        grad[9] = x[1] + x[3] -x[7] + x[8] - 1.0;
        grad[10] = x[0] + x[4];
        grad[11] = 0.0;
        grad[12] = x[0] + x[4];
        grad[13] = x[1] + x[3] -x[7] + x[8] - 1.0;
        grad[14] = 0.0;
        grad[15] = 0.0;
        grad[16] = -x[0] -x[4];
        grad[17] = x[0] + x[4];
        grad[18] = 0.0;
        grad[19] = -1.0;
        grad[20] = 0.0;
        grad[21] = -1.0;
        grad[22] = 0.0;
        grad[23] = 1.0;
        grad[24] = 1.0;
        grad[25] = 2.0;
        grad[26] = -1.0;
        grad[27] = 0.0;
        grad[28] = 0.0;
        grad[29] = 0.0;
        grad[30] = 0.0;
        grad[31] = 0.0;
        grad[32] = -1.0;
        grad[33] = 0.0;
        grad[34] = 0.0;
        grad[35] = 0.0;
        grad[36] = 0.0;
        grad[37] = 0.0;
        grad[38] = 0.0;
        grad[39] = 0.0;
        grad[40] = 0.0;
        grad[41] = 0.0;
        grad[42] = -1.0;
        grad[43] = 0.0;
        grad[44] = 0.0;
        grad[45] = 0.0;
        grad[46] = 0.0;
        grad[47] = 0.0;
        grad[48] = 0.0;
        grad[49] = 0.0;
        grad[50] = 0.0;
        grad[51] = 0.0;
        grad[52] = -1.0;
        grad[53] = 0.0;
        grad[54] = x[2];
        grad[55] = x[4];
        grad[56] = x[0] - 1.0;
        grad[57] = x[4];
        grad[58] = x[1] + x[3] -x[7] + x[8] - 1.0;
        grad[59] = 0.0;
        grad[60] = 0.0;
        grad[61] = -x[4];
        grad[62] = x[4];
        grad[63] = -x[2];
        grad[64] = -x[4];
        grad[65] = -x[0];
        grad[66] = -x[4];
        grad[67] = -x[1] -x[3] + x[7] -x[8] + 1.0;
        grad[68] = 0.0;
        grad[69] = 0.0;
        grad[70] = x[4];
        grad[71] = -x[4];
        grad[72] = 0.0;
        grad[73] = 0.0;
        grad[74] = 1.0;
        grad[75] = 1.0;
        grad[76] = 0.0;
        grad[77] = 0.0;
        grad[78] = 0.0;
        grad[79] = 0.0;
        grad[80] = 1.0;
        grad[81] = 0.0;
        grad[82] = 0.0;
        grad[83] = 0.0;
        grad[84] = -1.0;
        grad[85] = 0.0;
        grad[86] = 0.0;
        grad[87] = 0.0;
        grad[88] = 0.0;
        grad[89] = 0.0;
        grad[90] = 0.0;
        grad[91] = 0.0;
        grad[92] = 0.0;
        grad[93] = 0.0;
        grad[94] = 0.0;
        grad[95] = 0.0;
        grad[96] = 0.0;
        grad[97] = 0.0;
        grad[98] = -1.0;
        grad[99] = 0.0;
        grad[100] = 0.50;
        grad[101] = 0.0;
        grad[102] = 0.0;
        grad[103] = 0.0;
        grad[104] = 0.0;
        grad[105] = 0.0;
        grad[106] = 0.0;
        grad[107] = 0.0;
        grad[108] = 0.0;
        grad[109] = -0.50;
        grad[110] = 0.0;
        grad[111] = 0.0;
        grad[112] = 0.0;
        grad[113] = 0.0;
        grad[114] = 0.0;
        grad[115] = 0.0;
        grad[116] = 0.0;
    }

    return;
};
/**
    Inequality constraints for ilm
*/
void ilm_igad_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( 0.5*x[0]*x[1] - 0.5*x[0] - 0.5*x[2]);
    result[1] = ( -0.5*x[0] + 0.5*x[3]);
    result[2] = ( x[0] - 1.0);
    result[3] = ( -0.5*x[0]*x[1] + 0.5*x[2] - 0.5*x[3]);
    result[4] = ( 0.5*x[0]*x[1] - 0.5*x[0] + 0.5*x[2]);
    result[5] = ( -0.5*x[0] - 0.5*x[3]);
    result[6] = ( x[0] - 1.0);
    result[7] = ( -0.5*x[0]*x[1] - 0.5*x[2] + 0.5*x[3]);

    if (grad) {
        grad[0] = 0.5*x[1] - 0.5;
        grad[1] = 0.5*x[0];
        grad[2] = -0.500;
        grad[3] = 0.0;
        grad[4] = -0.500;
        grad[5] = 0.0;
        grad[6] = 0.0;
        grad[7] = 0.500;
        grad[8] = 1.00;
        grad[9] = 0.0;
        grad[10] = 0.0;
        grad[11] = 0.0;
        grad[12] = -0.5*x[1];
        grad[13] = -0.5*x[0];
        grad[14] = 0.500;
        grad[15] = -0.500;
        grad[16] = 0.5*x[1] - 0.5;
        grad[17] = 0.5*x[0];
        grad[18] = 0.500;
        grad[19] = 0.0;
        grad[20] = -0.500;
        grad[21] = 0.0;
        grad[22] = 0.0;
        grad[23] = -0.500;
        grad[24] = 1.00;
        grad[25] = 0.0;
        grad[26] = 0.0;
        grad[27] = 0.0;
        grad[28] = -0.5*x[1];
        grad[29] = -0.5*x[0];
        grad[30] = -0.500;
        grad[31] = 0.500;
    }

    return;
};

/**
    Inequality constraints for nph
*/
void nph_igad_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( -0.25*x[0]*x[1] - 0.75*x[1]*x[4] + x[1] - 0.25*x[2] + x[4] - 1.0);
    result[1] = ( 0.25*x[0]*x[1] + 0.75*x[1]*x[4] - x[1] + 0.25*x[2]);
    result[2] = ( -x[4]);
    result[3] = ( -0.25*x[0]*x[1] + x[0] - 0.75*x[1]*x[4] + x[1] + 0.75*x[2] - 1.0);
    result[4] = ( 0.25*x[0]*x[1] + 0.75*x[1]*x[4] - x[1] - 0.75*x[2]);
    result[5] = ( -x[0]);
    result[6] = ( 0.25*x[0] + x[3] - 0.75*x[4] - 1.0);
    result[7] = ( -0.25*x[0] + 0.75*x[4]);
    result[8] = ( -x[3]);

    if (grad) {
        grad[0] = -0.25*x[1];
        grad[1] = -0.25*x[0] - 0.75*x[4] + 1.0;
        grad[2] = -0.250;
        grad[3] = 0.0;
        grad[4] = 1.0 - 0.75*x[1];
        grad[5] = 0.25*x[1];
        grad[6] = 0.25*x[0] + 0.75*x[4] - 1.0;
        grad[7] = 0.250;
        grad[8] = 0.0;
        grad[9] = 0.75*x[1];
        grad[10] = 0.0;
        grad[11] = 0.0;
        grad[12] = 0.0;
        grad[13] = 0.0;
        grad[14] = -1.00;
        grad[15] = 1.0 - 0.25*x[1];
        grad[16] = -0.25*x[0] - 0.75*x[4] + 1.0;
        grad[17] = 0.750;
        grad[18] = 0.0;
        grad[19] = -0.75*x[1];
        grad[20] = 0.25*x[1];
        grad[21] = 0.25*x[0] + 0.75*x[4] - 1.0;
        grad[22] = -0.750;
        grad[23] = 0.0;
        grad[24] = 0.75*x[1];
        grad[25] = -1.00;
        grad[26] = 0.0;
        grad[27] = 0.0;
        grad[28] = 0.0;
        grad[29] = 0.0;
        grad[30] = 0.250;
        grad[31] = 0.0;
        grad[32] = 0.0;
        grad[33] = 1.00;
        grad[34] = -0.750;
        grad[35] = -0.250;
        grad[36] = 0.0;
        grad[37] = 0.0;
        grad[38] = 0.0;
        grad[39] = 0.750;
        grad[40] = 0.0;
        grad[41] = 0.0;
        grad[42] = 0.0;
        grad[43] = -1.00;
        grad[44] = 0.0;
    }

    return;
};

/**
    Inequality constraints for lct
*/
void lct_igad_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( -x[0]);
    result[1] = ( x[0] - 1.0);

    if (grad) {
        grad[0] = -1.00;
        grad[1] = 1.00;
    }

    return;
};

/**
    Inequality constraints for kals
*/
void kals_igad_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( -x[0]);
    result[1] = ( x[0] - 1.0);

    if (grad) {
        grad[0] = -1.00;
        grad[1] = 1.00;
    }

    return;
};

/**
    Inequality constraints for mel
*/
void mel_igad_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( -x[1]);
    result[1] = ( x[1] - 1.0);
    result[2] = ( -x[0]*x[2] - x[0]*x[3] + x[0] + x[2] + x[3] - 1.0);
    result[3] = ( x[0]*x[2] + x[0]*x[3] - x[0]);
    result[4] = ( -x[2]);
    result[5] = ( -x[3]);
    result[6] = ( x[1] - 0.5*x[2] - 0.5*x[3]);
    result[7] = ( -x[1] + 0.5*x[2] + 0.5*x[3] - 1.0);

    if (grad) {
        grad[0] = 0.0;
        grad[1] = -1.00;
        grad[2] = 0.0;
        grad[3] = 0.0;
        grad[4] = 0.0;
        grad[5] = 1.00;
        grad[6] = 0.0;
        grad[7] = 0.0;
        grad[8] = -x[2] - x[3] + 1.0;
        grad[9] = 0.0;
        grad[10] = 1.0 - x[0];
        grad[11] = 1.0 - x[0];
        grad[12] = x[2] + x[3] - 1.0;
        grad[13] = 0.0;
        grad[14] = x[0];
        grad[15] = x[0];
        grad[16] = 0.0;
        grad[17] = 0.0;
        grad[18] = -1.00;
        grad[19] = 0.0;
        grad[20] = 0.0;
        grad[21] = 0.0;
        grad[22] = 0.0;
        grad[23] = -1.00;
        grad[24] = 0.0;
        grad[25] = 1.00;
        grad[26] = -0.500;
        grad[27] = -0.500;
        grad[28] = 0.0;
        grad[29] = -1.00;
        grad[30] = 0.500;
        grad[31] = 0.500;
    }

    return;
};

//---------------------------------------------------------------------------
//---------------------------------Evans&Frost,2021--------------------------
//---------------------------------------------------------------------------


/**
    Inequality constraints for fluid
*/
void fluid_um_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( -x[0]);
    result[1] = ( x[0] - 1.0);

    if (grad) {
        grad[0] = -1.0;
        grad[1] = 1.0;
    }

    return;
};

/**
    Inequality constraints for ol
*/
void ol_um_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( x[0] - 1.0);
    result[1] = ( -x[0]);

    if (grad) {
        grad[0] = 1.0;
        grad[1] = -1.0;
    }

    return;
};

/**
    Inequality constraints for br
*/
void br_um_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( x[0] - 1.0);
    result[1] = ( -x[0]);

    if (grad) {
        grad[0] = 1.0;
        grad[1] = -1.0;
    }

    return;
};

/**
    Inequality constraints for ch
*/
void ch_um_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( x[0] - 1.0);
    result[1] = ( -x[0]);

    if (grad) {
        grad[0] = 1.0;
        grad[1] = -1.0;
    }

    return;
};

/**
    Inequality constraints for atg
*/
void atg_um_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( -x[0]*x[1] -x[0]*x[2] + x[0] + x[1]*x[3] + x[1] + x[2]*x[3] + x[2] -x[3] - 1.0);
    result[1] = ( x[0]*x[1] + x[0]*x[2] -x[0] -x[1]*x[3] -x[2]*x[3] + x[3]);
    result[2] = ( -x[2]);
    result[3] = ( -x[1]);
    result[4] = ( x[0] - 0.5*x[1]*x[3] - 0.5*x[2]*x[3] + 0.5*x[3] - 1.0);
    result[5] = ( -x[0] + 0.5*x[1]*x[3] + 0.5*x[2]*x[3] - 0.5*x[3]);
    result[6] = ( 0.5*x[1] + 0.5*x[2] - 1.0);
    result[7] = ( -0.5*x[1] - 0.5*x[2]);

    if (grad) {
        grad[0] = -x[1] -x[2] + 1.0;
        grad[1] = -x[0] + x[3] + 1.0;
        grad[2] = -x[0] + x[3] + 1.0;
        grad[3] = x[1] + x[2] - 1.0;
        grad[4] = x[1] + x[2] - 1.0;
        grad[5] = x[0] -x[3];
        grad[6] = x[0] -x[3];
        grad[7] = -x[1] -x[2] + 1.0;
        grad[8] = 0.0;
        grad[9] = 0.0;
        grad[10] = -1.0;
        grad[11] = 0.0;
        grad[12] = 0.0;
        grad[13] = -1.0;
        grad[14] = 0.0;
        grad[15] = 0.0;
        grad[16] = 1.0;
        grad[17] = -0.5*x[3];
        grad[18] = -0.5*x[3];
        grad[19] = -0.5*x[1] - 0.5*x[2] + 0.5;
        grad[20] = -1.0;
        grad[21] = 0.5*x[3];
        grad[22] = 0.5*x[3];
        grad[23] = 0.5*x[1] + 0.5*x[2] - 0.5;
        grad[24] = 0.0;
        grad[25] = 0.50;
        grad[26] = 0.50;
        grad[27] = 0.0;
        grad[28] = 0.0;
        grad[29] = -0.50;
        grad[30] = -0.50;
        grad[31] = 0.0;
    }

    return;
};

/**
    Inequality constraints for g
*/
void g_um_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( x[0] - 1.0);
    result[1] = ( -x[0]);

    if (grad) {
        grad[0] = 1.0;
        grad[1] = -1.0;
    }

    return;
};

/**
    Inequality constraints for ta
*/
void ta_um_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( -x[0]*x[3] + x[0] + x[3]*x[4] + x[3] -x[4] - 1.0);
    result[1] = ( x[0]*x[3] -x[0] -x[3]*x[4] + x[4]);
    result[2] = ( -x[3]);
    result[3] = ( -x[0]*x[1] -x[0]*x[2] + x[0] + x[1] + x[2] - 0.5*x[3]*x[4] + 0.5*x[4] - 1.0);
    result[4] = ( x[0]*x[1] + x[0]*x[2] -x[0] + 0.5*x[3]*x[4] - 0.5*x[4]);
    result[5] = ( -x[2]);
    result[6] = ( -x[1]);
    result[7] = ( x[1] + x[2] -x[3] - 1.0);
    result[8] = ( -x[1] -x[2] + x[3]);

    if (grad) {
        grad[0] = 1.0 -x[3];
        grad[1] = 0.0;
        grad[2] = 0.0;
        grad[3] = -x[0] + x[4] + 1.0;
        grad[4] = x[3] - 1.0;
        grad[5] = x[3] - 1.0;
        grad[6] = 0.0;
        grad[7] = 0.0;
        grad[8] = x[0] -x[4];
        grad[9] = 1.0 -x[3];
        grad[10] = 0.0;
        grad[11] = 0.0;
        grad[12] = 0.0;
        grad[13] = -1.0;
        grad[14] = 0.0;
        grad[15] = -x[1] -x[2] + 1.0;
        grad[16] = 1.0 -x[0];
        grad[17] = 1.0 -x[0];
        grad[18] = -0.5*x[4];
        grad[19] = 0.5 - 0.5*x[3];
        grad[20] = x[1] + x[2] - 1.0;
        grad[21] = x[0];
        grad[22] = x[0];
        grad[23] = 0.5*x[4];
        grad[24] = 0.5*x[3] - 0.5;
        grad[25] = 0.0;
        grad[26] = 0.0;
        grad[27] = -1.0;
        grad[28] = 0.0;
        grad[29] = 0.0;
        grad[30] = 0.0;
        grad[31] = -1.0;
        grad[32] = 0.0;
        grad[33] = 0.0;
        grad[34] = 0.0;
        grad[35] = 0.0;
        grad[36] = 1.0;
        grad[37] = 1.0;
        grad[38] = -1.0;
        grad[39] = 0.0;
        grad[40] = 0.0;
        grad[41] = -1.0;
        grad[42] = -1.0;
        grad[43] = 1.0;
        grad[44] = 0.0;
    }

    return;
};

/**
    Inequality constraints for chl
*/
void chl_um_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( -x[0]*x[1] + x[0]*x[3] + x[0] + x[1]*x[4] + x[1] -x[3]*x[4] -x[3] -x[4] - 1.0);
    result[1] = ( x[0]*x[1] -x[0]*x[3] -x[0] -x[1]*x[4] + x[3]*x[4] + x[4]);
    result[2] = ( -x[1] + x[3]);
    result[3] = ( x[0] - 0.25*x[1]*x[4] - 0.25*x[1]*x[5] - 0.25*x[2]*x[5] + 0.25*x[3]*x[4] - 0.25*x[3]*x[5] + 0.25*x[4] + 0.25*x[5] - 1.0);
    result[4] = ( -x[0] + 0.25*x[1]*x[4] + 0.25*x[1]*x[5] + 0.25*x[2]*x[5] - 0.25*x[3]*x[4] + 0.25*x[3]*x[5] - 0.25*x[4] - 0.25*x[5]);
    result[5] = ( -x[0]*x[1] -x[0]*x[2] -x[0]*x[3] + x[0] + x[1]*x[5] + x[1] + x[2]*x[5] + x[2] + x[3]*x[5] + x[3] -x[5] - 1.0);
    result[6] = ( x[0]*x[1] + x[0]*x[2] + x[0]*x[3] -x[0] -x[1]*x[5] -x[2]*x[5] -x[3]*x[5] + x[5]);
    result[7] = ( -x[2]);
    result[8] = ( -x[1] -x[3]);
    result[9] = ( x[1] + 0.5*x[2] - 1.0);
    result[10] = ( -x[1] - 0.5*x[2]);

    if (grad) {
        grad[0] = -x[1] + x[3] + 1.0;
        grad[1] = -x[0] + x[4] + 1.0;
        grad[2] = 0.0;
        grad[3] = x[0] -x[4] - 1.0;
        grad[4] = x[1] -x[3] - 1.0;
        grad[5] = 0.0;
        grad[6] = x[1] -x[3] - 1.0;
        grad[7] = x[0] -x[4];
        grad[8] = 0.0;
        grad[9] = -x[0] + x[4];
        grad[10] = -x[1] + x[3] + 1.0;
        grad[11] = 0.0;
        grad[12] = 0.0;
        grad[13] = -1.0;
        grad[14] = 0.0;
        grad[15] = 1.0;
        grad[16] = 0.0;
        grad[17] = 0.0;
        grad[18] = 1.0;
        grad[19] = -0.25*x[4] - 0.25*x[5];
        grad[20] = -0.25*x[5];
        grad[21] = 0.25*x[4] - 0.25*x[5];
        grad[22] = -0.25*x[1] + 0.25*x[3] + 0.25;
        grad[23] = -0.25*x[1] - 0.25*x[2] - 0.25*x[3] + 0.25;
        grad[24] = -1.0;
        grad[25] = 0.25*x[4] + 0.25*x[5];
        grad[26] = 0.25*x[5];
        grad[27] = -0.25*x[4] + 0.25*x[5];
        grad[28] = 0.25*x[1] - 0.25*x[3] - 0.25;
        grad[29] = 0.25*x[1] + 0.25*x[2] + 0.25*x[3] - 0.25;
        grad[30] = -x[1] -x[2] -x[3] + 1.0;
        grad[31] = -x[0] + x[5] + 1.0;
        grad[32] = -x[0] + x[5] + 1.0;
        grad[33] = -x[0] + x[5] + 1.0;
        grad[34] = 0.0;
        grad[35] = x[1] + x[2] + x[3] - 1.0;
        grad[36] = x[1] + x[2] + x[3] - 1.0;
        grad[37] = x[0] -x[5];
        grad[38] = x[0] -x[5];
        grad[39] = x[0] -x[5];
        grad[40] = 0.0;
        grad[41] = -x[1] -x[2] -x[3] + 1.0;
        grad[42] = 0.0;
        grad[43] = 0.0;
        grad[44] = -1.0;
        grad[45] = 0.0;
        grad[46] = 0.0;
        grad[47] = 0.0;
        grad[48] = 0.0;
        grad[49] = -1.0;
        grad[50] = 0.0;
        grad[51] = -1.0;
        grad[52] = 0.0;
        grad[53] = 0.0;
        grad[54] = 0.0;
        grad[55] = 1.0;
        grad[56] = 0.50;
        grad[57] = 0.0;
        grad[58] = 0.0;
        grad[59] = 0.0;
        grad[60] = 0.0;
        grad[61] = -1.0;
        grad[62] = -0.50;
        grad[63] = 0.0;
        grad[64] = 0.0;
        grad[65] = 0.0;
    }

    return;
};

/**
    Inequality constraints for anth
*/
void anth_um_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( x[0] - x[1]*x[3] + 1.5*x[2] + x[3] - 1.0);
    result[1] = ( -x[0] + x[1]*x[3] - 1.5*x[2] - x[3]);
    result[2] = ( x[0] - x[2] - 1.0);
    result[3] = ( -x[0] + x[2]);
    result[4] = ( -x[1]);
    result[5] = ( -x[0]*x[1] + x[0] + x[1]*x[3] + x[1] - x[3] - 1.0);
    result[6] = ( x[0]*x[1] - x[0] - x[1]*x[3] + x[3]);
    result[7] = ( -0.5*x[1]);
    result[8] = ( 0.5*x[1] - 1.0);

    if (grad) {
        grad[0] = 1.0;
        grad[1] = -x[3];
        grad[2] = 1.5;
        grad[3] = 1.0 - x[1];
        grad[4] = -1.0;
        grad[5] = x[3];
        grad[6] = -1.5;
        grad[7] = x[1] - 1.0;
        grad[8] = 1.0;
        grad[9] = 0.0;
        grad[10] = -1.0;
        grad[11] = 0.0;
        grad[12] = -1.0;
        grad[13] = 0.0;
        grad[14] = 1.0;
        grad[15] = 0.0;
        grad[16] = 0.0;
        grad[17] = -1.0;
        grad[18] = 0.0;
        grad[19] = 0.0;
        grad[20] = 1.0 - x[1];
        grad[21] = -x[0] + x[3] + 1.0;
        grad[22] = 0.0;
        grad[23] = x[1] - 1.0;
        grad[24] = x[1] - 1.0;
        grad[25] = x[0] - x[3];
        grad[26] = 0.0;
        grad[27] = 1.0 - x[1];
        grad[28] = 0.0;
        grad[29] = -0.50;
        grad[30] = 0.0;
        grad[31] = 0.0;
        grad[32] = 0.0;
        grad[33] = 0.50;
        grad[34] = 0.0;
        grad[35] = 0.0;
    }

    return;
};


/**
    Inequality constraints for spi
*/
void spi_um_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( -x[1]);
    result[1] = ( x[1] - 1.0);
    result[2] = ( x[0] - 1.0);
    result[3] = ( -x[0]);

    if (grad) {
        grad[0] = 0.0;
        grad[1] = -1.0;
        grad[2] = 0.0;
        grad[3] = 1.0;
        grad[4] = 1.0;
        grad[5] = 0.0;
        grad[6] = -1.0;
        grad[7] = 0.0;
    }

    return;
};

/**
    Inequality constraints for opx
*/
void opx_um_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( -x[0]*x[1] -x[0]*x[2] + x[0] + x[1] + x[2] - 0.5*x[3] - 1.0);
    result[1] = ( x[0]*x[1] + x[0]*x[2] -x[0] + 0.5*x[3]);
    result[2] = ( -x[2]);
    result[3] = ( -x[1]);
    result[4] = ( x[0] + 0.5*x[3] - 1.0);
    result[5] = ( -x[0] - 0.5*x[3]);
    result[6] = ( -0.5*x[1] - 0.5*x[2]);
    result[7] = ( 0.5*x[1] + 0.5*x[2] - 1.0);

    if (grad) {
        grad[0] = -x[1] -x[2] + 1.0;
        grad[1] = 1.0 -x[0];
        grad[2] = 1.0 -x[0];
        grad[3] = -0.50;
        grad[4] = x[1] + x[2] - 1.0;
        grad[5] = x[0];
        grad[6] = x[0];
        grad[7] = 0.50;
        grad[8] = 0.0;
        grad[9] = 0.0;
        grad[10] = -1.0;
        grad[11] = 0.0;
        grad[12] = 0.0;
        grad[13] = -1.0;
        grad[14] = 0.0;
        grad[15] = 0.0;
        grad[16] = 1.0;
        grad[17] = 0.0;
        grad[18] = 0.0;
        grad[19] = 0.50;
        grad[20] = -1.0;
        grad[21] = 0.0;
        grad[22] = 0.0;
        grad[23] = -0.50;
        grad[24] = 0.0;
        grad[25] = -0.50;
        grad[26] = -0.50;
        grad[27] = 0.0;
        grad[28] = 0.0;
        grad[29] = 0.50;
        grad[30] = 0.50;
        grad[31] = 0.0;
    }

    return;
};

/**
    Inequality constraints for po
*/
void po_um_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( x[0] - 1.0);
    result[1] = ( -x[0]);

    if (grad) {
        grad[0] = 1.0;
        grad[1] = -1.0;
    }

    return;
};

/**
    Inequality constraints for pl4tr
*/
void pl4tr_ume_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = (x[0] - 1.0);
    result[1] = (-x[0]);
    result[2] = (-0.25*x[0] - 0.25);
    result[3] = (0.25*x[0] - 0.75);

    if (grad) {
        grad[0] = 1.0;
        grad[1] = -1.0;
        grad[2] = -0.25;
        grad[3] = 0.25;
    }

    return;
};

/**
    Inequality constraints for hb
*/
void hb_ume_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = (x[3] - 1.0);
    result[1] = (-x[3]);
    result[2] = (x[0] -x[6] - 1.0);
    result[3] = (-x[0] + x[6]);
    result[4] = (-x[0]*x[1] -x[0]*x[5] + x[0] + x[1]*x[7] + x[1] + x[5]*x[7] + x[5] -x[7] - 1.0);
    result[5] = (x[0]*x[1] + x[0]*x[5] -x[0] -x[1]*x[7] -x[5]*x[7] + x[7]);
    result[6] = (-x[1]);
    result[7] = (-x[5]);
    result[8] = (-x[4]);
    result[9] = (-x[0]*x[2] -x[0]*x[4] + x[0] -x[1]*x[7] + x[2] + x[4] -x[5]*x[7] + 1.5*x[6] + x[7] - 1.0);
    result[10] = (x[0]*x[2] + x[0]*x[4] -x[0] + x[1]*x[7] + x[5]*x[7] - 1.5*x[6] -x[7]);
    result[11] = (-x[2]);
    result[12] = (0.5*x[1] - 0.5*x[2] + 0.25*x[3] + 0.5*x[5] - 1.0);
    result[13] = (-0.5*x[1] + 0.5*x[2] - 0.25*x[3] - 0.5*x[5]);

    if (grad) {
        grad[0] = 0.0;
        grad[1] = 0.0;
        grad[2] = 0.0;
        grad[3] = 1.0;
        grad[4] = 0.0;
        grad[5] = 0.0;
        grad[6] = 0.0;
        grad[7] = 0.0;
        grad[8] = 0.0;
        grad[9] = 0.0;
        grad[10] = 0.0;
        grad[11] = -1.0;
        grad[12] = 0.0;
        grad[13] = 0.0;
        grad[14] = 0.0;
        grad[15] = 0.0;
        grad[16] = 1.0;
        grad[17] = 0.0;
        grad[18] = 0.0;
        grad[19] = 0.0;
        grad[20] = 0.0;
        grad[21] = 0.0;
        grad[22] = -1.0;
        grad[23] = 0.0;
        grad[24] = -1.0;
        grad[25] = 0.0;
        grad[26] = 0.0;
        grad[27] = 0.0;
        grad[28] = 0.0;
        grad[29] = 0.0;
        grad[30] = 1.0;
        grad[31] = 0.0;
        grad[32] = -x[1] -x[5] + 1.0;
        grad[33] = -x[0] + x[7] + 1.0;
        grad[34] = 0.0;
        grad[35] = 0.0;
        grad[36] = 0.0;
        grad[37] = -x[0] + x[7] + 1.0;
        grad[38] = 0.0;
        grad[39] = x[1] + x[5] - 1.0;
        grad[40] = x[1] + x[5] - 1.0;
        grad[41] = x[0] -x[7];
        grad[42] = 0.0;
        grad[43] = 0.0;
        grad[44] = 0.0;
        grad[45] = x[0] -x[7];
        grad[46] = 0.0;
        grad[47] = -x[1] -x[5] + 1.0;
        grad[48] = 0.0;
        grad[49] = -1.0;
        grad[50] = 0.0;
        grad[51] = 0.0;
        grad[52] = 0.0;
        grad[53] = 0.0;
        grad[54] = 0.0;
        grad[55] = 0.0;
        grad[56] = 0.0;
        grad[57] = 0.0;
        grad[58] = 0.0;
        grad[59] = 0.0;
        grad[60] = 0.0;
        grad[61] = -1.0;
        grad[62] = 0.0;
        grad[63] = 0.0;
        grad[64] = 0.0;
        grad[65] = 0.0;
        grad[66] = 0.0;
        grad[67] = 0.0;
        grad[68] = -1.0;
        grad[69] = 0.0;
        grad[70] = 0.0;
        grad[71] = 0.0;
        grad[72] = -x[2] -x[4] + 1.0;
        grad[73] = -x[7];
        grad[74] = 1.0 -x[0];
        grad[75] = 0.0;
        grad[76] = 1.0 -x[0];
        grad[77] = -x[7];
        grad[78] = 1.5;
        grad[79] = -x[1] -x[5] + 1.0;
        grad[80] = x[2] + x[4] - 1.0;
        grad[81] = x[7];
        grad[82] = x[0];
        grad[83] = 0.0;
        grad[84] = x[0];
        grad[85] = x[7];
        grad[86] = -1.5;
        grad[87] = x[1] + x[5] - 1.0;
        grad[88] = 0.0;
        grad[89] = 0.0;
        grad[90] = -1.0;
        grad[91] = 0.0;
        grad[92] = 0.0;
        grad[93] = 0.0;
        grad[94] = 0.0;
        grad[95] = 0.0;
        grad[96] = 0.0;
        grad[97] = 0.50;
        grad[98] = -0.50;
        grad[99] = 0.25;
        grad[100] = 0.0;
        grad[101] = 0.50;
        grad[102] = 0.0;
        grad[103] = 0.0;
        grad[104] = 0.0;
        grad[105] = -0.50;
        grad[106] = 0.50;
        grad[107] = -0.25;
        grad[108] = 0.0;
        grad[109] = -0.50;
        grad[110] = 0.0;
        grad[111] = 0.0;
    }

    return;
};

/**
    Inequality constraints for aug
*/
void aug_ume_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = (-x[0]*x[1] -x[0]*x[4] + x[0] + x[1] + 0.5*x[3]*x[5] + 0.5*x[4]*x[5] + x[4] - 0.5*x[5] - 1.0);
    result[1] = (x[0]*x[1] + x[0]*x[4] -x[0] - 0.5*x[3]*x[5] - 0.5*x[4]*x[5] + 0.5*x[5]);
    result[2] = (-x[1] + x[2] -x[4]);
    result[3] = (-x[2]);
    result[4] = (-x[0]*x[3] -x[0]*x[4] + x[0] - 0.5*x[3]*x[5] + x[3] - 0.5*x[4]*x[5] + x[4] + 0.5*x[5] - 1.0);
    result[5] = (x[0]*x[3] + x[0]*x[4] -x[0] + 0.5*x[3]*x[5] + 0.5*x[4]*x[5] - 0.5*x[5]);
    result[6] = (-x[3]);
    result[7] = (-x[4]);
    result[8] = (0.5*x[1] - 0.5*x[6] - 1.0);
    result[9] = (-0.5*x[1] + 0.5*x[6]);
    result[10] = (0.5*x[1] + 0.5*x[6] - 1.0);
    result[11] = (-0.5*x[1] - 0.5*x[6]);

    if (grad) {
        grad[0] = -x[1] -x[4] + 1.0;
        grad[1] = 1.0 -x[0];
        grad[2] = 0.0;
        grad[3] = 0.5*x[5];
        grad[4] = -x[0] + 0.5*x[5] + 1.0;
        grad[5] = 0.5*x[3] + 0.5*x[4] - 0.5;
        grad[6] = 0.0;
        grad[7] = x[1] + x[4] - 1.0;
        grad[8] = x[0];
        grad[9] = 0.0;
        grad[10] = -0.5*x[5];
        grad[11] = x[0] - 0.5*x[5];
        grad[12] = -0.5*x[3] - 0.5*x[4] + 0.5;
        grad[13] = 0.0;
        grad[14] = 0.0;
        grad[15] = -1.0;
        grad[16] = 1.0;
        grad[17] = 0.0;
        grad[18] = -1.0;
        grad[19] = 0.0;
        grad[20] = 0.0;
        grad[21] = 0.0;
        grad[22] = 0.0;
        grad[23] = -1.0;
        grad[24] = 0.0;
        grad[25] = 0.0;
        grad[26] = 0.0;
        grad[27] = 0.0;
        grad[28] = -x[3] -x[4] + 1.0;
        grad[29] = 0.0;
        grad[30] = 0.0;
        grad[31] = -x[0] - 0.5*x[5] + 1.0;
        grad[32] = -x[0] - 0.5*x[5] + 1.0;
        grad[33] = -0.5*x[3] - 0.5*x[4] + 0.5;
        grad[34] = 0.0;
        grad[35] = x[3] + x[4] - 1.0;
        grad[36] = 0.0;
        grad[37] = 0.0;
        grad[38] = x[0] + 0.5*x[5];
        grad[39] = x[0] + 0.5*x[5];
        grad[40] = 0.5*x[3] + 0.5*x[4] - 0.5;
        grad[41] = 0.0;
        grad[42] = 0.0;
        grad[43] = 0.0;
        grad[44] = 0.0;
        grad[45] = -1.0;
        grad[46] = 0.0;
        grad[47] = 0.0;
        grad[48] = 0.0;
        grad[49] = 0.0;
        grad[50] = 0.0;
        grad[51] = 0.0;
        grad[52] = 0.0;
        grad[53] = -1.0;
        grad[54] = 0.0;
        grad[55] = 0.0;
        grad[56] = 0.0;
        grad[57] = 0.50;
        grad[58] = 0.0;
        grad[59] = 0.0;
        grad[60] = 0.0;
        grad[61] = 0.0;
        grad[62] = -0.50;
        grad[63] = 0.0;
        grad[64] = -0.50;
        grad[65] = 0.0;
        grad[66] = 0.0;
        grad[67] = 0.0;
        grad[68] = 0.0;
        grad[69] = 0.50;
        grad[70] = 0.0;
        grad[71] = 0.50;
        grad[72] = 0.0;
        grad[73] = 0.0;
        grad[74] = 0.0;
        grad[75] = 0.0;
        grad[76] = 0.50;
        grad[77] = 0.0;
        grad[78] = -0.50;
        grad[79] = 0.0;
        grad[80] = 0.0;
        grad[81] = 0.0;
        grad[82] = 0.0;
        grad[83] = -0.50;
    }

    return;
};


SS_ref NLopt_opt_mp_st_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_mp_st, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, st_mp_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};

SS_ref NLopt_opt_mp_sp_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_mp_sp, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, sp_mp_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};

SS_ref NLopt_opt_mp_sa_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_mp_sa, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, sa_mp_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};

SS_ref NLopt_opt_mp_fsp_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_mp_fsp, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, fsp_mp_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};

SS_ref NLopt_opt_mp_opx_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_mp_opx, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, opx_mp_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};

SS_ref NLopt_opt_mp_mu_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_mp_mu, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, mu_mp_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};

SS_ref NLopt_opt_mp_mt_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_mp_mt, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, mt_mp_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};

SS_ref NLopt_opt_mp_ma_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_mp_ma, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, ma_mp_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};
SS_ref NLopt_opt_mp_ilm_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_mp_ilm, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, ilm_mp_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    
    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};
SS_ref NLopt_opt_mp_ilmm_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_mp_ilmm, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, ilmm_mp_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};

SS_ref NLopt_opt_mp_g_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_mp_g, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, g_mp_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};

SS_ref NLopt_opt_mp_ep_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_mp_ep, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, ep_mp_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};

SS_ref NLopt_opt_mp_ctd_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_mp_ctd, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, ctd_mp_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};

SS_ref NLopt_opt_mp_chl_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_mp_chl, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, chl_mp_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};

SS_ref NLopt_opt_mp_cd_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_mp_cd, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, cd_mp_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};

SS_ref NLopt_opt_mp_bi_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_mp_bi, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, bi_mp_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};


SS_ref NLopt_opt_mp_liq_function(global_variable gv, SS_ref SS_ref_db){

    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_mp_liq, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, liq_mp_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};

SS_ref NLopt_opt_ig_fper_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_CCSAQ, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_ig_fper, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, fper_ig_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    
    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};

SS_ref NLopt_opt_ig_bi_function(global_variable gv, SS_ref SS_ref_db){

   int    n_em     = SS_ref_db.n_em;
   unsigned int n  = SS_ref_db.n_xeos;
   unsigned int m  = SS_ref_db.n_sf;
    
   double *x  = SS_ref_db.iguess; 

   for (int i = 0; i < (SS_ref_db.n_xeos); i++){
      SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
      SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
   }
   
   SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n));
   nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
   nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
   nlopt_set_min_objective(SS_ref_db.opt, obj_ig_bi, &SS_ref_db);
   nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, bi_ig_c, NULL, NULL);
   nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
   nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
   
   double minf;
   SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

   /* Send back needed local solution parameters */
   for (int i = 0; i < SS_ref_db.n_xeos; i++){
      SS_ref_db.xeos[i] = x[i];
   }
   SS_ref_db.df   = minf;
   nlopt_destroy(SS_ref_db.opt);

   return SS_ref_db;
};

SS_ref NLopt_opt_ig_cd_function(global_variable gv, SS_ref SS_ref_db){

	int    n_em     = SS_ref_db.n_em;
 	unsigned int n  = SS_ref_db.n_xeos;	
 	unsigned int m  = SS_ref_db.n_sf;	
   
	double *x  = SS_ref_db.iguess;   

	for (int i = 0; i < (SS_ref_db.n_xeos); i++){
		SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
		SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];		
	}
	
   SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
   nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
   nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
   nlopt_set_min_objective(SS_ref_db.opt, obj_ig_cd, &SS_ref_db);
   nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, cd_ig_c, NULL, NULL);
   nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
   nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
   
   double minf;
   SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

   /* Send back needed local solution parameters */
   for (int i = 0; i < SS_ref_db.n_xeos; i++){
      SS_ref_db.xeos[i] = x[i];
   }
 
   SS_ref_db.df   = minf;
   nlopt_destroy(SS_ref_db.opt);

   return SS_ref_db;
};

SS_ref NLopt_opt_ig_cpx_function(global_variable gv, SS_ref SS_ref_db){

   int    n_em     = SS_ref_db.n_em;
   unsigned int n  = SS_ref_db.n_xeos;
   unsigned int m  = SS_ref_db.n_sf;
   
   double *x  = SS_ref_db.iguess; 
   

   for (int i = 0; i < (SS_ref_db.n_xeos); i++){
      SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
      SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
   }
   
   SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
   nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
   nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
   nlopt_set_min_objective(SS_ref_db.opt, obj_ig_cpx, &SS_ref_db);
   nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, cpx_ig_c, NULL, NULL);
   nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
   nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
   
   double minf;
   SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

   /* Send back needed local solution parameters */
   for (int i = 0; i < SS_ref_db.n_xeos; i++){
      SS_ref_db.xeos[i] = x[i];
   }
 
   SS_ref_db.df   = minf;
   nlopt_destroy(SS_ref_db.opt);

   return SS_ref_db;
};

SS_ref NLopt_opt_ig_ep_function(global_variable gv, SS_ref SS_ref_db){

   int    n_em     = SS_ref_db.n_em;
   unsigned int n  = SS_ref_db.n_xeos;
   unsigned int m  = SS_ref_db.n_sf;

   double *x  = SS_ref_db.iguess; 

   for (int i = 0; i < (SS_ref_db.n_xeos); i++){
      SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
      SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
   }
   
   SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
   nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
   nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
   nlopt_set_min_objective(SS_ref_db.opt, obj_ig_ep, &SS_ref_db);
   nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, ep_ig_c, NULL, NULL);
   nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
   nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
   
   double minf;
   SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

   /* Send back needed local solution parameters */
   for (int i = 0; i < SS_ref_db.n_xeos; i++){
      SS_ref_db.xeos[i] = x[i];
   }
 
   SS_ref_db.df   = minf;
   nlopt_destroy(SS_ref_db.opt);

   return SS_ref_db;
};

SS_ref NLopt_opt_ig_fl_function(global_variable gv, SS_ref SS_ref_db){

   int    n_em     = SS_ref_db.n_em;
   unsigned int n  = SS_ref_db.n_xeos;
   unsigned int m  = SS_ref_db.n_sf;

   double *x  = SS_ref_db.iguess; 

   for (int i = 0; i < (SS_ref_db.n_xeos); i++){
      SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
      SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
   }
   
   SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
   nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
   nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
   nlopt_set_min_objective(SS_ref_db.opt, obj_ig_fl, &SS_ref_db);
   nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, fl_ig_c, NULL, NULL);
   nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
   nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
   
   double minf;
   SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

   /* Send back needed local solution parameters */
   for (int i = 0; i < SS_ref_db.n_xeos; i++){
      SS_ref_db.xeos[i] = x[i];
   }
 
   SS_ref_db.df   = minf;
   nlopt_destroy(SS_ref_db.opt);

   return SS_ref_db;
};

SS_ref NLopt_opt_ig_g_function(global_variable gv, SS_ref SS_ref_db){
    
   int    n_em     = SS_ref_db.n_em;
   unsigned int n  = SS_ref_db.n_xeos;
   unsigned int m  = SS_ref_db.n_sf;

   double *x  = SS_ref_db.iguess; 

   for (int i = 0; i < (SS_ref_db.n_xeos); i++){
      SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
      SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
   }
   
   SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
   nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
   nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
   nlopt_set_min_objective(SS_ref_db.opt, obj_ig_g, &SS_ref_db);
   nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, g_ig_c, NULL, NULL);
   nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
   nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
   
   double minf;
   SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);
   
   /* Send back needed local solution parameters */
   for (int i = 0; i < SS_ref_db.n_xeos; i++){
      SS_ref_db.xeos[i] = x[i];
   }
 
   SS_ref_db.df   = minf;
   nlopt_destroy(SS_ref_db.opt);
	
   return SS_ref_db;
};

SS_ref NLopt_opt_ig_hb_function(global_variable gv, SS_ref SS_ref_db){

   int    n_em     = SS_ref_db.n_em;
   unsigned int n  = SS_ref_db.n_xeos;
   unsigned int m  = SS_ref_db.n_sf;

   double *x  = SS_ref_db.iguess; 

   for (int i = 0; i < (SS_ref_db.n_xeos); i++){
      SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
      SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
   }
   
   SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
   nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
   nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
   nlopt_set_min_objective(SS_ref_db.opt, obj_ig_hb, &SS_ref_db);
   nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, hb_ig_c, NULL, NULL);
   nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
   nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);
    
   double minf;
   SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

   /* Send back needed local solution parameters */
   for (int i = 0; i < SS_ref_db.n_xeos; i++){
      SS_ref_db.xeos[i] = x[i];
   }
 
   SS_ref_db.df   = minf;
   nlopt_destroy(SS_ref_db.opt);

   return SS_ref_db;
};

SS_ref NLopt_opt_ig_ilm_function(global_variable gv, SS_ref SS_ref_db){

   int    n_em     = SS_ref_db.n_em;
   unsigned int n  = SS_ref_db.n_xeos;
   unsigned int m  = SS_ref_db.n_sf;

   double *x  = SS_ref_db.iguess; 

   for (int i = 0; i < (SS_ref_db.n_xeos); i++){
      SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
      SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
   }
   
   SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
   nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
   nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
   nlopt_set_min_objective(SS_ref_db.opt, obj_ig_ilm, &SS_ref_db);
   nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, ilm_ig_c, NULL, NULL);
   nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
   nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
   
   double minf;
   SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

   /* Send back needed local solution parameters */
   for (int i = 0; i < SS_ref_db.n_xeos; i++){
      SS_ref_db.xeos[i] = x[i];
   }
 
   SS_ref_db.df   = minf;
   nlopt_destroy(SS_ref_db.opt);

   return SS_ref_db;
};

SS_ref NLopt_opt_ig_liq_function(global_variable gv, SS_ref SS_ref_db){

   int    n_em     = SS_ref_db.n_em;
   unsigned int n  = SS_ref_db.n_xeos;
   unsigned int m  = SS_ref_db.n_sf;
   
   double *x  = SS_ref_db.iguess; 
   
   for (int i = 0; i < (SS_ref_db.n_xeos); i++){
      SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
      SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
   }

   SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
   nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
   nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
   nlopt_set_min_objective(SS_ref_db.opt, obj_ig_liq, &SS_ref_db);
   nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, liq_ig_c, NULL, NULL);
   nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
   nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
   
   double minf;
   SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

   /* Send back needed local solution parameters */
   for (int i = 0; i < SS_ref_db.n_xeos; i++){
      SS_ref_db.xeos[i] = x[i];
   }
 
   SS_ref_db.df   = minf;
   nlopt_destroy(SS_ref_db.opt);

   return SS_ref_db;
};

SS_ref NLopt_opt_ig_mu_function(global_variable gv, SS_ref SS_ref_db){

   int    n_em     = SS_ref_db.n_em;
   unsigned int n  = SS_ref_db.n_xeos;
   unsigned int m  = SS_ref_db.n_sf;

   double *x  = SS_ref_db.iguess; 

   for (int i = 0; i < (SS_ref_db.n_xeos); i++){
      SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
      SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
   }
   
   SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
   nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
   nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
   nlopt_set_min_objective(SS_ref_db.opt, obj_ig_mu, &SS_ref_db);
   nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, mu_ig_c, NULL, NULL);
   nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
   nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
   
   double minf;
   SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

   /* Send back needed local solution parameters */
   for (int i = 0; i < SS_ref_db.n_xeos; i++){
      SS_ref_db.xeos[i] = x[i];
   }
 
   SS_ref_db.df   = minf;
   nlopt_destroy(SS_ref_db.opt);

   return SS_ref_db;
};

SS_ref NLopt_opt_ig_ol_function(global_variable gv, SS_ref SS_ref_db){
   
   int    n_em     = SS_ref_db.n_em;
   unsigned int n  = SS_ref_db.n_xeos;
   unsigned int m  = SS_ref_db.n_sf;

   double *x  = SS_ref_db.iguess; 

   for (int i = 0; i < (SS_ref_db.n_xeos); i++){
      SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
      SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
   }
   
   SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
   nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
   nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
   nlopt_set_min_objective(SS_ref_db.opt, obj_ig_ol, &SS_ref_db);
   nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, ol_ig_c, NULL, NULL);
   nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
   nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
   
   double minf;
   SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

   /* Send back needed local solution parameters */
   for (int i = 0; i < SS_ref_db.n_xeos; i++){
      SS_ref_db.xeos[i] = x[i];
   }
 
   SS_ref_db.df   = minf;
   nlopt_destroy(SS_ref_db.opt);
  
   return SS_ref_db;
};

SS_ref NLopt_opt_ig_opx_function(global_variable gv, SS_ref SS_ref_db){
    
   int    n_em     = SS_ref_db.n_em;
   unsigned int n  = SS_ref_db.n_xeos;
   unsigned int m  = SS_ref_db.n_sf;

   double *x  = SS_ref_db.iguess; 

   for (int i = 0; i < (SS_ref_db.n_xeos); i++){
      SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
      SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
   }
   
   SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
   nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
   nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
   nlopt_set_min_objective(SS_ref_db.opt, obj_ig_opx, &SS_ref_db);
   nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, opx_ig_c, NULL, NULL);
   nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
   nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
   
   double minf;
   SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

   /* Send back needed local solution parameters */
   for (int i = 0; i < SS_ref_db.n_xeos; i++){
      SS_ref_db.xeos[i] = x[i];
   }
 
   SS_ref_db.df   = minf;
   nlopt_destroy(SS_ref_db.opt);
  
   return SS_ref_db;
};

SS_ref NLopt_opt_ig_fsp_function(global_variable gv, SS_ref SS_ref_db){
   
   int    n_em     = SS_ref_db.n_em;
   unsigned int n  = SS_ref_db.n_xeos;
   unsigned int m  = SS_ref_db.n_sf;

   double *x  = SS_ref_db.iguess; 

   for (int i = 0; i < (SS_ref_db.n_xeos); i++){
      SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
      SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
   }
   
   SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
   nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
   nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
   nlopt_set_min_objective(SS_ref_db.opt, obj_ig_fsp, &SS_ref_db);
   nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, fsp_ig_c, NULL, NULL);
   nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
   nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
   
   double minf;
   SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

   /* Send back needed local solution parameters */
   for (int i = 0; i < SS_ref_db.n_xeos; i++){
      SS_ref_db.xeos[i] = x[i];
   }
 
   SS_ref_db.df   = minf;
   nlopt_destroy(SS_ref_db.opt);
  
   return SS_ref_db;
};

SS_ref NLopt_opt_ig_spl_function(global_variable gv, SS_ref SS_ref_db){

   int    n_em     = SS_ref_db.n_em;
   unsigned int n  = SS_ref_db.n_xeos;
   unsigned int m  = SS_ref_db.n_sf;

   double *x  = SS_ref_db.iguess; 

   for (int i = 0; i < (SS_ref_db.n_xeos); i++){
      SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
      SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
   }
   
   SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n));
   nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
   nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
   nlopt_set_min_objective(SS_ref_db.opt, obj_ig_spl, &SS_ref_db);
   nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, spl_ig_c, &SS_ref_db, NULL);
   nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
   nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);

   double minf;
   SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

   /* Send back needed local solution parameters */
   for (int i = 0; i < SS_ref_db.n_xeos; i++){
      SS_ref_db.xeos[i]   = x[i];
   }
 
   SS_ref_db.df   = minf;
   nlopt_destroy(SS_ref_db.opt);
  	
   return SS_ref_db;
};

//---------------------------------------------------------------------------
//------------------------Alkaline dry Weller et al., 2024-------------------
//---------------------------------------------------------------------------
SS_ref NLopt_opt_igad_liq_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }

    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 

    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_igad_liq, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, liq_igad_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }

    SS_ref_db.df   = minf;

    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};

SS_ref NLopt_opt_igad_fsp_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_igad_fsp, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, fsp_igad_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);
 
    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};
SS_ref NLopt_opt_igad_spl_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_igad_spl, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, spl_igad_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};
SS_ref NLopt_opt_igad_g_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_igad_g, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, g_igad_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};
SS_ref NLopt_opt_igad_ol_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_igad_ol, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, ol_igad_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};
SS_ref NLopt_opt_igad_opx_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_igad_opx, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, opx_igad_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};
SS_ref NLopt_opt_igad_cpx_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_igad_cpx, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, cpx_igad_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};
SS_ref NLopt_opt_igad_ilm_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_igad_ilm, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, ilm_igad_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};
SS_ref NLopt_opt_igad_nph_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_igad_nph, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, nph_igad_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};
SS_ref NLopt_opt_igad_lct_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_igad_lct, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, lct_igad_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};
SS_ref NLopt_opt_igad_kals_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_igad_kals, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, kals_igad_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};
SS_ref NLopt_opt_igad_mel_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_igad_mel, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, mel_igad_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};

//---------------------------------------------------------------------------
//---------------------------------Evans&Frost,2021--------------------------
//---------------------------------------------------------------------------
SS_ref NLopt_opt_um_fluid_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_um_fluid, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, fluid_um_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};

SS_ref NLopt_opt_um_ol_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_um_ol, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, ol_um_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};

SS_ref NLopt_opt_um_br_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_um_br, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, br_um_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};

SS_ref NLopt_opt_um_ch_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_um_ch, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, ch_um_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};

SS_ref NLopt_opt_um_atg_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_um_atg, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, atg_um_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};

SS_ref NLopt_opt_um_g_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_um_g, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, g_um_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};

SS_ref NLopt_opt_um_ta_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_um_ta, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, ta_um_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);
    
    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};

SS_ref NLopt_opt_um_chl_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_um_chl, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, chl_um_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};

SS_ref NLopt_opt_um_anth_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_um_anth, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, anth_um_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};

SS_ref NLopt_opt_um_spi_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_um_spi, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, spi_um_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};

SS_ref NLopt_opt_um_opx_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_um_opx, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, opx_um_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};

SS_ref NLopt_opt_um_po_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_um_po, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, po_um_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};
SS_ref NLopt_opt_ume_pl4tr_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_ume_pl4tr, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, pl4tr_ume_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
     nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};
SS_ref NLopt_opt_ume_hb_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_ume_hb, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, hb_ume_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};
SS_ref NLopt_opt_ume_aug_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_ume_aug, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, aug_ume_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};

/**
    Equality constraints aq17
*/
void aq17_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d         = (SS_ref *) SS_ref_db;

    int n_em          = d->n_em;
    double *charge    = d->mat_phi;

    result[0] = - 1.0;
    result[1] =   0.0;
    for (int i = 0; i < n_em; i++){
        result[0] += x[i];
        result[1] += charge[i]*x[i];
    }

    int j = 0;
    if (grad) {
        for (int i = 0; i < n_em; i++){
            grad[j] = 1.0;
            j += 1;
        }
        for (int i = 0; i < n_em; i++){
            grad[j] = charge[i];
            j += 1;
        }
    }

    return;
};
/* NLopt function to minimize for aqueous fluid */
SS_ref NLopt_opt_aq17_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = 2;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_aq17, &SS_ref_db);
    nlopt_add_equality_mconstraint(SS_ref_db.opt, m, aq17_c, &SS_ref_db, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    // nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }

    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};


/**************************************************************************************/
/**************************************************************************************/
/**************************************************************************************/
/**************************************************************************************/

/**
    Inequality constraints for g
*/
void g_mtl_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = (-x[0]*x[1] - 1.0/3.0*x[0]*x[4] + x[0] + x[1]*x[3] + x[1] + 1.0/3.0*x[3]*x[4] -x[3] + 1.0/3.0*x[4] - 1.0);
    result[1] = (x[0]*x[1] + 1.0/3.0*x[0]*x[4] -x[0] -x[1]*x[3] - 1.0/3.0*x[3]*x[4] + x[3]);
    result[2] = (-x[1]);
    result[3] = (-1.0/3.0*x[4]);
    result[4] = (x[2] + 0.5*x[4] - 1.0);
    result[5] = (0.5*x[0]*x[2] - 1.5*x[1]*x[3] - 0.5*x[2] - 0.5*x[3]*x[4] + 1.5*x[3]);
    result[6] = (-0.5*x[0]*x[2] + 1.5*x[1]*x[3] + 0.5*x[3]*x[4] - 1.5*x[3]);
    result[7] = (-0.5*x[2] - 0.5*x[4]);

    if (grad) {
        grad[0] = -x[1] - 1.0/3.0*x[4] + 1.0;
        grad[1] = -x[0] + x[3] + 1.0;
        grad[2] = 0.0;
        grad[3] = x[1] + 1.0/3.0*x[4] - 1.0;
        grad[4] = -1.0/3.0*x[0] + 1.0/3.0*x[3] + 1.0/3.0;
        grad[5] = x[1] + 1.0/3.0*x[4] - 1.0;
        grad[6] = x[0] -x[3];
        grad[7] = 0.0;
        grad[8] = -x[1] - 1.0/3.0*x[4] + 1.0;
        grad[9] = 1.0/3.0*x[0] - 1.0/3.0*x[3];
        grad[10] = 0.0;
        grad[11] = -1.0;
        grad[12] = 0.0;
        grad[13] = 0.0;
        grad[14] = 0.0;
        grad[15] = 0.0;
        grad[16] = 0.0;
        grad[17] = 0.0;
        grad[18] = 0.0;
        grad[19] = -1.0/3.0;
        grad[20] = 0.0;
        grad[21] = 0.0;
        grad[22] = 1.0;
        grad[23] = 0.0;
        grad[24] = 0.50;
        grad[25] = 0.5*x[2];
        grad[26] = -1.5*x[3];
        grad[27] = 0.5*x[0] - 0.5;
        grad[28] = -1.5*x[1] - 0.5*x[4] + 1.5;
        grad[29] = -0.5*x[3];
        grad[30] = -0.5*x[2];
        grad[31] = 1.5*x[3];
        grad[32] = -0.5*x[0];
        grad[33] = 1.5*x[1] + 0.5*x[4] - 1.5;
        grad[34] = 0.5*x[3];
        grad[35] = 0.0;
        grad[36] = 0.0;
        grad[37] = -0.50;
        grad[38] = 0.0;
        grad[39] = -0.50;
    }

    return;
};

/**
    Inequality constraints for fp
*/
void fp_mtl_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = (x[0] - 1.0);
    result[1] = (-x[0]);

    if (grad) {
        grad[0] = 1.0;
        grad[1] = -1.0;
    }

    return;
};

/**
    Inequality constraints for mpv
*/
void mpv_mtl_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = (-x[2]);
    result[1] = (-x[0]*x[1] -x[0]*x[2] -x[0]*x[3] + x[0] + x[1] + x[2] + x[3] - 1.0);
    result[2] = (x[0]*x[1] + x[0]*x[2] + x[0]*x[3] -x[0]);
    result[3] = (-0.5*x[3]);
    result[4] = (-x[1] - 0.5*x[3]);
    result[5] = (-x[1]);
    result[6] = (x[1] - 1.0);

    if (grad) {
        grad[0] = 0.0;
        grad[1] = 0.0;
        grad[2] = -1.0;
        grad[3] = 0.0;
        grad[4] = -x[1] -x[2] -x[3] + 1.0;
        grad[5] = 1.0 -x[0];
        grad[6] = 1.0 -x[0];
        grad[7] = 1.0 -x[0];
        grad[8] = x[1] + x[2] + x[3] - 1.0;
        grad[9] = x[0];
        grad[10] = x[0];
        grad[11] = x[0];
        grad[12] = 0.0;
        grad[13] = 0.0;
        grad[14] = 0.0;
        grad[15] = -0.50;
        grad[16] = 0.0;
        grad[17] = -1.0;
        grad[18] = 0.0;
        grad[19] = -0.50;
        grad[20] = 0.0;
        grad[21] = -1.0;
        grad[22] = 0.0;
        grad[23] = 0.0;
        grad[24] = 0.0;
        grad[25] = 1.0;
        grad[26] = 0.0;
        grad[27] = 0.0;
    }

    return;
};

/**
    Inequality constraints for cpv
*/
void cpv_mtl_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = (-x[2]);
    result[1] = (-x[0]*x[1] -x[0]*x[2] -x[0]*x[3] + x[0] + x[1] + x[2] + x[3] - 1.0);
    result[2] = (x[0]*x[1] + x[0]*x[2] + x[0]*x[3] -x[0]);
    result[3] = (-0.5*x[3]);
    result[4] = (-x[1] - 0.5*x[3]);
    result[5] = (-x[1]);
    result[6] = (x[1] - 1.0);

    if (grad) {
        grad[0] = 0.0;
        grad[1] = 0.0;
        grad[2] = -1.0;
        grad[3] = 0.0;
        grad[4] = -x[1] -x[2] -x[3] + 1.0;
        grad[5] = 1.0 -x[0];
        grad[6] = 1.0 -x[0];
        grad[7] = 1.0 -x[0];
        grad[8] = x[1] + x[2] + x[3] - 1.0;
        grad[9] = x[0];
        grad[10] = x[0];
        grad[11] = x[0];
        grad[12] = 0.0;
        grad[13] = 0.0;
        grad[14] = 0.0;
        grad[15] = -0.50;
        grad[16] = 0.0;
        grad[17] = -1.0;
        grad[18] = 0.0;
        grad[19] = -0.50;
        grad[20] = 0.0;
        grad[21] = -1.0;
        grad[22] = 0.0;
        grad[23] = 0.0;
        grad[24] = 0.0;
        grad[25] = 1.0;
        grad[26] = 0.0;
        grad[27] = 0.0;
    }

    return;
};
/**
    Inequality constraints for crn
*/
void crn_mtl_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = (-x[0]*x[1] + x[0] + x[1] - 1.0);
    result[1] = (x[0]*x[1] -x[0]);
    result[2] = (-x[1]);
    result[3] = (-x[1]);
    result[4] = (x[1] - 1.0);

    if (grad) {
        grad[0] = 1.0 -x[1];
        grad[1] = 1.0 -x[0];
        grad[2] = x[1] - 1.0;
        grad[3] = x[0];
        grad[4] = 0.0;
        grad[5] = -1.0;
        grad[6] = 0.0;
        grad[7] = -1.0;
        grad[8] = 0.0;
        grad[9] = 1.0;
    }

    return;
};

/**
    Inequality constraints for cf
*/
void cf_mtl_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = (-x[3]);
    result[1] = (-x[1]*x[3] -x[1]*x[4] + x[1] -x[2]*x[3] -x[2]*x[4] + x[2] + x[3] + x[4] - 1.0);
    result[2] = (x[1]*x[3] + x[1]*x[4] -x[1] + x[2]*x[3] + x[2]*x[4] -x[2]);
    result[3] = (-x[4]);
    result[4] = (0.5*x[0]*x[1] - 0.5*x[0] + 0.5*x[2]*x[3] + 0.5*x[2]*x[4] - 0.5*x[2]);
    result[5] = (-0.5*x[0]*x[1] - 0.5*x[2]*x[3] - 0.5*x[2]*x[4] + 0.5*x[2]);
    result[6] = (x[0] + 0.5*x[4] - 1.0);
    result[7] = (-0.5*x[0] - 0.5*x[4]);

    if (grad) {
        grad[0] = 0.0;
        grad[1] = 0.0;
        grad[2] = 0.0;
        grad[3] = -1.0;
        grad[4] = 0.0;
        grad[5] = 0.0;
        grad[6] = -x[3] -x[4] + 1.0;
        grad[7] = -x[3] -x[4] + 1.0;
        grad[8] = -x[1] -x[2] + 1.0;
        grad[9] = -x[1] -x[2] + 1.0;
        grad[10] = 0.0;
        grad[11] = x[3] + x[4] - 1.0;
        grad[12] = x[3] + x[4] - 1.0;
        grad[13] = x[1] + x[2];
        grad[14] = x[1] + x[2];
        grad[15] = 0.0;
        grad[16] = 0.0;
        grad[17] = 0.0;
        grad[18] = 0.0;
        grad[19] = -1.0;
        grad[20] = 0.5*x[1] - 0.5;
        grad[21] = 0.5*x[0];
        grad[22] = 0.5*x[3] + 0.5*x[4] - 0.5;
        grad[23] = 0.5*x[2];
        grad[24] = 0.5*x[2];
        grad[25] = -0.5*x[1];
        grad[26] = -0.5*x[0];
        grad[27] = -0.5*x[3] - 0.5*x[4] + 0.5;
        grad[28] = -0.5*x[2];
        grad[29] = -0.5*x[2];
        grad[30] = 1.0;
        grad[31] = 0.0;
        grad[32] = 0.0;
        grad[33] = 0.0;
        grad[34] = 0.50;
        grad[35] = -0.50;
        grad[36] = 0.0;
        grad[37] = 0.0;
        grad[38] = 0.0;
        grad[39] = -0.50;
    }

    return;
};

/**
    Inequality constraints for nal
*/
void nal_mtl_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = (-x[4]);
    result[1] = (-x[1]*x[4] -x[1]*x[5] + x[1] + x[2]*x[4] + x[2]*x[5] -x[2] + x[4] + x[5] - 1.0);
    result[2] = (x[1]*x[4] + x[1]*x[5] -x[1] -x[2]*x[4] -x[2]*x[5] + x[2]);
    result[3] = (-x[5]);
    result[4] = (x[1] -x[3] - 1.0);
    result[5] = (-x[1] + x[3]);
    result[6] = (-0.5*x[0]*x[1] + 0.5*x[0] - 0.0833333333333333*x[1]*x[5] + 0.5*x[1] - 0.166666666666667*x[2]*x[4] - 0.166666666666667*x[2]*x[5] + 0.166666666666667*x[2] + 1.0/3.0*x[3] + 0.0833333333333333*x[5] - 0.5);
    result[7] = (0.5*x[0]*x[1] + 0.0833333333333333*x[1]*x[5] - 0.5*x[1] + 0.166666666666667*x[2]*x[4] + 0.166666666666667*x[2]*x[5] - 0.166666666666667*x[2] - 1.0/3.0*x[3]);
    result[8] = (-x[0]);
    result[9] = (0.5*x[0] - 0.0833333333333333*x[5] - 0.5);

    if (grad) {
        grad[0] = 0.0;
        grad[1] = 0.0;
        grad[2] = 0.0;
        grad[3] = 0.0;
        grad[4] = -1.0;
        grad[5] = 0.0;
        grad[6] = 0.0;
        grad[7] = -x[4] -x[5] + 1.0;
        grad[8] = x[4] + x[5] - 1.0;
        grad[9] = 0.0;
        grad[10] = -x[1] + x[2] + 1.0;
        grad[11] = -x[1] + x[2] + 1.0;
        grad[12] = 0.0;
        grad[13] = x[4] + x[5] - 1.0;
        grad[14] = -x[4] -x[5] + 1.0;
        grad[15] = 0.0;
        grad[16] = x[1] -x[2];
        grad[17] = x[1] -x[2];
        grad[18] = 0.0;
        grad[19] = 0.0;
        grad[20] = 0.0;
        grad[21] = 0.0;
        grad[22] = 0.0;
        grad[23] = -1.0;
        grad[24] = 0.0;
        grad[25] = 1.0;
        grad[26] = 0.0;
        grad[27] = -1.0;
        grad[28] = 0.0;
        grad[29] = 0.0;
        grad[30] = 0.0;
        grad[31] = -1.0;
        grad[32] = 0.0;
        grad[33] = 1.0;
        grad[34] = 0.0;
        grad[35] = 0.0;
        grad[36] = 0.5 - 0.5*x[1];
        grad[37] = -0.5*x[0] - 0.0833333333333333*x[5] + 0.5;
        grad[38] = -0.166666666666667*x[4] - 0.166666666666667*x[5] + 0.166666666666667;
        grad[39] = 1.0/3.0;
        grad[40] = -0.166666666666667*x[2];
        grad[41] = -0.0833333333333333*x[1] - 0.166666666666667*x[2] + 0.0833333333333333;
        grad[42] = 0.5*x[1];
        grad[43] = 0.5*x[0] + 0.0833333333333333*x[5] - 0.5;
        grad[44] = 0.166666666666667*x[4] + 0.166666666666667*x[5] - 0.166666666666667;
        grad[45] = -1.0/3.0;
        grad[46] = 0.166666666666667*x[2];
        grad[47] = 0.0833333333333333*x[1] + 0.166666666666667*x[2];
        grad[48] = -1.0;
        grad[49] = 0.0;
        grad[50] = 0.0;
        grad[51] = 0.0;
        grad[52] = 0.0;
        grad[53] = 0.0;
        grad[54] = 0.50;
        grad[55] = 0.0;
        grad[56] = 0.0;
        grad[57] = 0.0;
        grad[58] = 0.0;
        grad[59] = -0.0833333333333333;
    }

    return;
};

/**
    Inequality constraints for aki
*/
void aki_mtl_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = (-x[1]);
    result[1] = (-x[0]*x[1] + x[0] + x[1] - 1.0);
    result[2] = (x[0]*x[1] -x[0]);
    result[3] = (-x[1]);
    result[4] = (x[1] - 1.0);

    if (grad) {
        grad[0] = 0.0;
        grad[1] = -1.0;
        grad[2] = 1.0 -x[1];
        grad[3] = 1.0 -x[0];
        grad[4] = x[1] - 1.0;
        grad[5] = x[0];
        grad[6] = 0.0;
        grad[7] = -1.0;
        grad[8] = 0.0;
        grad[9] = 1.0;
    }

    return;
};

/**
    Inequality constraints for ol
*/
void ol_mtl_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = (x[0] - 1.0);
    result[1] = (-x[0]);

    if (grad) {
        grad[0] = 1.0;
        grad[1] = -1.0;
    }

    return;
};

/**
    Inequality constraints for wad
*/
void wad_mtl_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = (x[0] - 1.0);
    result[1] = (-x[0]);

    if (grad) {
        grad[0] = 1.0;
        grad[1] = -1.0;
    }

    return;
};

/**
    Inequality constraints for ring
*/
void ring_mtl_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = (x[0] - 1.0);
    result[1] = (-x[0]);

    if (grad) {
        grad[0] = 1.0;
        grad[1] = -1.0;
    }

    return;
};

/**
    Inequality constraints for cpx
*/
void cpx_mtl_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = (-x[0]*x[1] -x[0]*x[3] + x[0] -x[1]*x[4] + x[1] -x[3]*x[4] + x[3] + x[4] - 1.0);
    result[1] = (x[0]*x[1] + x[0]*x[3] -x[0] + x[1]*x[4] + x[3]*x[4] -x[4]);
    result[2] = (-x[1] -x[3]);
    result[3] = (x[0]*x[2] + x[1]*x[4] -x[2] + x[3]*x[4] -x[4]);
    result[4] = (-x[0]*x[2] -x[1]*x[4] -x[3]*x[4] + x[4]);
    result[5] = (x[2] + x[3] - 1.0);
    result[6] = (-x[3]);
    result[7] = (0.5*x[1] - 1.0);
    result[8] = (-0.5*x[1]);

    if (grad) {
        grad[0] = -x[1] -x[3] + 1.0;
        grad[1] = -x[0] -x[4] + 1.0;
        grad[2] = 0.0;
        grad[3] = -x[0] -x[4] + 1.0;
        grad[4] = -x[1] -x[3] + 1.0;
        grad[5] = x[1] + x[3] - 1.0;
        grad[6] = x[0] + x[4];
        grad[7] = 0.0;
        grad[8] = x[0] + x[4];
        grad[9] = x[1] + x[3] - 1.0;
        grad[10] = 0.0;
        grad[11] = -1.0;
        grad[12] = 0.0;
        grad[13] = -1.0;
        grad[14] = 0.0;
        grad[15] = x[2];
        grad[16] = x[4];
        grad[17] = x[0] - 1.0;
        grad[18] = x[4];
        grad[19] = x[1] + x[3] - 1.0;
        grad[20] = -x[2];
        grad[21] = -x[4];
        grad[22] = -x[0];
        grad[23] = -x[4];
        grad[24] = -x[1] -x[3] + 1.0;
        grad[25] = 0.0;
        grad[26] = 0.0;
        grad[27] = 1.0;
        grad[28] = 1.0;
        grad[29] = 0.0;
        grad[30] = 0.0;
        grad[31] = 0.0;
        grad[32] = 0.0;
        grad[33] = -1.0;
        grad[34] = 0.0;
        grad[35] = 0.0;
        grad[36] = 0.50;
        grad[37] = 0.0;
        grad[38] = 0.0;
        grad[39] = 0.0;
        grad[40] = 0.0;
        grad[41] = -0.50;
        grad[42] = 0.0;
        grad[43] = 0.0;
        grad[44] = 0.0;
    }

    return;
};

/**
    Inequality constraints for opx
*/
void opx_mtl_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = (-x[0]*x[1] -x[0]*x[2] + x[0] + x[1] - 0.5*x[3] - 1.0);
    result[1] = (x[0]*x[1] + x[0]*x[2] -x[0] + 0.5*x[3]);
    result[2] = (-x[1]);
    result[3] = (-x[2]);
    result[4] = (x[0] + x[2] + 0.5*x[3] - 1.0);
    result[5] = (-x[0] - 0.5*x[3]);
    result[6] = (0.5*x[1] - 1.0);
    result[7] = (-0.5*x[1]);

    if (grad) {
        grad[0] = -x[1] -x[2] + 1.0;
        grad[1] = 1.0 -x[0];
        grad[2] = -x[0];
        grad[3] = -0.50;
        grad[4] = x[1] + x[2] - 1.0;
        grad[5] = x[0];
        grad[6] = x[0];
        grad[7] = 0.50;
        grad[8] = 0.0;
        grad[9] = -1.0;
        grad[10] = 0.0;
        grad[11] = 0.0;
        grad[12] = 0.0;
        grad[13] = 0.0;
        grad[14] = -1.0;
        grad[15] = 0.0;
        grad[16] = 1.0;
        grad[17] = 0.0;
        grad[18] = 1.0;
        grad[19] = 0.50;
        grad[20] = -1.0;
        grad[21] = 0.0;
        grad[22] = 0.0;
        grad[23] = -0.50;
        grad[24] = 0.0;
        grad[25] = 0.50;
        grad[26] = 0.0;
        grad[27] = 0.0;
        grad[28] = 0.0;
        grad[29] = -0.50;
        grad[30] = 0.0;
        grad[31] = 0.0;
    }

    return;
};

/**
    Inequality constraints for hpx
*/
void hpx_mtl_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = (-x[0]*x[1] -x[0]*x[2] + x[0] + x[1] - 0.5*x[3] - 1.0);
    result[1] = (x[0]*x[1] + x[0]*x[2] -x[0] + 0.5*x[3]);
    result[2] = (-x[1]);
    result[3] = (-x[2]);
    result[4] = (x[0] + x[2] + 0.5*x[3] - 1.0);
    result[5] = (-x[0] - 0.5*x[3]);
    result[6] = (0.5*x[1] - 1.0);
    result[7] = (-0.5*x[1]);

    if (grad) {
        grad[0] = -x[1] -x[2] + 1.0;
        grad[1] = 1.0 -x[0];
        grad[2] = -x[0];
        grad[3] = -0.50;
        grad[4] = x[1] + x[2] - 1.0;
        grad[5] = x[0];
        grad[6] = x[0];
        grad[7] = 0.50;
        grad[8] = 0.0;
        grad[9] = -1.0;
        grad[10] = 0.0;
        grad[11] = 0.0;
        grad[12] = 0.0;
        grad[13] = 0.0;
        grad[14] = -1.0;
        grad[15] = 0.0;
        grad[16] = 1.0;
        grad[17] = 0.0;
        grad[18] = 1.0;
        grad[19] = 0.50;
        grad[20] = -1.0;
        grad[21] = 0.0;
        grad[22] = 0.0;
        grad[23] = -0.50;
        grad[24] = 0.0;
        grad[25] = 0.50;
        grad[26] = 0.0;
        grad[27] = 0.0;
        grad[28] = 0.0;
        grad[29] = -0.50;
        grad[30] = 0.0;
        grad[31] = 0.0;
    }

    return;
};

SS_ref NLopt_opt_mtl_g_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_mtl_g, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, g_mtl_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    
    double minf;
    if (gv.maxeval==1){  
       // we are only interested in evaluating the objective function  
       minf = obj_mtl_g(n, x, NULL, &SS_ref_db);
    }
    else{
      // do optimization
      SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);
    }
    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};
SS_ref NLopt_opt_mtl_fp_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_mtl_fp, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, fp_mtl_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    
    double minf;
    if (gv.maxeval==1){  
       // we are only interested in evaluating the objective function  
       minf = obj_mtl_fp(n, x, NULL, &SS_ref_db);
    }
    else{
      // do optimization
      SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);
    }
    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};
SS_ref NLopt_opt_mtl_mpv_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_mtl_mpv, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, mpv_mtl_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    
    double minf;
    if (gv.maxeval==1){  
       // we are only interested in evaluating the objective function  
       minf = obj_mtl_mpv(n, x, NULL, &SS_ref_db);
    }
    else{
      // do optimization
      SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);
    }
    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};
SS_ref NLopt_opt_mtl_cpv_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_mtl_cpv, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, cpv_mtl_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    
    double minf;
    if (gv.maxeval==1){  
       // we are only interested in evaluating the objective function  
       minf = obj_mtl_cpv(n, x, NULL, &SS_ref_db);
    }
    else{
      // do optimization
      SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);
    }
    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};
SS_ref NLopt_opt_mtl_crn_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_mtl_crn, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, crn_mtl_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    
    double minf;
    if (gv.maxeval==1){  
       // we are only interested in evaluating the objective function  
       minf = obj_mtl_crn(n, x, NULL, &SS_ref_db);
    }
    else{
      // do optimization
      SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);
    }
    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};
SS_ref NLopt_opt_mtl_cf_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_mtl_cf, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, cf_mtl_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    
    double minf;
    if (gv.maxeval==1){  
       // we are only interested in evaluating the objective function  
       minf = obj_mtl_cf(n, x, NULL, &SS_ref_db);
    }
    else{
      // do optimization
      SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);
    }
    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};
SS_ref NLopt_opt_mtl_nal_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_mtl_nal, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, nal_mtl_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    
    double minf;
    if (gv.maxeval==1){  
       // we are only interested in evaluating the objective function  
       minf = obj_mtl_nal(n, x, NULL, &SS_ref_db);
    }
    else{
      // do optimization
      SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);
    }
    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};
SS_ref NLopt_opt_mtl_aki_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_mtl_aki, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, aki_mtl_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    
    double minf;
    if (gv.maxeval==1){  
       // we are only interested in evaluating the objective function  
       minf = obj_mtl_aki(n, x, NULL, &SS_ref_db);
    }
    else{
      // do optimization
      SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);
    }
    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};
SS_ref NLopt_opt_mtl_ol_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_mtl_ol, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, ol_mtl_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    
    double minf;
    if (gv.maxeval==1){  
       // we are only interested in evaluating the objective function  
       minf = obj_mtl_ol(n, x, NULL, &SS_ref_db);
    }
    else{
      // do optimization
      SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);
    }
    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};
SS_ref NLopt_opt_mtl_wad_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_mtl_wad, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, wad_mtl_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    
    double minf;
    if (gv.maxeval==1){  
       // we are only interested in evaluating the objective function  
       minf = obj_mtl_wad(n, x, NULL, &SS_ref_db);
    }
    else{
      // do optimization
      SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);
    }
    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};
SS_ref NLopt_opt_mtl_ring_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_mtl_ring, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, ring_mtl_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    
    double minf;
    if (gv.maxeval==1){  
       // we are only interested in evaluating the objective function  
       minf = obj_mtl_ring(n, x, NULL, &SS_ref_db);
    }
    else{
      // do optimization
      SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);
    }
    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};
SS_ref NLopt_opt_mtl_cpx_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_mtl_cpx, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, cpx_mtl_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    
    double minf;
    if (gv.maxeval==1){  
       // we are only interested in evaluating the objective function  
       minf = obj_mtl_cpx(n, x, NULL, &SS_ref_db);
    }
    else{
      // do optimization
      SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);
    }
    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};
SS_ref NLopt_opt_mtl_opx_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_mtl_opx, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, opx_mtl_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    
    double minf;
    if (gv.maxeval==1){  
       // we are only interested in evaluating the objective function  
       minf = obj_mtl_opx(n, x, NULL, &SS_ref_db);
    }
    else{
      // do optimization
      SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);
    }
    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};
SS_ref NLopt_opt_mtl_hpx_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_mtl_hpx, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, hpx_mtl_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    
    double minf;
    if (gv.maxeval==1){  
       // we are only interested in evaluating the objective function  
       minf = obj_mtl_hpx(n, x, NULL, &SS_ref_db);
    }
    else{
      // do optimization
      SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);
    }
    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};





/**************************************************************************************/
/**************************************************************************************/
/*   Metapelite ext DB (White et al., 2014; Green et al., 2016; Evans & Forst, 2021)  */
/**************************************************************************************/
/**************************************************************************************/


/**
    Inequality constraints for liq_mp
*/
void liq_mpe_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( x[6] - 1.0);
    result[1] = ( -x[0]);
    result[2] = ( -x[1]*x[2]);
    result[3] = ( -x[1]*(1.0 - x[2]));
    result[4] = ( -x[3]);
    result[5] = ( x[3] + x[1] + x[6] + x[4] + x[0] - 1.0);
    result[6] = ( -x[4]);
    result[7] = ( -x[5]);
    result[8] = ( x[5] - 1.0);
    result[9] = ( -x[6]);

    if (grad) {
        grad[0] = 0.0;
        grad[1] = 0.0;
        grad[2] = 0.0;
        grad[3] = 0.0;
        grad[4] = 0.0;
        grad[5] = 0.0;
        grad[6] = 1.0;
        grad[7] = -1.0;
        grad[8] = 0.0;
        grad[9] = 0.0;
        grad[10] = 0.0;
        grad[11] = 0.0;
        grad[12] = 0.0;
        grad[13] = 0.0;
        grad[14] = 0.0;
        grad[15] = -x[2];
        grad[16] = -x[1];
        grad[17] = 0.0;
        grad[18] = 0.0;
        grad[19] = 0.0;
        grad[20] = 0.0;
        grad[21] = 0.0;
        grad[22] = x[2] - 1.0;
        grad[23] = x[1];
        grad[24] = 0.0;
        grad[25] = 0.0;
        grad[26] = 0.0;
        grad[27] = 0.0;
        grad[28] = 0.0;
        grad[29] = 0.0;
        grad[30] = 0.0;
        grad[31] = -1.0;
        grad[32] = 0.0;
        grad[33] = 0.0;
        grad[34] = 0.0;
        grad[35] = 1.0;
        grad[36] = 1.0;
        grad[37] = 0.0;
        grad[38] = 1.0;
        grad[39] = 1.0;
        grad[40] = 0.0;
        grad[41] = 1.0;
        grad[42] = 0.0;
        grad[43] = 0.0;
        grad[44] = 0.0;
        grad[45] = 0.0;
        grad[46] = -1.0;
        grad[47] = 0.0;
        grad[48] = 0.0;
        grad[49] = 0.0;
        grad[50] = 0.0;
        grad[51] = 0.0;
        grad[52] = 0.0;
        grad[53] = 0.0;
        grad[54] = -1.0;
        grad[55] = 0.0;
        grad[56] = 0.0;
        grad[57] = 0.0;
        grad[58] = 0.0;
        grad[59] = 0.0;
        grad[60] = 0.0;
        grad[61] = 1.0;
        grad[62] = 0.0;
        grad[63] = 0.0;
        grad[64] = 0.0;
        grad[65] = 0.0;
        grad[66] = 0.0;
        grad[67] = 0.0;
        grad[68] = 0.0;
        grad[69] = -1.0;
    }

    return;
};

/**
    Inequality constraints for st_mp
*/
void st_mpe_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( -x[1]*x[0] + x[1] + x[0] - 1.0);
    result[1] = ( x[1]*x[0] - x[0]);
    result[2] = ( -x[1]);
    result[3] = ( x[2] + 4.0/3.0*x[3] - 1.0);
    result[4] = ( -x[2]);
    result[5] = ( -x[3]);
    result[6] = ( -1./3.*x[3]);

    if (grad) {
        grad[0] = 1.0 - x[1];
        grad[1] = 1.0 - x[0];
        grad[2] = 0.0;
        grad[3] = 0.0;
        grad[4] = x[1] - 1.0;
        grad[5] = x[0];
        grad[6] = 0.0;
        grad[7] = 0.0;
        grad[8] = 0.0;
        grad[9] = -1.0;
        grad[10] = 0.0;
        grad[11] = 0.0;
        grad[12] = 0.0;
        grad[13] = 0.0;
        grad[14] = 1.0;
        grad[15] = 4.0/3.0;
        grad[16] = 0.0;
        grad[17] = 0.0;
        grad[18] = -1.0;
        grad[19] = 0.0;
        grad[20] = 0.0;
        grad[21] = 0.0;
        grad[22] = 0.0;
        grad[23] = -1.0;
        grad[24] = 0.0;
        grad[25] = 0.0;
        grad[26] = 0.0;
        grad[27] = -1./3.;
    }

    return;
};

/**
    Inequality constraints for sp_mp
*/
void sp_mpe_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( -x[1]);
    result[1] = ( x[1] + x[2] - 1.0);
    result[2] = ( -x[2]);
    result[3] = ( x[0] - 1.0);
    result[4] = ( -x[0]);

    if (grad) {
        grad[0] = 0.0;
        grad[1] = -1.0;
        grad[2] = 0.0;
        grad[3] = 0.0;
        grad[4] = 1.0;
        grad[5] = 1.0;
        grad[6] = 0.0;
        grad[7] = 0.0;
        grad[8] = -1.0;
        grad[9] = 1.0;
        grad[10] = 0.0;
        grad[11] = 0.0;
        grad[12] = -1.0;
        grad[13] = 0.0;
        grad[14] = 0.0;
    }

    return;
};

/**
    Inequality constraints for sa_mp
*/
void sa_mpe_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( -x[2]*x[0] + x[2] - 0.75*x[3] - x[0]*x[1] + x[0] + x[1] - 1.0);
    result[1] = ( x[2]*x[0] + 0.75*x[3] + x[0]*x[1] - x[0]);
    result[2] = ( -x[2]);
    result[3] = ( -x[1]);
    result[4] = ( 0.25*x[3] + x[0] - 1.0);
    result[5] = ( -0.25*x[3] - x[0]);
    result[6] = ( x[2] + x[1] - 1.0);
    result[7] = ( -x[2] - x[1]);

    if (grad) {
        grad[0] = -x[2] - x[1] + 1.0;
        grad[1] = 1.0 - x[0];
        grad[2] = 1.0 - x[0];
        grad[3] = -0.75;
        grad[4] = x[2] + x[1] - 1.0;
        grad[5] = x[0];
        grad[6] = x[0];
        grad[7] = 0.75;
        grad[8] = 0.0;
        grad[9] = 0.0;
        grad[10] = -1.0;
        grad[11] = 0.0;
        grad[12] = 0.0;
        grad[13] = -1.0;
        grad[14] = 0.0;
        grad[15] = 0.0;
        grad[16] = 1.0;
        grad[17] = 0.0;
        grad[18] = 0.0;
        grad[19] = 0.25;
        grad[20] = -1.0;
        grad[21] = 0.0;
        grad[22] = 0.0;
        grad[23] = -0.25;
        grad[24] = 0.0;
        grad[25] = 1.0;
        grad[26] = 1.0;
        grad[27] = 0.0;
        grad[28] = 0.0;
        grad[29] = -1.0;
        grad[30] = -1.0;
        grad[31] = 0.0;
    }

    return;
};

/**
    Inequality constraints for fsp_mp
*/
void fsp_mpe_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( x[0] + x[1] - 1.0);
    result[1] = ( -x[0]);
    result[2] = ( -x[1]);
    result[3] = ( -0.25*x[0] - 0.25);
    result[4] = ( 0.25*x[0] - 0.75);

    if (grad) {
        grad[0] = 1.0;
        grad[1] = 1.0;
        grad[2] = -1.0;
        grad[3] = 0.0;
        grad[4] = 0.0;
        grad[5] = -1.0;
        grad[6] = -0.25;
        grad[7] = 0.0;
        grad[8] = 0.25;
        grad[9] = 0.0;
    }

    return;
};

/**
    Inequality constraints for opx_mp
*/
void opx_mpe_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( 0.5*x[4]*x[5] - x[3]*x[0] + x[3] + 0.5*x[1]*x[5] - x[1]*x[0] + x[1] - 0.5*x[5] - x[0]*x[2] + x[0] + x[2] - 1.0);
    result[1] = ( -0.5*x[4]*x[5] + x[3]*x[0] - 0.5*x[1]*x[5] + x[1]*x[0] + 0.5*x[5] + x[0]*x[2] - x[0]);
    result[2] = ( -x[1]);
    result[3] = ( -x[3]);
    result[4] = ( -x[2]);
    result[5] = ( -0.5*x[4]*x[5] - x[4]*x[0] + x[4] - 0.5*x[1]*x[5] - x[1]*x[0] + x[1] + 0.5*x[5] + x[0] - 1.0);
    result[6] = ( 0.5*x[4]*x[5] + x[4]*x[0] + 0.5*x[1]*x[5] + x[1]*x[0] - 0.5*x[5] - x[0]);
    result[7] = ( -x[1]);
    result[8] = ( -x[4]);
    result[9] = ( -0.5*x[3] - 0.5*x[2]);
    result[10] = ( 0.5*x[3] + 0.5*x[2] - 1.0);

    if (grad) {
        grad[0] = -x[3] - x[1] - x[2] + 1.0;
        grad[1] = 0.5*x[5] - x[0] + 1.0;
        grad[2] = 1.0 - x[0];
        grad[3] = 1.0 - x[0];
        grad[4] = 0.5*x[5];
        grad[5] = 0.5*x[4] + 0.5*x[1] - 0.5;
        grad[6] = x[3] + x[1] + x[2] - 1.0;
        grad[7] = -0.5*x[5] + x[0];
        grad[8] = x[0];
        grad[9] = x[0];
        grad[10] = -0.5*x[5];
        grad[11] = -0.5*x[4] - 0.5*x[1] + 0.5;
        grad[12] = 0.0;
        grad[13] = -1.0;
        grad[14] = 0.0;
        grad[15] = 0.0;
        grad[16] = 0.0;
        grad[17] = 0.0;
        grad[18] = 0.0;
        grad[19] = 0.0;
        grad[20] = 0.0;
        grad[21] = -1.0;
        grad[22] = 0.0;
        grad[23] = 0.0;
        grad[24] = 0.0;
        grad[25] = 0.0;
        grad[26] = -1.0;
        grad[27] = 0.0;
        grad[28] = 0.0;
        grad[29] = 0.0;
        grad[30] = -x[4] - x[1] + 1.0;
        grad[31] = -0.5*x[5] - x[0] + 1.0;
        grad[32] = 0.0;
        grad[33] = 0.0;
        grad[34] = -0.5*x[5] - x[0] + 1.0;
        grad[35] = -0.5*x[4] - 0.5*x[1] + 0.5;
        grad[36] = x[4] + x[1] - 1.0;
        grad[37] = 0.5*x[5] + x[0];
        grad[38] = 0.0;
        grad[39] = 0.0;
        grad[40] = 0.5*x[5] + x[0];
        grad[41] = 0.5*x[4] + 0.5*x[1] - 0.5;
        grad[42] = 0.0;
        grad[43] = -1.0;
        grad[44] = 0.0;
        grad[45] = 0.0;
        grad[46] = 0.0;
        grad[47] = 0.0;
        grad[48] = 0.0;
        grad[49] = 0.0;
        grad[50] = 0.0;
        grad[51] = 0.0;
        grad[52] = -1.0;
        grad[53] = 0.0;
        grad[54] = 0.0;
        grad[55] = 0.0;
        grad[56] = -0.50;
        grad[57] = -0.50;
        grad[58] = 0.0;
        grad[59] = 0.0;
        grad[60] = 0.0;
        grad[61] = 0.0;
        grad[62] = 0.50;
        grad[63] = 0.50;
        grad[64] = 0.0;
        grad[65] = 0.0;
    }

    return;
};

/**
    Inequality constraints for mu_mp
*/
void mu_mpe_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( x[4] + x[3] - 1.0);
    result[1] = ( -x[3]);
    result[2] = ( -x[4]);
    result[3] = ( -x[0]*x[1] + x[0] + x[1] - 1.0);
    result[4] = ( x[0]*x[1] - x[0]);
    result[5] = ( -x[1]);
    result[6] = ( x[2] - 1.0);
    result[7] = ( -x[2]);
    result[8] = ( 0.5*x[4] + 0.5*x[1] - 1.0);
    result[9] = ( -0.5*x[4] - 0.5*x[1]);

    if (grad) {
        grad[0] = 0.0;
        grad[1] = 0.0;
        grad[2] = 0.0;
        grad[3] = 1.0;
        grad[4] = 1.0;
        grad[5] = 0.0;
        grad[6] = 0.0;
        grad[7] = 0.0;
        grad[8] = -1.0;
        grad[9] = 0.0;
        grad[10] = 0.0;
        grad[11] = 0.0;
        grad[12] = 0.0;
        grad[13] = 0.0;
        grad[14] = -1.0;
        grad[15] = 1.0 - x[1];
        grad[16] = 1.0 - x[0];
        grad[17] = 0.0;
        grad[18] = 0.0;
        grad[19] = 0.0;
        grad[20] = x[1] - 1.0;
        grad[21] = x[0];
        grad[22] = 0.0;
        grad[23] = 0.0;
        grad[24] = 0.0;
        grad[25] = 0.0;
        grad[26] = -1.0;
        grad[27] = 0.0;
        grad[28] = 0.0;
        grad[29] = 0.0;
        grad[30] = 0.0;
        grad[31] = 0.0;
        grad[32] = 1.0;
        grad[33] = 0.0;
        grad[34] = 0.0;
        grad[35] = 0.0;
        grad[36] = 0.0;
        grad[37] = -1.0;
        grad[38] = 0.0;
        grad[39] = 0.0;
        grad[40] = 0.0;
        grad[41] = 0.50;
        grad[42] = 0.0;
        grad[43] = 0.0;
        grad[44] = 0.50;
        grad[45] = 0.0;
        grad[46] = -0.50;
        grad[47] = 0.0;
        grad[48] = 0.0;
        grad[49] = -0.50;
    }

    return;
};

/**
    Inequality constraints for mt_mp
*/
void mt_mpe_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( 0.5*x[0] - 0.5);
    result[1] = ( 0.5*x[1] - x[0]);
    result[2] = ( -0.5*x[1] + 0.5*x[0] - 0.5);
    result[3] = ( -x[1]);
    result[4] = ( x[1] - 1.0);

    if (grad) {
        grad[0] = 0.50;
        grad[1] = 0.0;
        grad[2] = -1.0;
        grad[3] = 0.50;
        grad[4] = 0.50;
        grad[5] = -0.50;
        grad[6] = 0.0;
        grad[7] = -1.0;
        grad[8] = 0.0;
        grad[9] = 1.0;
    }

    return;
};

/**
    Inequality constraints for ma_mp
*/
void ma_mpe_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( x[4] + x[3] - 1.0);
    result[1] = ( -x[3]);
    result[2] = ( -x[4]);
    result[3] = ( -x[0]*x[1] + x[0] + x[1] - 1.0);
    result[4] = ( x[0]*x[1] - x[0]);
    result[5] = ( -x[1]);
    result[6] = ( x[2] - 1.0);
    result[7] = ( -x[2]);
    result[8] = ( 0.5*x[4] + 0.5*x[1] - 1.0);
    result[9] = ( -0.5*x[4] - 0.5*x[1]);

    if (grad) {
        grad[0] = 0.0;
        grad[1] = 0.0;
        grad[2] = 0.0;
        grad[3] = 1.0;
        grad[4] = 1.0;
        grad[5] = 0.0;
        grad[6] = 0.0;
        grad[7] = 0.0;
        grad[8] = -1.0;
        grad[9] = 0.0;
        grad[10] = 0.0;
        grad[11] = 0.0;
        grad[12] = 0.0;
        grad[13] = 0.0;
        grad[14] = -1.0;
        grad[15] = 1.0 - x[1];
        grad[16] = 1.0 - x[0];
        grad[17] = 0.0;
        grad[18] = 0.0;
        grad[19] = 0.0;
        grad[20] = x[1] - 1.0;
        grad[21] = x[0];
        grad[22] = 0.0;
        grad[23] = 0.0;
        grad[24] = 0.0;
        grad[25] = 0.0;
        grad[26] = -1.0;
        grad[27] = 0.0;
        grad[28] = 0.0;
        grad[29] = 0.0;
        grad[30] = 0.0;
        grad[31] = 0.0;
        grad[32] = 1.0;
        grad[33] = 0.0;
        grad[34] = 0.0;
        grad[35] = 0.0;
        grad[36] = 0.0;
        grad[37] = -1.0;
        grad[38] = 0.0;
        grad[39] = 0.0;
        grad[40] = 0.0;
        grad[41] = 0.50;
        grad[42] = 0.0;
        grad[43] = 0.0;
        grad[44] = 0.50;
        grad[45] = 0.0;
        grad[46] = -0.50;
        grad[47] = 0.0;
        grad[48] = 0.0;
        grad[49] = -0.50;
    }

    return;
};

/**
    Inequality constraints for ilm
*/
void ilm_mpe_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( -0.5*x[0] - 0.5*x[1]);
    result[1] = ( -0.5*x[0] + 0.5*x[1]);
    result[2] = ( x[0] - 1.0);
    result[3] = ( -0.5*x[0] + 0.5*x[1]);
    result[4] = ( -0.5*x[0] - 0.5*x[1]);
    result[5] = ( x[0] - 1.0);

    if (grad) {
        grad[0] = -0.500;
        grad[1] = -0.500;
        grad[2] = -0.500;
        grad[3] = 0.500;
        grad[4] = 1.00;
        grad[5] = 0.0;
        grad[6] = -0.500;
        grad[7] = 0.500;
        grad[8] = -0.500;
        grad[9] = -0.500;
        grad[10] = 1.00;
        grad[11] = 0.0;
    }

    return;
};
/**
    Inequality constraints for ilmm_mp
*/
void ilmm_mpe_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( -0.5*x[0] + 0.5*x[1] + 0.5*x[2] - 0.5*x[3]);
    result[1] = ( -0.5*x[0] + 0.5*x[1] + 0.5*x[2] + 0.5*x[3]);
    result[2] = ( - x[1]);
    result[3] = ( - x[2]);
    result[4] = (  x[0] - 1.0);
    result[5] = ( -0.5*x[0] + 0.5*x[1] + 0.5*x[2] + 0.5*x[3]);
    result[6] = ( -0.5*x[0] - 0.5*x[1] - 0.5*x[2] - 0.5*x[3]);

    if (grad) {
        grad[0] = -0.50;
        grad[1] = 0.50;
        grad[2] = 0.50;
        grad[3] = -0.50;
        grad[4] = -0.50;
        grad[5] = 0.50;
        grad[6] = 0.50;
        grad[7] = 0.50;
        grad[8] = 0.0;
        grad[9] = -1.0;
        grad[10] = 0.0;
        grad[11] = 0.0;
        grad[12] = 0.0;
        grad[13] = 0.0;
        grad[14] = -1.0;
        grad[15] = 0.0;
        grad[16] = 1.0;
        grad[17] = 0.0;
        grad[18] = 0.0;
        grad[19] = 0.0;
        grad[20] = -0.50;
        grad[21] = 0.50;
        grad[22] = 0.50;
        grad[23] = 0.50;
        grad[24] = -0.50;
        grad[25] = -0.50;
        grad[26] = -0.50;
        grad[27] = -0.50;
    }

    return;
};

/**
    Inequality constraints for g_mp
*/
void g_mpe_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( - x[2]*x[0] + x[2] - x[0]*x[1] + x[0] + x[1] - 1.0);
    result[1] = ( x[2]*x[0] + x[0]*x[1] - x[0]);
    result[2] = ( - x[2]);
    result[3] = ( - x[1]);
    result[4] = ( - 1.0 + x[3]);
    result[5] = ( - x[3]);

    if (grad) {
        grad[0] = -x[2] - x[1] + 1.0;
        grad[1] = 1.0 - x[0];
        grad[2] = 1.0 - x[0];
        grad[3] = 0.0;
        grad[4] = x[2] + x[1] - 1.0;
        grad[5] = x[0];
        grad[6] = x[0];
        grad[7] = 0.0;
        grad[8] = 0.0;
        grad[9] = 0.0;
        grad[10] = -1.0;
        grad[11] = 0.0;
        grad[12] = 0.0;
        grad[13] = -1.0;
        grad[14] = 0.0;
        grad[15] = 0.0;
        grad[16] = 0.0;
        grad[17] = 0.0;
        grad[18] = 0.0;
        grad[19] = 1.0;
        grad[20] = 0.0;
        grad[21] = 0.0;
        grad[22] = 0.0;
        grad[23] = -1.0;
    }

    return;
};

/**
    Inequality constraints for ep_mp
*/
void ep_mpe_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( -x[0] + x[1]);
    result[1] = ( x[0] - x[1] - 1.0);
    result[2] = ( -x[0] - x[1]);
    result[3] = ( x[0] + x[1] - 1.0);

    if (grad) {
        grad[0] = -1.0;
        grad[1] = 1.0;
        grad[2] = 1.0;
        grad[3] = -1.0;
        grad[4] = -1.0;
        grad[5] = -1.0;
        grad[6] = 1.0;
        grad[7] = 1.0;
    }

    return;
};

/**
    Inequality constraints for ctd_mp
*/
void ctd_mpe_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( x[2] - 1.0);
    result[1] = ( -x[2]);
    result[2] = ( x[1]*x[0] - x[0]);
    result[3] = ( -x[1]*x[0] + x[1] + x[0] - 1.0);
    result[4] = ( -x[1]);

    if (grad) {
        grad[0] = 0.0;
        grad[1] = 0.0;
        grad[2] = 1.0;
        grad[3] = 0.0;
        grad[4] = 0.0;
        grad[5] = -1.0;
        grad[6] = x[1] - 1.0;
        grad[7] = x[0];
        grad[8] = 0.0;
        grad[9] = 1.0 - x[1];
        grad[10] = 1.0 - x[0];
        grad[11] = 0.0;
        grad[12] = 0.0;
        grad[13] = -1.0;
        grad[14] = 0.0;
    }

    return;
}

/**
    Inequality constraints for chl_mp
*/
void chl_mpe_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( x[3]*x[5] - x[3]*x[0] + x[3] - x[5]*x[4] + x[5]*x[1] - x[5] + x[4]*x[0] - x[4] - x[0]*x[1] + x[0] + x[1] - 1.0);
    result[1] = ( -x[3]*x[5] + x[3]*x[0] + x[5]*x[4] - x[5]*x[1] + x[5] - x[4]*x[0] + x[0]*x[1] - x[0]);
    result[2] = ( x[4] - x[1]);
    result[3] = ( -0.25*x[2]*x[6] - 0.25*x[3]*x[5] - x[3]*x[0] + x[3] + 0.25*x[5]*x[4] - 0.25*x[5]*x[1] + 0.25*x[5] - 0.25*x[6]*x[4] - 0.25*x[6]*x[1] + 0.25*x[6] + x[0] - 1.0);
    result[4] = ( -x[3]);
    result[5] = ( 0.25*x[2]*x[6] + 0.25*x[3]*x[5] + x[3]*x[0] - 0.25*x[5]*x[4] + 0.25*x[5]*x[1] - 0.25*x[5] + 0.25*x[6]*x[4] + 0.25*x[6]*x[1] - 0.25*x[6] - x[0]);
    result[6] = ( x[2]*x[6] - x[2]*x[0] + x[2] + x[6]*x[4] + x[6]*x[1] - x[6] - x[4]*x[0] + x[4] - x[0]*x[1] + x[0] + x[1] - 1.0);
    result[7] = ( -x[2]*x[6] + x[2]*x[0] - x[6]*x[4] - x[6]*x[1] + x[6] + x[4]*x[0] + x[0]*x[1] - x[0]);
    result[8] = ( -x[2]);
    result[9] = ( -x[4] - x[1]);
    result[10] = ( 0.5*x[2] + x[1] - 1.0);
    result[11] = ( -0.5*x[2] - x[1]);

    if (grad) {
        grad[0] = -x[3] + x[4] - x[1] + 1.0;
        grad[1] = x[5] - x[0] + 1.0;
        grad[2] = 0.0;
        grad[3] = x[5] - x[0] + 1.0;
        grad[4] = -x[5] + x[0] - 1.0;
        grad[5] = x[3] - x[4] + x[1] - 1.0;
        grad[6] = 0.0;
        grad[7] = x[3] - x[4] + x[1] - 1.0;
        grad[8] = -x[5] + x[0];
        grad[9] = 0.0;
        grad[10] = -x[5] + x[0];
        grad[11] = x[5] - x[0];
        grad[12] = -x[3] + x[4] - x[1] + 1.0;
        grad[13] = 0.0;
        grad[14] = 0.0;
        grad[15] = -1.0;
        grad[16] = 0.0;
        grad[17] = 0.0;
        grad[18] = 1.0;
        grad[19] = 0.0;
        grad[20] = 0.0;
        grad[21] = 1.0 - x[3];
        grad[22] = -0.25*x[5] - 0.25*x[6];
        grad[23] = -0.25*x[6];
        grad[24] = -0.25*x[5] - x[0] + 1.0;
        grad[25] = 0.25*x[5] - 0.25*x[6];
        grad[26] = -0.25*x[3] + 0.25*x[4] - 0.25*x[1] + 0.25;
        grad[27] = -0.25*x[2] - 0.25*x[4] - 0.25*x[1] + 0.25;
        grad[28] = 0.0;
        grad[29] = 0.0;
        grad[30] = 0.0;
        grad[31] = -1.0;
        grad[32] = 0.0;
        grad[33] = 0.0;
        grad[34] = 0.0;
        grad[35] = x[3] - 1.0;
        grad[36] = 0.25*x[5] + 0.25*x[6];
        grad[37] = 0.25*x[6];
        grad[38] = 0.25*x[5] + x[0];
        grad[39] = -0.25*x[5] + 0.25*x[6];
        grad[40] = 0.25*x[3] - 0.25*x[4] + 0.25*x[1] - 0.25;
        grad[41] = 0.25*x[2] + 0.25*x[4] + 0.25*x[1] - 0.25;
        grad[42] = -x[2] - x[4] - x[1] + 1.0;
        grad[43] = x[6] - x[0] + 1.0;
        grad[44] = x[6] - x[0] + 1.0;
        grad[45] = 0.0;
        grad[46] = x[6] - x[0] + 1.0;
        grad[47] = 0.0;
        grad[48] = x[2] + x[4] + x[1] - 1.0;
        grad[49] = x[2] + x[4] + x[1] - 1.0;
        grad[50] = -x[6] + x[0];
        grad[51] = -x[6] + x[0];
        grad[52] = 0.0;
        grad[53] = -x[6] + x[0];
        grad[54] = 0.0;
        grad[55] = -x[2] - x[4] - x[1] + 1.0;
        grad[56] = 0.0;
        grad[57] = 0.0;
        grad[58] = -1.0;
        grad[59] = 0.0;
        grad[60] = 0.0;
        grad[61] = 0.0;
        grad[62] = 0.0;
        grad[63] = 0.0;
        grad[64] = -1.0;
        grad[65] = 0.0;
        grad[66] = 0.0;
        grad[67] = -1.0;
        grad[68] = 0.0;
        grad[69] = 0.0;
        grad[70] = 0.0;
        grad[71] = 1.0;
        grad[72] = 0.50;
        grad[73] = 0.0;
        grad[74] = 0.0;
        grad[75] = 0.0;
        grad[76] = 0.0;
        grad[77] = 0.0;
        grad[78] = -1.0;
        grad[79] = -0.50;
        grad[80] = 0.0;
        grad[81] = 0.0;
        grad[82] = 0.0;
        grad[83] = 0.0;
    }

    return;
};

/**
    Inequality constraints for cd_mp
*/
void cd_mpe_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( x[0]*x[1] - x[0]);
    result[1] = ( -x[0]*x[1] + x[0] + x[1] - 1.0);
    result[2] = ( -x[1]);
    result[3] = ( -x[2]);
    result[4] = ( x[2] - 1.0);

    if (grad) {
        grad[0] = x[1] - 1.0;
        grad[1] = x[0];
        grad[2] = 0.0;
        grad[3] = 1.0 - x[1];
        grad[4] = 1.0 - x[0];
        grad[5] = 0.0;
        grad[6] = 0.0;
        grad[7] = -1.0;
        grad[8] = 0.0;
        grad[9] = 0.0;
        grad[10] = 0.0;
        grad[11] = -1.0;
        grad[12] = 0.0;
        grad[13] = 0.0;
        grad[14] = 1.0;
    }

    return;
};

/**
    Inequality constraints for bi_mp
*/
void bi_mpe_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = (  - x[3]*x[0] + x[3] - 3.0*x[1]*x[0] + x[1] + 2.0/3.0*x[5] - x[4]*x[0] + x[4] - x[0]*x[2] + x[0] + x[2] - 1.0);
    result[1] = (  - x[1]);
    result[2] = (  + x[3]*x[0] + 3.0*x[1]*x[0] - 2.0/3.0*x[5] + x[4]*x[0] + x[0]*x[2] - x[0]);
    result[3] = (  - x[3]);
    result[4] = (  - x[4]);
    result[5] = (  - x[2]);
    result[6] = (  + x[1] - 1.0/3.0*x[5] + x[0] - 1.0);
    result[7] = (  - x[1]);
    result[8] = (  + 1.0/3.0*x[5] - x[0]);
    result[9] = (  + 0.5*x[3] + 0.5*x[2] - 0.5);
    result[10] = ( - 0.5*x[3] - 0.5*x[2] - 0.5);
    result[11] = ( x[4] - 1.0);
    result[12] = ( - x[4]);

    if (grad) {
        grad[0] = -x[3] - 3.0*x[1] - x[4] - x[2] + 1.0;
        grad[1] = 1.0 - 3.0*x[0];
        grad[2] = 1.0 - x[0];
        grad[3] = 1.0 - x[0];
        grad[4] = 1.0 - x[0];
        grad[5] = 2.0/3.0;
        grad[6] = 0.0;
        grad[7] = -1.0;
        grad[8] = 0.0;
        grad[9] = 0.0;
        grad[10] = 0.0;
        grad[11] = 0.0;
        grad[12] = x[3] + 3.0*x[1] + x[4] + x[2] - 1.0;
        grad[13] = 3.0*x[0];
        grad[14] = x[0];
        grad[15] = x[0];
        grad[16] = x[0];
        grad[17] = -2.0/3.0;
        grad[18] = 0.0;
        grad[19] = 0.0;
        grad[20] = 0.0;
        grad[21] = -1.0;
        grad[22] = 0.0;
        grad[23] = 0.0;
        grad[24] = 0.0;
        grad[25] = 0.0;
        grad[26] = 0.0;
        grad[27] = 0.0;
        grad[28] = -1.0;
        grad[29] = 0.0;
        grad[30] = 0.0;
        grad[31] = 0.0;
        grad[32] = -1.0;
        grad[33] = 0.0;
        grad[34] = 0.0;
        grad[35] = 0.0;
        grad[36] = 1.0;
        grad[37] = 1.0;
        grad[38] = 0.0;
        grad[39] = 0.0;
        grad[40] = 0.0;
        grad[41] = -1.0/3.0;
        grad[42] = 0.0;
        grad[43] = -1.0;
        grad[44] = 0.0;
        grad[45] = 0.0;
        grad[46] = 0.0;
        grad[47] = 0.0;
        grad[48] = -1.0;
        grad[49] = 0.0;
        grad[50] = 0.0;
        grad[51] = 0.0;
        grad[52] = 0.0;
        grad[53] = 1.0/3.0;
        grad[54] = 0.0;
        grad[55] = 0.0;
        grad[56] = 0.50;
        grad[57] = 0.50;
        grad[58] = 0.0;
        grad[59] = 0.0;
        grad[60] = 0.0;
        grad[61] = 0.0;
        grad[62] = -0.50;
        grad[63] = -0.50;
        grad[64] = 0.0;
        grad[65] = 0.0;
        grad[66] = 0.0;
        grad[67] = 0.0;
        grad[68] = 0.0;
        grad[69] = 0.0;
        grad[70] = 1.0;
        grad[71] = 0.0;
        grad[72] = 0.0;
        grad[73] = 0.0;
        grad[74] = 0.0;
        grad[75] = 0.0;
        grad[76] = -1.0;
        grad[77] = 0.0;
    }

    return;
};


/**
    Inequality constraints for fl
*/
void fl_mpe_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = (x[0] - 1.0);
    result[1] = (-x[0]);

    if (grad) {
        grad[0] = 1.0;
        grad[1] = -1.0;
    }

    return;
};

/**
    Inequality constraints for occm
*/
void occm_mpe_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = (x[0] + x[1] - 0.5*x[2] - 1.0);
    result[1] = (-x[0] + 0.5*x[2] - 0.25*x[3]);
    result[2] = (-x[1] + 0.25*x[3]);
    result[3] = (x[0] + x[1] + 0.5*x[2] - 1.0);
    result[4] = (-x[0] - 0.5*x[2] + 0.75*x[3]);
    result[5] = (-x[1] - 0.75*x[3]);
    result[6] = (x[0] + x[1] + 0.5*x[2] - 1.0);
    result[7] = (-x[0] - 0.5*x[2] - 0.25*x[3]);
    result[8] = (-x[1] + 0.25*x[3]);

    if (grad) {
        grad[0] = 1.0;
        grad[1] = 1.0;
        grad[2] = -0.50;
        grad[3] = 0.0;
        grad[4] = -1.0;
        grad[5] = 0.0;
        grad[6] = 0.50;
        grad[7] = -0.25;
        grad[8] = 0.0;
        grad[9] = -1.0;
        grad[10] = 0.0;
        grad[11] = 0.25;
        grad[12] = 1.0;
        grad[13] = 1.0;
        grad[14] = 0.50;
        grad[15] = 0.0;
        grad[16] = -1.0;
        grad[17] = 0.0;
        grad[18] = -0.50;
        grad[19] = 0.75;
        grad[20] = 0.0;
        grad[21] = -1.0;
        grad[22] = 0.0;
        grad[23] = -0.75;
        grad[24] = 1.0;
        grad[25] = 1.0;
        grad[26] = 0.50;
        grad[27] = 0.0;
        grad[28] = -1.0;
        grad[29] = 0.0;
        grad[30] = -0.50;
        grad[31] = -0.25;
        grad[32] = 0.0;
        grad[33] = -1.0;
        grad[34] = 0.0;
        grad[35] = 0.25;
    }

    return;
};

/**
    Inequality constraints for po
*/
void po_mpe_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = (x[0] - 1.0);
    result[1] = (-x[0]);

    if (grad) {
        grad[0] = 1.0;
        grad[1] = -1.0;
    }

    return;
};

/**
    Inequality constraints for hb
*/
void hb_mpe_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = (x[3] - 1.0);
    result[1] = (x[3]*x[4] - x[3]);
    result[2] = (-x[3]*x[4]);
    result[3] = (x[0] - x[8] - 1.0);
    result[4] = (-x[0] + x[8]);
    result[5] = (-x[0]*x[1] - x[0]*x[6] - x[0]*x[7] + x[0] + x[1]*x[9] + x[1] + x[6]*x[9] + x[6] + x[7]*x[9] + x[7] - x[9] - 1.0);
    result[6] = (x[0]*x[1] + x[0]*x[6] + x[0]*x[7] - x[0] - x[1]*x[9] - x[6]*x[9] - x[7]*x[9] + x[9]);
    result[7] = (-x[1]);
    result[8] = (-x[6]);
    result[9] = (-x[7]);
    result[10] = (-x[5]);
    result[11] = (-x[0]*x[2] - x[0]*x[5] + x[0] - x[1]*x[9] + x[2] + x[5] - x[6]*x[9] - x[7]*x[9] + 1.5*x[8] + x[9] - 1.0);
    result[12] = (x[0]*x[2] + x[0]*x[5] - x[0] + x[1]*x[9] + x[6]*x[9] + x[7]*x[9] - 1.5*x[8] - x[9]);
    result[13] = (-x[2]);
    result[14] = (0.5*x[1] - 0.5*x[2] + 0.25*x[3] + 0.5*x[6] + 0.5*x[7] - 1.0);
    result[15] = (-0.5*x[1] + 0.5*x[2] - 0.25*x[3] - 0.5*x[6] - 0.5*x[7]);
    result[16] = (x[7] - 1.0);
    result[17] = (-x[7]);

    if (grad) {
        grad[0] = 0.0;
        grad[1] = 0.0;
        grad[2] = 0.0;
        grad[3] = 1.0;
        grad[4] = 0.0;
        grad[5] = 0.0;
        grad[6] = 0.0;
        grad[7] = 0.0;
        grad[8] = 0.0;
        grad[9] = 0.0;
        grad[10] = 0.0;
        grad[11] = 0.0;
        grad[12] = 0.0;
        grad[13] = x[4] - 1.0;
        grad[14] = x[3];
        grad[15] = 0.0;
        grad[16] = 0.0;
        grad[17] = 0.0;
        grad[18] = 0.0;
        grad[19] = 0.0;
        grad[20] = 0.0;
        grad[21] = 0.0;
        grad[22] = 0.0;
        grad[23] = -x[4];
        grad[24] = -x[3];
        grad[25] = 0.0;
        grad[26] = 0.0;
        grad[27] = 0.0;
        grad[28] = 0.0;
        grad[29] = 0.0;
        grad[30] = 1.0;
        grad[31] = 0.0;
        grad[32] = 0.0;
        grad[33] = 0.0;
        grad[34] = 0.0;
        grad[35] = 0.0;
        grad[36] = 0.0;
        grad[37] = 0.0;
        grad[38] = -1.0;
        grad[39] = 0.0;
        grad[40] = -1.0;
        grad[41] = 0.0;
        grad[42] = 0.0;
        grad[43] = 0.0;
        grad[44] = 0.0;
        grad[45] = 0.0;
        grad[46] = 0.0;
        grad[47] = 0.0;
        grad[48] = 1.0;
        grad[49] = 0.0;
        grad[50] = -x[1] - x[6] - x[7] + 1.0;
        grad[51] = -x[0] + x[9] + 1.0;
        grad[52] = 0.0;
        grad[53] = 0.0;
        grad[54] = 0.0;
        grad[55] = 0.0;
        grad[56] = -x[0] + x[9] + 1.0;
        grad[57] = -x[0] + x[9] + 1.0;
        grad[58] = 0.0;
        grad[59] = x[1] + x[6] + x[7] - 1.0;
        grad[60] = x[1] + x[6] + x[7] - 1.0;
        grad[61] = x[0] - x[9];
        grad[62] = 0.0;
        grad[63] = 0.0;
        grad[64] = 0.0;
        grad[65] = 0.0;
        grad[66] = x[0] - x[9];
        grad[67] = x[0] - x[9];
        grad[68] = 0.0;
        grad[69] = -x[1] - x[6] - x[7] + 1.0;
        grad[70] = 0.0;
        grad[71] = -1.0;
        grad[72] = 0.0;
        grad[73] = 0.0;
        grad[74] = 0.0;
        grad[75] = 0.0;
        grad[76] = 0.0;
        grad[77] = 0.0;
        grad[78] = 0.0;
        grad[79] = 0.0;
        grad[80] = 0.0;
        grad[81] = 0.0;
        grad[82] = 0.0;
        grad[83] = 0.0;
        grad[84] = 0.0;
        grad[85] = 0.0;
        grad[86] = -1.0;
        grad[87] = 0.0;
        grad[88] = 0.0;
        grad[89] = 0.0;
        grad[90] = 0.0;
        grad[91] = 0.0;
        grad[92] = 0.0;
        grad[93] = 0.0;
        grad[94] = 0.0;
        grad[95] = 0.0;
        grad[96] = 0.0;
        grad[97] = -1.0;
        grad[98] = 0.0;
        grad[99] = 0.0;
        grad[100] = 0.0;
        grad[101] = 0.0;
        grad[102] = 0.0;
        grad[103] = 0.0;
        grad[104] = 0.0;
        grad[105] = -1.0;
        grad[106] = 0.0;
        grad[107] = 0.0;
        grad[108] = 0.0;
        grad[109] = 0.0;
        grad[110] = -x[2] - x[5] + 1.0;
        grad[111] = -x[9];
        grad[112] = 1.0 - x[0];
        grad[113] = 0.0;
        grad[114] = 0.0;
        grad[115] = 1.0 - x[0];
        grad[116] = -x[9];
        grad[117] = -x[9];
        grad[118] = 1.5;
        grad[119] = -x[1] - x[6] - x[7] + 1.0;
        grad[120] = x[2] + x[5] - 1.0;
        grad[121] = x[9];
        grad[122] = x[0];
        grad[123] = 0.0;
        grad[124] = 0.0;
        grad[125] = x[0];
        grad[126] = x[9];
        grad[127] = x[9];
        grad[128] = -1.5;
        grad[129] = x[1] + x[6] + x[7] - 1.0;
        grad[130] = 0.0;
        grad[131] = 0.0;
        grad[132] = -1.0;
        grad[133] = 0.0;
        grad[134] = 0.0;
        grad[135] = 0.0;
        grad[136] = 0.0;
        grad[137] = 0.0;
        grad[138] = 0.0;
        grad[139] = 0.0;
        grad[140] = 0.0;
        grad[141] = 0.50;
        grad[142] = -0.50;
        grad[143] = 0.25;
        grad[144] = 0.0;
        grad[145] = 0.0;
        grad[146] = 0.50;
        grad[147] = 0.50;
        grad[148] = 0.0;
        grad[149] = 0.0;
        grad[150] = 0.0;
        grad[151] = -0.50;
        grad[152] = 0.50;
        grad[153] = -0.25;
        grad[154] = 0.0;
        grad[155] = 0.0;
        grad[156] = -0.50;
        grad[157] = -0.50;
        grad[158] = 0.0;
        grad[159] = 0.0;
        grad[160] = 0.0;
        grad[161] = 0.0;
        grad[162] = 0.0;
        grad[163] = 0.0;
        grad[164] = 0.0;
        grad[165] = 0.0;
        grad[166] = 0.0;
        grad[167] = 1.0;
        grad[168] = 0.0;
        grad[169] = 0.0;
        grad[170] = 0.0;
        grad[171] = 0.0;
        grad[172] = 0.0;
        grad[173] = 0.0;
        grad[174] = 0.0;
        grad[175] = 0.0;
        grad[176] = 0.0;
        grad[177] = -1.0;
        grad[178] = 0.0;
        grad[179] = 0.0;
    }

    return;
};

/**
    Inequality constraints for aug
*/
void aug_mpe_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = (-x[0]*x[1] - x[0]*x[4] + x[0] + x[1] + 0.5*x[3]*x[5] + 0.5*x[4]*x[5] + x[4] - 0.5*x[5] - 1.0);
    result[1] = (x[0]*x[1] + x[0]*x[4] - x[0] - 0.5*x[3]*x[5] - 0.5*x[4]*x[5] + 0.5*x[5]);
    result[2] = (-x[1] + x[2] - x[4]);
    result[3] = (-x[2]);
    result[4] = (-x[0]*x[3] - x[0]*x[4] + x[0] - 0.5*x[3]*x[5] + x[3] - 0.5*x[4]*x[5] + x[4] + 0.5*x[5] - 1.0);
    result[5] = (x[0]*x[3] + x[0]*x[4] - x[0] + 0.5*x[3]*x[5] + 0.5*x[4]*x[5] - 0.5*x[5]);
    result[6] = (-x[3]);
    result[7] = (-x[4]);
    result[8] = (0.5*x[1] - 0.5*x[6] - 1.0);
    result[9] = (-0.5*x[1] + 0.5*x[6]);
    result[10] = (0.5*x[1] + 0.5*x[6] - 1.0);
    result[11] = (-0.5*x[1] - 0.5*x[6]);

    if (grad) {
        grad[0] = -x[1] - x[4] + 1.0;
        grad[1] = 1.0 - x[0];
        grad[2] = 0.0;
        grad[3] = 0.5*x[5];
        grad[4] = -x[0] + 0.5*x[5] + 1.0;
        grad[5] = 0.5*x[3] + 0.5*x[4] - 0.5;
        grad[6] = 0.0;
        grad[7] = x[1] + x[4] - 1.0;
        grad[8] = x[0];
        grad[9] = 0.0;
        grad[10] = -0.5*x[5];
        grad[11] = x[0] - 0.5*x[5];
        grad[12] = -0.5*x[3] - 0.5*x[4] + 0.5;
        grad[13] = 0.0;
        grad[14] = 0.0;
        grad[15] = -1.0;
        grad[16] = 1.0;
        grad[17] = 0.0;
        grad[18] = -1.0;
        grad[19] = 0.0;
        grad[20] = 0.0;
        grad[21] = 0.0;
        grad[22] = 0.0;
        grad[23] = -1.0;
        grad[24] = 0.0;
        grad[25] = 0.0;
        grad[26] = 0.0;
        grad[27] = 0.0;
        grad[28] = -x[3] - x[4] + 1.0;
        grad[29] = 0.0;
        grad[30] = 0.0;
        grad[31] = -x[0] - 0.5*x[5] + 1.0;
        grad[32] = -x[0] - 0.5*x[5] + 1.0;
        grad[33] = -0.5*x[3] - 0.5*x[4] + 0.5;
        grad[34] = 0.0;
        grad[35] = x[3] + x[4] - 1.0;
        grad[36] = 0.0;
        grad[37] = 0.0;
        grad[38] = x[0] + 0.5*x[5];
        grad[39] = x[0] + 0.5*x[5];
        grad[40] = 0.5*x[3] + 0.5*x[4] - 0.5;
        grad[41] = 0.0;
        grad[42] = 0.0;
        grad[43] = 0.0;
        grad[44] = 0.0;
        grad[45] = -1.0;
        grad[46] = 0.0;
        grad[47] = 0.0;
        grad[48] = 0.0;
        grad[49] = 0.0;
        grad[50] = 0.0;
        grad[51] = 0.0;
        grad[52] = 0.0;
        grad[53] = -1.0;
        grad[54] = 0.0;
        grad[55] = 0.0;
        grad[56] = 0.0;
        grad[57] = 0.50;
        grad[58] = 0.0;
        grad[59] = 0.0;
        grad[60] = 0.0;
        grad[61] = 0.0;
        grad[62] = -0.50;
        grad[63] = 0.0;
        grad[64] = -0.50;
        grad[65] = 0.0;
        grad[66] = 0.0;
        grad[67] = 0.0;
        grad[68] = 0.0;
        grad[69] = 0.50;
        grad[70] = 0.0;
        grad[71] = 0.50;
        grad[72] = 0.0;
        grad[73] = 0.0;
        grad[74] = 0.0;
        grad[75] = 0.0;
        grad[76] = 0.50;
        grad[77] = 0.0;
        grad[78] = -0.50;
        grad[79] = 0.0;
        grad[80] = 0.0;
        grad[81] = 0.0;
        grad[82] = 0.0;
        grad[83] = -0.50;
    }

    return;
};

/**
    Inequality constraints for dio
*/
void dio_mpe_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = (-x[0]*x[1] + x[0]*x[3] + x[0] + x[1]*x[5] + x[1] + x[3]*x[5] - x[3] - x[5] - 1.0);
    result[1] = (x[0]*x[1] - x[0]*x[3] - x[0] - x[1]*x[5] - x[3]*x[5] + x[5]);
    result[2] = (-x[1]*x[2] + x[4]);
    result[3] = (x[1]*x[2] - x[1] + x[3] - x[4]);
    result[4] = (-x[0]*x[1] - x[0]*x[3] + x[0] - x[1]*x[5] + x[1] - x[3]*x[5] + x[3] + x[5] - 1.0);
    result[5] = (x[0]*x[1] + x[0]*x[3] - x[0] + x[1]*x[5] + x[3]*x[5] - x[5]);
    result[6] = (-x[1]*x[2] - x[4]);
    result[7] = (x[1]*x[2] - x[1] - x[3] + x[4]);
    result[8] = (-x[1] + x[3]);
    result[9] = (x[1] - x[3] - 1.0);
    result[10] = (-x[1] - x[3]);
    result[11] = (x[1] + x[3] - 1.0);

    if (grad) {
        grad[0] = -x[1] + x[3] + 1.0;
        grad[1] = -x[0] + x[5] + 1.0;
        grad[2] = 0.0;
        grad[3] = x[0] + x[5] - 1.0;
        grad[4] = 0.0;
        grad[5] = x[1] + x[3] - 1.0;
        grad[6] = x[1] - x[3] - 1.0;
        grad[7] = x[0] - x[5];
        grad[8] = 0.0;
        grad[9] = -x[0] - x[5];
        grad[10] = 0.0;
        grad[11] = -x[1] - x[3] + 1.0;
        grad[12] = 0.0;
        grad[13] = -x[2];
        grad[14] = -x[1];
        grad[15] = 0.0;
        grad[16] = 1.0;
        grad[17] = 0.0;
        grad[18] = 0.0;
        grad[19] = x[2] - 1.0;
        grad[20] = x[1];
        grad[21] = 1.0;
        grad[22] = -1.0;
        grad[23] = 0.0;
        grad[24] = -x[1] - x[3] + 1.0;
        grad[25] = -x[0] - x[5] + 1.0;
        grad[26] = 0.0;
        grad[27] = -x[0] - x[5] + 1.0;
        grad[28] = 0.0;
        grad[29] = -x[1] - x[3] + 1.0;
        grad[30] = x[1] + x[3] - 1.0;
        grad[31] = x[0] + x[5];
        grad[32] = 0.0;
        grad[33] = x[0] + x[5];
        grad[34] = 0.0;
        grad[35] = x[1] + x[3] - 1.0;
        grad[36] = 0.0;
        grad[37] = -x[2];
        grad[38] = -x[1];
        grad[39] = 0.0;
        grad[40] = -1.0;
        grad[41] = 0.0;
        grad[42] = 0.0;
        grad[43] = x[2] - 1.0;
        grad[44] = x[1];
        grad[45] = -1.0;
        grad[46] = 1.0;
        grad[47] = 0.0;
        grad[48] = 0.0;
        grad[49] = -1.0;
        grad[50] = 0.0;
        grad[51] = 1.0;
        grad[52] = 0.0;
        grad[53] = 0.0;
        grad[54] = 0.0;
        grad[55] = 1.0;
        grad[56] = 0.0;
        grad[57] = -1.0;
        grad[58] = 0.0;
        grad[59] = 0.0;
        grad[60] = 0.0;
        grad[61] = -1.0;
        grad[62] = 0.0;
        grad[63] = -1.0;
        grad[64] = 0.0;
        grad[65] = 0.0;
        grad[66] = 0.0;
        grad[67] = 1.0;
        grad[68] = 0.0;
        grad[69] = 1.0;
        grad[70] = 0.0;
        grad[71] = 0.0;
    }

    return;
};


SS_ref NLopt_opt_mpe_st_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_mpe_st, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, st_mpe_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};

SS_ref NLopt_opt_mpe_sp_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_mpe_sp, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, sp_mpe_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};

SS_ref NLopt_opt_mpe_sa_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_mpe_sa, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, sa_mpe_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};

SS_ref NLopt_opt_mpe_fsp_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_mpe_fsp, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, fsp_mpe_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};

SS_ref NLopt_opt_mpe_opx_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_mpe_opx, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, opx_mpe_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};

SS_ref NLopt_opt_mpe_mu_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_mpe_mu, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, mu_mpe_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};

SS_ref NLopt_opt_mpe_mt_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_mpe_mt, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, mt_mpe_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};

SS_ref NLopt_opt_mpe_ma_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_mpe_ma, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, ma_mpe_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};
SS_ref NLopt_opt_mpe_ilm_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_mpe_ilm, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, ilm_mpe_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    
    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};
SS_ref NLopt_opt_mpe_ilmm_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_mpe_ilmm, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, ilmm_mpe_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};

SS_ref NLopt_opt_mpe_g_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_mpe_g, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, g_mpe_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};

SS_ref NLopt_opt_mpe_ep_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_mpe_ep, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, ep_mpe_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};

SS_ref NLopt_opt_mpe_ctd_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_mpe_ctd, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, ctd_mpe_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};

SS_ref NLopt_opt_mpe_chl_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_mpe_chl, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, chl_mpe_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};

SS_ref NLopt_opt_mpe_cd_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_mpe_cd, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, cd_mpe_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};

SS_ref NLopt_opt_mpe_bi_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_mpe_bi, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, bi_mpe_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};


SS_ref NLopt_opt_mpe_liq_function(global_variable gv, SS_ref SS_ref_db){

    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_mpe_liq, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, liq_mpe_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};

SS_ref NLopt_opt_mpe_fl_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_mpe_fl, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, fl_mpe_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};

SS_ref NLopt_opt_mpe_occm_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_mpe_occm, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, occm_mpe_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};

SS_ref NLopt_opt_mpe_po_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_mpe_po, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, po_mpe_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};
SS_ref NLopt_opt_mpe_hb_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_mpe_hb, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, hb_mpe_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};

SS_ref NLopt_opt_mpe_aug_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_mpe_aug, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, aug_mpe_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};

SS_ref NLopt_opt_mpe_dio_function(global_variable gv, SS_ref SS_ref_db){
    
    int    n_em     = SS_ref_db.n_em;
    unsigned int n  = SS_ref_db.n_xeos;
    unsigned int m  = SS_ref_db.n_sf;
    
    double *x  = SS_ref_db.iguess; 
    
    for (int i = 0; i < (SS_ref_db.n_xeos); i++){
       SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
       SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_mpe_dio, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, dio_mpe_c, NULL, NULL);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
       SS_ref_db.xeos[i] = x[i];
    }
    
    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);
    
    return SS_ref_db;
};

/**************************************************************************************/
/**************************************************************************************/
/**************************************************************************************/
/**************************************************************************************/

/** 
  attributes the right solution phase to the solution phase array and calculates xi
*/
void TC_mp_NLopt_opt_init(	        NLopt_type 			*NLopt_opt,
									global_variable 	 gv				){	
						 
	for (int iss = 0; iss < gv.len_ss; iss++){

		if      (strcmp( gv.SS_list[iss], "liq")   == 0 ){
			NLopt_opt[iss]  = NLopt_opt_mp_liq_function; 		}
		else if (strcmp( gv.SS_list[iss], "fsp") == 0){
			NLopt_opt[iss]  = NLopt_opt_mp_fsp_function; 		}
		else if (strcmp( gv.SS_list[iss], "bi")    == 0){
			NLopt_opt[iss]  = NLopt_opt_mp_bi_function; 		}
		else if (strcmp( gv.SS_list[iss], "g")     == 0){
			NLopt_opt[iss]  = NLopt_opt_mp_g_function; 			}
		else if (strcmp( gv.SS_list[iss], "ep")    == 0){
			NLopt_opt[iss]  = NLopt_opt_mp_ep_function; 		}
		else if (strcmp( gv.SS_list[iss], "ma")    == 0){
			NLopt_opt[iss]  = NLopt_opt_mp_ma_function; 		}
		else if (strcmp( gv.SS_list[iss], "mu")    == 0){
			NLopt_opt[iss]  = NLopt_opt_mp_mu_function; 		}
		else if (strcmp( gv.SS_list[iss], "opx")   == 0){
			NLopt_opt[iss]  = NLopt_opt_mp_opx_function; 		}
		else if (strcmp( gv.SS_list[iss], "sa")    == 0){
			NLopt_opt[iss]  = NLopt_opt_mp_sa_function; 		}
		else if (strcmp( gv.SS_list[iss], "cd")    == 0){
			NLopt_opt[iss]  = NLopt_opt_mp_cd_function; 		}
		else if (strcmp( gv.SS_list[iss], "st")    == 0){
			NLopt_opt[iss]  = NLopt_opt_mp_st_function; 		}
		else if (strcmp( gv.SS_list[iss], "chl")   == 0){
			NLopt_opt[iss]  = NLopt_opt_mp_chl_function; 		}
		else if (strcmp( gv.SS_list[iss], "ctd")   == 0){
			NLopt_opt[iss]  = NLopt_opt_mp_ctd_function; 		}
		else if (strcmp( gv.SS_list[iss], "sp")    == 0){
			NLopt_opt[iss]  = NLopt_opt_mp_sp_function; 		}
		else if (strcmp( gv.SS_list[iss], "ilm")   == 0){
			NLopt_opt[iss]  = NLopt_opt_mp_ilm_function; 		}
		else if (strcmp( gv.SS_list[iss], "ilmm")   == 0){
			NLopt_opt[iss]  = NLopt_opt_mp_ilmm_function; 		}
		else if (strcmp( gv.SS_list[iss], "mt")    == 0){
			NLopt_opt[iss]  = NLopt_opt_mp_mt_function; 		}
		else if (strcmp( gv.SS_list[iss], "aq17")  == 0){
			NLopt_opt[iss]  = NLopt_opt_aq17_function; 		    }
		else{
			printf("\nsolid solution '%s' is not in the database, cannot be initiated\n", gv.SS_list[iss]);	
		}	
	};			
}

void TC_mb_NLopt_opt_init(	        NLopt_type 			*NLopt_opt,
									global_variable 	 gv				){	
						 
	for (int iss = 0; iss < gv.len_ss; iss++){
        if (strcmp( gv.SS_list[iss], "liq")  == 0){
            NLopt_opt[iss]  = NLopt_opt_mb_liq_function;        }
        else if (strcmp( gv.SS_list[iss], "hb")  == 0){
            NLopt_opt[iss]  = NLopt_opt_mb_hb_function;         }
        else if (strcmp( gv.SS_list[iss], "aug")  == 0){
            NLopt_opt[iss]  = NLopt_opt_mb_aug_function;        }
        else if (strcmp( gv.SS_list[iss], "dio")  == 0){
            NLopt_opt[iss]  = NLopt_opt_mb_dio_function;        }
        else if (strcmp( gv.SS_list[iss], "opx")  == 0){
            NLopt_opt[iss]  = NLopt_opt_mb_opx_function;        }
        else if (strcmp( gv.SS_list[iss], "g")  == 0){
            NLopt_opt[iss]  = NLopt_opt_mb_g_function;          }
        else if (strcmp( gv.SS_list[iss], "ol")  == 0){
            NLopt_opt[iss]  = NLopt_opt_mb_ol_function;         }
        else if (strcmp( gv.SS_list[iss], "fsp")  == 0){
            NLopt_opt[iss]  = NLopt_opt_mb_fsp_function;        }
        else if (strcmp( gv.SS_list[iss], "abc")  == 0){
            NLopt_opt[iss]  = NLopt_opt_mb_abc_function;        }
        else if (strcmp( gv.SS_list[iss], "k4tr")  == 0){
            NLopt_opt[iss]  = NLopt_opt_mb_k4tr_function;       }
        else if (strcmp( gv.SS_list[iss], "sp")  == 0){
            NLopt_opt[iss]  = NLopt_opt_mb_sp_function;         }
        else if (strcmp( gv.SS_list[iss], "spl")  == 0){
            NLopt_opt[iss]  = NLopt_opt_mb_spl_function;        }
        else if (strcmp( gv.SS_list[iss], "ilm")  == 0){
            NLopt_opt[iss]  = NLopt_opt_mb_ilm_function;        }
        else if (strcmp( gv.SS_list[iss], "ilmm")  == 0){
            NLopt_opt[iss]  = NLopt_opt_mb_ilmm_function;       }
        else if (strcmp( gv.SS_list[iss], "ep")  == 0){
            NLopt_opt[iss]  = NLopt_opt_mb_ep_function;         }
        else if (strcmp( gv.SS_list[iss], "bi")  == 0){
            NLopt_opt[iss]  = NLopt_opt_mb_bi_function;         }
        else if (strcmp( gv.SS_list[iss], "mu")  == 0){
            NLopt_opt[iss]  = NLopt_opt_mb_mu_function;         }
        else if (strcmp( gv.SS_list[iss], "chl")  == 0){
            NLopt_opt[iss]  = NLopt_opt_mb_chl_function;        }
		else{
			printf("\nsolid solution '%s' is not in the database, cannot be initiated\n", gv.SS_list[iss]);	
		}	
	};				
}

void TC_ig_NLopt_opt_init(	        NLopt_type 			*NLopt_opt,
									global_variable 	 gv				){	
						 
	for (int iss = 0; iss < gv.len_ss; iss++){

		if      (strcmp( gv.SS_list[iss], "bi")  == 0 ){
			NLopt_opt[iss]  = NLopt_opt_ig_bi_function; 		}
		else if (strcmp( gv.SS_list[iss], "fper")  == 0){
			NLopt_opt[iss]  = NLopt_opt_ig_fper_function; 		}
		else if (strcmp( gv.SS_list[iss], "cd")  == 0){
			NLopt_opt[iss]  = NLopt_opt_ig_cd_function; 		}
		else if (strcmp( gv.SS_list[iss], "cpx") == 0){
			NLopt_opt[iss]  = NLopt_opt_ig_cpx_function; 		}
		else if (strcmp( gv.SS_list[iss], "ep")  == 0){
			NLopt_opt[iss]  = NLopt_opt_ig_ep_function; 		}
		else if (strcmp( gv.SS_list[iss], "fl")  == 0){
			NLopt_opt[iss]  = NLopt_opt_ig_fl_function; 		}
		else if (strcmp( gv.SS_list[iss], "g")   == 0){
			NLopt_opt[iss]  = NLopt_opt_ig_g_function; 		    }
		else if (strcmp( gv.SS_list[iss], "hb")  == 0){
			NLopt_opt[iss]  = NLopt_opt_ig_hb_function; 		}
		else if (strcmp( gv.SS_list[iss], "ilm") == 0){
			NLopt_opt[iss]  = NLopt_opt_ig_ilm_function; 		}
		else if (strcmp( gv.SS_list[iss], "liq") == 0){
			NLopt_opt[iss]  = NLopt_opt_ig_liq_function; 		}
		else if (strcmp( gv.SS_list[iss], "mu")  == 0){
			NLopt_opt[iss]  = NLopt_opt_ig_mu_function; 		}
		else if (strcmp( gv.SS_list[iss], "ol")  == 0){
			NLopt_opt[iss]  = NLopt_opt_ig_ol_function; 		}
		else if (strcmp( gv.SS_list[iss], "opx") == 0){
			NLopt_opt[iss]  = NLopt_opt_ig_opx_function; 		}
		else if (strcmp( gv.SS_list[iss], "fsp") == 0){
			NLopt_opt[iss]  = NLopt_opt_ig_fsp_function; 		}
		else if (strcmp( gv.SS_list[iss], "spl") == 0){
			NLopt_opt[iss]  = NLopt_opt_ig_spl_function; 		}
		else{
			printf("\nsolid solution '%s' is not in the database, cannot be initiated\n", gv.SS_list[iss]);	
		}	
	};				
}


void TC_igad_NLopt_opt_init(	        NLopt_type 			*NLopt_opt,
									global_variable 	 gv				){	
						 
	for (int iss = 0; iss < gv.len_ss; iss++){

		if (strcmp( gv.SS_list[iss], "cpx") == 0){
			NLopt_opt[iss]  = NLopt_opt_igad_cpx_function; 		}
		else if (strcmp( gv.SS_list[iss], "g")   == 0){
			NLopt_opt[iss]  = NLopt_opt_igad_g_function; 		   }
		else if (strcmp( gv.SS_list[iss], "ilm") == 0){
			NLopt_opt[iss]  = NLopt_opt_igad_ilm_function; 		}
		else if (strcmp( gv.SS_list[iss], "liq") == 0){
			NLopt_opt[iss]  = NLopt_opt_igad_liq_function; 		}
		else if (strcmp( gv.SS_list[iss], "ol")  == 0){
			NLopt_opt[iss]  = NLopt_opt_igad_ol_function; 		}
		else if (strcmp( gv.SS_list[iss], "opx") == 0){
			NLopt_opt[iss]  = NLopt_opt_igad_opx_function; 		}
		else if (strcmp( gv.SS_list[iss], "fsp") == 0){
			NLopt_opt[iss]  = NLopt_opt_igad_fsp_function; 		}
		else if (strcmp( gv.SS_list[iss], "spl") == 0){
			NLopt_opt[iss]  = NLopt_opt_igad_spl_function; 		}
		else if (strcmp( gv.SS_list[iss], "nph") == 0){
			NLopt_opt[iss]  = NLopt_opt_igad_nph_function; 	}
		else if (strcmp( gv.SS_list[iss], "lct") == 0){
			NLopt_opt[iss]  = NLopt_opt_igad_lct_function; 		}
		else if (strcmp( gv.SS_list[iss], "kals") == 0){
			NLopt_opt[iss]  = NLopt_opt_igad_kals_function; 	}
		else if (strcmp( gv.SS_list[iss], "mel") == 0){
			NLopt_opt[iss]  = NLopt_opt_igad_mel_function; 		}
		else{
			printf("\nsolid solution '%s' is not in the database, cannot be initiated\n", gv.SS_list[iss]);	
		}	
	};				
}

void TC_um_NLopt_opt_init(	        NLopt_type 			*NLopt_opt,
									global_variable 	 gv				){	
						 
	for (int iss = 0; iss < gv.len_ss; iss++){

		if      (strcmp( gv.SS_list[iss], "fl")  == 0 ){
			NLopt_opt[iss]  = NLopt_opt_um_fluid_function; 		}
		else if (strcmp( gv.SS_list[iss], "ol")  == 0){
			NLopt_opt[iss]  = NLopt_opt_um_ol_function; 		}
		else if (strcmp( gv.SS_list[iss], "br") == 0){
			NLopt_opt[iss]  = NLopt_opt_um_br_function; 		}
		else if (strcmp( gv.SS_list[iss], "ch")  == 0){
			NLopt_opt[iss]  = NLopt_opt_um_ch_function; 		}
		else if (strcmp( gv.SS_list[iss], "atg")  == 0){
			NLopt_opt[iss]  = NLopt_opt_um_atg_function; 		}
		else if (strcmp( gv.SS_list[iss], "g")   == 0){
			NLopt_opt[iss]  = NLopt_opt_um_g_function; 		    }
		else if (strcmp( gv.SS_list[iss], "ta")  == 0){
			NLopt_opt[iss]  = NLopt_opt_um_ta_function; 		}
		else if (strcmp( gv.SS_list[iss], "chl") == 0){
			NLopt_opt[iss]  = NLopt_opt_um_chl_function; 		}
		else if (strcmp( gv.SS_list[iss], "anth") == 0){
			NLopt_opt[iss]  = NLopt_opt_um_anth_function; 		}
		else if (strcmp( gv.SS_list[iss], "spi")  == 0){
			NLopt_opt[iss]  = NLopt_opt_um_spi_function; 		}
		else if (strcmp( gv.SS_list[iss], "opx") == 0){
			NLopt_opt[iss]  = NLopt_opt_um_opx_function; 		}
		else if (strcmp( gv.SS_list[iss], "po") == 0){
			NLopt_opt[iss]  = NLopt_opt_um_po_function; 		}
		else{
			printf("\nsolid solution '%s' is not in the database, cannot be initiated\n", gv.SS_list[iss]);	
		}	
	};						
}


void TC_um_ext_NLopt_opt_init(	    NLopt_type 			*NLopt_opt,
									global_variable 	 gv				){	
						 
	for (int iss = 0; iss < gv.len_ss; iss++){

		if      (strcmp( gv.SS_list[iss], "fl")  == 0 ){
			NLopt_opt[iss]  = NLopt_opt_um_fluid_function; 		}
		else if (strcmp( gv.SS_list[iss], "ol")  == 0){
			NLopt_opt[iss]  = NLopt_opt_um_ol_function; 		}
		else if (strcmp( gv.SS_list[iss], "br") == 0){
			NLopt_opt[iss]  = NLopt_opt_um_br_function; 		}
		else if (strcmp( gv.SS_list[iss], "ch")  == 0){
			NLopt_opt[iss]  = NLopt_opt_um_ch_function; 		}
		else if (strcmp( gv.SS_list[iss], "atg")  == 0){
			NLopt_opt[iss]  = NLopt_opt_um_atg_function; 		}
		else if (strcmp( gv.SS_list[iss], "g")   == 0){
			NLopt_opt[iss]  = NLopt_opt_um_g_function; 		    }
		else if (strcmp( gv.SS_list[iss], "ta")  == 0){
			NLopt_opt[iss]  = NLopt_opt_um_ta_function; 		}
		else if (strcmp( gv.SS_list[iss], "chl") == 0){
			NLopt_opt[iss]  = NLopt_opt_um_chl_function; 		}
		else if (strcmp( gv.SS_list[iss], "anth") == 0){
			NLopt_opt[iss]  = NLopt_opt_um_anth_function; 		}
		else if (strcmp( gv.SS_list[iss], "spi")  == 0){
			NLopt_opt[iss]  = NLopt_opt_um_spi_function; 		}
		else if (strcmp( gv.SS_list[iss], "opx") == 0){
			NLopt_opt[iss]  = NLopt_opt_um_opx_function; 		}
		else if (strcmp( gv.SS_list[iss], "po") == 0){
			NLopt_opt[iss]  = NLopt_opt_um_po_function; 		}
		else if (strcmp( gv.SS_list[iss], "pl4tr")  == 0){
			NLopt_opt[iss]  = NLopt_opt_ume_pl4tr_function; 	}
		else if (strcmp( gv.SS_list[iss], "hb") == 0){
			NLopt_opt[iss]  = NLopt_opt_ume_hb_function; 		}
		else if (strcmp( gv.SS_list[iss], "aug") == 0){
			NLopt_opt[iss]  = NLopt_opt_ume_aug_function; 		}
		else{
			printf("\nsolid solution '%s' is not in the database, cannot be initiated\n", gv.SS_list[iss]);	
		}	
	};						
}


void TC_mtl_NLopt_opt_init(	        NLopt_type 			*NLopt_opt,
									global_variable 	 gv				){	
						 
	for (int iss = 0; iss < gv.len_ss; iss++){

		if      (strcmp( gv.SS_list[iss], "g")  == 0 ){
			NLopt_opt[iss]  = NLopt_opt_mtl_g_function; 		}
		else if (strcmp( gv.SS_list[iss], "fp")  == 0){
			NLopt_opt[iss]  = NLopt_opt_mtl_fp_function; 		}
		else if (strcmp( gv.SS_list[iss], "mpv") == 0){
			NLopt_opt[iss]  = NLopt_opt_mtl_mpv_function; 		}
		else if (strcmp( gv.SS_list[iss], "cpv") == 0){
			NLopt_opt[iss]  = NLopt_opt_mtl_cpv_function; 		}
		else if (strcmp( gv.SS_list[iss], "crn")  == 0){
			NLopt_opt[iss]  = NLopt_opt_mtl_crn_function; 		}
		else if (strcmp( gv.SS_list[iss], "cf")  == 0){
			NLopt_opt[iss]  = NLopt_opt_mtl_cf_function; 		}
		else if (strcmp( gv.SS_list[iss], "nal")   == 0){
			NLopt_opt[iss]  = NLopt_opt_mtl_nal_function; 		    }
		else if (strcmp( gv.SS_list[iss], "aki")  == 0){
			NLopt_opt[iss]  = NLopt_opt_mtl_aki_function; 		}
		else if (strcmp( gv.SS_list[iss], "ol") == 0){
			NLopt_opt[iss]  = NLopt_opt_mtl_ol_function; 		}
		else if (strcmp( gv.SS_list[iss], "wad") == 0){
			NLopt_opt[iss]  = NLopt_opt_mtl_wad_function; 		}
		else if (strcmp( gv.SS_list[iss], "ring")  == 0){
			NLopt_opt[iss]  = NLopt_opt_mtl_ring_function; 		}
		else if (strcmp( gv.SS_list[iss], "cpx") == 0){
			NLopt_opt[iss]  = NLopt_opt_mtl_cpx_function; 		}
		else if (strcmp( gv.SS_list[iss], "opx") == 0){
			NLopt_opt[iss]  = NLopt_opt_mtl_opx_function; 		}
		else if (strcmp( gv.SS_list[iss], "hpx")  == 0){
			NLopt_opt[iss]  = NLopt_opt_mtl_hpx_function; 	}
		else{
			printf("\nsolid solution '%s' is not in the database, cannot be initiated\n", gv.SS_list[iss]);	
		}	
	};						
}


/** 
  attributes the right solution phase to the solution phase array and calculates xi
*/
void TC_mpe_NLopt_opt_init(	        NLopt_type 			*NLopt_opt,
									global_variable 	 gv				){	
						 
	for (int iss = 0; iss < gv.len_ss; iss++){

		if      (strcmp( gv.SS_list[iss], "liq")   == 0 ){
			NLopt_opt[iss]  = NLopt_opt_mpe_liq_function; 		}
		else if (strcmp( gv.SS_list[iss], "fsp") == 0){
			NLopt_opt[iss]  = NLopt_opt_mpe_fsp_function; 		}
		else if (strcmp( gv.SS_list[iss], "bi")    == 0){
			NLopt_opt[iss]  = NLopt_opt_mpe_bi_function; 		}
		else if (strcmp( gv.SS_list[iss], "g")     == 0){
			NLopt_opt[iss]  = NLopt_opt_mpe_g_function; 			}
		else if (strcmp( gv.SS_list[iss], "ep")    == 0){
			NLopt_opt[iss]  = NLopt_opt_mpe_ep_function; 		}
		else if (strcmp( gv.SS_list[iss], "ma")    == 0){
			NLopt_opt[iss]  = NLopt_opt_mpe_ma_function; 		}
		else if (strcmp( gv.SS_list[iss], "mu")    == 0){
			NLopt_opt[iss]  = NLopt_opt_mpe_mu_function; 		}
		else if (strcmp( gv.SS_list[iss], "opx")   == 0){
			NLopt_opt[iss]  = NLopt_opt_mpe_opx_function; 		}
		else if (strcmp( gv.SS_list[iss], "sa")    == 0){
			NLopt_opt[iss]  = NLopt_opt_mpe_sa_function; 		}
		else if (strcmp( gv.SS_list[iss], "cd")    == 0){
			NLopt_opt[iss]  = NLopt_opt_mpe_cd_function; 		}
		else if (strcmp( gv.SS_list[iss], "st")    == 0){
			NLopt_opt[iss]  = NLopt_opt_mpe_st_function; 		}
		else if (strcmp( gv.SS_list[iss], "chl")   == 0){
			NLopt_opt[iss]  = NLopt_opt_mpe_chl_function; 		}
		else if (strcmp( gv.SS_list[iss], "ctd")   == 0){
			NLopt_opt[iss]  = NLopt_opt_mpe_ctd_function; 		}
		else if (strcmp( gv.SS_list[iss], "sp")    == 0){
			NLopt_opt[iss]  = NLopt_opt_mpe_sp_function; 		}
		else if (strcmp( gv.SS_list[iss], "ilm")   == 0){
			NLopt_opt[iss]  = NLopt_opt_mpe_ilm_function; 		}
		else if (strcmp( gv.SS_list[iss], "ilmm")   == 0){
			NLopt_opt[iss]  = NLopt_opt_mpe_ilmm_function; 		}
		else if (strcmp( gv.SS_list[iss], "mt")    == 0){
			NLopt_opt[iss]  = NLopt_opt_mpe_mt_function; 		}
		else if (strcmp( gv.SS_list[iss], "fl")   == 0){
			NLopt_opt[iss]  = NLopt_opt_mpe_fl_function; 		}
		else if (strcmp( gv.SS_list[iss], "occm")   == 0){
			NLopt_opt[iss]  = NLopt_opt_mpe_occm_function; 		}
		else if (strcmp( gv.SS_list[iss], "aug")    == 0){
			NLopt_opt[iss]  = NLopt_opt_mpe_aug_function; 		}
		else if (strcmp( gv.SS_list[iss], "dio")   == 0){
			NLopt_opt[iss]  = NLopt_opt_mpe_dio_function; 		}
		else if (strcmp( gv.SS_list[iss], "hb")   == 0){
			NLopt_opt[iss]  = NLopt_opt_mpe_hb_function; 		}
		else if (strcmp( gv.SS_list[iss], "po")    == 0){
			NLopt_opt[iss]  = NLopt_opt_mpe_po_function; 		}
		else{
			printf("\nsolid solution '%s' is not in the database, cannot be initiated\n", gv.SS_list[iss]);	
		}	
	};			
}


void TC_NLopt_opt_init(	        	NLopt_type 			*NLopt_opt,
									global_variable 	 gv				){


	if (gv.EM_database == 0){				// metapelite database //
		TC_mp_NLopt_opt_init(	 				NLopt_opt,
												gv							);
	}
	if (gv.EM_database == 1){				// metabasite database //
		TC_mb_NLopt_opt_init(	 				NLopt_opt,
												gv							);
	}
	else if (gv.EM_database == 2){			// igneous database //
		TC_ig_NLopt_opt_init(	 				NLopt_opt,
												gv							);
	}
	else if (gv.EM_database == 3){			// igneous alkali dry database //
		TC_igad_NLopt_opt_init(	 				NLopt_opt,
												gv							);
	}
	else if (gv.EM_database == 4){			// ultramafic database //
		TC_um_NLopt_opt_init(	 				NLopt_opt,
												gv							);
	}
	else if (gv.EM_database == 5){			// ultramafic database //
		TC_um_ext_NLopt_opt_init(	 			NLopt_opt,
												gv							);
	}
	else if (gv.EM_database == 6){			// ultramafic database //
		TC_mtl_NLopt_opt_init(	 			NLopt_opt,
												gv							);
	}
	else if (gv.EM_database == 7){			// ultramafic database //
		TC_mpe_NLopt_opt_init(	 			NLopt_opt,
												gv							);
	}
}