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
#include "MAGEMin.h"
#include "gss_function.h"			// order of header file declaration is important
#include "objective_functions.h"
#include "NLopt_opt_function.h"
#include "toolkit.h"

#define nEl 11						// max number of non-zeros compoenents
#define eps_sf -1e-10				// eps to shift site fraction from zero


/** 
  local minimization for clinopyroxene
*/
void cpx_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( eps_sf + -x[8]*x[4] - x[8]*x[0] + x[8] - x[3]*x[4] - x[3]*x[0] + x[3] + x[4]*x[7] - x[4]*x[1] + x[4] + x[7]*x[0] - x[7] - x[0]*x[1] + x[0] + x[1] - 1.0);
    result[1] = ( eps_sf + x[8]*x[4] + x[8]*x[0] + x[3]*x[4] + x[3]*x[0] - x[4]*x[7] + x[4]*x[1] - x[4] - x[7]*x[0] + x[0]*x[1] - x[0]);
    result[2] = ( eps_sf + x[6] + x[5] - x[8] - x[3] + 2.0*x[7] - x[1]);
    result[3] = ( eps_sf + -x[5]);
    result[4] = ( eps_sf + -x[6]);
    result[5] = ( eps_sf + -x[7]);
    result[6] = ( eps_sf + x[8]*x[4] + x[3]*x[4] + x[2]*x[0] - x[2] - x[4]*x[7] + x[4]*x[1] - x[4]);
    result[7] = ( eps_sf + -x[8]*x[4] - x[3]*x[4] - x[2]*x[0] + x[4]*x[7] - x[4]*x[1] + x[4]);
    result[8] = ( eps_sf + x[8] + x[3] + x[2] - 1.0);
    result[9] = ( eps_sf + -x[3]);
    result[10] = ( eps_sf + -x[8]);
    result[11] = ( eps_sf + 0.5*x[1] - 1.0);
    result[12] = ( eps_sf + -0.5*x[1]);

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
void ep_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( eps_sf + -x[0] + x[1]);
    result[1] = ( eps_sf + x[0] - x[1] - 1.0);
    result[2] = ( eps_sf + -x[0] - x[1]);
    result[3] = ( eps_sf + x[0] + x[1] - 1.0);

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
void fl_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( eps_sf + x[6] + x[3] + x[2] + x[9] + x[5] + x[4] + x[8] + x[1] + x[7] + x[0] - 1.0);
    result[1] = ( eps_sf + -x[1]);
    result[2] = ( eps_sf + -x[0]);
    result[3] = ( eps_sf + -x[2]);
    result[4] = ( eps_sf + -x[3]);
    result[5] = ( eps_sf + -x[4]);
    result[6] = ( eps_sf + -x[5]);
    result[7] = ( eps_sf + -x[6]);
    result[8] = ( eps_sf + -x[7]);
    result[9] = ( eps_sf + -x[8]);
    result[10] = ( eps_sf + -x[9]);
    result[11] = ( eps_sf + x[9] - 1.0);

    if (grad) {
        grad[0] = 1.0;
        grad[1] = 1.0;
        grad[2] = 1.0;
        grad[3] = 1.0;
        grad[4] = 1.0;
        grad[5] = 1.0;
        grad[6] = 1.0;
        grad[7] = 1.0;
        grad[8] = 1.0;
        grad[9] = 1.0;
        grad[10] = 0.0;
        grad[11] = -1.0;
        grad[12] = 0.0;
        grad[13] = 0.0;
        grad[14] = 0.0;
        grad[15] = 0.0;
        grad[16] = 0.0;
        grad[17] = 0.0;
        grad[18] = 0.0;
        grad[19] = 0.0;
        grad[20] = -1.0;
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
        grad[52] = 0.0;
        grad[53] = 0.0;
        grad[54] = -1.0;
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
        grad[65] = -1.0;
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
        grad[76] = -1.0;
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
        grad[87] = -1.0;
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
        grad[100] = 0.0;
        grad[101] = 0.0;
        grad[102] = 0.0;
        grad[103] = 0.0;
        grad[104] = 0.0;
        grad[105] = 0.0;
        grad[106] = 0.0;
        grad[107] = 0.0;
        grad[108] = 0.0;
        grad[109] = -1.0;
        grad[110] = 0.0;
        grad[111] = 0.0;
        grad[112] = 0.0;
        grad[113] = 0.0;
        grad[114] = 0.0;
        grad[115] = 0.0;
        grad[116] = 0.0;
        grad[117] = 0.0;
        grad[118] = 0.0;
        grad[119] = 1.0;
    }

    return;
};


/** 
  local minimization for garnet
*/
void g_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( eps_sf + -x[1]*x[0] + x[1] + x[0] - 1.0);
    result[1] = ( eps_sf + x[1]*x[0] - x[0]);
    result[2] = ( eps_sf + -x[1]);
    result[3] = ( eps_sf + x[3] + x[2] + 2.0*x[4] - 1.0);
    result[4] = ( eps_sf + -x[3]);
    result[5] = ( eps_sf + -x[2]);
    result[6] = ( eps_sf + -x[4]);

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
    }

    return;
};


/** 
  local minimization for hornblende
*/
void hb_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( eps_sf + x[3] - 1.0);
    result[1] = ( eps_sf + x[3]*x[4] - x[3]);
    result[2] = ( eps_sf + -x[3]*x[4]);
    result[3] = ( eps_sf + -x[8] + x[0] - 1.0);
    result[4] = ( eps_sf + x[8] - x[0]);
    result[5] = ( eps_sf + x[9]*x[6] + x[9]*x[7] + x[9]*x[1] - x[9] - x[6]*x[0] + x[6] - x[7]*x[0] + x[7] - x[0]*x[1] + x[0] + x[1] - 1.0);
    result[6] = ( eps_sf + -x[9]*x[6] - x[9]*x[7] - x[9]*x[1] + x[9] + x[6]*x[0] + x[7]*x[0] + x[0]*x[1] - x[0]);
    result[7] = ( eps_sf + -x[1]);
    result[8] = ( eps_sf + -x[6]);
    result[9] = ( eps_sf + -x[7]);
    result[10] = ( eps_sf + -x[5]);
    result[11] = ( eps_sf + 1.5*x[8] - x[9]*x[6] - x[9]*x[7] - x[9]*x[1] + x[9] - x[5]*x[0] + x[5] - x[0]*x[2] + x[0] + x[2] - 1.0);
    result[12] = ( eps_sf + -1.5*x[8] + x[9]*x[6] + x[9]*x[7] + x[9]*x[1] - x[9] + x[5]*x[0] + x[0]*x[2] - x[0]);
    result[13] = ( eps_sf + -x[2]);
    result[14] = ( eps_sf + 0.25*x[3] + 0.5*x[6] + 0.5*x[7] + 0.5*x[1] - 0.5*x[2] - 1.0);
    result[15] = ( eps_sf + -0.25*x[3] - 0.5*x[6] - 0.5*x[7] - 0.5*x[1] + 0.5*x[2]);
    result[16] = ( eps_sf + x[7] - 1.0);

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
        grad[50] = -x[6] - x[7] - x[1] + 1.0;
        grad[51] = x[9] - x[0] + 1.0;
        grad[52] = 0.0;
        grad[53] = 0.0;
        grad[54] = 0.0;
        grad[55] = 0.0;
        grad[56] = x[9] - x[0] + 1.0;
        grad[57] = x[9] - x[0] + 1.0;
        grad[58] = 0.0;
        grad[59] = x[6] + x[7] + x[1] - 1.0;
        grad[60] = x[6] + x[7] + x[1] - 1.0;
        grad[61] = -x[9] + x[0];
        grad[62] = 0.0;
        grad[63] = 0.0;
        grad[64] = 0.0;
        grad[65] = 0.0;
        grad[66] = -x[9] + x[0];
        grad[67] = -x[9] + x[0];
        grad[68] = 0.0;
        grad[69] = -x[6] - x[7] - x[1] + 1.0;
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
        grad[110] = -x[5] - x[2] + 1.0;
        grad[111] = -x[9];
        grad[112] = 1.0 - x[0];
        grad[113] = 0.0;
        grad[114] = 0.0;
        grad[115] = 1.0 - x[0];
        grad[116] = -x[9];
        grad[117] = -x[9];
        grad[118] = 1.50;
        grad[119] = -x[6] - x[7] - x[1] + 1.0;
        grad[120] = x[5] + x[2] - 1.0;
        grad[121] = x[9];
        grad[122] = x[0];
        grad[123] = 0.0;
        grad[124] = 0.0;
        grad[125] = x[0];
        grad[126] = x[9];
        grad[127] = x[9];
        grad[128] = -1.50;
        grad[129] = x[6] + x[7] + x[1] - 1.0;
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
        grad[143] = 0.250;
        grad[144] = 0.0;
        grad[145] = 0.0;
        grad[146] = 0.50;
        grad[147] = 0.50;
        grad[148] = 0.0;
        grad[149] = 0.0;
        grad[150] = 0.0;
        grad[151] = -0.50;
        grad[152] = 0.50;
        grad[153] = -0.250;
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
    }

    return;
};


/** 
  local minimization for ilmenite
*/
void ilm_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( eps_sf + -0.5*x[0] - 0.5*x[1]);
    result[1] = ( eps_sf + -0.5*x[0] + 0.5*x[1]);
    result[2] = ( eps_sf + x[0] - 1.0);
    result[3] = ( eps_sf + -0.5*x[0] + 0.5*x[1]);
    result[4] = ( eps_sf + -0.5*x[0] - 0.5*x[1]);
    result[5] = ( eps_sf + x[0] - 1.0);

    if (grad) {
        grad[0] = -0.50;
        grad[1] = -0.50;
        grad[2] = -0.50;
        grad[3] = 0.50;
        grad[4] = 1.0;
        grad[5] = 0.0;
        grad[6] = -0.50;
        grad[7] = 0.50;
        grad[8] = -0.50;
        grad[9] = -0.50;
        grad[10] = 1.0;
        grad[11] = 0.0;
    }

    return;
};


/** 
  local minimization for liquid (melt)
*/
void liq_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( eps_sf + x[6] + x[3] + x[2] + x[10] + x[5] + x[4] + x[8] + x[1] + x[7] + x[0] - 0.25*x[9]*(-3.0*x[6] - 3.0*x[3] - 3.0*x[2] - 3.0*x[10] - 3.0*x[5] - 3.0*x[4] - 3.0*x[8] - 3.0*x[1] - 3.0*x[7] - 3.0*x[0] + 4.0) - 1.0);
    result[1] = ( eps_sf + -0.75*x[1]*x[9] - x[1] + x[9]);
    result[2] = ( eps_sf + -0.75*x[0]*x[9] - x[0] + x[9]);
    result[3] = ( eps_sf + -0.75*x[4]*x[9] - x[4]);
    result[4] = ( eps_sf + -0.75*x[5]*x[9] - x[5]);
    result[5] = ( eps_sf + -0.75*x[6]*x[9] - x[6]);
    result[6] = ( eps_sf + -0.75*x[7]*x[9] - x[7]);
    result[7] = ( eps_sf + -0.75*x[8]*x[9] - x[8]);
    result[8] = ( eps_sf + -x[9]);
    result[9] = ( eps_sf + -x[3] - x[2] - 0.75*x[9]*(x[3] + x[2]));
    result[10] = ( eps_sf + 0.75*x[10]*x[9] + x[10] - 1.0);
    result[11] = ( eps_sf + -4.0*x[2]);
    result[12] = ( eps_sf + -4.0*x[3]);
    result[13] = ( eps_sf + -x[0]);
    result[14] = ( eps_sf + -x[1]);
    result[15] = ( eps_sf + -4.0*x[3] - 4.0*x[2] - x[1] - x[0]);
    result[16] = ( eps_sf + -x[10]);
    result[17] = ( eps_sf + x[10] - 1.0);

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
        grad[97] = -1.0;
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
        grad[123] = -4.0;
        grad[124] = 0.0;
        grad[125] = 0.0;
        grad[126] = 0.0;
        grad[127] = 0.0;
        grad[128] = 0.0;
        grad[129] = 0.0;
        grad[130] = 0.0;
        grad[131] = 0.0;
        grad[132] = 0.0;
        grad[133] = 0.0;
        grad[134] = 0.0;
        grad[135] = -4.0;
        grad[136] = 0.0;
        grad[137] = 0.0;
        grad[138] = 0.0;
        grad[139] = 0.0;
        grad[140] = 0.0;
        grad[141] = 0.0;
        grad[142] = 0.0;
        grad[143] = -1.0;
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
        grad[155] = -1.0;
        grad[156] = 0.0;
        grad[157] = 0.0;
        grad[158] = 0.0;
        grad[159] = 0.0;
        grad[160] = 0.0;
        grad[161] = 0.0;
        grad[162] = 0.0;
        grad[163] = 0.0;
        grad[164] = 0.0;
        grad[165] = -1.0;
        grad[166] = -1.0;
        grad[167] = -4.0;
        grad[168] = -4.0;
        grad[169] = 0.0;
        grad[170] = 0.0;
        grad[171] = 0.0;
        grad[172] = 0.0;
        grad[173] = 0.0;
        grad[174] = 0.0;
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
        grad[185] = 0.0;
        grad[186] = -1.0;
        grad[187] = 0.0;
        grad[188] = 0.0;
        grad[189] = 0.0;
        grad[190] = 0.0;
        grad[191] = 0.0;
        grad[192] = 0.0;
        grad[193] = 0.0;
        grad[194] = 0.0;
        grad[195] = 0.0;
        grad[196] = 0.0;
        grad[197] = 1.0;
    }

    return;
};


/** 
  local minimization for muscovite
*/
void mu_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( eps_sf + x[4] + x[3] - 1.0);
    result[1] = ( eps_sf + -x[3]);
    result[2] = ( eps_sf + -x[4]);
    result[3] = ( eps_sf + -x[0]*x[1] + x[0] + x[1] - 1.0);
    result[4] = ( eps_sf + x[0]*x[1] - x[0]);
    result[5] = ( eps_sf + -x[1]);
    result[6] = ( eps_sf + x[2] - 1.0);
    result[7] = ( eps_sf + -x[2]);
    result[8] = ( eps_sf + 0.5*x[4] + 0.5*x[1] - 1.0);
    result[9] = ( eps_sf + -0.5*x[4] - 0.5*x[1]);

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
void ol_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( eps_sf + -x[2] + x[0] - 1.0);
    result[1] = ( eps_sf + x[2] - x[0]);
    result[2] = ( eps_sf + -x[1]*x[0] + x[1] + x[2] + x[0] - 1.0);
    result[3] = ( eps_sf + x[1]*x[0] - x[2] - x[0]);
    result[4] = ( eps_sf + -x[1]);

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
void opx_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( eps_sf + -x[7]*x[3] - x[7]*x[0] + x[7] + x[3]*x[5] - x[3]*x[1] + x[3] + x[5]*x[0] - x[5] - x[0]*x[1] + x[0] + x[1] - 1.0);
    result[1] = ( eps_sf + x[7]*x[3] + x[7]*x[0] - x[3]*x[5] + x[3]*x[1] - x[3] - x[5]*x[0] + x[0]*x[1] - x[0]);
    result[2] = ( eps_sf + x[6] + x[4] - x[7] + 2.0*x[5] - x[1]);
    result[3] = ( eps_sf + -x[4]);
    result[4] = ( eps_sf + -x[6]);
    result[5] = ( eps_sf + -x[5]);
    result[6] = ( eps_sf + -x[2]*x[0] + x[2] + x[7]*x[3] - x[7]*x[0] + x[7] - x[3]*x[5] + x[3]*x[1] - x[3] + x[0] - 1.0);
    result[7] = ( eps_sf + x[2]*x[0] - x[7]*x[3] + x[7]*x[0] + x[3]*x[5] - x[3]*x[1] + x[3] - x[0]);
    result[8] = ( eps_sf + -x[2]);
    result[9] = ( eps_sf + -x[7]);
    result[10] = ( eps_sf + 0.5*x[1] - 1.0);
    result[11] = ( eps_sf + -0.5*x[1]);

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
void pl4T_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( eps_sf + x[0] + x[1] - 1.0);
    result[1] = ( eps_sf + -x[0]);
    result[2] = ( eps_sf + -x[1]);
    result[3] = ( eps_sf + -0.25*x[0] - 0.25);
    result[4] = ( eps_sf + 0.25*x[0] - 0.75);

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
void spn_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( eps_sf + -2.0/3.0*x[4] + 1.0/3.0*x[3]*x[0] - 1.0/3.0*x[3] + 1.0/3.0*x[0] - 1.0/3.0);
    result[1] = ( eps_sf + -2.0/3.0*x[5] - 1.0/3.0*x[3]*x[0] - 1.0/3.0*x[0]);
    result[2] = ( eps_sf + 2.0/3.0*x[4] + 2.0/3.0*x[5] + 2.0/3.0*x[6] - 2.0/3.0*x[2]*x[1] - 2.0/3.0*x[3]*x[1] + 1.0/3.0*x[3] + 2.0/3.0*x[1] - 2.0/3.0);
    result[3] = ( eps_sf + -2.0/3.0*x[6] + 2.0/3.0*x[2]*x[1] + 2.0/3.0*x[3]*x[1] - 2.0/3.0*x[1]);
    result[4] = ( eps_sf + 1.0/3.0*x[4] + 1.0/3.0*x[3]*x[0] - 1.0/3.0*x[3] + 1.0/3.0*x[0] - 1.0/3.0);
    result[5] = ( eps_sf + 1.0/3.0*x[5] - 1.0/3.0*x[3]*x[0] - 1.0/3.0*x[0]);
    result[6] = ( eps_sf + -1.0/3.0*x[4] - 1.0/3.0*x[5] - 1.0/3.0*x[6] - 2.0/3.0*x[2]*x[1] + x[2] - 2.0/3.0*x[3]*x[1] + 5.0/6.0*x[3] + 2.0/3.0*x[1] - 2.0/3.0);
    result[7] = ( eps_sf + 1.0/3.0*x[6] + 2.0/3.0*x[2]*x[1] + 2.0/3.0*x[3]*x[1] - 2.0/3.0*x[1]);
    result[8] = ( eps_sf + -x[2]);
    result[9] = ( eps_sf + -0.5*x[3]);

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
void bi_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( eps_sf + -x[2]*x[0] + x[2] + 2.0/3.0*x[4] - x[3]*x[0] + x[3] - x[0]*x[1] + x[0] + x[1] - 1.0);
    result[1] = ( eps_sf +  x[2]*x[0] - 2.0/3.0*x[4] + x[3]*x[0] + x[0]*x[1] - x[0]);
    result[2] = ( eps_sf + -x[2]);
    result[3] = ( eps_sf + -x[3]);
    result[4] = ( eps_sf + -x[1]);
    result[5] = ( eps_sf + -1.0/3.0*x[4] + x[0] - 1.0);
    result[6] = ( eps_sf + 1.0/3.0*x[4] - x[0]);
    result[7] = ( eps_sf + 0.5*x[2] + 0.5*x[1] - 0.5);
    result[8] = ( eps_sf + -0.5*x[2] - 0.5*x[1] - 0.5);
    result[9] = ( eps_sf + x[3] - 1.0);

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
        grad[36] = 0.50;
        grad[37] = 0.50;
        grad[38] = 0.0;
        grad[39] = 0.0;
        grad[40] = 0.0;
        grad[41] = -0.50;
        grad[42] = -0.50;
        grad[43] = 0.0;
        grad[44] = 0.0;
        grad[45] = 0.0;
        grad[46] = 0.0;
        grad[47] = 0.0;
        grad[48] = 1.0;
        grad[49] = 0.0;
    }

    return;
};

/** 
  local minimization for cordierite
*/
void cd_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( eps_sf + -x[0]);
    result[1] = ( eps_sf + x[0] - 1.0);
    result[2] = ( eps_sf + -x[1]);
    result[3] = ( eps_sf + x[1] - 1.0);

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

typedef struct global_min_datas {
	global_variable 	 gv; 
	struct bulk_info 	 z_b;
	obj_type 			*SS_objective;
	sf_type 			*SS_sf;
	PP_ref 				*PP_ref_db;
	SS_ref 				*SS_ref_db;
	csd_phase_set  		*cp;
	
} global_min_data;


/**
	associate the array of pointer with the right solution phase
*/
void SS_sf_init_function(	sf_type 			*SS_sf,
							global_variable 	 gv				){	

	for (int iss = 0; iss < gv.len_ss; iss++){

		if      (strcmp( gv.SS_list[iss], "bi")  == 0 ){
			SS_sf[iss]  = bi_c; 		}
		else if (strcmp( gv.SS_list[iss], "cd")  == 0){
			SS_sf[iss]  = cd_c; 		}
		else if (strcmp( gv.SS_list[iss], "cpx") == 0){
			SS_sf[iss]  = cpx_c; 		}
		else if (strcmp( gv.SS_list[iss], "ep")  == 0){
			SS_sf[iss]  = ep_c; 		}
		else if (strcmp( gv.SS_list[iss], "fl")  == 0){
			SS_sf[iss]  = fl_c; 		}
		else if (strcmp( gv.SS_list[iss], "g")   == 0){
			SS_sf[iss]  = g_c; 			}
		else if (strcmp( gv.SS_list[iss], "hb")  == 0){
			SS_sf[iss]  = hb_c; 		}
		else if (strcmp( gv.SS_list[iss], "ilm") == 0){
			SS_sf[iss]  = ilm_c; 		}
		else if (strcmp( gv.SS_list[iss], "liq") == 0){
			SS_sf[iss]  = liq_c; 		}
		else if (strcmp( gv.SS_list[iss], "mu")  == 0){
			SS_sf[iss]  = mu_c; 		}
		else if (strcmp( gv.SS_list[iss], "ol")  == 0){
			SS_sf[iss]  = ol_c; 		}
		else if (strcmp( gv.SS_list[iss], "opx") == 0){
			SS_sf[iss]  = opx_c; 		}
		else if (strcmp( gv.SS_list[iss], "pl4T") == 0){
			SS_sf[iss]  = pl4T_c; 		}
		else if (strcmp( gv.SS_list[iss], "spn") == 0){
			SS_sf[iss]  = spn_c; 		}
		else{
			printf("\nsolid solution '%s' is not in the database, cannot be initiated\n", gv.SS_list[iss]);	
		}	
	};			
}

/**
  Equality constraint for AUGLAG global minimization approach
*/
void GM_eq(unsigned l, double *result, unsigned n, const double *x, double *grad, void *GM_db){
	global_min_data *d  = (global_min_data *)GM_db;
	int ph, ss, sf, h, i, j, k, ix, iy, iz, sf_ok;
		
	double Gph, alpha;
	
	iy = 0;
	/** loop through every non-zero components of the bulk-rock composition */
	for (k = 0; k < d->z_b.nzEl_val; k++){
		ix = 0;		

		result[iy] = d->z_b.bulk_rock[d->z_b.nzEl_array[k]];
		
		/** loop through all the solution models considered in the assemblage */
		for (i = 0; i < d->gv.n_cp_phase; i++){
			ph 	  = d->gv.cp_id[i];
			ss    = d->cp[ph].id;
			
			alpha = x[ix];
			
			/* get initial guess */
			for (j = 0; j < d->cp[ph].n_xeos; j++){
				d->SS_ref_db[ss].iguess[j] = x[ix+1+j];
			}	
			
			/* retrieve phase data at given set of compositional variables */
			Gph   = (*d->SS_objective[ss])(		 		 d->cp[ph].n_xeos,
														 d->SS_ref_db[ss].iguess, 
														 d->SS_ref_db[ss].dfx, 
														&d->SS_ref_db[ss]				);
			
			for (j = 0; j < d->SS_ref_db[ss].n_em; j++){
				result[iy] -= alpha*d->SS_ref_db[ss].p[j]*d->SS_ref_db[ss].Comp[j][d->z_b.nzEl_array[k]]*d->SS_ref_db[ss].factor*d->SS_ref_db[ss].z_em[j];
			}
			
			ix += d->SS_ref_db[ss].n_em;
		}	
		/** loop though pure phase */
		for (i = 0; i < d->gv.n_pp_phase; i++){
			ph 	  = d->gv.pp_id[i];
			
			alpha = x[ix];
			
			result[iy] -= alpha*d->PP_ref_db[ph].Comp[d->z_b.nzEl_array[k]]*d->PP_ref_db[ph].factor;
			
			ix += 1;
			//Gph   = d->PP_ref_db[ph].gb_lvl;
		}
		iy += 1;
	}
 
	
	for (k = 0; k < d->z_b.nzEl_val; k++){
		printf(" %+10f",result[k]);
	}
	printf("\n");
    return;
};

/**
  Sets inequality constraints for the global minimization approach
*/
void GM_ineq(unsigned m, double *result, unsigned n, const double *x, double *grad, void *GM_db){
	global_min_data *d  = (global_min_data *)GM_db;
	int ph, ss, sf, i, j, k, ix, iy, iz, sf_ok;
		
	ix 		  = 0;
	iy		  = 0;
	iz		  = 0;
	for (i = 0; i < d->gv.n_cp_phase; i++){
		ph 	  = d->gv.cp_id[i];
		ss    = d->cp[ph].id;
		
		for (j = 0; j < d->SS_ref_db[ss].n_xeos; j++){
			//d->SS_ref_db[ss].iguess[j] = x[ix+1+j];
			d->SS_ref_db[ss].iguess[j] = x[ix+j];
		}	

		(*d->SS_sf[ss])(		d->SS_ref_db[ss].n_sf, 
								d->SS_ref_db[ss].sf, 
								d->SS_ref_db[ss].n_xeos, 
								d->SS_ref_db[ss].iguess, 
								d->SS_ref_db[ss].dsf, 
								NULL							);
								
		for (j = 0; j < d->SS_ref_db[ss].n_sf; j++){
			result[iy] = d->SS_ref_db[ss].sf[j];
			iy        += 1;
		}
		
		if (grad){
			sf    = 0;													//increment for site fraction index
			//iz 	 += 1;		
			for (j = 0; j < d->SS_ref_db[ss].n_sf; j++){
				for (k = 0; k < d->SS_ref_db[ss].n_xeos; k++){	
					grad[iz+k] = d->SS_ref_db[ss].dsf[sf];
					sf   += 1;
				}
				iz		 += n;
			}
			iz		 += d->SS_ref_db[ss].n_xeos;
		}
		
		//ix += d->SS_ref_db[ss].n_em;
		ix += d->SS_ref_db[ss].n_xeos;
	}

    return;
};

/** 
  objective function constructor for generalized iminization with Augmented Lagrangian approach
*/
double GM_obj(		unsigned  		 n, 
					const double 	*x, 
					double 			*grad, 
					void 			*GM_db		)
{			
	global_min_data *d  = (global_min_data *)GM_db;
	int ph, ss, i, j, ix, sf_ok;
	
	double Gsys = 0.0;
	double Gph, alpha;
	ix 		  = 0;
	for (i = 0; i < d->gv.n_cp_phase; i++){
		ph 	  = d->gv.cp_id[i];
		ss    = d->cp[ph].id;

		//alpha = x[ix];
		alpha = d->cp[ph].ss_n;
		
		for (j = 0; j < d->cp[ph].n_xeos; j++){
			//d->SS_ref_db[ss].iguess[j] = x[ix+1+j];
			d->SS_ref_db[ss].iguess[j] = x[ix+j];
		}	
		
		Gph   = (*d->SS_objective[ss])(		 		 d->cp[ph].n_xeos,
													 d->SS_ref_db[ss].iguess, 
													 d->SS_ref_db[ss].dfx, 
													&d->SS_ref_db[ss]				);
		Gsys += alpha*Gph;
		
		printf(" [%4s %+12.5f %+12.5f]",d->gv.SS_list[ss],Gph,alpha);
		
		sf_ok = 1;
		for (j = 0; j < d->cp[ph].n_sf; j++){
			if (d->SS_ref_db[ss].sf[j] < 0.0){ sf_ok = 0;};
		}
		printf("SFOK? %d |",sf_ok);
		for (j = 0; j < d->cp[ph].n_xeos; j++){
			printf(" %+12.5f",d->SS_ref_db[ss].iguess[j]);
		}
		printf("\n");
		
		if (grad){
			//grad[ix] = Gph;
			//ix 		+= 1;
			for (j = 0; j < d->cp[ph].n_xeos; j++){
				grad[ix] = alpha*d->SS_ref_db[ss].dfx[j];
				ix      += 1;
			}
		}
	}

	//for (i = 0; i < d->gv.n_pp_phase; i++){
		//ph 	  = d->gv.pp_id[i];
		
		//alpha = x[ix];
		
		//Gph   = d->PP_ref_db[ph].gb_lvl;
		//Gsys += Gph*alpha;
		//printf(" [%4s %+12.5f %+12.5f]\n",d->gv.PP_list[ph],Gph,x[ix]);
		//if (grad){
			//grad[ix] = Gph;
			//ix 		+= 1;
		//}
	//}
	
	//for (i = 0; i < n; i++){
		//printf(" %+10f",grad[i]);
	//}
	printf(" Gsys: %+12.5f\n",Gsys);

	return Gsys;
};

/** 
	global minimization given fixed phase assemblage
	n: number of variables in the objective function (phase fractions + compositional variables)
	m: number of inequality constraints (site-fractions)
	l: number of equality constraints (active number of components)
*/
global_variable NLopt_global_opt_function(	struct bulk_info 	z_b,
											global_variable 	 gv, 
											PP_ref 				*PP_ref_db,
											SS_ref 				*SS_ref_db,
	 										csd_phase_set  		*cp
){
	int ph, ss, i, j; 
    unsigned int n, m, l, ix;
    
	/**
	   create objective function array to be send to NLopt
	*/
	obj_type 		SS_objective[gv.len_ss];
	sf_type 		SS_sf[gv.len_ss];
	
	SS_objective_init_function(			SS_objective,
										gv				);

	SS_sf_init_function(				SS_sf,
										gv				);

	/** create and fill data structure to be send to the objective function */								
	global_min_data GM_db;
	
	GM_db.gv 			= gv;
	GM_db.z_b			= z_b;
	GM_db.SS_objective 	= SS_objective;
	GM_db.SS_sf  	 	= SS_sf;
	GM_db.PP_ref_db 	= PP_ref_db;
	GM_db.SS_ref_db 	= SS_ref_db;
	GM_db.cp 			= cp;
	
	/** get id of active phases */
	//gv = get_pp_id(			gv					);
	
	gv = get_ss_id(			gv,
							cp					);
	
	//n  = gv.n_pp_phase;
	n = 0;	
	for (i = 0; i < gv.n_cp_phase; i++){
		ph  = gv.cp_id[i]; 
		//n  += cp[ph].n_em;
		n  += cp[ph].n_xeos;
	}

	/** get the number of equality constraints */
	l = z_b.nzEl_val;
	
	/** declare solution array 2BM */    
	double x[n]; 
	double lb[n]; 
	double ub[n]; 

	/**	initialize x array */
	ix = 0;
	m  = 0;
	for (i = 0; i < gv.n_cp_phase; i++){
		ph     = gv.cp_id[i]; 
		ss	   = cp[ph].id;
		m 	  += cp[ph].n_sf;

		SS_ref_db[ss] = rotate_hyperplane(	gv, 
											SS_ref_db[ss]	);
											
		//x[ix]  = cp[ph].ss_n;
		//lb[ix] = 0.0;
		//ub[ix] = 1.0;
		
		//ix    += 1;
		for (j = 0; j < SS_ref_db[ss].n_xeos; j++){
			x[ix]  = cp[ph].xeos[j];
			
			lb[ix] = SS_ref_db[ss].bounds_ref[j][0];
			ub[ix] = SS_ref_db[ss].bounds_ref[j][1];
			ix    += 1;
		}
	}
	//for (i = 0; i < gv.n_pp_phase; i++){
		//ph     = gv.pp_id[i];
		//x[ix]  = gv.pp_n[ph];
		//lb[ix] = 0.0;
		//ub[ix] = 1.0;
		//ix    += 1;
	//}

	double tol_sf[m];
	for (i = 0; i < m; i++){
		tol_sf[i] = 0.0;
	}

	double tol_eq[l];
	for (i = 0; i < l; i++){
		tol_eq[i] = 1e-8;
	}
   
   
	//double minf;

	//nlopt_opt opt_auglag = nlopt_create(NLOPT_AUGLAG, n);
	//nlopt_set_lower_bounds(opt_auglag, lb);
	//nlopt_set_upper_bounds(opt_auglag, ub);
	//nlopt_set_min_objective(opt_auglag, GM_obj, &GM_db);
	//nlopt_add_equality_mconstraint(opt_auglag, l, GM_eq, &GM_db, tol_eq);
	//nlopt_add_inequality_mconstraint(opt_auglag, m, GM_ineq, &GM_db, tol_sf);
	//nlopt_set_ftol_rel(opt_auglag, gv.obj_tol);
	//nlopt_set_maxeval(opt_auglag, 0);

	//nlopt_opt opt = nlopt_create(NLOPT_LD_CCSAQ, n);
	//nlopt_set_lower_bounds(opt, lb);
	//nlopt_set_upper_bounds(opt, ub);
	//nlopt_set_min_objective(opt, GM_obj, &GM_db);

	//nlopt_set_ftol_rel(opt, gv.obj_tol);
	//nlopt_set_maxeval(opt, 0);

	//nlopt_set_local_optimizer(opt_auglag, opt);
	
	//int status = nlopt_optimize(opt_auglag, x, &minf);
	//printf("Solution of global minimization: %+10f, status %d\n",minf,status);
	//for (i = 0; i < n; i++){
		//printf(" %+10f",x[i]);
	//}
	//printf("\n");

	//nlopt_destroy(opt_auglag);
	//nlopt_destroy(opt);
	
	
	

	nlopt_opt opt = nlopt_create(NLOPT_LD_CCSAQ, n);
	nlopt_set_lower_bounds(opt, lb);
	nlopt_set_upper_bounds(opt, ub);
	nlopt_set_min_objective(opt, GM_obj, &GM_db);
	nlopt_add_inequality_mconstraint(opt, m, GM_ineq, &GM_db, tol_sf);
	//nlopt_add_equality_mconstraint(opt, l, GM_eq,   &GM_db, tol_eq);
	nlopt_set_ftol_rel(opt, gv.obj_tol);
	nlopt_set_maxeval(opt, 0);

	//nlopt_set_maxeval(opt, gv.maxeval);

	double minf;
	int status = nlopt_optimize(opt, x, &minf);
	
	printf("Solution of global minimization: %+10f, status %d\n",minf,status);
	for (i = 0; i < n; i++){
		printf(" %+10f",x[i]);
	}
	printf("\n");

	nlopt_destroy(opt);


	return gv;
}

SS_ref NLopt_opt_bi_function(global_variable gv, SS_ref SS_ref_db){

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
   nlopt_set_min_objective(SS_ref_db.opt, obj_bi, &SS_ref_db);
   nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, bi_c, NULL, SS_ref_db.tol_sf);
   nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
   nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
   
   double minf;
   if (gv.maxeval==1){  
     // we are only interested in evaluating the objective function  
     minf = obj_bi(n, x, NULL, &SS_ref_db); 
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

SS_ref NLopt_opt_cd_function(global_variable gv, SS_ref SS_ref_db){

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
	nlopt_set_min_objective(SS_ref_db.opt, obj_cd, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, cd_c, NULL, SS_ref_db.tol_sf);
	nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
   
	double minf;
	if (gv.maxeval==1){  
          // we are only interested in evaluating the objective function  
          minf = obj_cd(n, x, NULL, &SS_ref_db); 
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

SS_ref NLopt_opt_cpx_function(global_variable gv, SS_ref SS_ref_db){

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
   nlopt_set_min_objective(SS_ref_db.opt, obj_cpx, &SS_ref_db);
   nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, cpx_c, NULL, SS_ref_db.tol_sf);
   nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
   nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
   
   double minf;
   if (gv.maxeval==1){  
     // we are only interested in evaluating the objective function  
     minf = obj_cpx(n, x, NULL, &SS_ref_db); 
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

SS_ref NLopt_opt_ep_function(global_variable gv, SS_ref SS_ref_db){

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
   nlopt_set_min_objective(SS_ref_db.opt, obj_ep, &SS_ref_db);
   nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, ep_c, NULL, SS_ref_db.tol_sf);
   nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
   nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
   
   double minf;
   if (gv.maxeval==1){  
     // we are only interested in evaluating the objective function  
     minf = obj_ep(n, x, NULL, &SS_ref_db); 
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

SS_ref NLopt_opt_fl_function(global_variable gv, SS_ref SS_ref_db){

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
   nlopt_set_min_objective(SS_ref_db.opt, obj_fl, &SS_ref_db);
   nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, fl_c, NULL, SS_ref_db.tol_sf);
   nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
   nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
   
   double minf;
   if (gv.maxeval==1){  
     // we are only interested in evaluating the objective function  
     minf = obj_fl(n, x, NULL, &SS_ref_db); 
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

SS_ref NLopt_opt_g_function(global_variable gv, SS_ref SS_ref_db){
    
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
   nlopt_set_min_objective(SS_ref_db.opt, obj_g, &SS_ref_db);
   nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, g_c, NULL, SS_ref_db.tol_sf);
   nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
   nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
   
   double minf;
   if (gv.maxeval==1){  
     // we are only interested in evaluating the objective function  
     minf = obj_g(n, x, NULL, &SS_ref_db); 
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

SS_ref NLopt_opt_hb_function(global_variable gv, SS_ref SS_ref_db){

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
   nlopt_set_min_objective(SS_ref_db.opt, obj_hb, &SS_ref_db);
   nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, hb_c, NULL, SS_ref_db.tol_sf);
   nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
   nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
   
   double minf;
   if (gv.maxeval==1){  
     // we are only interested in evaluating the objective function  
     minf = obj_hb(n, x, NULL, &SS_ref_db); 
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

SS_ref NLopt_opt_ilm_function(global_variable gv, SS_ref SS_ref_db){

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
   nlopt_set_min_objective(SS_ref_db.opt, obj_ilm, &SS_ref_db);
   nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, ilm_c, NULL, SS_ref_db.tol_sf);
   nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
   nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
   
   double minf;
   if (gv.maxeval==1){  
     // we are only interested in evaluating the objective function  
     minf = obj_ilm(n, x, NULL, &SS_ref_db); 
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

SS_ref NLopt_opt_liq_function(global_variable gv, SS_ref SS_ref_db){

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
   nlopt_set_min_objective(SS_ref_db.opt, obj_liq, &SS_ref_db);
   nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, liq_c, NULL, SS_ref_db.tol_sf);
   nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
   nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
   
   double minf;
   if (gv.maxeval==1){  
     // we are only interested in evaluating the objective function  
     minf = obj_liq(n, x, NULL, &SS_ref_db); 
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

SS_ref NLopt_opt_mu_function(global_variable gv, SS_ref SS_ref_db){

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
   nlopt_set_min_objective(SS_ref_db.opt, obj_mu, &SS_ref_db);
   nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, mu_c, NULL, SS_ref_db.tol_sf);
   nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
   nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
   
   double minf;
   if (gv.maxeval==1){  
     // we are only interested in evaluating the objective function  
     minf = obj_mu(n, x, NULL, &SS_ref_db); 
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

SS_ref NLopt_opt_ol_function(global_variable gv, SS_ref SS_ref_db){
   
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
   nlopt_set_min_objective(SS_ref_db.opt, obj_ol, &SS_ref_db);
   nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, ol_c, NULL, SS_ref_db.tol_sf);
   nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
   nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
   
   double minf;
   if (gv.maxeval==1){  
     // we are only interested in evaluating the objective function  
     minf = obj_ol(n, x, NULL, &SS_ref_db); 
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

SS_ref NLopt_opt_opx_function(global_variable gv, SS_ref SS_ref_db){
    
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
   nlopt_set_min_objective(SS_ref_db.opt, obj_opx, &SS_ref_db);
   nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, opx_c, NULL, SS_ref_db.tol_sf);
   nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
   nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
   
   double minf;
   if (gv.maxeval==1){  
     // we are only interested in evaluating the objective function  
     minf = obj_opx(n, x, NULL, &SS_ref_db); 
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

SS_ref NLopt_opt_pl4T_function(global_variable gv, SS_ref SS_ref_db){
   
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
   nlopt_set_min_objective(SS_ref_db.opt, obj_pl4T, &SS_ref_db);
   nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, pl4T_c, NULL, SS_ref_db.tol_sf);
   nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
   nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
   
   double minf;
   if (gv.maxeval==1){  
     // we are only interested in evaluating the objective function  
     minf = obj_pl4T(n, x, NULL, &SS_ref_db); 
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

SS_ref NLopt_opt_spn_function(global_variable gv, SS_ref SS_ref_db){

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
   nlopt_set_min_objective(SS_ref_db.opt, obj_spn, &SS_ref_db);
   nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, spn_c, NULL, SS_ref_db.tol_sf);
   nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
   nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
   
   double minf;

   if (gv.maxeval==1){  
     // we are only interested in evaluating the objective function  
     minf = obj_spn(n, x, NULL, &SS_ref_db); 
   }
   else{
     // do optimization
     SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);
   }
   
   /* Send back needed local solution parameters */
   for (int i = 0; i < SS_ref_db.n_xeos; i++){
      SS_ref_db.xeos[i]   = x[i];
   }
 
   SS_ref_db.df   = minf;
   nlopt_destroy(SS_ref_db.opt);
  	
   return SS_ref_db;
};
/** 
  attributes the right solution phase to the solution phase array and calculates xi
*/
SS_ref NLopt_opt_function(		global_variable gv,
								SS_ref 			SS_ref_db, 
								int     		index			){
								
	clock_t t; 
	t = clock();

	/* Associate the right solid-solution data */
	if 		(strcmp( gv.SS_list[index], "bi") == 0 ){
		SS_ref_db  = NLopt_opt_bi_function( gv, SS_ref_db);	}
	else if (strcmp( gv.SS_list[index], "cd")  == 0){
		SS_ref_db  = NLopt_opt_cd_function( gv, SS_ref_db);	}
	else if (strcmp( gv.SS_list[index], "cpx") == 0){
		SS_ref_db  = NLopt_opt_cpx_function( gv, SS_ref_db);}	
	else if (strcmp( gv.SS_list[index], "ep")  == 0){
		SS_ref_db  = NLopt_opt_ep_function( gv, SS_ref_db);	}
	else if (strcmp( gv.SS_list[index], "fl")  == 0){
		SS_ref_db  = NLopt_opt_fl_function( gv, SS_ref_db);	}		
	else if (strcmp( gv.SS_list[index], "g")   == 0){
		SS_ref_db  = NLopt_opt_g_function(  gv, SS_ref_db);	}
	else if (strcmp( gv.SS_list[index], "hb")  == 0){
		SS_ref_db  = NLopt_opt_hb_function( gv, SS_ref_db);	}	
	else if (strcmp( gv.SS_list[index], "ilm") == 0){
		SS_ref_db  = NLopt_opt_ilm_function( gv, SS_ref_db);}
	else if (strcmp( gv.SS_list[index], "liq") == 0){
		SS_ref_db  = NLopt_opt_liq_function( gv, SS_ref_db);}
	else if (strcmp( gv.SS_list[index], "mu")  == 0){
		SS_ref_db  = NLopt_opt_mu_function( gv, SS_ref_db);	}	
	else if (strcmp( gv.SS_list[index], "ol")  == 0){
		SS_ref_db  = NLopt_opt_ol_function( gv, SS_ref_db);	}
	else if (strcmp( gv.SS_list[index], "opx") == 0){
		SS_ref_db  = NLopt_opt_opx_function( gv, SS_ref_db);}	
	else if (strcmp( gv.SS_list[index], "pl4T")  == 0){
		SS_ref_db  = NLopt_opt_pl4T_function( gv, SS_ref_db);	}	
	else if (strcmp( gv.SS_list[index], "spn") == 0){
		SS_ref_db  = NLopt_opt_spn_function( gv, SS_ref_db);	
		}
	else{
		printf("\nsolid solution '%s index %d' is not in the database\n",gv.SS_list[index], index);	}	
		
   t = clock() - t; 
   SS_ref_db.LM_time = ((double)t)/CLOCKS_PER_SEC*1000; // in seconds 
   
	return SS_ref_db;
};


