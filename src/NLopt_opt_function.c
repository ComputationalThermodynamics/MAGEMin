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
#include "gss_function.h"		    	  // order of header file declaration is important
#include "objective_functions.h"
#include "NLopt_opt_function.h"
#include "toolkit.h"

#define eps_sf -1e-10				      // eps to shift site fraction from zero

/**************************************************************************************/
/**************************************************************************************/
/********************IGNEOUS ALKALINE DATABASE (Weller et al., 2023)*******************/
/**************************************************************************************/
/**************************************************************************************/

/**
    Inequality constraints for liq
*/
void liq_alk_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( eps_sf + x[0] - x[10]*(-x[0] - x[1] - x[2] - x[3] - x[4] - x[5] - x[6] - x[7] - x[8] - x[9] + 1.0) - 1.0/6.0*x[11]*(-7.0*x[0] - 7.0*x[1] - 7.0*x[2] - 7.0*x[3] - 7.0*x[4] - 7.0*x[5] - 7.0*x[6] - 7.0*x[7] - 7.0*x[8] - 7.0*x[9] + 3.0) - 1.0/6.0*x[12]*(x[0] + x[1] + x[2] + x[3] + x[4] + x[5] + x[6] + x[7] + x[8] + x[9] - 3.0) - 1.0/6.0*x[13]*(-7.0*x[0] - 7.0*x[1] - 7.0*x[2] - 7.0*x[3] - 7.0*x[4] - 7.0*x[5] - 7.0*x[6] - 7.0*x[7] - 7.0*x[8] - 7.0*x[9] + 3.0) + x[1] + x[2] + x[3] + x[4] + x[5] + x[6] + x[7] + x[8] + x[9] - 1.0);
    result[1] = ( eps_sf + -x[10]*x[1] + x[10] - 7.0/6.0*x[11]*x[1] + 0.5*x[11] + 1.0/6.0*x[12]*x[1] - 7.0/6.0*x[13]*x[1] + 0.5*x[13] - x[1]);
    result[2] = ( eps_sf + -x[0]*x[10] - 7.0/6.0*x[0]*x[11] + 1.0/6.0*x[0]*x[12] - 7.0/6.0*x[0]*x[13] - x[0] + x[10]);
    result[3] = ( eps_sf + -x[10]*x[4] - 7.0/6.0*x[11]*x[4] + x[11] + 1.0/6.0*x[12]*x[4] - 7.0/6.0*x[13]*x[4] - x[4]);
    result[4] = ( eps_sf + -x[10]*x[5] - 7.0/6.0*x[11]*x[5] + 1.0/6.0*x[12]*x[5] - 7.0/6.0*x[13]*x[5] - x[5]);
    result[5] = ( eps_sf + -x[10]*x[6] - 7.0/6.0*x[11]*x[6] + 1.0/6.0*x[12]*x[6] - 7.0/6.0*x[13]*x[6] - x[6]);
    result[6] = ( eps_sf + -x[10]*x[7] - 7.0/6.0*x[11]*x[7] + 1.0/6.0*x[12]*x[7] - 7.0/6.0*x[13]*x[7] - x[7]);
    result[7] = ( eps_sf + -x[10]*x[8] - 7.0/6.0*x[11]*x[8] + 1.0/6.0*x[12]*x[8] - 7.0/6.0*x[13]*x[8] + x[13] - x[8]);
    result[8] = ( eps_sf + -x[11]);
    result[9] = ( eps_sf + -x[10]);
    result[10] = ( eps_sf + -x[12]);
    result[11] = ( eps_sf + -x[13]);
    result[12] = ( eps_sf + -x[10]*(x[2] + x[3]) - 7.0/6.0*x[11]*(x[2] + x[3]) + 1.0/6.0*x[12]*(x[2] + x[3]) + 0.5*x[12] - 7.0/6.0*x[13]*(x[2] + x[3]) - x[2] - x[3]);
    result[13] = ( eps_sf + x[10]*x[9] + 7.0/6.0*x[11]*x[9] - 1.0/6.0*x[12]*x[9] + 7.0/6.0*x[13]*x[9] + x[9] - 1.0);
    result[14] = ( eps_sf + -x[10]*x[9] - 7.0/6.0*x[11]*x[9] + 1.0/6.0*x[12]*x[9] - 7.0/6.0*x[13]*x[9] - x[9]);
    result[15] = ( eps_sf + x[10]*x[9] + 7.0/6.0*x[11]*x[9] - 1.0/6.0*x[12]*x[9] + 7.0/6.0*x[13]*x[9] + x[9] - 1.0);
    result[16] = ( eps_sf + 2.0*x[12] - 4.0*x[2]*(x[10] + 7.0/6.0*x[11] - 1.0/6.0*x[12] + 7.0/6.0*x[13] + 1.0));
    result[17] = ( eps_sf + -4.0*x[3]*(x[10] + 7.0/6.0*x[11] - 1.0/6.0*x[12] + 7.0/6.0*x[13] + 1.0));
    result[18] = ( eps_sf + -x[0]*(x[10] + 7.0/6.0*x[11] - 1.0/6.0*x[12] + 7.0/6.0*x[13] + 1.0) + x[10]);
    result[19] = ( eps_sf + x[10] + 0.5*x[11] + 0.5*x[13] - x[1]*(x[10] + 7.0/6.0*x[11] - 1.0/6.0*x[12] + 7.0/6.0*x[13] + 1.0));
    result[20] = ( eps_sf + 2.0*x[10] + 0.5*x[11] + 2.0*x[12] + 0.5*x[13] - (x[0] + x[1] + 4.0*x[2] + 4.0*x[3])*(x[10] + 7.0/6.0*x[11] - 1.0/6.0*x[12] + 7.0/6.0*x[13] + 1.0));

    if (grad) {
        grad[0] = x[10] + 7.0/6.0*x[11] - 1.0/6.0*x[12] + 7.0/6.0*x[13] + 1.0;
        grad[1] = x[10] + 7.0/6.0*x[11] - 1.0/6.0*x[12] + 7.0/6.0*x[13] + 1.0;
        grad[2] = x[10] + 7.0/6.0*x[11] - 1.0/6.0*x[12] + 7.0/6.0*x[13] + 1.0;
        grad[3] = x[10] + 7.0/6.0*x[11] - 1.0/6.0*x[12] + 7.0/6.0*x[13] + 1.0;
        grad[4] = x[10] + 7.0/6.0*x[11] - 1.0/6.0*x[12] + 7.0/6.0*x[13] + 1.0;
        grad[5] = x[10] + 7.0/6.0*x[11] - 1.0/6.0*x[12] + 7.0/6.0*x[13] + 1.0;
        grad[6] = x[10] + 7.0/6.0*x[11] - 1.0/6.0*x[12] + 7.0/6.0*x[13] + 1.0;
        grad[7] = x[10] + 7.0/6.0*x[11] - 1.0/6.0*x[12] + 7.0/6.0*x[13] + 1.0;
        grad[8] = x[10] + 7.0/6.0*x[11] - 1.0/6.0*x[12] + 7.0/6.0*x[13] + 1.0;
        grad[9] = x[10] + 7.0/6.0*x[11] - 1.0/6.0*x[12] + 7.0/6.0*x[13] + 1.0;
        grad[10] = x[0] + x[1] + x[2] + x[3] + x[4] + x[5] + x[6] + x[7] + x[8] + x[9] - 1.0;
        grad[11] = 7.0/6.0*x[0] + 7.0/6.0*x[1] + 7.0/6.0*x[2] + 7.0/6.0*x[3] + 7.0/6.0*x[4] + 7.0/6.0*x[5] + 7.0/6.0*x[6] + 7.0/6.0*x[7] + 7.0/6.0*x[8] + 7.0/6.0*x[9] - 0.5;
        grad[12] = -1.0/6.0*x[0] - 1.0/6.0*x[1] - 1.0/6.0*x[2] - 1.0/6.0*x[3] - 1.0/6.0*x[4] - 1.0/6.0*x[5] - 1.0/6.0*x[6] - 1.0/6.0*x[7] - 1.0/6.0*x[8] - 1.0/6.0*x[9] + 0.5;
        grad[13] = 7.0/6.0*x[0] + 7.0/6.0*x[1] + 7.0/6.0*x[2] + 7.0/6.0*x[3] + 7.0/6.0*x[4] + 7.0/6.0*x[5] + 7.0/6.0*x[6] + 7.0/6.0*x[7] + 7.0/6.0*x[8] + 7.0/6.0*x[9] - 0.5;
        grad[14] = 0.0;
        grad[15] = -x[10] - 7.0/6.0*x[11] + 1.0/6.0*x[12] - 7.0/6.0*x[13] - 1.0;
        grad[16] = 0.0;
        grad[17] = 0.0;
        grad[18] = 0.0;
        grad[19] = 0.0;
        grad[20] = 0.0;
        grad[21] = 0.0;
        grad[22] = 0.0;
        grad[23] = 0.0;
        grad[24] = 1.0 - x[1];
        grad[25] = 0.5 - 7.0/6.0*x[1];
        grad[26] = 1.0/6.0*x[1];
        grad[27] = 0.5 - 7.0/6.0*x[1];
        grad[28] = -x[10] - 7.0/6.0*x[11] + 1.0/6.0*x[12] - 7.0/6.0*x[13] - 1.0;
        grad[29] = 0.0;
        grad[30] = 0.0;
        grad[31] = 0.0;
        grad[32] = 0.0;
        grad[33] = 0.0;
        grad[34] = 0.0;
        grad[35] = 0.0;
        grad[36] = 0.0;
        grad[37] = 0.0;
        grad[38] = 1.0 - x[0];
        grad[39] = -7.0/6.0*x[0];
        grad[40] = 1.0/6.0*x[0];
        grad[41] = -7.0/6.0*x[0];
        grad[42] = 0.0;
        grad[43] = 0.0;
        grad[44] = 0.0;
        grad[45] = 0.0;
        grad[46] = -x[10] - 7.0/6.0*x[11] + 1.0/6.0*x[12] - 7.0/6.0*x[13] - 1.0;
        grad[47] = 0.0;
        grad[48] = 0.0;
        grad[49] = 0.0;
        grad[50] = 0.0;
        grad[51] = 0.0;
        grad[52] = -x[4];
        grad[53] = 1.0 - 7.0/6.0*x[4];
        grad[54] = 1.0/6.0*x[4];
        grad[55] = -7.0/6.0*x[4];
        grad[56] = 0.0;
        grad[57] = 0.0;
        grad[58] = 0.0;
        grad[59] = 0.0;
        grad[60] = 0.0;
        grad[61] = -x[10] - 7.0/6.0*x[11] + 1.0/6.0*x[12] - 7.0/6.0*x[13] - 1.0;
        grad[62] = 0.0;
        grad[63] = 0.0;
        grad[64] = 0.0;
        grad[65] = 0.0;
        grad[66] = -x[5];
        grad[67] = -7.0/6.0*x[5];
        grad[68] = 1.0/6.0*x[5];
        grad[69] = -7.0/6.0*x[5];
        grad[70] = 0.0;
        grad[71] = 0.0;
        grad[72] = 0.0;
        grad[73] = 0.0;
        grad[74] = 0.0;
        grad[75] = 0.0;
        grad[76] = -x[10] - 7.0/6.0*x[11] + 1.0/6.0*x[12] - 7.0/6.0*x[13] - 1.0;
        grad[77] = 0.0;
        grad[78] = 0.0;
        grad[79] = 0.0;
        grad[80] = -x[6];
        grad[81] = -7.0/6.0*x[6];
        grad[82] = 1.0/6.0*x[6];
        grad[83] = -7.0/6.0*x[6];
        grad[84] = 0.0;
        grad[85] = 0.0;
        grad[86] = 0.0;
        grad[87] = 0.0;
        grad[88] = 0.0;
        grad[89] = 0.0;
        grad[90] = 0.0;
        grad[91] = -x[10] - 7.0/6.0*x[11] + 1.0/6.0*x[12] - 7.0/6.0*x[13] - 1.0;
        grad[92] = 0.0;
        grad[93] = 0.0;
        grad[94] = -x[7];
        grad[95] = -7.0/6.0*x[7];
        grad[96] = 1.0/6.0*x[7];
        grad[97] = -7.0/6.0*x[7];
        grad[98] = 0.0;
        grad[99] = 0.0;
        grad[100] = 0.0;
        grad[101] = 0.0;
        grad[102] = 0.0;
        grad[103] = 0.0;
        grad[104] = 0.0;
        grad[105] = 0.0;
        grad[106] = -x[10] - 7.0/6.0*x[11] + 1.0/6.0*x[12] - 7.0/6.0*x[13] - 1.0;
        grad[107] = 0.0;
        grad[108] = -x[8];
        grad[109] = -7.0/6.0*x[8];
        grad[110] = 1.0/6.0*x[8];
        grad[111] = 1.0 - 7.0/6.0*x[8];
        grad[112] = 0.0;
        grad[113] = 0.0;
        grad[114] = 0.0;
        grad[115] = 0.0;
        grad[116] = 0.0;
        grad[117] = 0.0;
        grad[118] = 0.0;
        grad[119] = 0.0;
        grad[120] = 0.0;
        grad[121] = 0.0;
        grad[122] = 0.0;
        grad[123] = -1.00;
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
        grad[135] = 0.0;
        grad[136] = -1.00;
        grad[137] = 0.0;
        grad[138] = 0.0;
        grad[139] = 0.0;
        grad[140] = 0.0;
        grad[141] = 0.0;
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
        grad[152] = -1.00;
        grad[153] = 0.0;
        grad[154] = 0.0;
        grad[155] = 0.0;
        grad[156] = 0.0;
        grad[157] = 0.0;
        grad[158] = 0.0;
        grad[159] = 0.0;
        grad[160] = 0.0;
        grad[161] = 0.0;
        grad[162] = 0.0;
        grad[163] = 0.0;
        grad[164] = 0.0;
        grad[165] = 0.0;
        grad[166] = 0.0;
        grad[167] = -1.00;
        grad[168] = 0.0;
        grad[169] = 0.0;
        grad[170] = -x[10] - 7.0/6.0*x[11] + 1.0/6.0*x[12] - 7.0/6.0*x[13] - 1.0;
        grad[171] = -x[10] - 7.0/6.0*x[11] + 1.0/6.0*x[12] - 7.0/6.0*x[13] - 1.0;
        grad[172] = 0.0;
        grad[173] = 0.0;
        grad[174] = 0.0;
        grad[175] = 0.0;
        grad[176] = 0.0;
        grad[177] = 0.0;
        grad[178] = -x[2] - x[3];
        grad[179] = -7.0/6.0*x[2] - 7.0/6.0*x[3];
        grad[180] = 1.0/6.0*x[2] + 1.0/6.0*x[3] + 0.5;
        grad[181] = -7.0/6.0*x[2] - 7.0/6.0*x[3];
        grad[182] = 0.0;
        grad[183] = 0.0;
        grad[184] = 0.0;
        grad[185] = 0.0;
        grad[186] = 0.0;
        grad[187] = 0.0;
        grad[188] = 0.0;
        grad[189] = 0.0;
        grad[190] = 0.0;
        grad[191] = x[10] + 7.0/6.0*x[11] - 1.0/6.0*x[12] + 7.0/6.0*x[13] + 1.0;
        grad[192] = x[9];
        grad[193] = 7.0/6.0*x[9];
        grad[194] = -1.0/6.0*x[9];
        grad[195] = 7.0/6.0*x[9];
        grad[196] = 0.0;
        grad[197] = 0.0;
        grad[198] = 0.0;
        grad[199] = 0.0;
        grad[200] = 0.0;
        grad[201] = 0.0;
        grad[202] = 0.0;
        grad[203] = 0.0;
        grad[204] = 0.0;
        grad[205] = -x[10] - 7.0/6.0*x[11] + 1.0/6.0*x[12] - 7.0/6.0*x[13] - 1.0;
        grad[206] = -x[9];
        grad[207] = -7.0/6.0*x[9];
        grad[208] = 1.0/6.0*x[9];
        grad[209] = -7.0/6.0*x[9];
        grad[210] = 0.0;
        grad[211] = 0.0;
        grad[212] = 0.0;
        grad[213] = 0.0;
        grad[214] = 0.0;
        grad[215] = 0.0;
        grad[216] = 0.0;
        grad[217] = 0.0;
        grad[218] = 0.0;
        grad[219] = x[10] + 7.0/6.0*x[11] - 1.0/6.0*x[12] + 7.0/6.0*x[13] + 1.0;
        grad[220] = x[9];
        grad[221] = 7.0/6.0*x[9];
        grad[222] = -1.0/6.0*x[9];
        grad[223] = 7.0/6.0*x[9];
        grad[224] = 0.0;
        grad[225] = 0.0;
        grad[226] = -4.0*x[10] - 4.66666666666667*x[11] + 2.0/3.0*x[12] - 4.66666666666667*x[13] - 4.0;
        grad[227] = 0.0;
        grad[228] = 0.0;
        grad[229] = 0.0;
        grad[230] = 0.0;
        grad[231] = 0.0;
        grad[232] = 0.0;
        grad[233] = 0.0;
        grad[234] = -4.0*x[2];
        grad[235] = -4.66666666666667*x[2];
        grad[236] = 2.0/3.0*x[2] + 2.0;
        grad[237] = -4.66666666666667*x[2];
        grad[238] = 0.0;
        grad[239] = 0.0;
        grad[240] = 0.0;
        grad[241] = -4.0*x[10] - 4.66666666666667*x[11] + 2.0/3.0*x[12] - 4.66666666666667*x[13] - 4.0;
        grad[242] = 0.0;
        grad[243] = 0.0;
        grad[244] = 0.0;
        grad[245] = 0.0;
        grad[246] = 0.0;
        grad[247] = 0.0;
        grad[248] = -4.0*x[3];
        grad[249] = -4.66666666666667*x[3];
        grad[250] = 2.0/3.0*x[3];
        grad[251] = -4.66666666666667*x[3];
        grad[252] = -x[10] - 7.0/6.0*x[11] + 1.0/6.0*x[12] - 7.0/6.0*x[13] - 1.0;
        grad[253] = 0.0;
        grad[254] = 0.0;
        grad[255] = 0.0;
        grad[256] = 0.0;
        grad[257] = 0.0;
        grad[258] = 0.0;
        grad[259] = 0.0;
        grad[260] = 0.0;
        grad[261] = 0.0;
        grad[262] = 1.0 - x[0];
        grad[263] = -7.0/6.0*x[0];
        grad[264] = 1.0/6.0*x[0];
        grad[265] = -7.0/6.0*x[0];
        grad[266] = 0.0;
        grad[267] = -x[10] - 7.0/6.0*x[11] + 1.0/6.0*x[12] - 7.0/6.0*x[13] - 1.0;
        grad[268] = 0.0;
        grad[269] = 0.0;
        grad[270] = 0.0;
        grad[271] = 0.0;
        grad[272] = 0.0;
        grad[273] = 0.0;
        grad[274] = 0.0;
        grad[275] = 0.0;
        grad[276] = 1.0 - x[1];
        grad[277] = 0.5 - 7.0/6.0*x[1];
        grad[278] = 1.0/6.0*x[1];
        grad[279] = 0.5 - 7.0/6.0*x[1];
        grad[280] = -x[10] - 7.0/6.0*x[11] + 1.0/6.0*x[12] - 7.0/6.0*x[13] - 1.0;
        grad[281] = -x[10] - 7.0/6.0*x[11] + 1.0/6.0*x[12] - 7.0/6.0*x[13] - 1.0;
        grad[282] = -4.0*x[10] - 4.66666666666667*x[11] + 2.0/3.0*x[12] - 4.66666666666667*x[13] - 4.0;
        grad[283] = -4.0*x[10] - 4.66666666666667*x[11] + 2.0/3.0*x[12] - 4.66666666666667*x[13] - 4.0;
        grad[284] = 0.0;
        grad[285] = 0.0;
        grad[286] = 0.0;
        grad[287] = 0.0;
        grad[288] = 0.0;
        grad[289] = 0.0;
        grad[290] = -x[0] - x[1] - 4.0*x[2] - 4.0*x[3] + 2.0;
        grad[291] = -7.0/6.0*x[0] - 7.0/6.0*x[1] - 4.66666666666667*x[2] - 4.66666666666667*x[3] + 0.5;
        grad[292] = 1.0/6.0*x[0] + 1.0/6.0*x[1] + 2.0/3.0*x[2] + 2.0/3.0*x[3] + 2.0;
        grad[293] = -7.0/6.0*x[0] - 7.0/6.0*x[1] - 4.66666666666667*x[2] - 4.66666666666667*x[3] + 0.5;
    }

    return;
};

/**
    Inequality constraints for fl
*/
void fl_alk_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( eps_sf + x[0] + x[1] + x[2] - 1.0);
    result[1] = ( eps_sf + -x[0]);
    result[2] = ( eps_sf + -x[1]);
    result[3] = ( eps_sf + -x[2]);
    result[4] = ( eps_sf + x[2] - 1.0);

    if (grad) {
        grad[0] = 1.00;
        grad[1] = 1.00;
        grad[2] = 1.00;
        grad[3] = -1.00;
        grad[4] = 0.0;
        grad[5] = 0.0;
        grad[6] = 0.0;
        grad[7] = -1.00;
        grad[8] = 0.0;
        grad[9] = 0.0;
        grad[10] = 0.0;
        grad[11] = -1.00;
        grad[12] = 0.0;
        grad[13] = 0.0;
        grad[14] = 1.00;
    }

    return;
};

/**
    Inequality constraints for fsp
*/
void fsp_alk_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( eps_sf + x[0] + x[1] - 1.0);
    result[1] = ( eps_sf + -x[0]);
    result[2] = ( eps_sf + -x[1]);
    result[3] = ( eps_sf + -0.25*x[0] - 0.25);
    result[4] = ( eps_sf + 0.25*x[0] - 0.75);

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
    Inequality constraints for spn
*/
void spn_alk_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( eps_sf + 1.0/3.0*x[0]*x[3] + 1.0/3.0*x[0] - 1.0/3.0*x[3] - 2.0/3.0*x[4] - 1.0/3.0);
    result[1] = ( eps_sf + -1.0/3.0*x[0]*x[3] - 1.0/3.0*x[0] - 2.0/3.0*x[5]);
    result[2] = ( eps_sf + -2.0/3.0*x[1]*x[2] - 2.0/3.0*x[1]*x[3] + 2.0/3.0*x[1] + 1.0/3.0*x[3] + 2.0/3.0*x[4] + 2.0/3.0*x[5] + 2.0/3.0*x[6] - 2.0/3.0);
    result[3] = ( eps_sf + 2.0/3.0*x[1]*x[2] + 2.0/3.0*x[1]*x[3] - 2.0/3.0*x[1] - 2.0/3.0*x[6]);
    result[4] = ( eps_sf + 1.0/3.0*x[0]*x[3] + 1.0/3.0*x[0] - 1.0/3.0*x[3] + 1.0/3.0*x[4] - 1.0/3.0);
    result[5] = ( eps_sf + -1.0/3.0*x[0]*x[3] - 1.0/3.0*x[0] + 1.0/3.0*x[5]);
    result[6] = ( eps_sf + -2.0/3.0*x[1]*x[2] - 2.0/3.0*x[1]*x[3] + 2.0/3.0*x[1] + x[2] + 0.833333333333333*x[3] - 1.0/3.0*x[4] - 1.0/3.0*x[5] - 1.0/3.0*x[6] - 2.0/3.0);
    result[7] = ( eps_sf + 2.0/3.0*x[1]*x[2] + 2.0/3.0*x[1]*x[3] - 2.0/3.0*x[1] + 1.0/3.0*x[6]);
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
        grad[45] = 0.833333333333333 - 2.0/3.0*x[1];
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
void g_alk_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( eps_sf + -x[0]*x[1] + x[0] + x[1] - 1.0);
    result[1] = ( eps_sf + x[0]*x[1] - x[0]);
    result[2] = ( eps_sf + -x[1]);
    result[3] = ( eps_sf + x[2] + x[3] + 2.0*x[4] - 1.0);
    result[4] = ( eps_sf + -x[3]);
    result[5] = ( eps_sf + -x[2]);
    result[6] = ( eps_sf + -x[4]);
    result[7] = ( eps_sf + -x[4]);

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
void ol_alk_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( eps_sf + x[0] - x[2] - 1.0);
    result[1] = ( eps_sf + -x[0] + x[2]);
    result[2] = ( eps_sf + -x[0]*x[1] + x[0] + x[1] + x[2] - 1.0);
    result[3] = ( eps_sf + x[0]*x[1] - x[0] - x[2]);
    result[4] = ( eps_sf + -x[1]);

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
void opx_alk_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( eps_sf + -x[0]*x[1] + x[0]*x[5] - x[0]*x[7] + x[0] - x[1]*x[3] + x[1] + x[3]*x[5] - x[3]*x[7] + x[3] - x[5] + x[7] - 1.0);
    result[1] = ( eps_sf + x[0]*x[1] - x[0]*x[5] + x[0]*x[7] - x[0] + x[1]*x[3] - x[3]*x[5] + x[3]*x[7] - x[3]);
    result[2] = ( eps_sf + -x[1] + x[4] + 2.0*x[5] + x[6] - x[7]);
    result[3] = ( eps_sf + -x[4]);
    result[4] = ( eps_sf + -x[6]);
    result[5] = ( eps_sf + -x[5]);
    result[6] = ( eps_sf + -x[0]*x[2] - x[0]*x[7] + x[0] + x[1]*x[3] + x[2] - x[3]*x[5] + x[3]*x[7] - x[3] + x[7] - 1.0);
    result[7] = ( eps_sf + x[0]*x[2] + x[0]*x[7] - x[0] - x[1]*x[3] + x[3]*x[5] - x[3]*x[7] + x[3]);
    result[8] = ( eps_sf + -x[2]);
    result[9] = ( eps_sf + -x[7]);
    result[10] = ( eps_sf + 0.5*x[1] - 1.0);
    result[11] = ( eps_sf + -0.5*x[1]);

    if (grad) {
        grad[0] = -x[1] + x[5] - x[7] + 1.0;
        grad[1] = -x[0] - x[3] + 1.0;
        grad[2] = 0.0;
        grad[3] = -x[1] + x[5] - x[7] + 1.0;
        grad[4] = 0.0;
        grad[5] = x[0] + x[3] - 1.0;
        grad[6] = 0.0;
        grad[7] = -x[0] - x[3] + 1.0;
        grad[8] = x[1] - x[5] + x[7] - 1.0;
        grad[9] = x[0] + x[3];
        grad[10] = 0.0;
        grad[11] = x[1] - x[5] + x[7] - 1.0;
        grad[12] = 0.0;
        grad[13] = -x[0] - x[3];
        grad[14] = 0.0;
        grad[15] = x[0] + x[3];
        grad[16] = 0.0;
        grad[17] = -1.00;
        grad[18] = 0.0;
        grad[19] = 0.0;
        grad[20] = 1.00;
        grad[21] = 2.00;
        grad[22] = 1.00;
        grad[23] = -1.00;
        grad[24] = 0.0;
        grad[25] = 0.0;
        grad[26] = 0.0;
        grad[27] = 0.0;
        grad[28] = -1.00;
        grad[29] = 0.0;
        grad[30] = 0.0;
        grad[31] = 0.0;
        grad[32] = 0.0;
        grad[33] = 0.0;
        grad[34] = 0.0;
        grad[35] = 0.0;
        grad[36] = 0.0;
        grad[37] = 0.0;
        grad[38] = -1.00;
        grad[39] = 0.0;
        grad[40] = 0.0;
        grad[41] = 0.0;
        grad[42] = 0.0;
        grad[43] = 0.0;
        grad[44] = 0.0;
        grad[45] = -1.00;
        grad[46] = 0.0;
        grad[47] = 0.0;
        grad[48] = -x[2] - x[7] + 1.0;
        grad[49] = x[3];
        grad[50] = 1.0 - x[0];
        grad[51] = x[1] - x[5] + x[7] - 1.0;
        grad[52] = 0.0;
        grad[53] = -x[3];
        grad[54] = 0.0;
        grad[55] = -x[0] + x[3] + 1.0;
        grad[56] = x[2] + x[7] - 1.0;
        grad[57] = -x[3];
        grad[58] = x[0];
        grad[59] = -x[1] + x[5] - x[7] + 1.0;
        grad[60] = 0.0;
        grad[61] = x[3];
        grad[62] = 0.0;
        grad[63] = x[0] - x[3];
        grad[64] = 0.0;
        grad[65] = 0.0;
        grad[66] = -1.00;
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
        grad[79] = -1.00;
        grad[80] = 0.0;
        grad[81] = 0.500;
        grad[82] = 0.0;
        grad[83] = 0.0;
        grad[84] = 0.0;
        grad[85] = 0.0;
        grad[86] = 0.0;
        grad[87] = 0.0;
        grad[88] = 0.0;
        grad[89] = -0.500;
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
void cpx_alk_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( eps_sf + -x[0]*x[1] - x[0]*x[3] + x[0]*x[7] - x[0]*x[8] + x[0] - x[1]*x[4] + x[1] - x[3]*x[4] + x[3] + x[4]*x[7] - x[4]*x[8] + x[4] - x[7] + x[8] - 1.0);
    result[1] = ( eps_sf + x[0]*x[1] + x[0]*x[3] - x[0]*x[7] + x[0]*x[8] - x[0] + x[1]*x[4] + x[3]*x[4] - x[4]*x[7] + x[4]*x[8] - x[4]);
    result[2] = ( eps_sf + -x[1] - x[3] + x[5] + x[6] + 2.0*x[7] - x[8]);
    result[3] = ( eps_sf + -x[5]);
    result[4] = ( eps_sf + -x[6]);
    result[5] = ( eps_sf + -x[7]);
    result[6] = ( eps_sf + x[0]*x[2] + x[1]*x[4] - x[2] + x[3]*x[4] - x[4]*x[7] + x[4]*x[8] - x[4]);
    result[7] = ( eps_sf + -x[0]*x[2] - x[1]*x[4] - x[3]*x[4] + x[4]*x[7] - x[4]*x[8] + x[4]);
    result[8] = ( eps_sf + x[2] + x[3] + x[8] - 1.0);
    result[9] = ( eps_sf + -x[3]);
    result[10] = ( eps_sf + -x[8]);
    result[11] = ( eps_sf + 0.5*x[1] - 1.0);
    result[12] = ( eps_sf + -0.5*x[1]);

    if (grad) {
        grad[0] = -x[1] - x[3] + x[7] - x[8] + 1.0;
        grad[1] = -x[0] - x[4] + 1.0;
        grad[2] = 0.0;
        grad[3] = -x[0] - x[4] + 1.0;
        grad[4] = -x[1] - x[3] + x[7] - x[8] + 1.0;
        grad[5] = 0.0;
        grad[6] = 0.0;
        grad[7] = x[0] + x[4] - 1.0;
        grad[8] = -x[0] - x[4] + 1.0;
        grad[9] = x[1] + x[3] - x[7] + x[8] - 1.0;
        grad[10] = x[0] + x[4];
        grad[11] = 0.0;
        grad[12] = x[0] + x[4];
        grad[13] = x[1] + x[3] - x[7] + x[8] - 1.0;
        grad[14] = 0.0;
        grad[15] = 0.0;
        grad[16] = -x[0] - x[4];
        grad[17] = x[0] + x[4];
        grad[18] = 0.0;
        grad[19] = -1.00;
        grad[20] = 0.0;
        grad[21] = -1.00;
        grad[22] = 0.0;
        grad[23] = 1.00;
        grad[24] = 1.00;
        grad[25] = 2.00;
        grad[26] = -1.00;
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
        grad[42] = -1.00;
        grad[43] = 0.0;
        grad[44] = 0.0;
        grad[45] = 0.0;
        grad[46] = 0.0;
        grad[47] = 0.0;
        grad[48] = 0.0;
        grad[49] = 0.0;
        grad[50] = 0.0;
        grad[51] = 0.0;
        grad[52] = -1.00;
        grad[53] = 0.0;
        grad[54] = x[2];
        grad[55] = x[4];
        grad[56] = x[0] - 1.0;
        grad[57] = x[4];
        grad[58] = x[1] + x[3] - x[7] + x[8] - 1.0;
        grad[59] = 0.0;
        grad[60] = 0.0;
        grad[61] = -x[4];
        grad[62] = x[4];
        grad[63] = -x[2];
        grad[64] = -x[4];
        grad[65] = -x[0];
        grad[66] = -x[4];
        grad[67] = -x[1] - x[3] + x[7] - x[8] + 1.0;
        grad[68] = 0.0;
        grad[69] = 0.0;
        grad[70] = x[4];
        grad[71] = -x[4];
        grad[72] = 0.0;
        grad[73] = 0.0;
        grad[74] = 1.00;
        grad[75] = 1.00;
        grad[76] = 0.0;
        grad[77] = 0.0;
        grad[78] = 0.0;
        grad[79] = 0.0;
        grad[80] = 1.00;
        grad[81] = 0.0;
        grad[82] = 0.0;
        grad[83] = 0.0;
        grad[84] = -1.00;
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
        grad[98] = -1.00;
        grad[99] = 0.0;
        grad[100] = 0.500;
        grad[101] = 0.0;
        grad[102] = 0.0;
        grad[103] = 0.0;
        grad[104] = 0.0;
        grad[105] = 0.0;
        grad[106] = 0.0;
        grad[107] = 0.0;
        grad[108] = 0.0;
        grad[109] = -0.500;
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
void ilm_alk_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( eps_sf + 0.5*x[0]*x[1] - 0.5*x[0] - 0.5*x[2]);
    result[1] = ( eps_sf + -0.5*x[0] + 0.5*x[3]);
    result[2] = ( eps_sf + x[0] - 1.0);
    result[3] = ( eps_sf + -0.5*x[0]*x[1] + 0.5*x[2] - 0.5*x[3]);
    result[4] = ( eps_sf + 0.5*x[0]*x[1] - 0.5*x[0] + 0.5*x[2]);
    result[5] = ( eps_sf + -0.5*x[0] - 0.5*x[3]);
    result[6] = ( eps_sf + x[0] - 1.0);
    result[7] = ( eps_sf + -0.5*x[0]*x[1] - 0.5*x[2] + 0.5*x[3]);

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
    Inequality constraints for ness
*/
void ness_alk_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( eps_sf + -0.25*x[0]*x[1] - 0.75*x[1]*x[4] + x[1] - 0.25*x[2] + x[4] - 1.0);
    result[1] = ( eps_sf + 0.25*x[0]*x[1] + 0.75*x[1]*x[4] - x[1] + 0.25*x[2]);
    result[2] = ( eps_sf + -x[4]);
    result[3] = ( eps_sf + -0.25*x[0]*x[1] + x[0] - 0.75*x[1]*x[4] + x[1] + 0.75*x[2] - 1.0);
    result[4] = ( eps_sf + 0.25*x[0]*x[1] + 0.75*x[1]*x[4] - x[1] - 0.75*x[2]);
    result[5] = ( eps_sf + -x[0]);
    result[6] = ( eps_sf + 0.25*x[0] + x[3] - 0.75*x[4] - 1.0);
    result[7] = ( eps_sf + -0.25*x[0] + 0.75*x[4]);
    result[8] = ( eps_sf + -x[3]);

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
void lct_alk_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( eps_sf + -x[0]);
    result[1] = ( eps_sf + x[0] - 1.0);

    if (grad) {
        grad[0] = -1.00;
        grad[1] = 1.00;
    }

    return;
};

/**
    Inequality constraints for kals
*/
void kals_alk_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( eps_sf + -x[0]);
    result[1] = ( eps_sf + x[0] - 1.0);

    if (grad) {
        grad[0] = -1.00;
        grad[1] = 1.00;
    }

    return;
};

/**
    Inequality constraints for mel
*/
void mel_alk_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( eps_sf + -x[1]);
    result[1] = ( eps_sf + x[1] - 1.0);
    result[2] = ( eps_sf + -x[0]*x[2] - x[0]*x[3] + x[0] + x[2] + x[3] - 1.0);
    result[3] = ( eps_sf + x[0]*x[2] + x[0]*x[3] - x[0]);
    result[4] = ( eps_sf + -x[2]);
    result[5] = ( eps_sf + -x[3]);
    result[6] = ( eps_sf + x[1] - 0.5*x[2] - 0.5*x[3]);
    result[7] = ( eps_sf + -x[1] + 0.5*x[2] + 0.5*x[3] - 1.0);

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

/**
    Inequality constraints for hb
*/
void hb_alk_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( eps_sf + x[3] - 1.0);
    result[1] = ( eps_sf + x[3]*x[4] - x[3]);
    result[2] = ( eps_sf + -x[3]*x[4]);
    result[3] = ( eps_sf + x[0] - x[8] - 1.0);
    result[4] = ( eps_sf + -x[0] + x[8]);
    result[5] = ( eps_sf + -x[0]*x[1] - x[0]*x[6] - x[0]*x[7] + x[0] + x[1]*x[9] + x[1] + x[6]*x[9] + x[6] + x[7]*x[9] + x[7] - x[9] - 1.0);
    result[6] = ( eps_sf + x[0]*x[1] + x[0]*x[6] + x[0]*x[7] - x[0] - x[1]*x[9] - x[6]*x[9] - x[7]*x[9] + x[9]);
    result[7] = ( eps_sf + -x[1]);
    result[8] = ( eps_sf + -x[6]);
    result[9] = ( eps_sf + -x[7]);
    result[10] = ( eps_sf + -x[5]);
    result[11] = ( eps_sf + -x[0]*x[2] - x[0]*x[5] + x[0] - x[1]*x[9] + x[2] + x[5] - x[6]*x[9] - x[7]*x[9] + 1.5*x[8] + x[9] - 1.0);
    result[12] = ( eps_sf + x[0]*x[2] + x[0]*x[5] - x[0] + x[1]*x[9] + x[6]*x[9] + x[7]*x[9] - 1.5*x[8] - x[9]);
    result[13] = ( eps_sf + -x[2]);
    result[14] = ( eps_sf + 0.5*x[1] - 0.5*x[2] + 0.25*x[3] + 0.5*x[6] + 0.5*x[7] - 1.0);
    result[15] = ( eps_sf + -0.5*x[1] + 0.5*x[2] - 0.25*x[3] - 0.5*x[6] - 0.5*x[7]);
    result[16] = ( eps_sf + x[7] - 1.0);
    result[17] = ( eps_sf + -x[7]);

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
        grad[110] = -x[2] - x[5] + 1.0;
        grad[111] = -x[9];
        grad[112] = 1.0 - x[0];
        grad[113] = 0.0;
        grad[114] = 0.0;
        grad[115] = 1.0 - x[0];
        grad[116] = -x[9];
        grad[117] = -x[9];
        grad[118] = 1.50;
        grad[119] = -x[1] - x[6] - x[7] + 1.0;
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
    Inequality constraints for bi
*/
void bi_alk_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( eps_sf + -x[0]*x[1] - x[0]*x[2] - x[0]*x[3] + x[0] + x[1] + x[2] + x[3] + 2.0/3.0*x[4] - 1.0);
    result[1] = ( eps_sf + x[0]*x[1] + x[0]*x[2] + x[0]*x[3] - x[0] - 2.0/3.0*x[4]);
    result[2] = ( eps_sf + -x[2]);
    result[3] = ( eps_sf + -x[3]);
    result[4] = ( eps_sf + -x[1]);
    result[5] = ( eps_sf + x[0] - 1.0/3.0*x[4] - 1.0);
    result[6] = ( eps_sf + -x[0] + 1.0/3.0*x[4]);
    result[7] = ( eps_sf + 0.5*x[1] + 0.5*x[2] - 0.5);
    result[8] = ( eps_sf + -0.5*x[1] - 0.5*x[2] - 0.5);
    result[9] = ( eps_sf + x[3] - 1.0);
    result[10] = ( eps_sf + -x[3]);

    if (grad) {
        grad[0] = -x[1] - x[2] - x[3] + 1.0;
        grad[1] = 1.0 - x[0];
        grad[2] = 1.0 - x[0];
        grad[3] = 1.0 - x[0];
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
    Inequality constraints for ep
*/
void ep_alk_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( eps_sf + -x[0] + x[1]);
    result[1] = ( eps_sf + x[0] - x[1] - 1.0);
    result[2] = ( eps_sf + -x[0] - x[1]);
    result[3] = ( eps_sf + x[0] + x[1] - 1.0);

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
    Inequality constraints for cd
*/
void cd_alk_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( eps_sf + -x[0]);
    result[1] = ( eps_sf + x[0] - 1.0);
    result[2] = ( eps_sf + -x[1]);
    result[3] = ( eps_sf + x[1] - 1.0);

    if (grad) {
        grad[0] = -1.00;
        grad[1] = 0.0;
        grad[2] = 1.00;
        grad[3] = 0.0;
        grad[4] = 0.0;
        grad[5] = -1.00;
        grad[6] = 0.0;
        grad[7] = 1.00;
    }

    return;
};

/**
    Inequality constraints for liq
*/
void liq_igd_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( eps_sf + x[0] - x[10]*(-x[0] - x[1] - x[2] - x[3] - x[4] - x[5] - x[6] - x[7] - x[8] - x[9] + 1.0) - 2./3.*x[11]*(-x[0] - x[1] - x[2] - x[3] - x[4] - x[5] - x[6] - x[7] - x[8] - x[9]) - 1./6.*x[12]*(x[0] + x[1] + x[2] + x[3] + x[4] + x[5] + x[6] + x[7] + x[8] + x[9] - 3.0) - 1./3.*x[13]*(-x[0] - x[1] - x[2] - x[3] - x[4] - x[5] - x[6] - x[7] - x[8] - x[9]) + x[1] + x[2] + x[3] + x[4] + x[5] + x[6] + x[7] + x[8] + x[9] - 1.0);
    result[1] = ( eps_sf + -x[10]*x[1] + x[10] - 2./3.*x[11]*x[1] + 1./6.*x[12]*x[1] - 1./3.*x[13]*x[1] - x[1]);
    result[2] = ( eps_sf + -x[0]*x[10] - 2./3.*x[0]*x[11] + 1./6.*x[0]*x[12] - 1./3.*x[0]*x[13] - x[0] + x[10]);
    result[3] = ( eps_sf + -x[10]*x[4] - 2./3.*x[11]*x[4] + x[11] + 1./6.*x[12]*x[4] - 1./3.*x[13]*x[4] - x[4]);
    result[4] = ( eps_sf + -x[10]*x[5] - 2./3.*x[11]*x[5] + 1./6.*x[12]*x[5] - 1./3.*x[13]*x[5] - x[5]);
    result[5] = ( eps_sf + -x[10]*x[6] - 2./3.*x[11]*x[6] + 1./6.*x[12]*x[6] - 1./3.*x[13]*x[6] - x[6]);
    result[6] = ( eps_sf + -x[10]*x[7] - 2./3.*x[11]*x[7] + 1./6.*x[12]*x[7] - 1./3.*x[13]*x[7] - x[7]);
    result[7] = ( eps_sf + -x[10]*x[8] - 2./3.*x[11]*x[8] + 1./6.*x[12]*x[8] - 1./3.*x[13]*x[8] + x[13] - x[8]);
    result[8] = ( eps_sf + -x[11]);
    result[9] = ( eps_sf + -x[10]);
    result[10] = ( eps_sf + -x[12]);
    result[11] = ( eps_sf + -x[13]);
    result[12] = ( eps_sf + -x[10]*(x[2] + x[3]) - 2./3.*x[11]*(x[2] + x[3]) + 1./6.*x[12]*(x[2] + x[3]) + 0.5*x[12] - 1./3.*x[13]*(x[2] + x[3]) - x[2] - x[3]);
    result[13] = ( eps_sf + x[10]*x[9] + 2./3.*x[11]*x[9] - 1./6.*x[12]*x[9] + 1./3.*x[13]*x[9] + x[9] - 1.0);
    result[14] = ( eps_sf + -x[10]*x[9] - 2./3.*x[11]*x[9] + 1./6.*x[12]*x[9] - 1./3.*x[13]*x[9] - x[9]);
    result[15] = ( eps_sf + x[10]*x[9] + 2./3.*x[11]*x[9] - 1./6.*x[12]*x[9] + 1./3.*x[13]*x[9] + x[9] - 1.0);
    result[16] = ( eps_sf + 2.0*x[12] - 4.0*x[2]*(x[10] + 2./3.*x[11] - 1./6.*x[12] + 1./3.*x[13] + 1.0));
    result[17] = ( eps_sf + -4.0*x[3]*(x[10] + 2./3.*x[11] - 1./6.*x[12] + 1./3.*x[13] + 1.0));
    result[18] = ( eps_sf + -x[0]*(x[10] + 2./3.*x[11] - 1./6.*x[12] + 1./3.*x[13] + 1.0) + x[10]);
    result[19] = ( eps_sf + x[10] - x[1]*(x[10] + 2./3.*x[11] - 1./6.*x[12] + 1./3.*x[13] + 1.0));
    result[20] = ( eps_sf + 2.0*x[10] + 2.0*x[12] - (x[0] + x[1] + 4.0*x[2] + 4.0*x[3])*(x[10] + 2./3.*x[11] - 1./6.*x[12] + 1./3.*x[13] + 1.0));

    if (grad) {
        grad[0] = x[10] + 2./3.*x[11] - 1./6.*x[12] + 1./3.*x[13] + 1.0;
        grad[1] = x[10] + 2./3.*x[11] - 1./6.*x[12] + 1./3.*x[13] + 1.0;
        grad[2] = x[10] + 2./3.*x[11] - 1./6.*x[12] + 1./3.*x[13] + 1.0;
        grad[3] = x[10] + 2./3.*x[11] - 1./6.*x[12] + 1./3.*x[13] + 1.0;
        grad[4] = x[10] + 2./3.*x[11] - 1./6.*x[12] + 1./3.*x[13] + 1.0;
        grad[5] = x[10] + 2./3.*x[11] - 1./6.*x[12] + 1./3.*x[13] + 1.0;
        grad[6] = x[10] + 2./3.*x[11] - 1./6.*x[12] + 1./3.*x[13] + 1.0;
        grad[7] = x[10] + 2./3.*x[11] - 1./6.*x[12] + 1./3.*x[13] + 1.0;
        grad[8] = x[10] + 2./3.*x[11] - 1./6.*x[12] + 1./3.*x[13] + 1.0;
        grad[9] = x[10] + 2./3.*x[11] - 1./6.*x[12] + 1./3.*x[13] + 1.0;
        grad[10] = x[0] + x[1] + x[2] + x[3] + x[4] + x[5] + x[6] + x[7] + x[8] + x[9] - 1.0;
        grad[11] = 2./3.*x[0] + 2./3.*x[1] + 2./3.*x[2] + 2./3.*x[3] + 2./3.*x[4] + 2./3.*x[5] + 2./3.*x[6] + 2./3.*x[7] + 2./3.*x[8] + 2./3.*x[9];
        grad[12] = -1./6.*x[0] - 1./6.*x[1] - 1./6.*x[2] - 1./6.*x[3] - 1./6.*x[4] - 1./6.*x[5] - 1./6.*x[6] - 1./6.*x[7] - 1./6.*x[8] - 1./6.*x[9] + 0.5;
        grad[13] = 1./3.*x[0] + 1./3.*x[1] + 1./3.*x[2] + 1./3.*x[3] + 1./3.*x[4] + 1./3.*x[5] + 1./3.*x[6] + 1./3.*x[7] + 1./3.*x[8] + 1./3.*x[9];
        grad[14] = 0.0;
        grad[15] = -x[10] - 2./3.*x[11] + 1./6.*x[12] - 1./3.*x[13] - 1.0;
        grad[16] = 0.0;
        grad[17] = 0.0;
        grad[18] = 0.0;
        grad[19] = 0.0;
        grad[20] = 0.0;
        grad[21] = 0.0;
        grad[22] = 0.0;
        grad[23] = 0.0;
        grad[24] = 1.0 - x[1];
        grad[25] = -2./3.*x[1];
        grad[26] = 1./6.*x[1];
        grad[27] = -1./3.*x[1];
        grad[28] = -x[10] - 2./3.*x[11] + 1./6.*x[12] - 1./3.*x[13] - 1.0;
        grad[29] = 0.0;
        grad[30] = 0.0;
        grad[31] = 0.0;
        grad[32] = 0.0;
        grad[33] = 0.0;
        grad[34] = 0.0;
        grad[35] = 0.0;
        grad[36] = 0.0;
        grad[37] = 0.0;
        grad[38] = 1.0 - x[0];
        grad[39] = -2./3.*x[0];
        grad[40] = 1./6.*x[0];
        grad[41] = -1./3.*x[0];
        grad[42] = 0.0;
        grad[43] = 0.0;
        grad[44] = 0.0;
        grad[45] = 0.0;
        grad[46] = -x[10] - 2./3.*x[11] + 1./6.*x[12] - 1./3.*x[13] - 1.0;
        grad[47] = 0.0;
        grad[48] = 0.0;
        grad[49] = 0.0;
        grad[50] = 0.0;
        grad[51] = 0.0;
        grad[52] = -x[4];
        grad[53] = 1.0 - 2./3.*x[4];
        grad[54] = 1./6.*x[4];
        grad[55] = -1./3.*x[4];
        grad[56] = 0.0;
        grad[57] = 0.0;
        grad[58] = 0.0;
        grad[59] = 0.0;
        grad[60] = 0.0;
        grad[61] = -x[10] - 2./3.*x[11] + 1./6.*x[12] - 1./3.*x[13] - 1.0;
        grad[62] = 0.0;
        grad[63] = 0.0;
        grad[64] = 0.0;
        grad[65] = 0.0;
        grad[66] = -x[5];
        grad[67] = -2./3.*x[5];
        grad[68] = 1./6.*x[5];
        grad[69] = -1./3.*x[5];
        grad[70] = 0.0;
        grad[71] = 0.0;
        grad[72] = 0.0;
        grad[73] = 0.0;
        grad[74] = 0.0;
        grad[75] = 0.0;
        grad[76] = -x[10] - 2./3.*x[11] + 1./6.*x[12] - 1./3.*x[13] - 1.0;
        grad[77] = 0.0;
        grad[78] = 0.0;
        grad[79] = 0.0;
        grad[80] = -x[6];
        grad[81] = -2./3.*x[6];
        grad[82] = 1./6.*x[6];
        grad[83] = -1./3.*x[6];
        grad[84] = 0.0;
        grad[85] = 0.0;
        grad[86] = 0.0;
        grad[87] = 0.0;
        grad[88] = 0.0;
        grad[89] = 0.0;
        grad[90] = 0.0;
        grad[91] = -x[10] - 2./3.*x[11] + 1./6.*x[12] - 1./3.*x[13] - 1.0;
        grad[92] = 0.0;
        grad[93] = 0.0;
        grad[94] = -x[7];
        grad[95] = -2./3.*x[7];
        grad[96] = 1./6.*x[7];
        grad[97] = -1./3.*x[7];
        grad[98] = 0.0;
        grad[99] = 0.0;
        grad[100] = 0.0;
        grad[101] = 0.0;
        grad[102] = 0.0;
        grad[103] = 0.0;
        grad[104] = 0.0;
        grad[105] = 0.0;
        grad[106] = -x[10] - 2./3.*x[11] + 1./6.*x[12] - 1./3.*x[13] - 1.0;
        grad[107] = 0.0;
        grad[108] = -x[8];
        grad[109] = -2./3.*x[8];
        grad[110] = 1./6.*x[8];
        grad[111] = 1.0 - 1./3.*x[8];
        grad[112] = 0.0;
        grad[113] = 0.0;
        grad[114] = 0.0;
        grad[115] = 0.0;
        grad[116] = 0.0;
        grad[117] = 0.0;
        grad[118] = 0.0;
        grad[119] = 0.0;
        grad[120] = 0.0;
        grad[121] = 0.0;
        grad[122] = 0.0;
        grad[123] = -1.0;
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
        grad[135] = 0.0;
        grad[136] = -1.0;
        grad[137] = 0.0;
        grad[138] = 0.0;
        grad[139] = 0.0;
        grad[140] = 0.0;
        grad[141] = 0.0;
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
        grad[152] = -1.0;
        grad[153] = 0.0;
        grad[154] = 0.0;
        grad[155] = 0.0;
        grad[156] = 0.0;
        grad[157] = 0.0;
        grad[158] = 0.0;
        grad[159] = 0.0;
        grad[160] = 0.0;
        grad[161] = 0.0;
        grad[162] = 0.0;
        grad[163] = 0.0;
        grad[164] = 0.0;
        grad[165] = 0.0;
        grad[166] = 0.0;
        grad[167] = -1.0;
        grad[168] = 0.0;
        grad[169] = 0.0;
        grad[170] = -x[10] - 2./3.*x[11] + 1./6.*x[12] - 1./3.*x[13] - 1.0;
        grad[171] = -x[10] - 2./3.*x[11] + 1./6.*x[12] - 1./3.*x[13] - 1.0;
        grad[172] = 0.0;
        grad[173] = 0.0;
        grad[174] = 0.0;
        grad[175] = 0.0;
        grad[176] = 0.0;
        grad[177] = 0.0;
        grad[178] = -x[2] - x[3];
        grad[179] = -2./3.*x[2] - 2./3.*x[3];
        grad[180] = 1./6.*x[2] + 1./6.*x[3] + 0.5;
        grad[181] = -1./3.*x[2] - 1./3.*x[3];
        grad[182] = 0.0;
        grad[183] = 0.0;
        grad[184] = 0.0;
        grad[185] = 0.0;
        grad[186] = 0.0;
        grad[187] = 0.0;
        grad[188] = 0.0;
        grad[189] = 0.0;
        grad[190] = 0.0;
        grad[191] = x[10] + 2./3.*x[11] - 1./6.*x[12] + 1./3.*x[13] + 1.0;
        grad[192] = x[9];
        grad[193] = 2./3.*x[9];
        grad[194] = -1./6.*x[9];
        grad[195] = 1./3.*x[9];
        grad[196] = 0.0;
        grad[197] = 0.0;
        grad[198] = 0.0;
        grad[199] = 0.0;
        grad[200] = 0.0;
        grad[201] = 0.0;
        grad[202] = 0.0;
        grad[203] = 0.0;
        grad[204] = 0.0;
        grad[205] = -x[10] - 2./3.*x[11] + 1./6.*x[12] - 1./3.*x[13] - 1.0;
        grad[206] = -x[9];
        grad[207] = -2./3.*x[9];
        grad[208] = 1./6.*x[9];
        grad[209] = -1./3.*x[9];
        grad[210] = 0.0;
        grad[211] = 0.0;
        grad[212] = 0.0;
        grad[213] = 0.0;
        grad[214] = 0.0;
        grad[215] = 0.0;
        grad[216] = 0.0;
        grad[217] = 0.0;
        grad[218] = 0.0;
        grad[219] = x[10] + 2./3.*x[11] - 1./6.*x[12] + 1./3.*x[13] + 1.0;
        grad[220] = x[9];
        grad[221] = 2./3.*x[9];
        grad[222] = -1./6.*x[9];
        grad[223] = 1./3.*x[9];
        grad[224] = 0.0;
        grad[225] = 0.0;
        grad[226] = -4.0*x[10] - 2.66666666666667*x[11] + 2./3.*x[12] - 1.33333333333333*x[13] - 4.0;
        grad[227] = 0.0;
        grad[228] = 0.0;
        grad[229] = 0.0;
        grad[230] = 0.0;
        grad[231] = 0.0;
        grad[232] = 0.0;
        grad[233] = 0.0;
        grad[234] = -4.0*x[2];
        grad[235] = -2.66666666666667*x[2];
        grad[236] = 2./3.*x[2] + 2.0;
        grad[237] = -1.33333333333333*x[2];
        grad[238] = 0.0;
        grad[239] = 0.0;
        grad[240] = 0.0;
        grad[241] = -4.0*x[10] - 2.66666666666667*x[11] + 2./3.*x[12] - 1.33333333333333*x[13] - 4.0;
        grad[242] = 0.0;
        grad[243] = 0.0;
        grad[244] = 0.0;
        grad[245] = 0.0;
        grad[246] = 0.0;
        grad[247] = 0.0;
        grad[248] = -4.0*x[3];
        grad[249] = -2.66666666666667*x[3];
        grad[250] = 2./3.*x[3];
        grad[251] = -1.33333333333333*x[3];
        grad[252] = -x[10] - 2./3.*x[11] + 1./6.*x[12] - 1./3.*x[13] - 1.0;
        grad[253] = 0.0;
        grad[254] = 0.0;
        grad[255] = 0.0;
        grad[256] = 0.0;
        grad[257] = 0.0;
        grad[258] = 0.0;
        grad[259] = 0.0;
        grad[260] = 0.0;
        grad[261] = 0.0;
        grad[262] = 1.0 - x[0];
        grad[263] = -2./3.*x[0];
        grad[264] = 1./6.*x[0];
        grad[265] = -1./3.*x[0];
        grad[266] = 0.0;
        grad[267] = -x[10] - 2./3.*x[11] + 1./6.*x[12] - 1./3.*x[13] - 1.0;
        grad[268] = 0.0;
        grad[269] = 0.0;
        grad[270] = 0.0;
        grad[271] = 0.0;
        grad[272] = 0.0;
        grad[273] = 0.0;
        grad[274] = 0.0;
        grad[275] = 0.0;
        grad[276] = 1.0 - x[1];
        grad[277] = -2./3.*x[1];
        grad[278] = 1./6.*x[1];
        grad[279] = -1./3.*x[1];
        grad[280] = -x[10] - 2./3.*x[11] + 1./6.*x[12] - 1./3.*x[13] - 1.0;
        grad[281] = -x[10] - 2./3.*x[11] + 1./6.*x[12] - 1./3.*x[13] - 1.0;
        grad[282] = -4.0*x[10] - 2.66666666666667*x[11] + 2./3.*x[12] - 1.33333333333333*x[13] - 4.0;
        grad[283] = -4.0*x[10] - 2.66666666666667*x[11] + 2./3.*x[12] - 1.33333333333333*x[13] - 4.0;
        grad[284] = 0.0;
        grad[285] = 0.0;
        grad[286] = 0.0;
        grad[287] = 0.0;
        grad[288] = 0.0;
        grad[289] = 0.0;
        grad[290] = -x[0] - x[1] - 4.0*x[2] - 4.0*x[3] + 2.0;
        grad[291] = -2./3.*x[0] - 2./3.*x[1] - 2.66666666666667*x[2] - 2.66666666666667*x[3];
        grad[292] = 1./6.*x[0] + 1./6.*x[1] + 2./3.*x[2] + 2./3.*x[3] + 2.0;
        grad[293] = -1./3.*x[0] - 1./3.*x[1] - 1.33333333333333*x[2] - 1.33333333333333*x[3];
    }

    return;
};

/**
    Inequality constraints for fl
*/
void fl_igd_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( eps_sf + x[0] + x[1] + x[2] - 1.0);
    result[1] = ( eps_sf + -x[0]);
    result[2] = ( eps_sf + -x[1]);
    result[3] = ( eps_sf + -x[2]);
    result[4] = ( eps_sf + x[2] - 1.0);

    if (grad) {
        grad[0] = 1.0;
        grad[1] = 1.0;
        grad[2] = 1.0;
        grad[3] = -1.0;
        grad[4] = 0.0;
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
    Inequality constraints for fsp
*/
void fsp_igd_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
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
        grad[6] = -0.2500;
        grad[7] = 0.0;
        grad[8] = 0.2500;
        grad[9] = 0.0;
    }

    return;
};

/**
    Inequality constraints for spn
*/
void spn_igd_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( eps_sf + 1./3.*x[0]*x[3] + 1./3.*x[0] - 1./3.*x[3] - 2./3.*x[4] - 1./3.);
    result[1] = ( eps_sf + -1./3.*x[0]*x[3] - 1./3.*x[0] - 2./3.*x[5]);
    result[2] = ( eps_sf + -2./3.*x[1]*x[2] - 2./3.*x[1]*x[3] + 2./3.*x[1] + 1./3.*x[3] + 2./3.*x[4] + 2./3.*x[5] + 2./3.*x[6] - 2./3.);
    result[3] = ( eps_sf + 2./3.*x[1]*x[2] + 2./3.*x[1]*x[3] - 2./3.*x[1] - 2./3.*x[6]);
    result[4] = ( eps_sf + 1./3.*x[0]*x[3] + 1./3.*x[0] - 1./3.*x[3] + 1./3.*x[4] - 1./3.);
    result[5] = ( eps_sf + -1./3.*x[0]*x[3] - 1./3.*x[0] + 1./3.*x[5]);
    result[6] = ( eps_sf + -2./3.*x[1]*x[2] - 2./3.*x[1]*x[3] + 2./3.*x[1] + x[2] + 0.833333333333333*x[3] - 1./3.*x[4] - 1./3.*x[5] - 1./3.*x[6] - 2./3.);
    result[7] = ( eps_sf + 2./3.*x[1]*x[2] + 2./3.*x[1]*x[3] - 2./3.*x[1] + 1./3.*x[6]);
    result[8] = ( eps_sf + -x[2]);
    result[9] = ( eps_sf + -0.5*x[3]);

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
        grad[45] = 0.833333333333333 - 2./3.*x[1];
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
    Inequality constraints for g
*/
void g_igd_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( eps_sf + -x[0]*x[1] + x[0] + x[1] - 1.0);
    result[1] = ( eps_sf + x[0]*x[1] - x[0]);
    result[2] = ( eps_sf + -x[1]);
    result[3] = ( eps_sf + x[2] + x[3] + 2.0*x[4] - 1.0);
    result[4] = ( eps_sf + -x[3]);
    result[5] = ( eps_sf + -x[2]);
    result[6] = ( eps_sf + -x[4]);
    result[7] = ( eps_sf + -x[4]);

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
        grad[35] = 0.0;
        grad[36] = 0.0;
        grad[37] = 0.0;
        grad[38] = 0.0;
        grad[39] = -1.0;
    }

    return;
};

/**
    Inequality constraints for ol
*/
void ol_igd_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( eps_sf + x[0] - x[2] - 1.0);
    result[1] = ( eps_sf + -x[0] + x[2]);
    result[2] = ( eps_sf + -x[0]*x[1] + x[0] + x[1] + x[2] - 1.0);
    result[3] = ( eps_sf + x[0]*x[1] - x[0] - x[2]);
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
    Inequality constraints for opx
*/
void opx_igd_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( eps_sf + -x[0]*x[1] + x[0]*x[5] - x[0]*x[7] + x[0] - x[1]*x[3] + x[1] + x[3]*x[5] - x[3]*x[7] + x[3] - x[5] + x[7] - 1.0);
    result[1] = ( eps_sf + x[0]*x[1] - x[0]*x[5] + x[0]*x[7] - x[0] + x[1]*x[3] - x[3]*x[5] + x[3]*x[7] - x[3]);
    result[2] = ( eps_sf + -x[1] + x[4] + 2.0*x[5] + x[6] - x[7]);
    result[3] = ( eps_sf + -x[4]);
    result[4] = ( eps_sf + -x[6]);
    result[5] = ( eps_sf + -x[5]);
    result[6] = ( eps_sf + -x[0]*x[2] - x[0]*x[7] + x[0] + x[1]*x[3] + x[2] - x[3]*x[5] + x[3]*x[7] - x[3] + x[7] - 1.0);
    result[7] = ( eps_sf + x[0]*x[2] + x[0]*x[7] - x[0] - x[1]*x[3] + x[3]*x[5] - x[3]*x[7] + x[3]);
    result[8] = ( eps_sf + -x[2]);
    result[9] = ( eps_sf + -x[7]);
    result[10] = ( eps_sf + 0.5*x[1] - 1.0);
    result[11] = ( eps_sf + -0.5*x[1]);

    if (grad) {
        grad[0] = -x[1] + x[5] - x[7] + 1.0;
        grad[1] = -x[0] - x[3] + 1.0;
        grad[2] = 0.0;
        grad[3] = -x[1] + x[5] - x[7] + 1.0;
        grad[4] = 0.0;
        grad[5] = x[0] + x[3] - 1.0;
        grad[6] = 0.0;
        grad[7] = -x[0] - x[3] + 1.0;
        grad[8] = x[1] - x[5] + x[7] - 1.0;
        grad[9] = x[0] + x[3];
        grad[10] = 0.0;
        grad[11] = x[1] - x[5] + x[7] - 1.0;
        grad[12] = 0.0;
        grad[13] = -x[0] - x[3];
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
        grad[48] = -x[2] - x[7] + 1.0;
        grad[49] = x[3];
        grad[50] = 1.0 - x[0];
        grad[51] = x[1] - x[5] + x[7] - 1.0;
        grad[52] = 0.0;
        grad[53] = -x[3];
        grad[54] = 0.0;
        grad[55] = -x[0] + x[3] + 1.0;
        grad[56] = x[2] + x[7] - 1.0;
        grad[57] = -x[3];
        grad[58] = x[0];
        grad[59] = -x[1] + x[5] - x[7] + 1.0;
        grad[60] = 0.0;
        grad[61] = x[3];
        grad[62] = 0.0;
        grad[63] = x[0] - x[3];
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
void cpx_igd_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( eps_sf + -x[0]*x[1] - x[0]*x[3] + x[0]*x[7] - x[0]*x[8] + x[0] - x[1]*x[4] + x[1] - x[3]*x[4] + x[3] + x[4]*x[7] - x[4]*x[8] + x[4] - x[7] + x[8] - 1.0);
    result[1] = ( eps_sf + x[0]*x[1] + x[0]*x[3] - x[0]*x[7] + x[0]*x[8] - x[0] + x[1]*x[4] + x[3]*x[4] - x[4]*x[7] + x[4]*x[8] - x[4]);
    result[2] = ( eps_sf + -x[1] - x[3] + x[5] + x[6] + 2.0*x[7] - x[8]);
    result[3] = ( eps_sf + -x[5]);
    result[4] = ( eps_sf + -x[6]);
    result[5] = ( eps_sf + -x[7]);
    result[6] = ( eps_sf + x[0]*x[2] + x[1]*x[4] - x[2] + x[3]*x[4] - x[4]*x[7] + x[4]*x[8] - x[4]);
    result[7] = ( eps_sf + -x[0]*x[2] - x[1]*x[4] - x[3]*x[4] + x[4]*x[7] - x[4]*x[8] + x[4]);
    result[8] = ( eps_sf + x[2] + x[3] + x[8] - 1.0);
    result[9] = ( eps_sf + -x[3]);
    result[10] = ( eps_sf + -x[8]);
    result[11] = ( eps_sf + 0.5*x[1] - 1.0);
    result[12] = ( eps_sf + -0.5*x[1]);

    if (grad) {
        grad[0] = -x[1] - x[3] + x[7] - x[8] + 1.0;
        grad[1] = -x[0] - x[4] + 1.0;
        grad[2] = 0.0;
        grad[3] = -x[0] - x[4] + 1.0;
        grad[4] = -x[1] - x[3] + x[7] - x[8] + 1.0;
        grad[5] = 0.0;
        grad[6] = 0.0;
        grad[7] = x[0] + x[4] - 1.0;
        grad[8] = -x[0] - x[4] + 1.0;
        grad[9] = x[1] + x[3] - x[7] + x[8] - 1.0;
        grad[10] = x[0] + x[4];
        grad[11] = 0.0;
        grad[12] = x[0] + x[4];
        grad[13] = x[1] + x[3] - x[7] + x[8] - 1.0;
        grad[14] = 0.0;
        grad[15] = 0.0;
        grad[16] = -x[0] - x[4];
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
        grad[58] = x[1] + x[3] - x[7] + x[8] - 1.0;
        grad[59] = 0.0;
        grad[60] = 0.0;
        grad[61] = -x[4];
        grad[62] = x[4];
        grad[63] = -x[2];
        grad[64] = -x[4];
        grad[65] = -x[0];
        grad[66] = -x[4];
        grad[67] = -x[1] - x[3] + x[7] - x[8] + 1.0;
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
void ilm_igd_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( eps_sf + 0.5*x[0]*x[1] - 0.5*x[0] - 0.5*x[2]);
    result[1] = ( eps_sf + -0.5*x[0] + 0.5*x[3]);
    result[2] = ( eps_sf + x[0] - 1.0);
    result[3] = ( eps_sf + -0.5*x[0]*x[1] + 0.5*x[2] - 0.5*x[3]);
    result[4] = ( eps_sf + 0.5*x[0]*x[1] - 0.5*x[0] + 0.5*x[2]);
    result[5] = ( eps_sf + -0.5*x[0] - 0.5*x[3]);
    result[6] = ( eps_sf + x[0] - 1.0);
    result[7] = ( eps_sf + -0.5*x[0]*x[1] - 0.5*x[2] + 0.5*x[3]);

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
    Inequality constraints for hb
*/
void hb_igd_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( eps_sf + x[3] - 1.0);
    result[1] = ( eps_sf + x[3]*x[4] - x[3]);
    result[2] = ( eps_sf + -x[3]*x[4]);
    result[3] = ( eps_sf + x[0] - x[8] - 1.0);
    result[4] = ( eps_sf + -x[0] + x[8]);
    result[5] = ( eps_sf + -x[0]*x[1] - x[0]*x[6] - x[0]*x[7] + x[0] + x[1]*x[9] + x[1] + x[6]*x[9] + x[6] + x[7]*x[9] + x[7] - x[9] - 1.0);
    result[6] = ( eps_sf + x[0]*x[1] + x[0]*x[6] + x[0]*x[7] - x[0] - x[1]*x[9] - x[6]*x[9] - x[7]*x[9] + x[9]);
    result[7] = ( eps_sf + -x[1]);
    result[8] = ( eps_sf + -x[6]);
    result[9] = ( eps_sf + -x[7]);
    result[10] = ( eps_sf + -x[5]);
    result[11] = ( eps_sf + -x[0]*x[2] - x[0]*x[5] + x[0] - x[1]*x[9] + x[2] + x[5] - x[6]*x[9] - x[7]*x[9] + 1.5*x[8] + x[9] - 1.0);
    result[12] = ( eps_sf + x[0]*x[2] + x[0]*x[5] - x[0] + x[1]*x[9] + x[6]*x[9] + x[7]*x[9] - 1.5*x[8] - x[9]);
    result[13] = ( eps_sf + -x[2]);
    result[14] = ( eps_sf + 0.5*x[1] - 0.5*x[2] + 0.25*x[3] + 0.5*x[6] + 0.5*x[7] - 1.0);
    result[15] = ( eps_sf + -0.5*x[1] + 0.5*x[2] - 0.25*x[3] - 0.5*x[6] - 0.5*x[7]);
    result[16] = ( eps_sf + x[7] - 1.0);
    result[17] = ( eps_sf + -x[7]);

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
    Inequality constraints for bi
*/
void bi_igd_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( eps_sf + -x[0]*x[1] - x[0]*x[2] - x[0]*x[3] + x[0] + x[1] + x[2] + x[3] + 2./3.*x[4] - 1.0);
    result[1] = ( eps_sf + x[0]*x[1] + x[0]*x[2] + x[0]*x[3] - x[0] - 2./3.*x[4]);
    result[2] = ( eps_sf + -x[2]);
    result[3] = ( eps_sf + -x[3]);
    result[4] = ( eps_sf + -x[1]);
    result[5] = ( eps_sf + x[0] - 1./3.*x[4] - 1.0);
    result[6] = ( eps_sf + -x[0] + 1./3.*x[4]);
    result[7] = ( eps_sf + 0.5*x[1] + 0.5*x[2] - 0.5);
    result[8] = ( eps_sf + -0.5*x[1] - 0.5*x[2] - 0.5);
    result[9] = ( eps_sf + x[3] - 1.0);
    result[10] = ( eps_sf + -x[3]);

    if (grad) {
        grad[0] = -x[1] - x[2] - x[3] + 1.0;
        grad[1] = 1.0 - x[0];
        grad[2] = 1.0 - x[0];
        grad[3] = 1.0 - x[0];
        grad[4] = 2./3.;
        grad[5] = x[1] + x[2] + x[3] - 1.0;
        grad[6] = x[0];
        grad[7] = x[0];
        grad[8] = x[0];
        grad[9] = -2./3.;
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
        grad[29] = -1./3.;
        grad[30] = -1.0;
        grad[31] = 0.0;
        grad[32] = 0.0;
        grad[33] = 0.0;
        grad[34] = 1./3.;
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
        grad[50] = 0.0;
        grad[51] = 0.0;
        grad[52] = 0.0;
        grad[53] = -1.0;
        grad[54] = 0.0;
    }

    return;
};

/**
    Inequality constraints for ep
*/
void ep_igd_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
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
    Inequality constraints for cd
*/
void cd_igd_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
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


/**
    Inequality constraints for liq_mp
*/
void liq_mp_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( eps_sf + x[6] - 1.0);
    result[1] = ( eps_sf + -x[0]);
    result[2] = ( eps_sf + -x[1]*x[2]);
    result[3] = ( eps_sf + -x[1]*(1.0 - x[2]));
    result[4] = ( eps_sf + -x[3]);
    result[5] = ( eps_sf + x[3] + x[1] + x[6] + x[4] + x[0] - 1.0);
    result[6] = ( eps_sf + -x[4]);
    result[7] = ( eps_sf + -x[5]);
    result[8] = ( eps_sf + x[5] - 1.0);
    result[9] = ( eps_sf + -x[6]);

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
    result[0] = ( eps_sf + -x[1]*x[0] + x[1] + x[0] - 1.0);
    result[1] = ( eps_sf + x[1]*x[0] - x[0]);
    result[2] = ( eps_sf + -x[1]);
    result[3] = ( eps_sf + x[2] + 1.33333333333333*x[3] - 1.0);
    result[4] = ( eps_sf + -x[2]);
    result[5] = ( eps_sf + -x[3]);
    result[6] = ( eps_sf + -1./3.*x[3]);

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
        grad[15] = 1.33333333333333;
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
    result[0] = ( eps_sf + -x[1]);
    result[1] = ( eps_sf + x[1] + x[2] - 1.0);
    result[2] = ( eps_sf + -x[2]);
    result[3] = ( eps_sf + x[0] - 1.0);
    result[4] = ( eps_sf + -x[0]);

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
    result[0] = ( eps_sf + -x[2]*x[0] + x[2] - 0.75*x[3] - x[0]*x[1] + x[0] + x[1] - 1.0);
    result[1] = ( eps_sf + x[2]*x[0] + 0.75*x[3] + x[0]*x[1] - x[0]);
    result[2] = ( eps_sf + -x[2]);
    result[3] = ( eps_sf + -x[1]);
    result[4] = ( eps_sf + 0.25*x[3] + x[0] - 1.0);
    result[5] = ( eps_sf + -0.25*x[3] - x[0]);
    result[6] = ( eps_sf + x[2] + x[1] - 1.0);
    result[7] = ( eps_sf + -x[2] - x[1]);

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
    Inequality constraints for pl4tr_mp
*/
void pl4tr_mp_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
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
    result[0] = ( eps_sf + 0.5*x[4]*x[5] - x[3]*x[0] + x[3] + 0.5*x[1]*x[5] - x[1]*x[0] + x[1] - 0.5*x[5] - x[0]*x[2] + x[0] + x[2] - 1.0);
    result[1] = ( eps_sf + -0.5*x[4]*x[5] + x[3]*x[0] - 0.5*x[1]*x[5] + x[1]*x[0] + 0.5*x[5] + x[0]*x[2] - x[0]);
    result[2] = ( eps_sf + -x[1]);
    result[3] = ( eps_sf + -x[3]);
    result[4] = ( eps_sf + -x[2]);
    result[5] = ( eps_sf + -0.5*x[4]*x[5] - x[4]*x[0] + x[4] - 0.5*x[1]*x[5] - x[1]*x[0] + x[1] + 0.5*x[5] + x[0] - 1.0);
    result[6] = ( eps_sf + 0.5*x[4]*x[5] + x[4]*x[0] + 0.5*x[1]*x[5] + x[1]*x[0] - 0.5*x[5] - x[0]);
    result[7] = ( eps_sf + -x[1]);
    result[8] = ( eps_sf + -x[4]);
    result[9] = ( eps_sf + -0.5*x[3] - 0.5*x[2]);
    result[10] = ( eps_sf + 0.5*x[3] + 0.5*x[2] - 1.0);

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
    Inequality constraints for mt_mp
*/
void mt_mp_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( eps_sf + 0.5*x[0] - 0.5);
    result[1] = ( eps_sf + 0.5*x[1] - x[0]);
    result[2] = ( eps_sf + -0.5*x[1] + 0.5*x[0] - 0.5);
    result[3] = ( eps_sf + -x[1]);
    result[4] = ( eps_sf + x[1] - 1.0);

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
    Inequality constraints for ilm_mp
*/
void ilm_mp_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( eps_sf + 0.5*x[2] - 0.5*x[3] - 0.5*x[0] + 0.5*x[1]);
    result[1] = ( eps_sf + 0.5*x[2] + 0.5*x[3] - 0.5*x[0] + 0.5*x[1]);
    result[2] = ( eps_sf + -x[1]);
    result[3] = ( eps_sf + -x[2]);
    result[4] = ( eps_sf + x[0] - 1.0);
    result[5] = ( eps_sf + -0.5*x[2] - 0.5*x[3] - 0.5*x[0] - 0.5*x[1]);

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
        grad[21] = -0.50;
        grad[22] = -0.50;
        grad[23] = -0.50;
    }

    return;
};

/**
    Inequality constraints for g_mp
*/
void g_mp_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( eps_sf - x[2]*x[0] + x[2] - x[0]*x[1] + x[0] + x[1] - 1.0);
    result[1] = ( eps_sf + x[2]*x[0] + x[0]*x[1] - x[0]);
    result[2] = ( eps_sf - x[2]);
    result[3] = ( eps_sf - x[1]);
    result[4] = ( eps_sf - 1.0 + x[3]);
    result[5] = ( eps_sf - x[3]);

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
    Inequality constraints for ctd_mp
*/
void ctd_mp_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( eps_sf + x[2] - 1.0);
    result[1] = ( eps_sf + -x[2]);
    result[2] = ( eps_sf + x[1]*x[0] - x[0]);
    result[3] = ( eps_sf + -x[1]*x[0] + x[1] + x[0] - 1.0);
    result[4] = ( eps_sf + -x[1]);

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
    result[0] = ( eps_sf + x[3]*x[5] - x[3]*x[0] + x[3] - x[5]*x[4] + x[5]*x[1] - x[5] + x[4]*x[0] - x[4] - x[0]*x[1] + x[0] + x[1] - 1.0);
    result[1] = ( eps_sf + -x[3]*x[5] + x[3]*x[0] + x[5]*x[4] - x[5]*x[1] + x[5] - x[4]*x[0] + x[0]*x[1] - x[0]);
    result[2] = ( eps_sf + x[4] - x[1]);
    result[3] = ( eps_sf + -0.25*x[2]*x[6] - 0.25*x[3]*x[5] - x[3]*x[0] + x[3] + 0.25*x[5]*x[4] - 0.25*x[5]*x[1] + 0.25*x[5] - 0.25*x[6]*x[4] - 0.25*x[6]*x[1] + 0.25*x[6] + x[0] - 1.0);
    result[4] = ( eps_sf + -x[3]);
    result[5] = ( eps_sf + 0.25*x[2]*x[6] + 0.25*x[3]*x[5] + x[3]*x[0] - 0.25*x[5]*x[4] + 0.25*x[5]*x[1] - 0.25*x[5] + 0.25*x[6]*x[4] + 0.25*x[6]*x[1] - 0.25*x[6] - x[0]);
    result[6] = ( eps_sf + x[2]*x[6] - x[2]*x[0] + x[2] + x[6]*x[4] + x[6]*x[1] - x[6] - x[4]*x[0] + x[4] - x[0]*x[1] + x[0] + x[1] - 1.0);
    result[7] = ( eps_sf + -x[2]*x[6] + x[2]*x[0] - x[6]*x[4] - x[6]*x[1] + x[6] + x[4]*x[0] + x[0]*x[1] - x[0]);
    result[8] = ( eps_sf + -x[2]);
    result[9] = ( eps_sf + -x[4] - x[1]);
    result[10] = ( eps_sf + 0.5*x[2] + x[1] - 1.0);
    result[11] = ( eps_sf + -0.5*x[2] - x[1]);

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
    result[0] = ( eps_sf + -x[0]);
    result[1] = ( eps_sf + x[0] - 1.0);
    result[2] = ( eps_sf + -x[1]);
    result[3] = ( eps_sf + -x[2]);
    result[4] = ( eps_sf + x[2] - 1.0);

    if (grad) {
        grad[0] = -1.0;
        grad[1] = 0.0;
        grad[2] = 0.0;
        grad[3] = 1.0;
        grad[4] = 0.0;
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
    result[0] = ( eps_sf  - x[3]*x[0] + x[3] - 3.0*x[1]*x[0] + x[1] + 2.0/3.0*x[5] - x[4]*x[0] + x[4] - x[0]*x[2] + x[0] + x[2] - 1.0);
    result[1] = ( eps_sf  - x[1]);
    result[2] = ( eps_sf  + x[3]*x[0] + 3.0*x[1]*x[0] - 2.0/3.0*x[5] + x[4]*x[0] + x[0]*x[2] - x[0]);
    result[3] = ( eps_sf  - x[3]);
    result[4] = ( eps_sf  - x[4]);
    result[5] = ( eps_sf  - x[2]);
    result[6] = ( eps_sf  + x[1] - 1.0/3.0*x[5] + x[0] - 1.0);
    result[7] = ( eps_sf  - x[1]);
    result[8] = ( eps_sf  + 1.0/3.0*x[5] - x[0]);
    result[9] = ( eps_sf  + 0.5*x[3] + 0.5*x[2] - 0.5);
    result[10] = ( eps_sf - 0.5*x[3] - 0.5*x[2] - 0.5);
    result[11] = ( eps_sf + x[4] - 1.0);
    result[12] = ( eps_sf - x[4]);

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
  local minimization for clinopyroxene
*/
void cpx_ig_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
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
void ep_ig_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
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
void fl_ig_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( eps_sf + x[0] + x[1] + x[2] + x[3] + x[4] + x[5] + x[6] + x[7] + x[8] + x[9] - 1.0);
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
void hb_ig_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
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
    Inequality constraints for ilm
*/
void ilm_ig_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( eps_sf + 0.5*x[0]*x[1] - 0.5*x[0] - 0.5*x[2]);
    result[1] = ( eps_sf + -0.5*x[0] + 0.5*x[3]);
    result[2] = ( eps_sf + x[0] - 1.0);
    result[3] = ( eps_sf + -0.5*x[0]*x[1] + 0.5*x[2] - 0.5*x[3]);
    result[4] = ( eps_sf + 0.5*x[0]*x[1] - 0.5*x[0] + 0.5*x[2]);
    result[5] = ( eps_sf + -0.5*x[0] - 0.5*x[3]);
    result[6] = ( eps_sf + x[0] - 1.0);
    result[7] = ( eps_sf + -0.5*x[0]*x[1] - 0.5*x[2] + 0.5*x[3]);

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
    Inequality constraints for liqHw
*/
void liq_ig_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( eps_sf +x[6] +x[3] +x[2] +x[10] +x[5] +x[4] +x[8] +x[1] +x[7] +x[0] - 0.25*x[9]*(-3.0*x[6] - 3.0*x[3] - 3.0*x[2] - 3.0*x[10] - 3.0*x[5] - 3.0*x[4] - 3.0*x[8] - 3.0*x[1] - 3.0*x[7] - 3.0*x[0] + 4.0) - 1.0);
    result[1] = ( eps_sf + -0.75*x[1]*x[9] - x[1] +x[9]);
    result[2] = ( eps_sf + -0.75*x[0]*x[9] - x[0] +x[9]);
    result[3] = ( eps_sf + -0.75*x[4]*x[9] - x[4]);
    result[4] = ( eps_sf + -0.75*x[5]*x[9] - x[5]);
    result[5] = ( eps_sf + -0.75*x[6]*x[9] - x[6]);
    result[6] = ( eps_sf + -0.75*x[7]*x[9] - x[7]);
    result[7] = ( eps_sf + -0.75*x[8]*x[9] - x[8]);
    result[8] = ( eps_sf + -x[9]);
    result[9] = ( eps_sf + -x[3] - x[2] - 0.75*x[9]*(x[3] + x[2]));
    result[10] = ( eps_sf + 0.75*x[10]*x[9] +x[10] - 1.0);
    result[11] = ( eps_sf + -4.0*x[2]*(0.75*x[9] + 1.0));
    result[12] = ( eps_sf + -4.0*x[3]*(0.75*x[9] + 1.0));
    result[13] = ( eps_sf + -x[0]*(0.75*x[9] + 1.0) +x[9]);
    result[14] = ( eps_sf + -x[1]*(0.75*x[9] + 1.0) +x[9]);
    result[15] = ( eps_sf + 2.0*x[9] - (0.75*x[9] + 1.0)*(4.0*x[3] + 4.0*x[2] + x[1] + x[0]));
    result[16] = ( eps_sf + -x[10]*(0.75*x[9] + 1.0));
    result[17] = ( eps_sf + 0.75*x[10]*x[9] +x[10] - 1.0);

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
void ol_ig_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
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
void opx_ig_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
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
void pl4T_ig_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
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
void spn_ig_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *SS_ref_db){
    result[0] = ( eps_sf + 1./3.*x[0]*x[3] + 1./3.*x[0] - 1./3.*x[3] - 2./3.*x[4] - 1./3.);
    result[1] = ( eps_sf + -1./3.*x[0]*x[3] - 1./3.*x[0] - 2./3.*x[5]);
    result[2] = ( eps_sf + -2./3.*x[1]*x[2] - 2./3.*x[1]*x[3] + 2./3.*x[1] + 1./3.*x[3] + 2./3.*x[4] + 2./3.*x[5] + 2./3.*x[6] - 2./3.);
    result[3] = ( eps_sf + 2./3.*x[1]*x[2] + 2./3.*x[1]*x[3] - 2./3.*x[1] - 2./3.*x[6]);
    result[4] = ( eps_sf + 1./3.*x[0]*x[3] + 1./3.*x[0] - 1./3.*x[3] + 1./3.*x[4] - 1./3.);
    result[5] = ( eps_sf + -1./3.*x[0]*x[3] - 1./3.*x[0] + 1./3.*x[5]);
    result[6] = ( eps_sf + -2./3.*x[1]*x[2] - 2./3.*x[1]*x[3] + 2./3.*x[1] + x[2] + 0.833333333333333*x[3] - 1./3.*x[4] - 1./3.*x[5] - 1./3.*x[6] - 2./3.);
    result[7] = ( eps_sf + 2./3.*x[1]*x[2] + 2./3.*x[1]*x[3] - 2./3.*x[1] + 1./3.*x[6]);
    result[8] = ( eps_sf + -x[2]);
    result[9] = ( eps_sf + -0.5*x[3]);

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
        grad[45] = 0.833333333333333 - 2./3.*x[1];
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
    result[0] = ( eps_sf + -x[2]*x[0] + x[2] + 2.0/3.0*x[4] - x[3]*x[0] + x[3] - x[0]*x[1] + x[0] + x[1] - 1.0);
    result[1] = ( eps_sf + x[2]*x[0] - 2.0/3.0*x[4] + x[3]*x[0] + x[0]*x[1] - x[0]);
    result[2] = ( eps_sf + -x[2]);
    result[3] = ( eps_sf + -x[3]);
    result[4] = ( eps_sf + -x[1]);
    result[5] = ( eps_sf + -1.0/3.0*x[4] + x[0] - 1.0);
    result[6] = ( eps_sf + 1.0/3.0*x[4] - x[0]);
    result[7] = ( eps_sf + 0.5*x[2] + 0.5*x[1] - 0.5);
    result[8] = ( eps_sf + -0.5*x[2] - 0.5*x[1] - 0.5);
    result[9] = ( eps_sf + x[3] - 1.0);
    result[10] = ( eps_sf + -x[3]);

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



//---------------------------------------------------------------------------
//---------------------------------Evans&Frost,2021--------------------------
//---------------------------------------------------------------------------


/**
    Inequality constraints for fluid
*/
void fluid_um_c(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    result[0] = ( eps_sf + -x[0]);
    result[1] = ( eps_sf + x[0] - 1.0);

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
    result[0] = ( eps_sf + x[0] - 1.0);
    result[1] = ( eps_sf + -x[0]);

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
    result[0] = ( eps_sf + x[0] - 1.0);
    result[1] = ( eps_sf + -x[0]);

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
    result[0] = ( eps_sf + x[0] - 1.0);
    result[1] = ( eps_sf + -x[0]);

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
    result[0] = ( eps_sf + -x[0]*x[1] -x[0]*x[2] + x[0] + x[1]*x[3] + x[1] + x[2]*x[3] + x[2] -x[3] - 1.0);
    result[1] = ( eps_sf + x[0]*x[1] + x[0]*x[2] -x[0] -x[1]*x[3] -x[2]*x[3] + x[3]);
    result[2] = ( eps_sf + -x[2]);
    result[3] = ( eps_sf + -x[1]);
    result[4] = ( eps_sf + x[0] - 0.5*x[1]*x[3] - 0.5*x[2]*x[3] + 0.5*x[3] - 1.0);
    result[5] = ( eps_sf + -x[0] + 0.5*x[1]*x[3] + 0.5*x[2]*x[3] - 0.5*x[3]);
    result[6] = ( eps_sf + 0.5*x[1] + 0.5*x[2] - 1.0);
    result[7] = ( eps_sf + -0.5*x[1] - 0.5*x[2]);

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
    result[0] = ( eps_sf + x[0] - 1.0);
    result[1] = ( eps_sf + -x[0]);

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
    result[0] = ( eps_sf + -x[0]*x[3] + x[0] + x[3]*x[4] + x[3] -x[4] - 1.0);
    result[1] = ( eps_sf + x[0]*x[3] -x[0] -x[3]*x[4] + x[4]);
    result[2] = ( eps_sf + -x[3]);
    result[3] = ( eps_sf + -x[0]*x[1] -x[0]*x[2] + x[0] + x[1] + x[2] - 0.5*x[3]*x[4] + 0.5*x[4] - 1.0);
    result[4] = ( eps_sf + x[0]*x[1] + x[0]*x[2] -x[0] + 0.5*x[3]*x[4] - 0.5*x[4]);
    result[5] = ( eps_sf + -x[2]);
    result[6] = ( eps_sf + -x[1]);
    result[7] = ( eps_sf + x[1] + x[2] -x[3] - 1.0);
    result[8] = ( eps_sf + -x[1] -x[2] + x[3]);

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
    result[0] = ( eps_sf + -x[0]*x[1] + x[0]*x[3] + x[0] + x[1]*x[4] + x[1] -x[3]*x[4] -x[3] -x[4] - 1.0);
    result[1] = ( eps_sf + x[0]*x[1] -x[0]*x[3] -x[0] -x[1]*x[4] + x[3]*x[4] + x[4]);
    result[2] = ( eps_sf + -x[1] + x[3]);
    result[3] = ( eps_sf + x[0] - 0.25*x[1]*x[4] - 0.25*x[1]*x[5] - 0.25*x[2]*x[5] + 0.25*x[3]*x[4] - 0.25*x[3]*x[5] + 0.25*x[4] + 0.25*x[5] - 1.0);
    result[4] = ( eps_sf + -x[0] + 0.25*x[1]*x[4] + 0.25*x[1]*x[5] + 0.25*x[2]*x[5] - 0.25*x[3]*x[4] + 0.25*x[3]*x[5] - 0.25*x[4] - 0.25*x[5]);
    result[5] = ( eps_sf + -x[0]*x[1] -x[0]*x[2] -x[0]*x[3] + x[0] + x[1]*x[5] + x[1] + x[2]*x[5] + x[2] + x[3]*x[5] + x[3] -x[5] - 1.0);
    result[6] = ( eps_sf + x[0]*x[1] + x[0]*x[2] + x[0]*x[3] -x[0] -x[1]*x[5] -x[2]*x[5] -x[3]*x[5] + x[5]);
    result[7] = ( eps_sf + -x[2]);
    result[8] = ( eps_sf + -x[1] -x[3]);
    result[9] = ( eps_sf + x[1] + 0.5*x[2] - 1.0);
    result[10] = ( eps_sf + -x[1] - 0.5*x[2]);

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
    result[0] = ( eps_sf + x[0] - x[1]*x[3] + 1.5*x[2] + x[3] - 1.0);
    result[1] = ( eps_sf + -x[0] + x[1]*x[3] - 1.5*x[2] - x[3]);
    result[2] = ( eps_sf + x[0] - x[2] - 1.0);
    result[3] = ( eps_sf + -x[0] + x[2]);
    result[4] = ( eps_sf + -x[1]);
    result[5] = ( eps_sf + -x[0]*x[1] + x[0] + x[1]*x[3] + x[1] - x[3] - 1.0);
    result[6] = ( eps_sf + x[0]*x[1] - x[0] - x[1]*x[3] + x[3]);
    result[7] = ( eps_sf + -0.5*x[1]);
    result[8] = ( eps_sf + 0.5*x[1] - 1.0);

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
    result[0] = ( eps_sf + -x[1]);
    result[1] = ( eps_sf + x[1] - 1.0);
    result[2] = ( eps_sf + x[0] - 1.0);
    result[3] = ( eps_sf + -x[0]);

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
    result[0] = ( eps_sf + -x[0]*x[1] -x[0]*x[2] + x[0] + x[1] + x[2] - 0.5*x[3] - 1.0);
    result[1] = ( eps_sf + x[0]*x[1] + x[0]*x[2] -x[0] + 0.5*x[3]);
    result[2] = ( eps_sf + -x[2]);
    result[3] = ( eps_sf + -x[1]);
    result[4] = ( eps_sf + x[0] + 0.5*x[3] - 1.0);
    result[5] = ( eps_sf + -x[0] - 0.5*x[3]);
    result[6] = ( eps_sf + -0.5*x[1] - 0.5*x[2]);
    result[7] = ( eps_sf + 0.5*x[1] + 0.5*x[2] - 1.0);

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
    result[0] = ( eps_sf + x[0] - 1.0);
    result[1] = ( eps_sf + -x[0]);

    if (grad) {
        grad[0] = 1.0;
        grad[1] = -1.0;
    }

    return;
};

typedef struct global_min_datas {
	global_variable 	 gv; 
	bulk_info 	       z_b;
	obj_type 			    *SS_objective;
	sf_type 			    *SS_sf;
	PP_ref 				    *PP_ref_db;
	SS_ref 				    *SS_ref_db;
	csd_phase_set  		*cp;
	
} global_min_data;


SS_ref NLopt_opt_alk_liq_function(global_variable gv, SS_ref SS_ref_db){
    
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
    nlopt_set_min_objective(SS_ref_db.opt, obj_alk_liq, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, liq_alk_c, NULL, SS_ref_db.tol_sf);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    if (gv.maxeval==1){  
       // we are only interested in evaluating the objective function  
       minf = obj_alk_liq(n, x, NULL, &SS_ref_db);
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
SS_ref NLopt_opt_alk_fl_function(global_variable gv, SS_ref SS_ref_db){
    
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
    nlopt_set_min_objective(SS_ref_db.opt, obj_alk_fl, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, fl_alk_c, NULL, SS_ref_db.tol_sf);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    if (gv.maxeval==1){  
       // we are only interested in evaluating the objective function  
       minf = obj_alk_fl(n, x, NULL, &SS_ref_db);
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
SS_ref NLopt_opt_alk_fsp_function(global_variable gv, SS_ref SS_ref_db){
    
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
    nlopt_set_min_objective(SS_ref_db.opt, obj_alk_fsp, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, fsp_alk_c, NULL, SS_ref_db.tol_sf);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    if (gv.maxeval==1){  
       // we are only interested in evaluating the objective function  
       minf = obj_alk_fsp(n, x, NULL, &SS_ref_db);
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
SS_ref NLopt_opt_alk_spn_function(global_variable gv, SS_ref SS_ref_db){
    
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
    nlopt_set_min_objective(SS_ref_db.opt, obj_alk_spn, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, spn_alk_c, NULL, SS_ref_db.tol_sf);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    if (gv.maxeval==1){  
       // we are only interested in evaluating the objective function  
       minf = obj_alk_spn(n, x, NULL, &SS_ref_db);
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
SS_ref NLopt_opt_alk_g_function(global_variable gv, SS_ref SS_ref_db){
    
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
    nlopt_set_min_objective(SS_ref_db.opt, obj_alk_g, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, g_alk_c, NULL, SS_ref_db.tol_sf);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    if (gv.maxeval==1){  
       // we are only interested in evaluating the objective function  
       minf = obj_alk_g(n, x, NULL, &SS_ref_db);
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
SS_ref NLopt_opt_alk_ol_function(global_variable gv, SS_ref SS_ref_db){
    
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
    nlopt_set_min_objective(SS_ref_db.opt, obj_alk_ol, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, ol_alk_c, NULL, SS_ref_db.tol_sf);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    if (gv.maxeval==1){  
       // we are only interested in evaluating the objective function  
       minf = obj_alk_ol(n, x, NULL, &SS_ref_db);
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
SS_ref NLopt_opt_alk_opx_function(global_variable gv, SS_ref SS_ref_db){
    
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
    nlopt_set_min_objective(SS_ref_db.opt, obj_alk_opx, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, opx_alk_c, NULL, SS_ref_db.tol_sf);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    if (gv.maxeval==1){  
       // we are only interested in evaluating the objective function  
       minf = obj_alk_opx(n, x, NULL, &SS_ref_db);
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
SS_ref NLopt_opt_alk_cpx_function(global_variable gv, SS_ref SS_ref_db){
    
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
    nlopt_set_min_objective(SS_ref_db.opt, obj_alk_cpx, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, cpx_alk_c, NULL, SS_ref_db.tol_sf);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    if (gv.maxeval==1){  
       // we are only interested in evaluating the objective function  
       minf = obj_alk_cpx(n, x, NULL, &SS_ref_db);
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
SS_ref NLopt_opt_alk_ilm_function(global_variable gv, SS_ref SS_ref_db){
    
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
    nlopt_set_min_objective(SS_ref_db.opt, obj_alk_ilm, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, ilm_alk_c, NULL, SS_ref_db.tol_sf);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    if (gv.maxeval==1){  
       // we are only interested in evaluating the objective function  
       minf = obj_alk_ilm(n, x, NULL, &SS_ref_db);
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
SS_ref NLopt_opt_alk_ness_function(global_variable gv, SS_ref SS_ref_db){
    
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
    nlopt_set_min_objective(SS_ref_db.opt, obj_alk_ness, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, ness_alk_c, NULL, SS_ref_db.tol_sf);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    if (gv.maxeval==1){  
       // we are only interested in evaluating the objective function  
       minf = obj_alk_ness(n, x, NULL, &SS_ref_db);
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
SS_ref NLopt_opt_alk_lct_function(global_variable gv, SS_ref SS_ref_db){
    
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
    nlopt_set_min_objective(SS_ref_db.opt, obj_alk_lct, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, lct_alk_c, NULL, SS_ref_db.tol_sf);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    if (gv.maxeval==1){  
       // we are only interested in evaluating the objective function  
       minf = obj_alk_lct(n, x, NULL, &SS_ref_db);
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
SS_ref NLopt_opt_alk_kals_function(global_variable gv, SS_ref SS_ref_db){
    
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
    nlopt_set_min_objective(SS_ref_db.opt, obj_alk_kals, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, kals_alk_c, NULL, SS_ref_db.tol_sf);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    if (gv.maxeval==1){  
       // we are only interested in evaluating the objective function  
       minf = obj_alk_kals(n, x, NULL, &SS_ref_db);
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
SS_ref NLopt_opt_alk_mel_function(global_variable gv, SS_ref SS_ref_db){
    
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
    nlopt_set_min_objective(SS_ref_db.opt, obj_alk_mel, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, mel_alk_c, NULL, SS_ref_db.tol_sf);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    if (gv.maxeval==1){  
       // we are only interested in evaluating the objective function  
       minf = obj_alk_mel(n, x, NULL, &SS_ref_db);
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
SS_ref NLopt_opt_alk_hb_function(global_variable gv, SS_ref SS_ref_db){
    
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
    nlopt_set_min_objective(SS_ref_db.opt, obj_alk_hb, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, hb_alk_c, NULL, SS_ref_db.tol_sf);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    if (gv.maxeval==1){  
       // we are only interested in evaluating the objective function  
       minf = obj_alk_hb(n, x, NULL, &SS_ref_db);
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
SS_ref NLopt_opt_alk_bi_function(global_variable gv, SS_ref SS_ref_db){
    
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
    nlopt_set_min_objective(SS_ref_db.opt, obj_alk_bi, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, bi_alk_c, NULL, SS_ref_db.tol_sf);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    if (gv.maxeval==1){  
       // we are only interested in evaluating the objective function  
       minf = obj_alk_bi(n, x, NULL, &SS_ref_db);
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
SS_ref NLopt_opt_alk_ep_function(global_variable gv, SS_ref SS_ref_db){
    
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
    nlopt_set_min_objective(SS_ref_db.opt, obj_alk_ep, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, ep_alk_c, NULL, SS_ref_db.tol_sf);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    if (gv.maxeval==1){  
       // we are only interested in evaluating the objective function  
       minf = obj_alk_ep(n, x, NULL, &SS_ref_db);
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
SS_ref NLopt_opt_alk_cd_function(global_variable gv, SS_ref SS_ref_db){
    
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
    nlopt_set_min_objective(SS_ref_db.opt, obj_alk_cd, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, cd_alk_c, NULL, SS_ref_db.tol_sf);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    if (gv.maxeval==1){  
       // we are only interested in evaluating the objective function  
       minf = obj_alk_cd(n, x, NULL, &SS_ref_db);
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


SS_ref NLopt_opt_igd_liq_function(global_variable gv, SS_ref SS_ref_db){
    
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
    nlopt_set_min_objective(SS_ref_db.opt, obj_igd_liq, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, liq_igd_c, NULL, SS_ref_db.tol_sf);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    if (gv.maxeval==1){  
       // we are only interested in evaluating the objective function  
       minf = obj_igd_liq(n, x, NULL, &SS_ref_db);
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
SS_ref NLopt_opt_igd_fl_function(global_variable gv, SS_ref SS_ref_db){
    
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
    nlopt_set_min_objective(SS_ref_db.opt, obj_igd_fl, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, fl_igd_c, NULL, SS_ref_db.tol_sf);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    if (gv.maxeval==1){  
       // we are only interested in evaluating the objective function  
       minf = obj_igd_fl(n, x, NULL, &SS_ref_db);
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
SS_ref NLopt_opt_igd_fsp_function(global_variable gv, SS_ref SS_ref_db){
    
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
    nlopt_set_min_objective(SS_ref_db.opt, obj_igd_fsp, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, fsp_igd_c, NULL, SS_ref_db.tol_sf);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    if (gv.maxeval==1){  
       // we are only interested in evaluating the objective function  
       minf = obj_igd_fsp(n, x, NULL, &SS_ref_db);
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
SS_ref NLopt_opt_igd_spn_function(global_variable gv, SS_ref SS_ref_db){
    
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
    nlopt_set_min_objective(SS_ref_db.opt, obj_igd_spn, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, spn_igd_c, NULL, SS_ref_db.tol_sf);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    if (gv.maxeval==1){  
       // we are only interested in evaluating the objective function  
       minf = obj_igd_spn(n, x, NULL, &SS_ref_db);
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
SS_ref NLopt_opt_igd_g_function(global_variable gv, SS_ref SS_ref_db){
    
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
    nlopt_set_min_objective(SS_ref_db.opt, obj_igd_g, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, g_igd_c, NULL, SS_ref_db.tol_sf);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    if (gv.maxeval==1){  
       // we are only interested in evaluating the objective function  
       minf = obj_igd_g(n, x, NULL, &SS_ref_db);
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
SS_ref NLopt_opt_igd_ol_function(global_variable gv, SS_ref SS_ref_db){
    
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
    nlopt_set_min_objective(SS_ref_db.opt, obj_igd_ol, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, ol_igd_c, NULL, SS_ref_db.tol_sf);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    if (gv.maxeval==1){  
       // we are only interested in evaluating the objective function  
       minf = obj_igd_ol(n, x, NULL, &SS_ref_db);
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
SS_ref NLopt_opt_igd_opx_function(global_variable gv, SS_ref SS_ref_db){
    
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
    nlopt_set_min_objective(SS_ref_db.opt, obj_igd_opx, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, opx_igd_c, NULL, SS_ref_db.tol_sf);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    if (gv.maxeval==1){  
       // we are only interested in evaluating the objective function  
       minf = obj_igd_opx(n, x, NULL, &SS_ref_db);
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
SS_ref NLopt_opt_igd_cpx_function(global_variable gv, SS_ref SS_ref_db){
    
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
    nlopt_set_min_objective(SS_ref_db.opt, obj_igd_cpx, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, cpx_igd_c, NULL, SS_ref_db.tol_sf);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    if (gv.maxeval==1){  
       // we are only interested in evaluating the objective function  
       minf = obj_igd_cpx(n, x, NULL, &SS_ref_db);
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
SS_ref NLopt_opt_igd_ilm_function(global_variable gv, SS_ref SS_ref_db){
    
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
    nlopt_set_min_objective(SS_ref_db.opt, obj_igd_ilm, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, ilm_igd_c, NULL, SS_ref_db.tol_sf);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    if (gv.maxeval==1){  
       // we are only interested in evaluating the objective function  
       minf = obj_igd_ilm(n, x, NULL, &SS_ref_db);
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
SS_ref NLopt_opt_igd_hb_function(global_variable gv, SS_ref SS_ref_db){
    
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
    nlopt_set_min_objective(SS_ref_db.opt, obj_igd_hb, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, hb_igd_c, NULL, SS_ref_db.tol_sf);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    if (gv.maxeval==1){  
       // we are only interested in evaluating the objective function  
       minf = obj_igd_hb(n, x, NULL, &SS_ref_db);
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
SS_ref NLopt_opt_igd_bi_function(global_variable gv, SS_ref SS_ref_db){
    
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
    nlopt_set_min_objective(SS_ref_db.opt, obj_igd_bi, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, bi_igd_c, NULL, SS_ref_db.tol_sf);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    if (gv.maxeval==1){  
       // we are only interested in evaluating the objective function  
       minf = obj_igd_bi(n, x, NULL, &SS_ref_db);
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
SS_ref NLopt_opt_igd_ep_function(global_variable gv, SS_ref SS_ref_db){
    
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
    nlopt_set_min_objective(SS_ref_db.opt, obj_igd_ep, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, ep_igd_c, NULL, SS_ref_db.tol_sf);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    if (gv.maxeval==1){  
       // we are only interested in evaluating the objective function  
       minf = obj_igd_ep(n, x, NULL, &SS_ref_db);
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
SS_ref NLopt_opt_igd_cd_function(global_variable gv, SS_ref SS_ref_db){
    
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
    nlopt_set_min_objective(SS_ref_db.opt, obj_igd_cd, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, cd_igd_c, NULL, SS_ref_db.tol_sf);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    if (gv.maxeval==1){  
       // we are only interested in evaluating the objective function  
       minf = obj_igd_cd(n, x, NULL, &SS_ref_db);
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
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, st_mp_c, NULL, SS_ref_db.tol_sf);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    if (gv.maxeval==1){  
       // we are only interested in evaluating the objective function  
       minf = obj_mp_st(n, x, NULL, &SS_ref_db);
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
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, sp_mp_c, NULL, SS_ref_db.tol_sf);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    if (gv.maxeval==1){  
       // we are only interested in evaluating the objective function  
       minf = obj_mp_sp(n, x, NULL, &SS_ref_db);
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
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, sa_mp_c, NULL, SS_ref_db.tol_sf);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    if (gv.maxeval==1){  
       // we are only interested in evaluating the objective function  
       minf = obj_mp_sa(n, x, NULL, &SS_ref_db);
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

SS_ref NLopt_opt_mp_pl4tr_function(global_variable gv, SS_ref SS_ref_db){
    
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
    nlopt_set_min_objective(SS_ref_db.opt, obj_mp_pl4tr, &SS_ref_db);
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, pl4tr_mp_c, NULL, SS_ref_db.tol_sf);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    if (gv.maxeval==1){  
       // we are only interested in evaluating the objective function  
       minf = obj_mp_pl4tr(n, x, NULL, &SS_ref_db);
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
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, opx_mp_c, NULL, SS_ref_db.tol_sf);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    if (gv.maxeval==1){  
       // we are only interested in evaluating the objective function  
       minf = obj_mp_opx(n, x, NULL, &SS_ref_db);
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
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, mu_mp_c, NULL, SS_ref_db.tol_sf);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    if (gv.maxeval==1){  
       // we are only interested in evaluating the objective function  
       minf = obj_mp_mu(n, x, NULL, &SS_ref_db);
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
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, mt_mp_c, NULL, SS_ref_db.tol_sf);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    if (gv.maxeval==1){  
       // we are only interested in evaluating the objective function  
       minf = obj_mp_mt(n, x, NULL, &SS_ref_db);
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
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, ma_mp_c, NULL, SS_ref_db.tol_sf);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    if (gv.maxeval==1){  
       // we are only interested in evaluating the objective function  
       minf = obj_mp_ma(n, x, NULL, &SS_ref_db);
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
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, ilm_mp_c, NULL, SS_ref_db.tol_sf);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    if (gv.maxeval==1){  
       // we are only interested in evaluating the objective function  
       minf = obj_mp_ilm(n, x, NULL, &SS_ref_db);
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
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, g_mp_c, NULL, SS_ref_db.tol_sf);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    if (gv.maxeval==1){  
       // we are only interested in evaluating the objective function  
       minf = obj_mp_g(n, x, NULL, &SS_ref_db);
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
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, ep_mp_c, NULL, SS_ref_db.tol_sf);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    if (gv.maxeval==1){  
       // we are only interested in evaluating the objective function  
       minf = obj_mp_ep(n, x, NULL, &SS_ref_db);
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
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, ctd_mp_c, NULL, SS_ref_db.tol_sf);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    if (gv.maxeval==1){  
       // we are only interested in evaluating the objective function  
       minf = obj_mp_ctd(n, x, NULL, &SS_ref_db);
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
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, chl_mp_c, NULL, SS_ref_db.tol_sf);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    if (gv.maxeval==1){  
       // we are only interested in evaluating the objective function  
       minf = obj_mp_chl(n, x, NULL, &SS_ref_db);
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
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, cd_mp_c, NULL, SS_ref_db.tol_sf);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    if (gv.maxeval==1){  
       // we are only interested in evaluating the objective function  
       minf = obj_mp_cd(n, x, NULL, &SS_ref_db);
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
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, bi_mp_c, NULL, SS_ref_db.tol_sf);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    if (gv.maxeval==1){  
       // we are only interested in evaluating the objective function  
       minf = obj_mp_bi(n, x, NULL, &SS_ref_db);
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
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, liq_mp_c, NULL, SS_ref_db.tol_sf);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    if (gv.maxeval==1){  
       // we are only interested in evaluating the objective function  
       minf = obj_mp_liq(n, x, NULL, &SS_ref_db);
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
   nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, bi_ig_c, NULL, SS_ref_db.tol_sf);
   nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
   nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
   
   double minf;
   if (gv.maxeval==1){  
     // we are only interested in evaluating the objective function  
     minf = obj_ig_bi(n, x, NULL, &SS_ref_db); 
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
  nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, cd_ig_c, NULL, SS_ref_db.tol_sf);
	nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
  nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
   
	double minf;
	if (gv.maxeval==1){  
          // we are only interested in evaluating the objective function  
          minf = obj_ig_cd(n, x, NULL, &SS_ref_db); 
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
   nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, cpx_ig_c, NULL, SS_ref_db.tol_sf);
   nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
   nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
   
   double minf;
   if (gv.maxeval==1){  
     // we are only interested in evaluating the objective function  
     minf = obj_ig_cpx(n, x, NULL, &SS_ref_db); 
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
   nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, ep_ig_c, NULL, SS_ref_db.tol_sf);
   nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
   nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
   
   double minf;
   if (gv.maxeval==1){  
     // we are only interested in evaluating the objective function  
     minf = obj_ig_ep(n, x, NULL, &SS_ref_db); 
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
   nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, fl_ig_c, NULL, SS_ref_db.tol_sf);
   nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
   nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
   
   double minf;
   if (gv.maxeval==1){  
     // we are only interested in evaluating the objective function  
     minf = obj_ig_fl(n, x, NULL, &SS_ref_db); 
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
   nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, g_ig_c, NULL, SS_ref_db.tol_sf);
   nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
   nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
   
   double minf;
   if (gv.maxeval==1){  
     // we are only interested in evaluating the objective function  
     minf = obj_ig_g(n, x, NULL, &SS_ref_db); 
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
   nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, hb_ig_c, NULL, SS_ref_db.tol_sf);
   nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
   nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
   
   double minf;
   if (gv.maxeval==1){  
     // we are only interested in evaluating the objective function  
     minf = obj_ig_hb(n, x, NULL, &SS_ref_db); 
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
   nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, ilm_ig_c, NULL, SS_ref_db.tol_sf);
   nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
   nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
   
   double minf;
   if (gv.maxeval==1){  
     // we are only interested in evaluating the objective function  
     minf = obj_ig_ilm(n, x, NULL, &SS_ref_db); 
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
   nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, liq_ig_c, NULL, SS_ref_db.tol_sf);
   nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
   nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
   
   double minf;
   if (gv.maxeval==1){  
     // we are only interested in evaluating the objective function  
     minf = obj_ig_liq(n, x, NULL, &SS_ref_db); 
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
   nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, mu_ig_c, NULL, SS_ref_db.tol_sf);
   nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
   nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
   
   double minf;
   if (gv.maxeval==1){  
     // we are only interested in evaluating the objective function  
     minf = obj_ig_mu(n, x, NULL, &SS_ref_db); 
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
   nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, ol_ig_c, NULL, SS_ref_db.tol_sf);
   nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
   nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
   
   double minf;
   if (gv.maxeval==1){  
     // we are only interested in evaluating the objective function  
     minf = obj_ig_ol(n, x, NULL, &SS_ref_db); 
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
   nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, opx_ig_c, NULL, SS_ref_db.tol_sf);
   nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
   nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
   
   double minf;
   if (gv.maxeval==1){  
     // we are only interested in evaluating the objective function  
     minf = obj_ig_opx(n, x, NULL, &SS_ref_db); 
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

SS_ref NLopt_opt_ig_pl4T_function(global_variable gv, SS_ref SS_ref_db){
   
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
   nlopt_set_min_objective(SS_ref_db.opt, obj_ig_pl4T, &SS_ref_db);
   nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, pl4T_ig_c, NULL, SS_ref_db.tol_sf);
   nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
   nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
   
   double minf;
   if (gv.maxeval==1){  
     // we are only interested in evaluating the objective function  
     minf = obj_ig_pl4T(n, x, NULL, &SS_ref_db); 
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

SS_ref NLopt_opt_ig_spn_function(global_variable gv, SS_ref SS_ref_db){

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
   nlopt_set_min_objective(SS_ref_db.opt, obj_ig_spn, &SS_ref_db);
   nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, spn_ig_c, &SS_ref_db, SS_ref_db.tol_sf);
   nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
   nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
   
   double minf;

   if (gv.maxeval==1){  
     // we are only interested in evaluating the objective function  
     minf = obj_ig_spn(n, x, NULL, &SS_ref_db); 
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
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, fluid_um_c, NULL, SS_ref_db.tol_sf);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    if (gv.maxeval==1){  
       // we are only interested in evaluating the objective function  
       minf = obj_um_fluid(n, x, NULL, &SS_ref_db);
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
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, ol_um_c, NULL, SS_ref_db.tol_sf);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    if (gv.maxeval==1){  
       // we are only interested in evaluating the objective function  
       minf = obj_um_ol(n, x, NULL, &SS_ref_db);
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
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, br_um_c, NULL, SS_ref_db.tol_sf);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    if (gv.maxeval==1){  
       // we are only interested in evaluating the objective function  
       minf = obj_um_br(n, x, NULL, &SS_ref_db);
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
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, ch_um_c, NULL, SS_ref_db.tol_sf);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    if (gv.maxeval==1){  
       // we are only interested in evaluating the objective function  
       minf = obj_um_ch(n, x, NULL, &SS_ref_db);
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
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, atg_um_c, NULL, SS_ref_db.tol_sf);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    if (gv.maxeval==1){  
       // we are only interested in evaluating the objective function  
       minf = obj_um_atg(n, x, NULL, &SS_ref_db);
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
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, g_um_c, NULL, SS_ref_db.tol_sf);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    if (gv.maxeval==1){  
       // we are only interested in evaluating the objective function  
       minf = obj_um_g(n, x, NULL, &SS_ref_db);
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
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, ta_um_c, NULL, SS_ref_db.tol_sf);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    if (gv.maxeval==1){  
       // we are only interested in evaluating the objective function  
       minf = obj_um_ta(n, x, NULL, &SS_ref_db);
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
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, chl_um_c, NULL, SS_ref_db.tol_sf);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    if (gv.maxeval==1){  
       // we are only interested in evaluating the objective function  
       minf = obj_um_chl(n, x, NULL, &SS_ref_db);
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
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, anth_um_c, NULL, SS_ref_db.tol_sf);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    if (gv.maxeval==1){  
       // we are only interested in evaluating the objective function  
       minf = obj_um_anth(n, x, NULL, &SS_ref_db);
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
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, spi_um_c, NULL, SS_ref_db.tol_sf);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    if (gv.maxeval==1){  
       // we are only interested in evaluating the objective function  
       minf = obj_um_spi(n, x, NULL, &SS_ref_db);
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
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, opx_um_c, NULL, SS_ref_db.tol_sf);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    if (gv.maxeval==1){  
       // we are only interested in evaluating the objective function  
       minf = obj_um_opx(n, x, NULL, &SS_ref_db);
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
    nlopt_add_inequality_mconstraint(SS_ref_db.opt, m, po_um_c, NULL, SS_ref_db.tol_sf);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    nlopt_set_maxtime(SS_ref_db.opt, gv.maxgmTime);

    double minf;
    if (gv.maxeval==1){  
       // we are only interested in evaluating the objective function  
       minf = obj_um_po(n, x, NULL, &SS_ref_db);
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

/** 
  attributes the right solution phase to the solution phase array and calculates xi
*/
SS_ref NLopt_opt_function(		global_variable   gv,
								              SS_ref 			      SS_ref_db, 
								              int     		      index			){
								
	clock_t t; 
	t = clock();

	/* Associate the right solid-solution data */
  if(gv.EM_database == 0){ 
    if 		(strcmp( gv.SS_list[index], "liq") == 0 ){
      SS_ref_db  = NLopt_opt_mp_liq_function( gv, SS_ref_db);	}
    else if (strcmp( gv.SS_list[index], "pl4tr")  == 0){
      SS_ref_db  = NLopt_opt_mp_pl4tr_function( gv, SS_ref_db);	}
    else if (strcmp( gv.SS_list[index], "bi")  == 0){
      SS_ref_db  = NLopt_opt_mp_bi_function( gv, SS_ref_db);	}
    else if (strcmp( gv.SS_list[index], "g")  == 0){
      SS_ref_db  = NLopt_opt_mp_g_function( gv, SS_ref_db);	}
    else if (strcmp( gv.SS_list[index], "ep")  == 0){
      SS_ref_db  = NLopt_opt_mp_ep_function( gv, SS_ref_db);	}
    else if (strcmp( gv.SS_list[index], "ma")  == 0){
      SS_ref_db  = NLopt_opt_mp_ma_function( gv, SS_ref_db);	}
    else if (strcmp( gv.SS_list[index], "mu")  == 0){
      SS_ref_db  = NLopt_opt_mp_mu_function( gv, SS_ref_db);	}
    else if (strcmp( gv.SS_list[index], "opx")  == 0){
      SS_ref_db  = NLopt_opt_mp_opx_function( gv, SS_ref_db);	}
    else if (strcmp( gv.SS_list[index], "sa")  == 0){
      SS_ref_db  = NLopt_opt_mp_sa_function( gv, SS_ref_db);	}
    else if (strcmp( gv.SS_list[index], "cd")  == 0){
      SS_ref_db  = NLopt_opt_mp_cd_function( gv, SS_ref_db);	}
    else if (strcmp( gv.SS_list[index], "st")  == 0){
      SS_ref_db  = NLopt_opt_mp_st_function( gv, SS_ref_db);	}
    else if (strcmp( gv.SS_list[index], "chl")  == 0){
      SS_ref_db  = NLopt_opt_mp_chl_function( gv, SS_ref_db);	}
    else if (strcmp( gv.SS_list[index], "ctd")  == 0){
      SS_ref_db  = NLopt_opt_mp_ctd_function( gv, SS_ref_db);	}
    else if (strcmp( gv.SS_list[index], "sp")  == 0){
      SS_ref_db  = NLopt_opt_mp_sp_function( gv, SS_ref_db);	}
    else if (strcmp( gv.SS_list[index], "ilm")  == 0){
      SS_ref_db  = NLopt_opt_mp_ilm_function( gv, SS_ref_db);	}
    else if (strcmp( gv.SS_list[index], "mt")  == 0){
      SS_ref_db  = NLopt_opt_mp_mt_function( gv, SS_ref_db);	}
    else if (strcmp( gv.SS_list[index], "mt") == 0){
      SS_ref_db  = NLopt_opt_mp_mt_function( gv, SS_ref_db);	
      }
    else{
      printf("\nsolid solution '%s index %d' is not in the database\n",gv.SS_list[index], index);	
      }
  }	
  else if(gv.EM_database == 2){          // igneous
    if 		(strcmp( gv.SS_list[index], "bi") == 0 ){
      SS_ref_db  = NLopt_opt_ig_bi_function( gv, SS_ref_db);	}
    else if (strcmp( gv.SS_list[index], "cd")  == 0){
      SS_ref_db  = NLopt_opt_ig_cd_function( gv, SS_ref_db);	}
    else if (strcmp( gv.SS_list[index], "cpx") == 0){
      SS_ref_db  = NLopt_opt_ig_cpx_function( gv, SS_ref_db);}	
    else if (strcmp( gv.SS_list[index], "ep")  == 0){
      SS_ref_db  = NLopt_opt_ig_ep_function( gv, SS_ref_db);	}
    else if (strcmp( gv.SS_list[index], "fl")  == 0){
      SS_ref_db  = NLopt_opt_ig_fl_function( gv, SS_ref_db);	}		
    else if (strcmp( gv.SS_list[index], "g")   == 0){
      SS_ref_db  = NLopt_opt_ig_g_function(  gv, SS_ref_db);	}
    else if (strcmp( gv.SS_list[index], "hb")  == 0){
      SS_ref_db  = NLopt_opt_ig_hb_function( gv, SS_ref_db);	}	
    else if (strcmp( gv.SS_list[index], "ilm") == 0){
      SS_ref_db  = NLopt_opt_ig_ilm_function( gv, SS_ref_db);}
    else if (strcmp( gv.SS_list[index], "liq") == 0){
      SS_ref_db  = NLopt_opt_ig_liq_function( gv, SS_ref_db);}
    else if (strcmp( gv.SS_list[index], "mu")  == 0){
      SS_ref_db  = NLopt_opt_ig_mu_function( gv, SS_ref_db);	}	
    else if (strcmp( gv.SS_list[index], "ol")  == 0){
      SS_ref_db  = NLopt_opt_ig_ol_function( gv, SS_ref_db);	}
    else if (strcmp( gv.SS_list[index], "opx") == 0){
      SS_ref_db  = NLopt_opt_ig_opx_function( gv, SS_ref_db);}	
    else if (strcmp( gv.SS_list[index], "pl4T")  == 0){
      SS_ref_db  = NLopt_opt_ig_pl4T_function( gv, SS_ref_db);	}	
    else if (strcmp( gv.SS_list[index], "spn") == 0){
      SS_ref_db  = NLopt_opt_ig_spn_function( gv, SS_ref_db);	
      }
    else{
      printf("\nsolid solution '%s index %d' is not in the database\n",gv.SS_list[index], index);	
      }	
  }
  else if (gv.EM_database == 3){
    if (strcmp( gv.SS_list[index], "liq")  == 0){
        SS_ref_db  = NLopt_opt_igd_liq_function( gv, SS_ref_db);}
    else if (strcmp( gv.SS_list[index], "fl")  == 0){
        SS_ref_db  = NLopt_opt_igd_fl_function( gv, SS_ref_db);}
    else if (strcmp( gv.SS_list[index], "fsp")  == 0){
        SS_ref_db  = NLopt_opt_igd_fsp_function( gv, SS_ref_db);}
    else if (strcmp( gv.SS_list[index], "spn")  == 0){
        SS_ref_db  = NLopt_opt_igd_spn_function( gv, SS_ref_db);}
    else if (strcmp( gv.SS_list[index], "g")  == 0){
        SS_ref_db  = NLopt_opt_igd_g_function( gv, SS_ref_db);}
    else if (strcmp( gv.SS_list[index], "ol")  == 0){
        SS_ref_db  = NLopt_opt_igd_ol_function( gv, SS_ref_db);}
    else if (strcmp( gv.SS_list[index], "opx")  == 0){
        SS_ref_db  = NLopt_opt_igd_opx_function( gv, SS_ref_db);}
    else if (strcmp( gv.SS_list[index], "cpx")  == 0){
        SS_ref_db  = NLopt_opt_igd_cpx_function( gv, SS_ref_db);}
    else if (strcmp( gv.SS_list[index], "ilm")  == 0){
        SS_ref_db  = NLopt_opt_igd_ilm_function( gv, SS_ref_db);}
    else if (strcmp( gv.SS_list[index], "hb")  == 0){
        SS_ref_db  = NLopt_opt_igd_hb_function( gv, SS_ref_db);}
    else if (strcmp( gv.SS_list[index], "bi")  == 0){
        SS_ref_db  = NLopt_opt_igd_bi_function( gv, SS_ref_db);}
    else if (strcmp( gv.SS_list[index], "ep")  == 0){
        SS_ref_db  = NLopt_opt_igd_ep_function( gv, SS_ref_db);}
    else if (strcmp( gv.SS_list[index], "cd")  == 0){
        SS_ref_db  = NLopt_opt_igd_cd_function( gv, SS_ref_db);}
    else{
        printf("\nsolid solution '%s index %d' is not in the database\n",gv.SS_list[index], index);
    }
  }
  else if (gv.EM_database == 4){
      if (strcmp( gv.SS_list[index], "fluid")  == 0){
         SS_ref_db  = NLopt_opt_um_fluid_function( gv, SS_ref_db);}
      else if (strcmp( gv.SS_list[index], "ol")  == 0){
         SS_ref_db  = NLopt_opt_um_ol_function( gv, SS_ref_db);}
      else if (strcmp( gv.SS_list[index], "br")  == 0){
         SS_ref_db  = NLopt_opt_um_br_function( gv, SS_ref_db);}
      else if (strcmp( gv.SS_list[index], "ch")  == 0){
         SS_ref_db  = NLopt_opt_um_ch_function( gv, SS_ref_db);}
      else if (strcmp( gv.SS_list[index], "atg")  == 0){
         SS_ref_db  = NLopt_opt_um_atg_function( gv, SS_ref_db);}
      else if (strcmp( gv.SS_list[index], "g")  == 0){
         SS_ref_db  = NLopt_opt_um_g_function( gv, SS_ref_db);}
      else if (strcmp( gv.SS_list[index], "ta")  == 0){
         SS_ref_db  = NLopt_opt_um_ta_function( gv, SS_ref_db);}
      else if (strcmp( gv.SS_list[index], "chl")  == 0){
         SS_ref_db  = NLopt_opt_um_chl_function( gv, SS_ref_db);}
      else if (strcmp( gv.SS_list[index], "anth")  == 0){
         SS_ref_db  = NLopt_opt_um_anth_function( gv, SS_ref_db);}
      else if (strcmp( gv.SS_list[index], "spi")  == 0){
         SS_ref_db  = NLopt_opt_um_spi_function( gv, SS_ref_db);}
      else if (strcmp( gv.SS_list[index], "opx")  == 0){
         SS_ref_db  = NLopt_opt_um_opx_function( gv, SS_ref_db);}
      else if (strcmp( gv.SS_list[index], "po")  == 0){
         SS_ref_db  = NLopt_opt_um_po_function( gv, SS_ref_db);}
      else{
         printf("\nsolid solution '%s index %d' is not in the database\n",gv.SS_list[index], index);
      }
   }
   t = clock() - t; 
   SS_ref_db.LM_time = ((double)t)/CLOCKS_PER_SEC*1000; // in seconds 

	return SS_ref_db;
};


