/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **   Project      : MAGEMin
 **   License      : GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007
 **   Developers   : Nicolas Riel, Boris Kaus
 **   Contributors : Moccetti, N. B., Dominguez, H., Assunção J., Green E., Dolejš, D., Berlie N., and Rummel L.
 **   Organization : Institute of Geosciences, Johannes-Gutenberg University, Mainz
 **   Contact      : nriel[at]uni-mainz.de, kaus[at]uni-mainz.de
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/
#ifndef __BR_GSS_FUNCTION_H_
#define __BR_GSS_FUNCTION_H_

#include "../initialize.h"

SS_ref G_SS_br_ctd_function(SS_ref SS_ref_db, char *research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps);
SS_ref G_SS_br_car_function(SS_ref SS_ref_db, char *research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps);
SS_ref G_SS_br_chl_function(SS_ref SS_ref_db, char *research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps);
SS_ref G_SS_br_mica_function(SS_ref SS_ref_db, char *research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps);
SS_ref G_SS_br_talc_function(SS_ref SS_ref_db, char *research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps);
SS_ref G_SS_br_ilm_function(SS_ref SS_ref_db, char *research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps);
SS_ref G_SS_br_bt_function(SS_ref SS_ref_db, char *research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps);
SS_ref G_SS_br_ol_function(SS_ref SS_ref_db, char *research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps);
SS_ref G_SS_br_ep_function(SS_ref SS_ref_db, char *research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps);
SS_ref G_SS_br_opx_function(SS_ref SS_ref_db, char *research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps);
SS_ref G_SS_br_amph_function(SS_ref SS_ref_db, char *research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps);
SS_ref G_SS_br_spl_function(SS_ref SS_ref_db, char *research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps);
SS_ref G_SS_br_stau_function(SS_ref SS_ref_db, char *research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps);
SS_ref G_SS_br_crd_function(SS_ref SS_ref_db, char *research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps);
SS_ref G_SS_br_grt_function(SS_ref SS_ref_db, char *research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps);
SS_ref G_SS_br_omph_function(SS_ref SS_ref_db, char *research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps);
SS_ref G_SS_br_amphx_function(SS_ref SS_ref_db, char *research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps);
SS_ref G_SS_br_fsp_function(SS_ref SS_ref_db, char *research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps);

SS_ref G_SS_br_EM_function( global_variable gv, SS_ref SS_ref_db, int EM_dataset, bulk_info z_b, char *name);

#endif
