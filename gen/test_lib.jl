include("magemin_library.jl")


# Initialize 
gv = LibMAGEMin.global_variable_init();


test    = 0;
P       = 10.1
T       = 800.1
LibMAGEMin.get_bulk(bulk_rock, test, gv.len_ox)
LibMAGEMin.norm_array(		bulk_rock, gv.len_ox);	
z_b = LibMAGEMin.zeros_in_bulk(	bulk_rock, P, T);
