using MAGEMin_C

# Initialize database  - new way
dtb          =  "ig"
data        =   Initialize_MAGEMin(dtb, verbose=true);
test        =   0         #KLB1
data        =   use_predefined_bulk_rock(data, test);

# Call optimization routine for given P & T & bulk_rock
P           =   8.0
T           =   1400.0
out         =   point_wise_minimization(P,T, data);

elem_TE     = ["SiO2","TiO2","Al2O3","O","FeO","MnO","MgO","CaO","K2O","Na2O","H2O","Li","Be","B","Sc","V","Cr","Ni","Cu","Zn","Rb","Sr","Y","Zr","Nb","Cs","Ba","La","Ce","Pr","Nd","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","Pb","Th","U"]
bulk_TE     = [50.60777553,0.953497243,13.70435413,0.19,11.28130762,0.202560796,8.496024312,9.502380068,0.700881685,2.07434927,4,29.14258603,0.434482763,29.69200003,38.23663423,257.4346716,529.333066,208.2057375,88.87615683,91.7592182,16.60777308,163.4533209,20.74016207,66.90677472,3.808354064,1.529226981,122.8449739,6.938172601,16.04827796,2.253943183,10.18276823,3.3471043,0.915941652,3.28230146,1.417695298,3.851230952,0.914966282,2.20425,0.343734976,2.136202593,0.323405135,1.841502082,0.330971265,5.452969044,1.074692058,0.290233271]

KDs_dtb     = get_OL_KDs_database()


C0          = adjust_chemical_system(    KDs_dtb,bulk_TE,elem_TE)

out_TE      = TE_prediction(C0,KDs_dtb, "CB",out,dtb)


if ~isnothing(out_TE.bulk_cor_wt)      #then we can recompute the equilibrium after removing the SiO2 entering zircon
    sys_in      = "wt"    
    out         = single_point_minimization(P, T, data, X=out_TE.bulk_cor_wt, Xoxides=out.oxides, sys_in=sys_in)
    out_TE      = TE_prediction(C0,KDs_dtb, "CB",out,dtb)
end
