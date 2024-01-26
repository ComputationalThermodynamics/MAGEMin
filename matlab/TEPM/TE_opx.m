function [Di]= TE_opx(constant,T,P,ri_CN_8,ri_CN_6,comp,mc,fpe,ri_8,ri_6,order_flag,TE_name,fig,p,clr)
%% Constant and input 
melt_SiO2=comp(1,1);
melt_MgO=comp(1,6);
melt_Al2O3=comp(1,3);
melt_FeO=comp(1,4);
Mg_n=melt_MgO./(melt_MgO+melt_FeO);
wt_MgO_Opx=comp(2,6);
wt_FeO_Opx=comp(2,4);
Mg_n_opx=wt_MgO_Opx./(wt_MgO_Opx+wt_FeO_Opx);
% DMg= comp(2,5)/comp(1,5)
%% Structural formula 
% fudge factor 
Factor_opx= 12./(4.*mc(2,1)+4.*mc(2,2)+6.*mc(2,3)+2.*mc(2,4)+2.*mc(2,5)+2.*mc(2,6)+2.*mc(2,7)+2.*mc(2,8)+2.*mc(2,9));
% Structural formula   
XCa = Factor_opx*mc(2,7); % structural formula of opx
XFe = Factor_opx*mc(2,4);
XMg = Factor_opx*mc(2,6);
XSi = Factor_opx*mc(2,1);
XAl=2*Factor_opx*mc(2,3);
XNa=2*Factor_opx*mc(2,8);
XK= 2*Factor_opx*mc(2,9);
XAl4= 2-XSi;
XWo= XCa/(XCa+XFe+XMg);
XEn= XMg/(XCa+XFe+XMg);
XFs= XFe/(XCa+XFe+XMg);
X_Fe_Mg=1-(XCa+XNa);
XMg_M2=X_Fe_Mg*Mg_n_opx;
XFe_M2= 1-(XCa+XNa+XK+XMg_M2);
XMg_M1=XMg-XMg_M2;
XFe_M1=XFe-XFe_M2;

% DMg= XMg_M2/mc(1,6)
% if melt_SiO2>56
    DMg=1; 
% end

%% REE prediction + Sc
% r0 3+
% ro_3=(0.69+0.43.*XCa+0.23.*XMg);
ro_3=(0.69+0.43.*XCa+0.23.*XMg_M2);
% E3+
E_3= ((-1.37+1.85.*ro_3-0.53.*XCa).*10^3);
% D0 3+ = DLa 
Do=exp(-5.37+(38700./(constant(2).*T))+3.56.*XCa+3.54.*XAl4);
% Di 3+
Di_REE_Y= Do.*exp(fpe.*E_3.*(ro_3.*0.5.*((ri_CN_8{1,1}(1:15)-ro_3).^2)+1/3.*((ri_CN_8{1,1}(1:15)-ro_3).^3))./T);
D_Sc=exp(5.04-4.37*Mg_n_opx-0.41*log(wt_FeO_Opx)-0.318*log(melt_MgO));
Di_REE_Y_Sc=[Di_REE_Y D_Sc];
REE_plot=Do.*exp(fpe.*E_3.*(ro_3.*0.5.*((ri_8{1,1}-ro_3).^2)+1/3.*((ri_8{1,1}-ro_3).^3))./T);
%% LILE prediction
% 1+ cations
% D_Cs
D_Cs= exp(23.6-1.46.*P-0.26.*melt_SiO2-32.49.*Mg_n);
if D_Cs>1
D_Cs = exp(2.901-17.32*Mg_n);
end
% D_Rb
D_Rb= exp(6.995148 - 0.1015.*melt_SiO2 - 0.0571.*melt_Al2O3 + 0.4567.*log(melt_MgO) - 16.0222.*Mg_n);
% D_K
D_K=exp(-13.04737+0.9458.*P+0.4307.*XWo+0.13.*melt_SiO2-0.9374.*log(melt_FeO)+1.3376.*log(melt_MgO));

Di_LILE_1= [D_Cs D_Rb D_K];
% 2+ cations
ro_2= ro_3+0.12;
E_2=2/3*E_3;
Di_LILE_2= DMg.*exp(fpe.*E_2.*(ro_2.*0.5.*(0.89^2-ri_CN_8{1,3}.^2)-1/3.*(0.89^3-ri_CN_8{1,3}.^3))./T); %DMg=1 between basaltic melt and opx, Wood and Blundy,2013
LILE2_plot=DMg.*exp(fpe.*E_2.*(ro_2.*0.5.*(0.89^2-ri_8{1,3}.^2)-1/3.*(0.89^3-ri_8{1,3}.^3))./T);

%% HFSE prediction 
% 4+ cations
    %r0_4+
    ro_4=(0.618+0.032.*XCa+0.030.*XMg_M1);
    %E_4
    E_4= 2203;
    %D0_4+
    Do_4=exp(-4.825+(31780./(constant(2).*T))+4.17.*XAl4+8.55.*XCa.*XMg_M2-2.62.*XFe_M1);
    %D_HFSE    
    D_HFSE_4= Do_4.*exp(fpe.*E_4.*(ro_4.*0.5.*((ri_CN_6{1,4}-ro_4).^2)+1/3.*((ri_CN_6{1,4}-ro_4).^3))./T);
    HFSE4_plot=Do_4.*exp(fpe.*E_4.*(ro_4.*0.5.*((ri_6{1,4}-ro_4).^2)+1/3.*((ri_6{1,4}-ro_4).^3))./T);
    D_HFSE_4= D_HFSE_4(1:3) ;
% 5+ cations
    D_Ta=exp(-1.68+0.62.*XAl4);
    D_Nb=D_Ta.*(1.17-3.16.*XAl4);
    if D_Nb<0
       XAl4=0.369;
       D_Nb=D_Ta.*(1.17-3.16.*XAl4);
    end
    D_HFSE_5=[D_Ta D_Nb];    
%% Actinides
D_U= exp(3.967-11.668.*Mg_n_opx);
D_Th_1= -6.713332 + 19.48971.*XAl4 + 2.11.*log(melt_FeO) - 2.55.*log(melt_MgO); 
D_Th=exp(D_Th_1);
D_Act= [D_U D_Th];  
%% Concatenate the Di 
Di1= { Di_LILE_1 Di_LILE_2 Di_REE_Y_Sc D_HFSE_4 D_Act D_HFSE_5};
Di=cell2mat(Di1);
Di_plot={Di_LILE_1 LILE2_plot REE_plot HFSE4_plot D_Act};
% Di_plot=cell2mat(Di2);
%% plot
if fig==1
opxplot(T,P,XWo,XEn,XFs,Di,Di_plot,ri_CN_8,ri_CN_6,ri_8,ri_6,order_flag,TE_name,p,clr)
end