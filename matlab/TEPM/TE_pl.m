function [Di]= TE_pl(constant,fpe,T,P,ri_CN_8,comp,mc,ri_8,order_flag,TE_name,fig,p,clr)
%% Constant and input
% Reference Temperature, Dohmen and Blundy, 2014
Tr= 1565; % kelvin
SiO2_melt= comp(1,1);
MgO_melt= comp(1,6);
%% Structural formula and XAn  
% fudge factor 
Factor_pl= 16./(4.*mc(2,1)+4.*mc(2,2)+6.*mc(2,3)+2.*mc(2,4)+2.*mc(2,5)+2.*mc(2,6)+2.*mc(2,7)+2.*mc(2,8)+2.*mc(2,9));
% Molar fraction Na and Ca
XCa= Factor_pl.*mc(2,7); % structural formula of pl       
XNa= 2.*Factor_pl.*mc(2,8); 
XK= 2.*Factor_pl.*mc(2,9);    
XAn= XCa./(XCa+XNa+XK);
XAb= XNa./(XCa+XNa+XK);
XOr= XK./(XCa+XNa+XK);

DCa= exp((-19107 +467.*SiO2_melt)/(constant(2)*T));
DNa= exp((13621-18990.*XAn)/(constant(2)*T));
%% REE prediction + Transitive metals
% r0 3+
ro_3=(1.331-0.068.*XAn);
% E3+
E_3= (152-31.*XAn);
% D0 3+ = DLa 
DLa=((DCa.^2)./DNa).*exp((4400./(constant(2).*T))-(30.8./constant(2)));
% Di 3+
Di_REE_Y= DLa.*exp(fpe.*E_3.*(ro_3*0.5.*(ri_CN_8{1,1}(1).^2-ri_CN_8{1,1}(1:15).^2)-1/3.*(ri_CN_8{1,1}(1).^3-ri_CN_8{1,1}(1:15).^3))./T);
REE_plot= DLa.*exp(fpe.*E_3.*(ro_3.*0.5.*(ri_CN_8{1,1}(1).^2-ri_8{1,1}.^2)-1/3.*(ri_CN_8{1,1}(1).^3-ri_8{1,1}.^3))./T);
D_Sc= exp((-23000 - 7550*MgO_melt)/(constant(2)*T));
Di_REE_Y_Sc=[Di_REE_Y D_Sc];
%% LILE prediction
% 1+ cations
% r0 1+
ro_1=(1.24-0.017.*XAn);
% E1+
E_1=(49.05+17.16.*XAn);
% D0 1+ = DNa 
% Di 1+
Di_LILE_1=DNa.*exp(fpe.*E_1.*(ro_1.*0.5.*(1.18^2-ri_CN_8{1,2}.^2)-1/3.*(1.18.^3-ri_CN_8{1,2}.^3))./T);
LILE1_plot= DNa.*exp(fpe.*E_1.*(ro_1.*0.5.*(1.18^2-ri_8{1,2}.^2)-1/3.*(1.18^3-ri_8{1,2}.^3))./T);
% 2+ cations
% r0 2+  
ro_2= 1.2895+1.3e-4.*(T-Tr)+XAn.*(-0.0952-4e-5.*(T-Tr));
% E2+
E_2 = (120.03824-0.3686.*XAn);
% D0 2+ = DCa 
% Di 2+
Di_LILE_2=DCa.*exp(fpe.*E_2.*(ro_2.*0.5*(1.12^2-ri_CN_8{1,3}.^2)-1/3.*(1.12^3-ri_CN_8{1,3}.^3))./T);
LILE2_plot= DCa.*exp(fpe.*E_2.*(ro_2.*0.5.*(1.12^2-ri_8{1,3}.^2)-1/3.*(1.12^3-ri_8{1,3}.^3))./T);
%% HFSE prediction 
% 4+ cations
D_Ti = exp((-15400-28900.*(XAn))./(constant(2).*T));
% D_Zr1= exp((-15300-90400.*(XAn))./(constant(2).*T))
D_Zr= exp((-187270+2260.*SiO2_melt)/(constant(2).*T));
D_Hf= exp((-126860+1399.*SiO2_melt)/(constant(2).*T));
D_HFSE_4=[D_Ti D_Hf D_Zr];
% 5+ cations
% D_Ta= exp((-271.6+3.157.*SiO2_melt+8.926.*MgO_melt+26.501.*XAn)/(constant(2).*T));
% D_Nb=exp((-417.35+5.3.*SiO2_melt+10.34.*MgO_melt+40.88.*XAn)/(constant(2).*(T)));
D_Ta= exp((-271600+3157.*SiO2_melt+8926.*MgO_melt+26501.*XAn)/(constant(2).*T));
D_Nb=exp((-417350+5300.*SiO2_melt+10340.*MgO_melt+40880.*XAn)/(constant(2).*T));
D_HFSE_5=[D_Ta D_Nb];
%% Actinides
D_U= exp((-7450-48400.*(XAn))./(constant(2).*T));
D_Th=exp((-11646-31369.*(XAn))./(constant(2).*T));
D_Act= [D_U D_Th]; % add D_Th  

%% Concatenate the Di 
Di1= {Di_LILE_1 Di_LILE_2 Di_REE_Y_Sc D_HFSE_4 D_Act D_HFSE_5};
Di=cell2mat(Di1);

Di_plot={LILE1_plot LILE2_plot REE_plot D_HFSE_4 D_Act};
% Di_plot=cell2mat(Di2);
%% plot
if fig==1
plplot(T,P,XAn,XAb,XOr,Di,Di_plot,ri_CN_8,ri_8,order_flag,TE_name,p,clr)
end