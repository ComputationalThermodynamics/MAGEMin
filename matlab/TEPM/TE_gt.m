function [Di]= TE_gt(constant,T,P,ri_CN_8,ri_CN_6,mc,fpe,comp,ri_8,ri_6,order_flag,TE_name,fig,p,clr)
%% Constant and input  
XFe_melt=mc(1,4);
XMg_melt=mc(1,6);
XSi_melt=mc(1,1);
Melt_polymerisation= 0.008458.*comp(1,1)+0.1734;  
DFe=comp(1,4)/comp(2,4);
%% Structural formula     
% fudge factor 
Factor_gt= 24./(4.*mc(2,1)+4.*mc(2,2)+6.*mc(2,3)+2.*mc(2,4)+2.*mc(2,5)+2.*mc(2,6)+2.*mc(2,7)+2.*mc(2,8)+2.*mc(2,9));
% Structural formula   
XCa = Factor_gt*mc(2,7); 
XFe = Factor_gt*mc(2,4);
XMg = Factor_gt*mc(2,6);  
XAl = 2*Factor_gt*mc(2,3); 
XGr= XCa./(XCa+XMg+XFe);
XPy=XMg./(XCa+XMg+XFe);
XAlm=XFe./(XCa+XMg+XFe);
% DMg= mc(2,6)/mc(1,6); 
DMg= exp((25.820-0.1415*T+5.418*P)/(3*constant(2)*T))/exp((19.000*((XCa)^2))/(constant(2)*T));
%% REE prediction + Sc
% if XAlm<XPy
%     fprintf('Pyrope')
% end
 % r0 3+
ro_3=(0.780+0.155.*XCa);
% E3+
E_3= ((-1.62+2.29.*ro_3).*10^3);
% D0 3+ = DLa 
Do=exp(-2.75+((91700-91.35.*P.*(38-P))./(constant(2).*T))-1.42.*XCa)  ;
% else
%  % r0 3+
% ro_3_1=0.9302*XPy+0.993*XGr+0.916*XAlm-0.0044*(P-3)+0.000058*(T-1818)
% % E3+
% E_3= 2.826*(1.38+ro_3)^-3+12.4*P-0.072*T+237*(XAl)
% % D0 3+ = DLa 
% gamma_Fe=exp(19000*(XCa^2)/(constant(2)*T))
% Do1=exp((400.290+4.586*P-0.218*T)./(constant(2).*T))/((gamma_Fe*DFe)^2)   
% end
Di_REE_Y= Do.*exp(fpe.*E_3.*(ro_3.*0.5.*((ro_3-ri_CN_8{1,1}(1:15)).^2)-1/3.*((ro_3-ri_CN_8{1,1}(1:15)).^3))./T);
REE_plot= Do.*exp(fpe.*E_3.*(ro_3.*0.5.*((ro_3-ri_8{1,1}).^2)-1/3.*((ro_3-ri_8{1,1}).^3))./T);
DSc=5.79; %Adam and Green, 2006
Di_REE_Y_Sc=[Di_REE_Y DSc];
%% LILE prediction
% 1+ cations
% D_Cs
D_Cs= 0.001;
% D_Rb
D_Rb= exp(0.47029 -4.644.*XMg+0.656.*(XMg.^2)); %rsquared=0.57
% D_K
D_K=exp(4.72-1.4.*P-50.*XGr+8.09.*P.*XGr);%rsquared=0.57
Di_LILE_1= [D_Cs D_Rb D_K];
% 2+ cations
% r0 2+
ro_2= (ro_3+0.12);
% ro_2= 1.02;
% E2+
E_2 = 2/3.*E_3;
% D0 2+ = DMg
% Di 2+
Di_LILE_2= DMg.*exp(fpe.*E_2.*(ro_2.*0.5.*(0.89^2-ri_CN_8{1,3}.^2)-1/3.*(0.89^3-ri_CN_8{1,3}.^3))./T);
LILE2_plot= DMg.*exp(fpe.*E_2.*(ro_2.*0.5.*(0.89^2-ri_8{1,3}.^2)-1/3.*(0.89^3-ri_8{1,3}.^3))./T);
%% Actinides
% D_Th=((XFe_melt+XMg_melt)^4)*(XSi_melt.^2)*exp(2.*(11.46-(24200/T)...
%     +8.6*((1-(XFe_melt+XMg_melt))^2)-2.08*((1-XGr)^2))) % Salters and Longhi, 2002
D_Th=((XFe_melt+XMg_melt)^4)*(XSi_melt.^2)*((11.46-(24200/T)...
    +8.6*((1-(XFe_melt+XMg_melt))^2)-2.08*((1-XGr)^2))^2); % Salters and Longhi, 2002
D_U= 3.6.*D_Th+0.003; % Elkins et al, 2008
D_Act= [D_U D_Th];  
%% HFSE prediction 
% 4+ cations
D_Ti= 0.0037.*exp(7.386.*Melt_polymerisation);
%r0_4+
ro_X_4=ro_3-0.1;
ro_Y_4=0.67;
%E_4
E_X_4=1325;
E_Y_4= 15870.*exp(-0.1086.*(XGr.*100))+409.6.*exp(0.01191.*(XGr.*100));  
a= XGr>0.19 && XGr<0.4 ;  
% special=linspace(0.575,0.865,5000);
special1=linspace(0.575,0.725,2000);
special2=linspace(0.725,0.865,2000);
special=[special1; special2];
if a ==1 
    D_Hf=D_Th.*exp(fpe.*E_X_4.*(ro_X_4.*0.5.*(ri_CN_8{1,4}(1).^2-ri_CN_8{1,4}(2).^2)-1/3.*(ri_CN_8{1,4}(1).^3-ri_CN_8{1,4}(2).^3))./T)+...
         D_Ti.*exp(fpe.*E_Y_4.*(ro_Y_4.*0.5.*(ri_CN_6{1,4}(1).^2-ri_CN_6{1,4}(2).^2)-1/3.*(ri_CN_6{1,4}(1).^3-ri_CN_6{1,4}(2).^3))./T);
    D_Zr= D_Th.*exp(fpe.*E_X_4.*(ro_X_4.*0.5.*(ri_CN_8{1,4}(1).^2-ri_CN_8{1,4}(3).^2)-1/3.*(ri_CN_8{1,4}(1).^3-ri_CN_8{1,4}(3).^3))./T)+...
          D_Ti.*exp(fpe.*E_Y_4.*(ro_Y_4.*0.5.*(ri_CN_6{1,4}(1).^2-ri_CN_6{1,4}(3).^2)-1/3.*(ri_CN_6{1,4}(1).^3-ri_CN_6{1,4}(3).^3))./T);
    HFSE4_plot=D_Th.*exp(fpe.*E_X_4.*(ro_X_4.*0.5.*(ri_CN_8{1,4}(1).^2-special2.^2)-1/3.*(ri_CN_8{1,4}(1).^3-special2.^3))./T)+...
               D_Ti.*exp(fpe.*E_Y_4.*(ro_Y_4.*0.5.*(ri_CN_6{1,4}(1).^2-special1.^2)-1/3.*(ri_CN_6{1,4}(1).^3-special1.^3))./T);

else
    D_Hf=D_Ti.*exp(fpe.*E_Y_4.*(ro_Y_4.*0.5.*(ri_CN_6{1,4}(1).^2-ri_CN_6{1,4}(2).^2)-1/3.*(ri_CN_6{1,4}(1).^3-ri_CN_6{1,4}(2).^3))./T);
    D_Zr= D_Ti.*exp(fpe.*E_Y_4.*(ro_Y_4.*0.5.*(ri_CN_6{1,4}(1).^2-ri_CN_6{1,4}(3).^2)-1/3.*(ri_CN_6{1,4}(1).^3-ri_CN_6{1,4}(3).^3))./T);
    HFSE4_plot=D_Ti.*exp(fpe.*E_Y_4.*(ro_Y_4.*0.5.*(ri_CN_6{1,4}(1).^2-ri_6{1,4}.^2)-1/3.*(ri_CN_6{1,4}(1).^3-ri_6{1,4}.^3))./T);
end
D_HFSE_4= [D_Ti D_Hf D_Zr];
% 5+ cations
D_Ta=0.0146;
D_Nb=0.0290;
D_HFSE_5=[D_Ta D_Nb];    
   
%% Concatenate the Di 
Di1= {Di_LILE_1 Di_LILE_2 Di_REE_Y_Sc D_HFSE_4 D_Act D_HFSE_5};
Di=cell2mat(Di1);

Di_plot={LILE2_plot REE_plot HFSE4_plot D_Act Di_LILE_1};
% Di_plot=cell2mat(Di2);
%% plot
if fig==1
gtplot(T,P,XGr,XPy,XAlm,Di,Di_plot,ri_CN_8,ri_CN_6,ri_8,ri_6,order_flag,TE_name,special,p,clr)
end