function [Di]= TE_cpx(constant,T,P,ri_CN_8,ri_CN_6,comp,mc,fpe,ri_8,ri_6,order_flag,TE_name,fig,p,clr)
%% Constant and input  
%MgO WT% melt, SiO2 WT% melt
X_H2O_melt= mc(1,10);H2O_melt=comp(1,10);
X_Mg_melt=mc(1,6);
X_Fe_melt=mc(1,4);
Mg_n_melt=X_Mg_melt/(X_Mg_melt+X_Fe_melt);
z0=(mc(1,10)+mc(1,9)+mc(1,8)+2.*mc(1,7)+2.*mc(1,6)+2*mc(1,5)+2*mc(1,4)+4*mc(1,2))/15;
% z0=(mc(1,10)+mc(1,9)+mc(1,8)+2.*mc(1,7)+2.*mc(1,6)+2*mc(1,5)+2*mc(1,4)+3*mc(1,3)+4*mc(1,2)+4*mc(1,1))/22
% z0=1;

%% Structural formula 
% fudge factor 
Factor_cpx= 12./(4.*mc(2,1)+4.*mc(2,2)+6.*mc(2,3)+2.*mc(2,4)+2.*mc(2,5)+2.*mc(2,6)+2.*mc(2,7)+2.*mc(2,8)+2.*mc(2,9));

% Structural formula   
XMg = Factor_cpx*mc(2,6);
XFe=  Factor_cpx*mc(2,4);
XCa=  Factor_cpx*mc(2,7);
XAl = 2*Factor_cpx*mc(2,3);
XSi = Factor_cpx*mc(2,1);
XNa = 2*Factor_cpx*mc(2,8);
XK= 2*Factor_cpx*mc(2,9);
XAl4= 2-XSi;
XAl6= XAl-XAl4;
DCa= mc(2,7)/mc(1,7);
% DNa= mc(2,8)/mc(1,8);
% D_Ti=mc(2,2)/mc(1,2);
% if mc(2,2)==0
% fpe_Ti= exp((37800-5375.*P-660.*P^2+2280.*((4-z0).^2))./(constant(2).*T))% Problem with Z0=> value too small 
fpe_Ti=exp((46837-4257*P-9967*P^2)./(constant(2).*T));%first model of Hill et al, 2011
% ghgh=2280*(4-z0)^2
Y=((1-(XK+XNa))*(XAl4)^2+2*(XK+XNa)*(XAl4)*(XSi))+...
    ((XK+XNa)*(XSi)^2+2*(1-(XK+XNa))*XAl6*XSi*XAl4)*exp(-16500/(constant(2).*T))...
    +((1-(XNa+XK))*XSi^2)*exp(-4*16500./(constant(2).*T));
D_Ti = Y.*fpe_Ti;
% end
Mg_n_cpx=(XMg/(XMg+XFe));
% if mc(2,8)==0
    DNa= exp((10367+2100*P-165*P^2)/T-10.27+0.36*P-0.018*P^2);
% end
X_Fe_Mg=1-(XCa+XNa+XK);
XMg_M2=X_Fe_Mg*Mg_n_cpx;
XFe_M2= 1-(XCa+XNa+XK+XMg_M2);
XMg_M1=XMg-XMg_M2;
% XFe_M1=XFe-XFe_M2;
XWo= XCa/(XCa+XFe+XMg);
XEn= XMg/(XCa+XFe+XMg);
XFs= XFe/(XCa+XFe+XMg);
XCa_M2= 1-(XFe_M2+XNa+XMg_M2);
% M2= XCa+XFe_M2+XMg_M2+XNa
% XMg_M1=XMg-XMg_M2;
% XFe_M1=XFe-XFe_M2;
%% REE prediction 
% dG=-7.14+71900/(constant(2)*T)
% dG2=-7+79000/(constant(2)*T)
% disp(T)
% figure(1);
% plot(T,dG,'k.');hold on
% % plot(T,dG2,'r.');
% if comp(1,1)>58
% ro_M2_3=(1.066-0.08.*XAl6-0.212.*XMg_M2);
% E_3= ((2.27.*(ro_M2_3)-2).*10^3);
% Do=exp(-7.14+(71900./(constant(2).*T))+4.80.*XAl4+1.98.*XMg_M2-0.91.*X_H2O_melt);
% % Do=exp(-7.05+(76000./(constant(2).*T))+4.8.*XAl4+1.98.*XMg_M2-0.91.*X_H2O_melt);
% Di_REE_Y= Do.*exp(fpe.*E_3.*(ro_M2_3.*0.5.*((ro_M2_3-ri_CN_8{1,1}(1,1:15)).^2)-1/3.*((ro_M2_3-ri_CN_8{1,1}(1,1:15)).^3))./T);
% REE_plot= Do.*exp(fpe.*E_3.*(ro_M2_3.*0.5.*((ro_M2_3-ri_8{1,1}).^2)-1/3.*((ro_M2_3-ri_8{1,1}).^3))./T);
% else
ro_M2_3=(1.066-0.104.*XAl6-0.212.*XMg_M2);
E_3= ((2.27.*(ro_M2_3)-2).*10^3);
Do=exp(-7.14+(71900./(constant(2).*T))+4.37.*XAl4+1.98.*XMg_M2-0.91.*X_H2O_melt);
Di_REE_Y= Do.*exp(fpe.*E_3.*(ro_M2_3.*0.5.*((ro_M2_3-ri_CN_8{1,1}(1,1:15)).^2)-1/3.*((ro_M2_3-ri_CN_8{1,1}(1,1:15)).^3))./T);
REE_plot= Do.*exp(fpe.*E_3.*(ro_M2_3.*0.5.*((ro_M2_3-ri_8{1,1}).^2)-1/3.*((ro_M2_3-ri_8{1,1}).^3))./T);
% end
% ro_M2_3= 0.974+0.067*XCa_M2-0.051*XAl6
% E_3=318.6+6.9*P-0.036*T
% Do1=(Mg_n_melt/XMg_M1)*exp((88750-65.644*T+7050*P-770*P^2)/(constant(2)*T))
% gamma_melt=-0.32+0.0563/(1-X_H2O_melt)+1.84*(1-X_H2O_melt)-0.58*(1-X_H2O_melt)^2
% Do=Do1*gamma_melt*(1-X_H2O_melt)*100/(100-H2O_melt)
% Di_REE_Y= Do.*exp(fpe.*E_3.*(ro_M2_3.*0.5.*((ro_M2_3-ri_CN_8{1,1}(1,1:15)).^2)-1/3.*((ro_M2_3-ri_CN_8{1,1}(1,1:15)).^3))./T);
% REE_plot= Do.*exp(fpe.*E_3.*(ro_M2_3.*0.5.*((ro_M2_3-ri_8{1,1}).^2)-1/3.*((ro_M2_3-ri_8{1,1}).^3))./T);
%% Transitive metals Sc
%Z0 = charge of melt 
fpe_Sc= exp((255.646-0.149*T+4.233*(P^2)-2.280*((3-z0)^2))/(constant(2)*T));
% fpe_Sc= exp((255646-149*T+4233*P^2-2280*((3-z0)^2))/(constant(2)*T));
L= ((XK+XNa)*(XSi^2)+2*(1-(XK+XNa))*XAl6*XSi*XAl4)+(((XAl4^2)*(1-(XK+XNa)))...
    +2*(XK+XNa)*XAl4*XSi+(1-(XK+XNa))*XSi)*exp(-28.000/(constant(2)*T));
DSc=L*fpe_Sc;
Di_REE_Y_Sc=[Di_REE_Y DSc];
%% LILE prediction
% 1+ cations
% r0 1+
    ro_1=(1.066-0.104.*XAl6-0.212.*XMg_M2)+0.20;
%     ro_1= 0.974+0.067*XCa-0.051*XAl6+0.12
    % E1+
    E_1=1/3.*E_3;
    % D0 1+ = DNa 
    % Di 1+
    Di_LILE_1= DNa.*exp(fpe.*E_1.*(ro_1.*0.5.*(1.18^2-ri_CN_8{1,2}.^2)-1/3.*(1.18^3-ri_CN_8{1,2}.^3))./T);
    LILE1_plot= DNa.*exp(fpe.*E_1.*(ro_1.*0.5.*(1.18^2-ri_8{1,2}.^2)-1/3.*(1.18^3-ri_8{1,2}.^3))./T);
%     Obs_data=[0.0018 0.03; 0.0018 0.04];
% if T== 1663
% figure(1)
% semilogy(ri_CN_8{1,2},Di_LILE_1,'bo','markersize', 10)
% hold on
% semilogy(ri_8{1,2},LILE1_plot,'b-')
% semilogy(ri_CN_8{1,2}(2:3),Obs_data(2,:),'rsq','markersize', 10)
% 
% end
% if T== 1613
% figure(1)
% semilogy(ri_CN_8{1,2},Di_LILE_1,'ko','markersize', 10)
% hold on
% semilogy(ri_8{1,2},LILE1_plot,'k-')
% semilogy(ri_CN_8{1,2}(2:3),Obs_data(1,:),'gsq','markersize', 10)
% end
% se\f
% 2+ cations
    % r0 2+
    ro_2= (1.066-0.104.*XAl6-0.212.*XMg_M2)+0.06;
    % E2+
    E_2 = 2/3.*E_3;
    % D0 2+ = DCa 
    % Di 2+
    Di_LILE_2= DCa.*exp(fpe.*E_2.*(ro_2.*0.5.*(1.12^2-ri_CN_8{1,3}.^2)-1/3.*(1.12^3-ri_CN_8{1,3}.^3))./T);
    LILE2_plot= DCa.*exp(fpe.*E_2.*(ro_2.*0.5.*(1.12^2-ri_8{1,3}.^2)-1/3.*(1.12^3-ri_8{1,3}.^3))./T);

%% HFSE prediction 
% 4+ cations
%r0_4+
ro_M1_4=(0.64-0.008.*P+0.071.*XAl6)+0.004;
%E_4
E_4= (10473-5.09.*T-201.54.*P+14633.*XAl4)+331;

%D_HFSE    
D_Zr= D_Ti.*exp(fpe.*E_4.*(ro_M1_4.*0.5.*(ri_CN_6{1,4}(1).^2-ri_CN_6{1,4}(3).^2)-1/3.*(ri_CN_6{1,4}(1).^3-ri_CN_6{1,4}(3).^3))./T);
D_Hf= D_Ti.*exp(fpe.*E_4.*(ro_M1_4.*0.5.*(ri_CN_6{1,4}(1).^2-ri_CN_6{1,4}(2).^2)-1/3.*(ri_CN_6{1,4}(1).^3-ri_CN_6{1,4}(2).^3))./T);

D_HFSE_4=[D_Ti D_Hf D_Zr];
HFSE4_plot= D_Ti.*exp(fpe.*E_4.*(ro_M1_4.*0.5.*(ri_CN_6{1,4}(1).^2-ri_6{1,4}.^2)-1/3.*(ri_CN_6{1,4}(1).^3-ri_6{1,4}.^3))./T);

% 5+ cations
D_Ta= exp(-2.127+3.769.*XAl4);
D_Nb=0.003+0.292.*D_Ta;
D_HFSE_5=[D_Ta D_Nb];
%% Actinides
% r04+=r03+
% E4+
E_4_M2=4/3.*E_3;
gamma_Mg=exp(7500/(constant(2).*T));
gamma_Th=exp(fpe.*E_4_M2.*((ro_M2_3.*0.5.*((ri_CN_8{1,4}(5)-ro_M2_3).^2))+1/3.*((ri_CN_8{1,4}(5)-ro_M2_3).^3))./T);    
D_Th=X_Mg_melt./(gamma_Mg.*gamma_Th.*XMg_M1).*exp((214.79-0.757.*T+16.42.*P-1.5.*P.^2)./(constant(2).*T)) ;
D_U= D_Th.*exp(-910.17.*E_4_M2.*(ro_M2_3.*0.5.*(ri_CN_8{1,4}(5).^2-ri_CN_8{1,4}(4).^2)-1/3.*(ri_CN_8{1,4}(5).^3-ri_CN_8{1,4}(4).^3))./T);
Act_plot= D_Th.*exp(-910.17.*E_4_M2.*(ro_M2_3.*0.5.*(ri_CN_8{1,4}(5).^2-ri_8{1,4}.^2)-1/3.*(ri_CN_8{1,4}(5).^3-ri_8{1,4}.^3))./T);

D_Act= [D_U D_Th];
%% Concatenate the Di 
Di1= {Di_LILE_1 Di_LILE_2 Di_REE_Y_Sc D_HFSE_4 D_Act D_HFSE_5};
Di=cell2mat(Di1);
Di_plot={LILE1_plot LILE2_plot REE_plot HFSE4_plot Act_plot};
% Di_plot=cell2mat(Di2);
%% plot
if fig==1
cpxplot(T,P,XWo,XEn,XFs,Di,Di_plot,ri_CN_8,ri_CN_6,ri_8,ri_6,order_flag,TE_name,p,clr)
end