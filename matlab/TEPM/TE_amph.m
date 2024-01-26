function [Di]= TE_amph(P,T,mc,comp,ri_CN_8,ri_CN_6,fpe,ri_8,ri_6,order_flag,TE_name,fig,TE_name1,p,clr)
%% Constant and input  
%MgO WT% melt, SiO2 WT% melt
melt_SiO2= comp(1,1);
X=sum(mc(1,:),2); %sum of melt cations
Xnf=mc(1,8)+mc(1,9)+mc(1,1)+mc(1,9);
Melt_polymerisation= 0.008458.*melt_SiO2+0.1734; % rsquared= 0.86 % fit to convert SiO2 melt in melt_polymerisation made by the unpublished work of Pierre Bouilhol from Tiepolo et al, 2007 data.
% Melt_polymerisation= Xnf/X;
DCa= comp(2,7)./comp(1,7);
%% REE prediction 
% load coeff3.mat


La=0.0135.*exp(0.0478.*melt_SiO2); Ce=0.0194.*exp(0.0578.*melt_SiO2); Pr=0.0266.*exp(0.0596.*melt_SiO2); Nd=0.0388.*exp(0.0607.*melt_SiO2); Sm=0.0417.*exp(0.0607.*melt_SiO2); Eu=0.0445.*exp(0.0607.*melt_SiO2); Gd=0.0466.*exp(0.0606.*melt_SiO2); Tb=0.0482.*exp(0.0605.*melt_SiO2); 
Dy=0.0492.*exp(0.0603.*melt_SiO2); Y=0.0495.*exp(0.0601.*melt_SiO2); Ho=0.0496.*exp(0.06.*melt_SiO2); Er=0.0495.*exp(0.0597.*melt_SiO2); Tm=0.0491.*exp(0.0593.*melt_SiO2); Yb=0.0485.*exp(0.0589.*melt_SiO2); Lu=0.0478.*exp(0.0585.*melt_SiO2);
Di_REE_Y=[La Ce Pr Nd Sm Eu Gd Tb Dy Y Ho Er Tm Yb Lu];


% Di_REE_Y=coeff3(:,1).*exp(Melt_polymerisation.*coeff3(:,2));
a=linspace(0.013,0.050,101);b=linspace(0.045,0.07,101);
REE_plot=a.*exp(b.*melt_SiO2);


% La= exp(1.21*log(DCa)-2.5);Ce=exp(1.39*log(DCa)-1.93); Pr=0.0266.*exp(0.0596.*melt_SiO2); Nd=exp(1.57*log(DCa)-1.02); Sm=exp(1.60*log(DCa)-0.46); 
% Eu=exp(1.73*log(DCa)-0.45); Gd=exp(1.5*log(DCa)-0.03); Tb=exp(1.6*log(DCa)+0.02);Dy=exp(1.6*log(DCa)+0.06); Y=exp(1.79*log(DCa)-0.16); 
% Ho=exp(1.64*log(DCa)+0.04); Er=exp(1.67*log(DCa)-0.03); Tm=exp(1.58*log(DCa)-0.03); Yb=exp(1.7*log(DCa)-0.22); Lu=exp(1.8*log(DCa)-0.36); 
% D_Sc=exp(1.67*log(DCa)+1.17);
% Di_REE_Y=[La Ce Pr Nd Sm Eu Gd Tb Dy Y Ho Er Tm Yb Lu];

%% Transitive metals Sc for Tiepolo et al.,2007
D_Sc=0.0169.*exp(9.33.*Melt_polymerisation);
Di_REE_Y_Sc=[Di_REE_Y D_Sc];
%% LILE prediction
ri_12=linspace(1.45,1.95,101);
% 1+ cations
% D_Rb= exp(-2.9008+3.0454.*P-(0.1094+4.1559.*P).*Melt_polymerisation); % rsquared= 0.51
% Model for K
r_CN_12=[1.88,1.72,1.64];
ro_1=1.52;
% E_1=92;
E_1=86.75;
% Do_1=2.85;% ro, E, Do are average data from Dape and Baker, 2000
Do_1=2.99;
Di_LILE_1=Do_1.*exp(fpe.*E_1.*(ro_1.*0.5.*((r_CN_12-ro_1).^2)+1/3.*((r_CN_12-ro_1).^3))./T); 
LILE1_plot=Do_1.*exp(fpe.*E_1.*(ro_1.*0.5.*((ri_12-ro_1).^2)+1/3.*((ri_12-ro_1).^3))./T); 

% Di_LILE_1=[D_Cs D_Rb D_K];
% 2+ cations
r_Ba_CN_12=1.61;
ro_2_A=1.52;
E_2_A=315;
Do_2_A=2.6;
D_Ba= Do_2_A.*exp(fpe.*E_2_A.*(ro_2_A.*0.5.*((r_Ba_CN_12-ro_2_A).^2)+1/3.*((r_Ba_CN_12-ro_2_A).^3))./T); 
D_Ba_plot= Do_2_A.*exp(fpe.*E_2_A.*(ro_2_A.*0.5.*((ri_12-ro_2_A).^2)+1/3.*((ri_12-ro_2_A).^3))./T);

D_Sr=exp(0.251*log(DCa)-0.998);
Di_LILE_2=[D_Ba D_Sr];

%% HFSE prediction & Transitive metals Sc Nandedkar et al., 2016
% a1=[1.67;0.251;1.36;1.82;3.20;3.09;1.14;0.78];b1=[1.17;-0.998;-0.99;-1.94;-5.8;-5.56;-1.51;-1.39];
% for p=1:size(a1)
%     Nand(p)=exp(a1(p).*log(DCa)+b1(p));
% end
% 
% Di_LILE_2=[D_Ba Nand(2)];
% Di_REE_Y_Sc=[Di_REE_Y Nand(1)];
%  D_HFSE_4=[D_Ti Nand(3) Nand(4)];
%  D_HFSE_5=[Nand(6) Nand(5)];

%% HFSE Tiepolo
% 4+ cations
% E_4= 1810;
% ro_4=0.65;
% Do_4=1.73; 
E_4= 1393.2;
ro_4=0.65;
Do_4=1.27; 
if melt_SiO2<65
D_Ti=0.0273.*exp(6.9531.*Melt_polymerisation);
D_Hf= 0.0035.*exp(8.894.*Melt_polymerisation);
D_Zr=0.001493.*exp(9.352.*Melt_polymerisation);
HFSE4_plot=[D_Ti D_Hf D_Zr]; 
else 
D_Ti= Do_4.*exp(fpe.*E_4.*(ro_4.*0.5.*((ri_CN_6{1,4}(1)-ro_4).^2)+1/3.*((ri_CN_6{1,4}(1)-ro_4).^3))./T);
D_Hf= Do_4.*exp(fpe.*E_4.*(ro_4.*0.5.*((ri_CN_6{1,4}(2)-ro_4).^2)+1/3.*((ri_CN_6{1,4}(2)-ro_4).^3))./T);
D_Zr= Do_4.*exp(fpe.*E_4.*(ro_4.*0.5.*((ri_CN_6{1,4}(3)-ro_4).^2)+1/3.*((ri_CN_6{1,4}(3)-ro_4).^3))./T);
HFSE4_plot= Do_4.*exp(fpe.*E_4.*(ro_4.*0.5.*((ri_6{1,4}-ro_4).^2)+1/3.*((ri_6{1,4}-ro_4).^3))./T);
end
D_HFSE_4=[D_Ti D_Hf D_Zr];    

% 5+ cations
E_5= 1985;
ro_5=0.73;
Do_5=2.10; % ro, Do and E parameters are the average data from the work of Dalpe and Baker, 2000 on Calcic Amphibole
if melt_SiO2<56
D_Nb= 6.5e-6.*exp(18.46.*Melt_polymerisation);%rsquared= 0.40
D_Ta=D_Nb.*(3.48e-5.*exp(16.51.*Melt_polymerisation));%rsquared= 0.37
else 
D_Nb= Do_5.*exp(fpe.*E_5.*(ro_5.*0.5.*((ri_CN_6{1,5}(2)-ro_5).^2)+1/3.*((ri_CN_6{1,5}(2)-ro_5).^3))./T);
D_Ta= Do_5.*exp(fpe.*E_5.*(ro_5.*0.5.*((ri_CN_6{1,5}(1)-ro_5).^2)+1/3.*((ri_CN_6{1,5}(1)-ro_5).^3))./T);
end
D_HFSE_5=[D_Ta D_Nb];

%% Actinides
D_U=9e-6.*exp(12.15.*Melt_polymerisation); %rsquared= 0.66
D_Th=1e-5.*exp(12.24.*Melt_polymerisation);%rsquared= 0.62
D_Act= [D_U D_Th]; 
%% Concatenate the Di 
Di1= {Di_LILE_1 Di_LILE_2 Di_REE_Y_Sc D_HFSE_4 D_Act D_HFSE_5};
Di=cell2mat(Di1);

Di_plot={REE_plot LILE1_plot D_Ba_plot D_Sr HFSE4_plot D_Act};
% Di_plot=cell2mat(Di2);
%% plot
if fig==1
amphplot(T,P,Di,Di_plot,ri_CN_8,ri_CN_6,ri_8,ri_6,order_flag,TE_name,TE_name1,p,clr)
end