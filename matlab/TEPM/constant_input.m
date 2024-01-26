function [constant,ri_CN_8,ri_8,fpe,ri_CN_6,ri_6,Molar_mass,flag,order_flag,TE_name1,TE_name,clr]=constant_input
%% Constant 
pi = 3.14159265359;    
R = 8.3144598;
N_avo = 602.2140857;
fpe=-4*pi*N_avo/R;
constant=[pi R N_avo];

% Ionic radius for 8-fold coordinated lattices
REE_Sc_CN_8= [1.16,1.143,1.126,1.09,1.079,1.066,1.053,1.04,1.027,1.019,1.015,1.004,0.994,0.985,0.977,0.870];
LILE_1_CN_8= [1.74,1.61,1.51];
LILE_2_CN_8=[1.42,1.26];
HFSE_Act_4_CN_8=[0.7399,0.83,0.84,1,1.05];
HFSE_5_CN_8=[0.7295,0.74];
ri_CN_8= {REE_Sc_CN_8 LILE_1_CN_8 LILE_2_CN_8 HFSE_Act_4_CN_8 HFSE_5_CN_8};
% for plot
REE_8=linspace(0.7,1.6,101);
% REE_8=linspace(0.8,1.23,101);% for px plot
% REE_8=linspace(0.7,1.6,101);% for pl plot
LILE1_8=linspace(1.3,1.74,101);
LILE2_8=linspace(0.7,1.42,101);
HFSE4_8=linspace(0.7,1.3,101);
HFSE5_8=linspace(min(ri_CN_8{1,5}),max(ri_CN_8{1,5}),101);
ri_8= {REE_8 LILE1_8 LILE2_8 HFSE4_8 HFSE5_8};

% Ionic radius for 6-fold coordinated lattices
REE_Sc_CN_6= [1.032,1.01,0.99,0.983,0.958,0.947,0.938,0.923,0.912,0.9,0.901,0.890,0.88,0.868,0.861,0.745];
LILE_1_CN_6= [1.67,1.52,1.38];   
LILE_2_CN_6=[1.35,1.18];
HFSE_Act_4_CN_6=[0.605,0.71,0.72,0.89,0.94];
HFSE_5_CN_6=[0.6295,0.64]; 
ri_CN_6= {REE_Sc_CN_6 LILE_1_CN_6 LILE_2_CN_6 HFSE_Act_4_CN_6 HFSE_5_CN_6};
% for plot
REE_6=linspace(min(ri_CN_6{1,1}),max(ri_CN_6{1,1}),101);
LILE1_6=linspace(1.2,1.8,101);
LILE2_6=linspace(1.1,1.5,101);
HFSE4_6=linspace(min(ri_CN_6{1,4}),max(ri_CN_6{1,4}),101);
HFSE5_6=linspace(min(ri_CN_6{1,5}),max(ri_CN_6{1,5}),101);
ri_6= {REE_6 LILE1_6 LILE2_6 HFSE4_6 HFSE5_6};

Molar_name = {'SiO2','TiO2','Al2O3','FeO','MnO','MgO','CaO','Na2O','K2O','H2O'};
Molar_mass = [60.084,79.979,101.961,71.846,70.937,40.304,56.077,61.979,94.196,18.01528];

%% Rearranging trace elements for plot

TE_name1={' Cs ',' Rb ',' K ',' Ba ',' Sr ',' La ',' Ce ',' Pr ',' Nd ',' Sm ',' Eu ',' Gd ',' Tb '...
    ,' Dy ',' Y ',' Ho ',' Er ',' Tm ',' Yb ',' Lu ',' Sc ',' Ti ',' Hf ',' Zr ',' U ',' Th ',' Ta ',' Nb '};
flag=1:28;

TE_name={'Cs','Rb','Ba','Th','U','K','Ta','Nb','La','Ce','Pr','Sr','Nd',...
    'Sm','Hf','Zr','Eu','Ti','Gd','Tb','Dy','Y','Ho','Er','Tm','Yb','Lu','Sc'};
order_flag=[1 2 4 26 25 3 27 28 6 7 8 5 9 10 23 24 11 22 12 13 14 15 16 17 18 19 20 21];

%% Color plot
clr=[0 0 0; 1 0 0; 0.5 0 0.5; 0 0 1; 0.3 0.75 0.93; 0 1 0; 1 0.9 0; 0.8 0.8 0.8;...
    1 0.5 0; 0.8 0.5 0.2; 0.5 0.7 0.2; 0.8 0.6 0];
