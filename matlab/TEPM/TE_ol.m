function [Di]= TE_ol(T,P,ri_CN_8,comp,mc,fpe,ri_8,fig,clr)
%% Constant and input  
MgO= comp(1,6);
melt_MgO=MgO;
SiO2= comp(1,1);
melt_SiO2=SiO2;
Mg_n= 100*(comp(1,6)./(comp(1,6)+comp(1,4)));
    
%% Structural formula  
% fudge factor 
Factor_ol= 8./(4.*mc(2,1)+4.*mc(2,2)+6.*mc(2,3)+2.*mc(2,4)+2.*mc(2,5)+2.*mc(2,6)+2.*mc(2,7)+2.*mc(2,8)+2.*mc(2,9));
XFe = Factor_ol*mc(2,4);
XMg = Factor_ol*mc(2,6);
X_Fo= XMg./(XMg+XFe);

XAl = 1/117.3*(0.22+0.06*P-0.0047*X_Fo+1.54*mc(1,3)); 

DMg= comp(2,6)/comp(1,6);
%% REE prediction + Sc
% r0 3+
% ro_3=0.725;
ro_3=0.72;
% E3+
% E_3=442;
E_3=426;
% D0 3+  
% Do=exp(-0.67-0.17.*P+117.30.*XAl-0.0147.*X_Fo);
Do=exp(-0.45-0.11.*P+1.54*mc(1,3)-0.0194.*X_Fo);

% Di 3+
% Di_REE_Y_Sc= Do.*exp(fpe.*E_3.*(ro_3.*0.5.*((ri_CN_6{1,1}-ro_3).^2)+1/3.*((ri_CN_6{1,1}-ro_3).^3))./T) % use of ionic radius 8-fold coordinated following the Yao et al, 2012 recomendations
Di_REE_Y_Sc= Do.*exp(fpe.*E_3.*(ro_3.*0.5.*((ri_CN_8{1,1}-ro_3).^2)+1/3.*((ri_CN_8{1,1}-ro_3).^3))./T); % use of ionic radius 8-fold coordinated following the Yao et al, 2012 recomendations

%% LILE prediction
% 1+ cations
    % D_Cs
    if melt_MgO>1.11
    D_Cs= exp(0.055-3.05.*log(melt_MgO)) ;
    % D_Rb
    D_Rb= exp(0.055-3.05.*log(melt_MgO));
    % D_K
    D_K=exp(0.055-3.05.*log(melt_MgO));
    else
    D_Cs= 0.035;
    % D_Rb
    D_Rb= 0.035;
    % D_K
    D_K=0.035;  
    end
    Di_LILE_1= [D_Cs D_Rb D_K];
    
% 2+ cations
    % r0 2+
%     ro_2= 0.72;
    ro_2= 0.89;

    % E2+
%     E_2 = 2/3.*E_3;
    E_2=240;
    % D0 2+ = DMg
    % Di 2+
    Di_LILE_2= DMg.*exp(fpe.*E_2.*(ro_2.*0.5.*(ro_2^2-ri_CN_8{1,3}.^2)+1/3.*(ri_CN_8{1,3}.^3-ro_2^3))./T);
%     Di_LILE_2= DMg.*exp(fpe.*E_2.*(ro_2.*0.5.*(ro_2^2-ri_CN_6{1,3}.^2)+1/3.*(ri_CN_6{1,3}.^3-ro_2^3))./T)

%     lile=disp(Di_LILE_2)
        
%% HFSE prediction 
% 4+ cations
    if melt_SiO2>57.1
        D_Ti=exp(-9.72+0.11.*melt_SiO2);
        elseif melt_MgO<0.493
        D_Ti=0.29-0.52.*melt_MgO;
        elseif Mg_n<8.35
        D_Ti=0.32-0.034.*Mg_n; 
    else
        D_Ti= 0.0302; %%%%%%%%%%%%% PIcked up from sparsed data of Bedard, 2005 ===> not accurate
    end
    if melt_MgO>5.2
    D_Hf=exp(-2.98-0.76.*log(melt_MgO));
    D_Zr= exp(-0.25-1.88.*log(melt_MgO));       
    else
    D_Hf= 0.04;   
    D_Zr= 0.036; 
    end    
    D_HFSE_4= [D_Ti D_Hf D_Zr];
% 5+ cations
    if melt_MgO>3.79
    D_Ta=exp(-1.5*log(melt_MgO)); 
    else
    D_Ta=0.126;
    end
    if melt_MgO>4.4
    D_Nb=exp(0.28-3.29.*log(melt_MgO)); 
    else
    D_Nb=0.01025;
    end
    D_HFSE_5=[D_Ta D_Nb];       
%% Actinides
if melt_MgO>5.91
D_U=exp(7.08-5.67.*log(melt_MgO)) ;
else
D_U=0.048;
end
if melt_MgO>4.9
D_Th=exp(3.8-4.22.*log(melt_MgO)); 
else
D_Th=0.0542;
end
D_Act= [D_U D_Th];     
%% Concatenate the Di 
Di= { Di_LILE_1 Di_LILE_2 Di_REE_Y_Sc D_HFSE_4 D_Act D_HFSE_5};
Di=cell2mat(Di);
