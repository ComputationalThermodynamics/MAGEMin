clear all;
close all;
set(0,'defaultfigurecolor',[1 1 1])
%export_fig file.png -m2.5
dir   = '../OUTPUT/';


A_    = '__RESIDUAL_NORM.txt';  dir_A = strcat(dir,A_);
B_    = '__GAMMA.txt';          dir_B = strcat(dir,B_);
C_    = '__MIN_TIME.txt';       dir_C = strcat(dir,C_);
D_    = '__SUM_XI.txt';         dir_D = strcat(dir,D_);
E_    = '__GIBBS_SYSTEM.txt';   dir_E = strcat(dir,E_);
F_    = '__BR_RESIDUAL.txt';    dir_F = strcat(dir,F_);
G_    = '__PHASE_FRACTION.txt'; dir_G = strcat(dir,G_);
H_    = '__GRAD_MAX.txt';       dir_H = strcat(dir,H_);

A = importdata(dir_A,'\t');
B = importdata(dir_B,'\t');
C = importdata(dir_C,'\t');
D = importdata(dir_D,'\t');
E = importdata(dir_E,'\t');
F = importdata(dir_F,'\t');
% G = importdata(dir_G,'\t');
H = importdata(dir_H,'\t');

sum_ss_time = sum(C.data);

n_ite = size(A.data,1);
n_ss = size(C.data,2);
n_Ox  = 11-1; 

% FIGURE RESIDUALS
fighandle = figure('position', [100, 100, 1500, 800]);
subplot(232)
hold on
plot(A.data(2:n_ite),'k-')
plot(A.data(2:n_ite),'r.')
xlabel('# Iteration');
ylabel('Residual wt%');
title('\bf Residual evolution');
set(gca, 'YScale', 'log')

subplot(234)
m = ["spl","bi","cd","cpx","ep","g","amp","ilm","liq","mu","ol","opx","pl"];
b = [];
n = [];
for i=1:n_ss
    if sum_ss_time(i) > 0.0
        hold on
        a = plot(C.data(2:n_ite,i),'k.');
        a = plot(C.data(2:n_ite,i),'-');

        b = [b,a];
        n = [n,m(i)];
        xlabel('# Iteration');
        ylabel('Min time (ms)');
    end
end
ylim([0 1])
legend(b, n);
title('\bf Min time per SS');

subplot(235)
m = ["spl","bi","cd","cpx","ep","g","amp","ilm","liq","mu","ol","opx","pl"];
b = [];
n = [];
for i=1:n_ss
    if sum_ss_time(i) > 0.0
        hold on
        a = plot(D.data(2:n_ite,i),'k.');
        a = plot(D.data(2:n_ite,i),'-');

        b = [b,a];
        n = [n,m(i)];
        xlabel('# Iteration');
        ylabel('Sum(xi) (ms)');
    end
end
ylim([0 10])
legend(b, n);
title('\bf Sum (xi) per SS');

subplot(236)
m = ["spl","bi","cd","cpx","ep","g","amp","ilm","liq","mu","ol","opx","pl"];
b = [];
n = [];
for i=1:n_ss
    if sum_ss_time(i) > 0.0
        hold on
        a = plot(H.data(2:n_ite,i),'k.');
        a = plot(H.data(2:n_ite,i),'-');

        b = [b,a];
        n = [n,m(i)];
        xlabel('# Iteration');
        ylabel('Grad max');
    end
end
ylim([0 10])
legend(b, n);
title('\bf Grad max per SS');



subplot(231)
hold on
plot(E.data(2:n_ite),'k-')
plot(E.data(2:n_ite),'r.')
xlabel('# Iteration');
ylabel('Gibbs system');
title('\bf Gibbs system evolution');

subplot(233)
m = ["SiO2","Al2O3","CaO","MgO","FeO","K2O","Na2O","TiO2","O","Cr2O3","H2O"];
b = [];
n = [];
for i=1:n_Ox
    hold on
    a = plot(F.data(2:n_ite,i),'k.');
    a = plot(F.data(2:n_ite,i),'-');

    b = [b,a];
    n = [n,m(i)];
    xlabel('# Iteration');
    ylabel('Oxide residual');
end
legend(b, n);
title('\bf Evolution of Mass residual per Oxide');

% FIGURE GAMMA
fighandle = figure('position', [100, 100, 1000, 600]);
m = ["SiO2","Al2O3","CaO","MgO","FeO","K2O","Na2O","TiO2","O","Cr2O3","H2O"];
for i=1:n_Ox
    
    subplot(3,4,i)
    hold on
    a = plot(B.data(2:n_ite,i),'k-');
    a = plot(B.data(2:n_ite,i),'r.');
    legend(a, m(i));
    xlabel('# Iteration');
    ylabel('Gamma (G)');
end



