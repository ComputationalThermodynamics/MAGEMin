function Gamma = Gamma_estimation_NN1(T,P)
% This estimates Gamma, using a trained neural network.
% WARNING: it requires the deep learning toolbox to be installed, and 
% is mainly for debugging purpeses (we can port the results to a matlab 
% script later). 


fname = 'KLB1_UHR_nn1';

load(fname);

for i=1:10
   Gamma(i,:) = net{i}([T(:) P(:)]')';
end
Gamma(11,:) = 0;    % H2O
