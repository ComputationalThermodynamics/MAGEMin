% This is a routine to prepare data from the PseudoSectionData structure for
% feeding to the neural network


nPoints = length(PseudoSectionData.PhaseData);
TC      = zeros(nPoints,1);
for iPoint=1:nPoints
    Data = PseudoSectionData.PhaseData{iPoint};
    if isfield(Data,'TC_recompute')
        if Data.TC_recompute.success
            TC(iPoint) = 1;
        end
    end
end

Gamma =  PseudoSectionData.Gamma';
Gamma_SiO2  =   Gamma(:,1);
Gamma_Al2O3 =   Gamma(:,2);
Gamma_CaO   =   Gamma(:,3);
Gamma_MgO   =   Gamma(:,4);
Gamma_FeO   =   Gamma(:,5);
Gamma_K2O   =   Gamma(:,6);
Gamma_Na2O  =   Gamma(:,7);
Gamma_TiO2  =   Gamma(:,8);
Gamma_O     =   Gamma(:,9);
Gamma_Cr2O3 =   Gamma(:,10);
Gamma_H2O   =   Gamma(:,11);


TP_vec = PseudoSectionData.TP_vec;


indT     =  find(TC);
TPvec_TC =  PseudoSectionData.TP_vec(indT,:);
Gibbs_TC =  PseudoSectionData.Gibbs(indT);
Gamma_TC =  PseudoSectionData.Gamma(:,indT)';



Train_NN = true;
if Train_NN

    
%     for i=1:10
    for i=6:6
        
     


        % Train the shallow NN
        x = TP_vec';
%         t = Gamma(:,i)';
        
        % Test: first fit a plane through the data and only train the NN on
        % the deviation of this plane
        DM      =   [TP_vec, ones(size(TP_vec(:,1)))]; 
        Coeff   =   DM\Gamma(:,i);                                % coefficients of the plane that has Gamma_Lin = Coeff(1)*T + Coeff(2)*P + Coeff(3);
        
        % Subtract plane from data
        t       =   Gamma(:,i)-(TP_vec(:,1)*Coeff(1) + TP_vec(:,2)*Coeff(2) + Coeff(3));
        t       = t'; 
           
        
        PlaneData(:,i) = Coeff;
        
        
%         t = Gamma';
        
        
        
        % Choose a Training Function
        % For a list of all training functions type: help nntrain
        % 'trainlm' is usually fastest.
        % 'trainbr' takes longer but may be better for challenging problems.
        % 'trainscg' uses less memory. Suitable in low memory situations.
        trainFcn = 'trainlm';  % Bayesian Regularization backpropagation.
        
        % Create a Fitting Network
        hiddenLayerSize    = [30 30];
%         hiddenLayerSize    = [3 3 3 3 3 3];
        
        
        net{i}             = fitnet(hiddenLayerSize,trainFcn);
        
        net{i}.layers{2}.transferFcn = 'tansig';
        net{i}.layers{1}.transferFcn = 'tansig';
        
%         for j=1:length(net{i}.layers)
% %             net{i}.layers{j}.transferFcn = 'logsig';
%         end
%         
        % Choose Input and Output Pre/Post-Processing Functions
        % For a list of all processing functions type: help nnprocess
        net{i}.input.processFcns = {'removeconstantrows','mapminmax'};
        net{i}.output.processFcns = {'removeconstantrows','mapminmax'};
        
        % Setup Division of Data for Training, Validation, Testing
        % For a list of all data division functions type: help nndivision
        net{i}.divideFcn = 'dividerand';  % Divide data randomly
        net{i}.divideMode = 'sample';  % Divide up every sample
        net{i}.divideParam.trainRatio = 70/100;
        net{i}.divideParam.valRatio = 15/100;
        net{i}.divideParam.testRatio = 15/100;
        
%         net{i}.trainParam.epochs=2000;
        
        % Choose a Performance Function
        % For a list of all performance functions type: help nnperformance
        net{i}.performFcn = 'mse';  % Mean Squared Error
        
        % Choose Plot Functions
        % For a list of all plot functions type: help nnplot
        net{i}.plotFcns = {'plotperform','plottrainstate','ploterrhist', ...
            'plotregression', 'plotfit'};

        % Train the Network
        [net{i},tr] = train(net{i},x,t);
        
        % Test the Network
        y               = net{i}(x);
        e               = gsubtract(t,y);
        performance(i)  = perform(net{i},t,y)
        
      
        
        if 1==1
            figure(111)
            subplot(121)
            plot3(x(1,:), x(2,:), t,'.',x(1,:), x(2,:), y,'ro'),title('original data vs NN prediction')
            xlabel('Temp'), ylabel('Pres'), zlabel('Gamma K2O - linear trend')
            grid on
            
            subplot(122)
            plot3(x(1,:), x(2,:), t-y,'k.'), title('difference')
            xlabel('Temp'), ylabel('Pres'),
            grid on
        end
        
    end
    
    
    
    
    
end