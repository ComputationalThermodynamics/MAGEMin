function PhaseData = Update_Gamma_onPoints(PhaseData,TP_vec, elements, newPoints,PrPts,Computation)
% This updates the Gamma for new points based on previous calculations.
% Two different methods are available:
%       1) Average from coarser mesh, which is known if we employ AMR
%       2) Average from a certain # of closests surrounding points



% Method = {'Inherit from coarse mesh','NearestPoint','NearestNeighbour'};
% Method=Method{1};


Method = Computation.Gamma_Method;
% Method = 'Inherit from coarse mesh';
% Method = 'AMR interpolation';

switch Method
    case 'Inherit from coarse mesh'
        
        if ~isempty(PrPts)
            for i=1:length(newPoints)
                Prev = PrPts(i,:); Prev=Prev(Prev>0);
                Gamma = zeros(size(PhaseData{1}.Gamma));
                for j=1:length(Prev)
                    Gamma = Gamma + PhaseData{Prev(j)}.Gamma;
                end
                PhaseData{newPoints(i)}.Gamma = Gamma/j;
                
            end
        end
        
    case 'AMR interpolation'
        oldPoints = 1:length(PhaseData);
        rmPoints = [];
        Gamma_vec = zeros(length(PhaseData),11);
        for i=1:length(PhaseData)
            if ~isfield(PhaseData{i},'Gamma')
                rmPoints(end+1) = i;
            else
                Gamma_vec(i,:) = PhaseData{i}.Gamma;
            end
        end
        oldPoints(rmPoints)=[];
        
        
        
        TP_vec_old  = TP_vec(oldPoints,:);
        Ti          = TP_vec(newPoints,1);
        Pi          = TP_vec(newPoints,2);
        
        % update Gamma
        Gamma_new = Interpolate_AMR_grid(elements, TP_vec_old, Gamma_vec, Ti, Pi);
        Gamma_new(isnan(Gamma_new))=0;
        
        for i=1:length(newPoints)
            PhaseData{newPoints(i)}.Gamma = Gamma_new(i,:);
        end
        
        
        
        
    case 'Average surrounding points'
        
        TP_norm = max(TP_vec)-min(TP_vec);
        
        oldPoints = 1:length(PhaseData);
        rmPoints = [];
        for i=1:length(PhaseData)
            if ~isfield(PhaseData{i},'Gamma')
                rmPoints(end+1) = i;
            end
        end
        oldPoints(rmPoints)=[];
        
        
        % Set this to all newPoints
        for i=1:length(newPoints)
            id = newPoints(i);
            
            if isfield(Computation,'num_Neighbours')
                num_Neighbours  =   Computation.num_Neighbours;
            else
                num_Neighbours  =   1;
            end
            
            % Get closest points to this one
            [idx] = nearestneighbour([TP_vec(id,:)./TP_norm]', [TP_vec(oldPoints(:),:)./TP_norm]','n',num_Neighbours);
            
            %             if idx(1)==id
            %                 idx(1)  =   [];
            %             end
            %             if length(idx)==0
            %                 warning('something iffy')
            %             else
            %
            %
            for j=1:length(idx)
                Gamma_vec(:,j) = PhaseData{idx(j)}.Gamma;
            end
            Gamma = mean(Gamma_vec,2);      % average
            %             end
            PhaseData{id}.Gamma = Gamma;
            
        end
        
        
    case 'Neural Network 1'
        
        % Estimate Gamma from a NN
        Gamma = Gamma_estimation_NN1(TP_vec(newPoints,1), TP_vec(newPoints,2)); 
        
         for i=1:length(newPoints)
            id = newPoints(i);
            PhaseData{id}.Gamma = Gamma(:,i);
         end
        
    otherwise
        error('Method to compute Gamma not yet implemented')
end

