function [Line] = Extract_Reactions_MAGEMin(LineData, PhaseBound)
% Experimental routine to compute the reaction line
%

% % Extract boundary of domain from the points of the field.
% %  NOTE:this requires a continuous field, not two at different locations!
% ind       = find(PseudoSectionData.NumAssemblage==Data.NumAssemblage);
% TP_points = PseudoSectionData.TP_vec(ind,:);
% id        = boundary(TP_points(:,1),TP_points(:,2));
% TP_bound  = TP_points(id,:);
% 
% simplify = false;
% if simplify
%     TP_poly   = polyshape(TP_bound);
%     TP_bound  = TP_poly.Vertices;
% end

TP_points =  PhaseBound.Field1_TP;
id        = boundary(TP_points(:,1),TP_points(:,2));
TP_bound  = TP_points(id,:);


if ~LineData.CalcTatP
    dT          = diff(LineData.T_val)/(LineData.numPoints-1)
    T_line      = LineData.T_val(1):dT:LineData.T_val(2);
    idx         = nearestneighbour(T_line,TP_points(:,1)');
    P_line      = TP_points(idx,2)';
    dmaxA       = [0 -1];

else
    dP          = diff(LineData.P_val)/(LineData.numPoints-1)
    P_line      = LineData.P_val(1):dP:LineData.P_val(2);
    idx         = nearestneighbour(P_line,TP_points(:,2)');
    T_line      = TP_points(idx,1)';
    dmaxA       = [50 0];

end
TP_bound = [T_line(:), P_line(:)];
Tnorm    = polynorm(TP_bound);

Tol = [1 0.01];
for i=1:length(TP_bound)
    result = Bisection_T( TP_bound(i,:)-dmaxA/10, dmaxA, Tol, PhaseBound.Field1_StablePhases );
    
    if ~isnan(result(1))
        TP_bound(i,:) = result;
    else
        result = Bisection_T( TP_bound(i,:)-dmaxA/2, dmaxA, Tol, PhaseBound.Field1_StablePhases );
        TP_bound(i,:) = result;
    end
    
    
end

Line.T = TP_bound(:,1);
Line.P = TP_bound(:,2);



    % Bisection vs. T
    function [ x] = Bisection_T( A, dmaxA, Tol, Field1 )
        
        B    = A + dmaxA;
        dA   = (B-A)/20;
        
        % Scan initial start point to ensure that our starting point is
        % within Field1
        fA = -1;
        N  = 0
        while fA<1 & N<20
            fA = perform_calculation(A(2),A(1), Field1);
            N  = N+1
            if fA<1
                A = A + dA;
            end
        end
        
        if fA==1
            C = A;
            % Found a feasible starting point
            N  = 0;
            while any(abs(B-A) > Tol)
                C  = (A+B)/2 ;
                fC = perform_calculation(C(2),C(1), Field1);
                if fA*fC>0
                    A=C;
                else
                    B=C;
                end
                N = N+1;
            end
            x=C;
            
        else
            % Did not find a point that is in the required field
            x = NaN*A;
        end
        
    end




    % Calculation routine
    function [PhasePresent, Variance, StableSol] = perform_calculation(P,T, Field1 )
        
        
        % Perform computation
        command = ['./MAGEMin --test=0 --Pres=',num2str(P),' --Temp=',num2str(T)];
        system(command);
        
        % Read data
        [PhaseData] = ReadData_MAGEMin();
        
        % Stable Solutions
        StableSol = PhaseData{1}.StableSolutions;
        Variance  = 10 + 2 - length(StableSol);
        
        PhasePresent = -1;
        if isempty(setxor(StableSol, Field1))
            PhasePresent = 1;
        end
        
    end



    function [T] = polynorm(xy)
        if xy(end,:)~=xy(1,:)
           xy = [xy; xy(1,:)];
        end
        
        % First derivative
        dsxy = diff(xy);
        ds = sqrt( sum(dsxy.^2,2));
        T  = dsxy./repmat(ds,[1 2]);        % normal
        
        % The normal is defined inbetween two points; average it here
        T = [T(end,:); T; T(1,:)];
        T = (T(2:end,:) + T(1:end-1,:))/2;
        
        T = [T(:,2), -T(:,1)];  %v direction normal to point
        
    end


end






