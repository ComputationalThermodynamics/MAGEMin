function tf = ishg2()
%ISHG2 Determine if Handle Graphics version 2 (HG2) engine is supported.
%
% Usage:
%
%   TF = ISHG2()
%
% Outputs:
%
%   TF <1x1 logical>
%     - Indicates if HG2 is supported by version of Matlab currently in use
%     - HG2 was introduced in Matlab version 8.4 (R2014b)
%     - True if HG2 is supported, false otherwise
%
% References:
%
%   https://www.mathworks.com/matlabcentral/answers/136834-determine-if-using-hg2#answer_187810

    tf = [100, 1] * sscanf(version, '%d.', 2) >= 804;

end
