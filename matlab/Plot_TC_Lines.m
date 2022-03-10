function h = Plot_TC_Lines(varargin)
% This adds lines, computed with TC, to the phase diagram
linecolor= 'r';
linewidth= 2;

PseudoSectionData = varargin{1};
if nargin>1
    UIAxes      = varargin{2};
    linecolor   = varargin{3};
    linewidth   = varargin{4};
    
    hold(UIAxes,'on')
    
else
    UIAxes=[];
    hold on
end
h=[];

if isfield(PseudoSectionData,'LinesPlot')
   LinesPlot =  PseudoSectionData.LinesPlot;
   
   for i=1:length(LinesPlot)
       
       Line = LinesPlot{i};
       
       if ~isempty(Line)
           if isempty(UIAxes)
               h(i) =  plot(Line.T,Line.P,linecolor,'linewidth',linewidth);
           else
               if ~isempty(Line.P)
                   h(i) =  plot(UIAxes,Line.T,Line.P,linecolor,'linewidth',linewidth);
               end
           end
       end
   end
   
end