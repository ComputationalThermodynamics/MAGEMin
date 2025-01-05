function h = copyUIAxes(varargin)
% use COPYUIAXES to copy the content of a UI axis to a new figure.
% COPYUIAXES receives the handle to a UI axes and copies all of its
% children and most of its properties to a new axis. If the UI axis
% has a legend, the legend is copied too (but see tech notes below). 
% The handle to a colobar may also be included to copy the colorbar.
% Requires Matlab r2016a or later. 
%
% COPYUIAXES(uiax) creates a new figure and axis in default positions
% and copies the content of the UI axes onto the new axes. 'uiax' is
% required to be a UIAxis handle created by uiaxes or AppDesigner.
% If the axis contains a non-empty legend property, the legend is 
% copied too (but see tech notes below).
%
% COPYUIAXES(uiax,destination) copies the content of the UI axis to 
% destination object which can either be a figure or axes handle. If
% the handle is a figure the axes will be created within the figure. 
%
% h = COPYUIAXES(___) returns a scalar structure containing all graphics
% handles that were produced.  
%
% COPYUIAXES(___, 'legend', h) specifies the legend handle (h) to 
% copy to the new destination. If the lengend handle is provided 
% here, it overrides any legend detected from within the axis handle.
%
% COPYUIAXES(___, 'colorbar', h) specifies the colorbar handle (h)  
% to copy to the new destination (but see tech notes below). 
%
% COPYUIAXES(___', 'copyPosition', true) will set the destination figure 
% size equal to the source figure size and will create axes in the same 
% position as the source axes.  The legend and colorbar positions will
% also be copied (but see tech notes below). If 'destination' is an 
% axes handle, 'copyPosition' is ignored.
%
% COPYUIAXES(___, 'listIgnoredProps', TF) when TF is true, a table
% of graphic objects appears in the command window listing the properties
% of each object that were ignored. Some properties are not editable
% while others are intentionally ignored.  
%
% Examples (for r2016b or later)
%     fh = uifigure();
%     uiax = uiaxes(fh);
%     hold(uiax,'on')
%     grid(uiax,'on')
%     x = linspace(0,3*pi,200);
%     y = cos(x) + rand(1,200);
%     ph = scatter(uiax,x,y,25,linspace(1,10,200),'filled');
%     lh = plot(uiax, x, smooth(x,y), 'k-'); 
%     title(uiax,'copyUIaxes Demo','FontSize',18)
%     xlabel(uiax,'x axis')
%     ylabel(uiax,'y axis')
%     zlabel(uiax,'z axis')
%     cb = colorbar(uiax); 
%     cb.Ticks = -2:2:12;
%     caxis(uiax,[-2,12])
%     ylabel(cb,'cb ylabel')
%     lh = legend(uiax,[ph,lh],{'Raw','lowess'},'Location','NorthEast'); 
% 
%     % Ex 1: Copy axis to another figure
%     h = copyUIAxes(uiax);
%     
%     % Ex 2: specify subplot axis as destination
%     fNew = figure();
%     axh = subplot(2,2,1,'Parent',fNew);
%     h = copyUIAxes(uiax,axh);
% 
%     % Ex 3: Copy colorbar and specify legend handle inputs
%     h = copyUIAxes(uiax, 'colorbar', cb, 'legend', lh);
%     
%     % Ex 4: See which properties were not copied
%     h = copyUIAxes(uiax, 'colorbar', cb, 'listIgnoredProps', true);
%
%     % Ex 5: Copy position, too. 
%     h = copyUIAxes(uiax, 'colorbar', cb, 'copyPosition', true);
%    
% Source: <a href = "https://www.mathworks.com/matlabcentral/fileexchange/73103-copyuiaxes">copyUIAxes</a>
% Contact author: <a href = "https://www.mathworks.com/matlabcentral/profile/authors/3753776-adam-danz">Adam Danz</a> 
% Copyright (c) 2020  All rights reserved
% To follow discussion on matlab central see note [1].
%% Version history
% 191000  vs 1.0.0 first upload to FEX
% 200209  vs 1.1.0 return fig handle; improved error handling; improved listIgnoredProps text.
% 200319  vs 1.2.0 Now works with categoricalHistograms; if set() results in warnings, a final 
%       warning msg appears; added some commented-out troubleshooting code; new newFieldOrder();
%       Added class to addRequired(uiax) error msg; drawnow() updates figs before copy; axes
%       made visible at the end only if UIAxes are visible;
% 200529 vs 2.0.0 Updated for Matlab r2020a which contained several new properties for uiaxes.
%       Now works with tiled layout plots. 
% 200608 vs 3.0.0 Properties are set outside of the debug loop (accidentally left over in 2.0.0).
%       Added copyPosition option. New conditions to set fig visibility at the end, removed ax visibility 
%       conditions. All internally generated errors will have 'CUIAx' error ID. Now copying legend title.  
%       Using ancestor() to get figure handle instead of Parent prop. Colormap copied for <r2018a.
% 201113 vs 4.0.0 Added support to yyaxis() plots & subtitle(). PositionConstraint added to excludeList
%       to avoid tiledlayout warning & other errors (see [13]). Colorbar's Layout property not copied (see [25]).
% 210317 vs 4.0.1 Fixed Matlab version retrieval problem with r2021a; started using verLessThan.  
%% Tech notes
% Mimicking copyobj() for UIaxes is not trivial.  Many of the axes properties
% and sub-handles are hidden, obscured, read-only, or vary between Matlab releases 
% so some properties are not copied.  See listIgnoredProps to list those properties.  
% 
% Copying the legend is not currently supported so we have to create a new one.
% This is also not trivial since we do not have a 1:1 match between line 
% object handles and their DisplayName values from within the legend handle. 
% So we must re-create the legend and search for the DisplayName values in 
% the axes' children.  This may result in inconsistencies with the original 
% legend.  For example, if there is more than 1 object with the same DisplayName
% value then the legend will pair with the first one which may not be the same 
% object in the original legend.  In these cases it is better to delete the new
% legend and just recreate it after copying axes.  Similar problems exist with 
% copying the colorbar.  The colorbar cannot be copied so we must recreated it 
% and do our best to copy the properties. 
%
% Copying the axes position within the figure was a giant task.  For example, an 
% axes in AppDesigner could be a child of Tab > TabGroup > Panel > GridLayout > Figure
% and each of those objects' positions are relative to their parent, potentially
% with different units. Furthermore, when the legend location is 'none' or when the 
% colorbar location is 'manual', their positions are relative to the UIAxes in the 
% source figure but those object positions are relative to the figure in regular axes. 
% So a lot of work is done to reproduce the desintation positions within the figure.  
%
% See footnotes for additional details regarding certain parts of the code. 
%% input parser
persistent vs
p = inputParser();
p.FunctionName = mfilename;
addRequired(p, 'uiax', @(h)assert(isa(h,'matlab.ui.control.UIAxes'),'CUIAx:SourceType',...
    sprintf(['copyUIAxes() only copies from UIAxes. You''re attemping to copy a %s ',...
    'object. Use copyobj() to copy objects from other axis types.'], class(h))))
addOptional(p, 'destination', [], @(h)assert(isempty(h)||isgraphics(h,'axes')||isgraphics(h,'figure'),... % see notes [9,10]
    'CUIAx:DestType','Destination must be an axis handle, figure handle, or empty.'));  
addParameter(p, 'legend', [], @(h)assert(isempty(h)||isgraphics(h,'legend'),'CUIAx:Leg','It must be a valid legend handle or empty.'));
addParameter(p, 'colorbar', [], @(h)assert(isempty(h)||isgraphics(h,'colorbar'),'CUIAx:CB','It must be a valid colorbar handle or empty.'));
addParameter(p, 'listIgnoredProps', false, @(h)assert(islogical(h),'CUIAx:IgnPrp','It must be a logical true|false.'));
addParameter(p, 'copyPosition', false, @(h)assert(islogical(h),'CUIAx:CopyPos','It must be a logical true|false.'));
% Undocumented: yyNextAx is true when copying the 2nd yyaxis, detected internally.  Otherwise, false.
addParameter(p, 'yyNextAx', false, @(h)assert(islogical(h),'CUIAx:yyNextAx','It must be a logical true|false.')); % undocumented
parse(p,varargin{:})
% Store matlab version comparisons
if isempty(vs)
    vs.lessThan2018a = verLessThan('Matlab','9.4'); 
    vs.lessThan2020a = verLessThan('Matlab','9.8');
end
%% Produce figure and axes if needed
if isempty(p.Results.destination)
    h.figure = figure('Visible','off'); 
    h.axes = axes(h.figure); 
    destinationCreatedInternally = [true, true]; % [fig, ax]
elseif isgraphics(p.Results.destination, 'figure')
    h.figure = p.Results.destination; 
    h.axes = axes(h.figure); 
    destinationCreatedInternally = [false, true]; % [fig, ax] 
else % assume it's an axis handle; see notes [9,10]
    h.figure = ancestor(p.Results.destination,'figure'); 
    h.axes = p.Results.destination; 
    destinationCreatedInternally = [false, false]; % [fig, ax]
end
drawnow(); pause(0.05); % To avoid lagging graphics problems
% Get fig handle of source axes
sourceFig = ancestor(p.Results.uiax,'figure'); 
%% Adjust axes with CategoricalRulers [16]
% As far as I'm aware, this only affects histogram() plots with 
% categorical data along the x or y axis. I've tested pie() and
% scatter plots with categorical data and there are no problems. 
% histogram2() and other hist plots do not allow categorical inputs (as of r2019b). 
% This section must come before assigning data to the axis. 
% Get handle(s) to categorical histogram in UIAxes.
catHistHandle = findall(p.Results.uiax.Children, 'Type', 'categoricalhistogram'); 
% Detect if X and Y axis are CategoricalRuler
if ~isempty(catHistHandle) && strcmpi(class(p.Results.uiax.XAxis), 'matlab.graphics.axis.decorator.CategoricalRuler')
    matlab.graphics.internal.configureAxes(h.axes, catHistHandle(1).Data, 1); % set xaxis to categorical [16]
end
if ~isempty(catHistHandle) && strcmpi(class(p.Results.uiax.YAxis), 'matlab.graphics.axis.decorator.CategoricalRuler')
    matlab.graphics.internal.configureAxes(h.axes, 1, catHistHandle(1).Data); % set yaxis to categorical [16]
end
%% Detect yyaxis plots
% yyaxis plots are detected by the number of yaxis handles and are sent through this function
% for each axis. YAxisLocation is read-only in yyaxis plots and is added to the exclusion list. 
isyyaxis = numel(p.Results.uiax.YAxis) == 2; 
if isyyaxis
    yyaxis(h.axes, p.Results.uiax.YAxisLocation)
    yyexcludeList = 'YAxisLocation'; 
else
    yyexcludeList = ''; 
end
%% Set destination axes location
% This is only done if 1) requested by user, 2) axes were produced internally.
posCopied = setPosition(p, h, destinationCreatedInternally, true, [], sourceFig);
%% Copy axis children and (most) properties
% Copy all children from UIAxes to new axis
copyobj(p.Results.uiax.Children, h.axes)
% Anonymous func used to move selected fields to the end of a structure 
% INPUTS: 
%   S: input structure
%   propList: nx1 or 1xn cell array of chars listing fields of S to put at the end. If field is missing, ignored.
% OUTPUT: The structure S with the reordered fields. 
% EXAMPLE; Snew = newFieldOrder(uiaxGoodParams, {'XLim','YLim','ZLim'})
newFieldOrder = @(S, propList)orderfields(S,[find(~ismember(fields(S),propList)); find(ismember(fields(S),propList))]); 
% Choose selected properties to copy
copyLaterList = {'Title';'Subtitle';'XLabel';'YLabel';'ZLabel'}; % see note [13]
excludeList = [{'Parent';'Children';'XAxis';'YAxis';'ZAxis';'Position';'OuterPosition';'InnerPosition';...
    'PositionConstraint';'Units';'Toolbar';'Layout'};yyexcludeList;copyLaterList]; 
[uiaxGoodParams, uiaxbadProps] = getGoodParams(p.Results.uiax, h.axes, excludeList, 'axis');  % see note [2]
uiaxbadProps(ismember(uiaxbadProps.axis,copyLaterList),:) = []; %rm the fields that we'll copy later.
% move axis limit fields last or they may not be set properly; see note [19]
uiaxGoodParams = newFieldOrder(uiaxGoodParams, {'XLim','YLim','ZLim'}); 
% clear the last warning in order to detect new warnings.
[lastWarnMsg, lastWarnID] = lastwarn(); 
lastwarn('') % Clear last warning so we can detect if one appears. 
% set properties 
set(h.axes, uiaxGoodParams)
% Set colormap prior to r2018a (see note [24])
if vs.lessThan2018a && isprop(p.Results.uiax, 'Colormap')
    colormap(h.axes, p.Results.uiax.Colormap); 
    uiaxbadProps(strcmpi(uiaxbadProps.axis,'Colormap'),:) = []; 
end
% For trouble shooting, loop through each property rather then setting them all at once.
% To remove fields for testing purposes, uiaxGoodParams = rmfield(uiaxGoodParams, {'XTick', 'XLim'}); 
%     warning('TROUBLESHOOTING MODE')
%     axProps = fieldnames(uiaxGoodParams);
%     for i = 1:numel(axProps)
%         fprintf('Property #%d: %s (%d/%d)\n', i, axProps{i}, i, numel(axProps));
%         set(h.axes, axProps{i}, uiaxGoodParams.(axProps{i}))
%     end
% Copy properties of special-cases listed in copyLaterList
% The order of commands below overcome hysteresis problems that occur when set() comes too soon after get(). 
h.axesTitle = title(h.axes, p.Results.uiax.Title.String);
if ~vs.lessThan2020a % subtitles added in r2020b
    h.axesSubtitle = subtitle(h.axes, p.Results.uiax.Subtitle.String);
    [axSubtGoodParams, axSubtbadProps] = getGoodParams(p.Results.uiax.Subtitle, h.axesSubtitle, {'Parent'; 'Position'}, 'axesSubtitle');
elseif p.Results.listIgnoredProps
    axSubtbadProps = table(); 
end
h.axesXLabel = xlabel(h.axes, p.Results.uiax.XLabel.String);
h.axesYLabel = ylabel(h.axes, p.Results.uiax.YLabel.String);
h.axesZLabel = zlabel(h.axes, p.Results.uiax.ZLabel.String);
[axTtlGoodParams, axTtlbadProps] = getGoodParams(p.Results.uiax.Title, h.axesTitle, {'Parent'; 'Position'}, 'axesTitle');
[axXLGoodParams, axXLbadProps] = getGoodParams(p.Results.uiax.XLabel, h.axesXLabel, {'Parent'; 'Position'}, 'axesXLabel');
[axYLGoodParams, axYLbadProps] = getGoodParams(p.Results.uiax.YLabel, h.axesYLabel, {'Parent'; 'Position'}, 'axesYLabel');
[axZLGoodParams, axZLbadProps] = getGoodParams(p.Results.uiax.ZLabel, h.axesZLabel, {'Parent'; 'Position'}, 'axesZLabel');
set(h.axesTitle, axTtlGoodParams)
set(h.axesXLabel, axXLGoodParams)
set(h.axesYLabel, axYLGoodParams)
set(h.axesZLabel, axZLGoodParams)
if ~vs.lessThan2020a % subtitles added in r2020b
    set(h.axesSubtitle, axSubtGoodParams)
end
%% If this is a yyaxis, send the other axis through
if isyyaxis
     % Switch YAxisLocation side
    if strcmpi(p.Results.uiax.YAxisLocation, 'left')
        nextside = 'right';
    else
        nextside = 'left';
    end
    yyaxis(p.Results.uiax, nextside)
    yyaxis(h.axes, nextside)
    if ~p.Results.yyNextAx
        newChildrenYY = p.Results.uiax.Children; 
        copyUIAxes(p.Results.uiax, h.axes, 'yyNextAx', true);
    else
        % Done copying 2nd yyaxis; return to the first one again.
        return
    end
else 
    newChildrenYY = []; 
end
%% Detect legend and copy if one exists (see note [14])
if (any(strcmpi(properties(p.Results.uiax),'Legend')) && ~isempty(p.Results.uiax.Legend)) ... % see notes [6,8]
        || ~isempty(p.Results.legend) 
    % if Legend was provided, use that handle, otherwise use the detected one. 
    if ~isempty(p.Results.legend)
        legHand = p.Results.legend;
    else
        legHand = p.Results.uiax.Legend; 
    end
    % Search for objects in new axes that have matching displayNames values as legend strings (see note [11])
    newChildren = [h.axes.Children; newChildrenYY]; 
    hasDisplayName = isprop(newChildren,'DisplayName'); 
    displames = get(newChildren(hasDisplayName),'DisplayName'); 
    if isempty(displames)
        displames = ''; 
    end
    [~,legIdx] = ismember(legHand.String, displames); 
    legObjHands = newChildren(legIdx); 
    
    % Create new legend and copy selected properties
    h.legend = legend(h.axes,legObjHands,legHand.String); % see note [7]
    excludeList = {'String'; 'Parent'; 'Children'; 'Position'; 'Units'; 'UIContextMenu'; 'ContextMenu'}; 
    if strcmpi(legHand.Location, 'none') 
        excludeList{end+1} = 'Location'; % see note [20]
    end
    [legGoodParams, legbadProps] = getGoodParams(legHand, h.legend, excludeList, 'legend');  % see note [3]
    set(h.legend, legGoodParams)
    % For trouble shooting, loop through each property rather then setting them all at once.
    %     legProps = fields(legGoodParams); 
    %     for i = 1:numel(legProps)
    %         fprintf('Property #%d: %s\n', i, legProps{i});
    %         set(h.legend, legProps{i}, legGoodParams.(legProps{i}));
    %     end
    % Copy legend title; see note [13]
    h.legendTitle = title(h.legend, legHand.Title.String);
    [legTtlGoodParams, legTtlbadProps] = getGoodParams(legHand.Title, h.legendTitle, {'Parent'; 'Children'}, 'legendTitle');
    set(h.legendTitle, legTtlGoodParams)
    
    % Reset position of axes i& set position of legend if requirements are met. 
    posCopied = setPosition(p, h, destinationCreatedInternally, false, legHand, sourceFig); % see note [20]
else
    legbadProps = table(); 
    legTtlbadProps = table(); 
    legHand = []; 
end
%% Detect colorbar and copy if one exists (see note [14, 24])
% Note, as of r2019b, there is no way I know of to detect colorbar or get its handle from ui axes (email me if you know how).
if  ~isempty(p.Results.colorbar) 
    %Copy colorbar & selected properties
    h.colorbar = colorbar(h.axes); % see note [3,25]
    excludeList = {'Parent'; 'Children'; 'Position'; 'Units'; 'UIContextMenu'; 'ContextMenu';'Layout'};
    if strcmpi(p.Results.colorbar.Location, 'manual')
        excludeList{end+1} = 'Location'; % see note [20]
    end
    [cbGoodParams, cbbadProps] = getGoodParams(p.Results.colorbar, h.colorbar, excludeList, 'colorbar');
    set(h.colorbar, cbGoodParams)
    
    % % For trouble shooting, loop through each property rather then setting them all at once.
    % % To remove fields for testing purposes, cbTtlGoodParams = rmfield(cbTtlGoodParams, {'XTick', 'XLim'});
    % warning('TROUBLESHOOTING MODE')
    % cbProps = fieldnames(cbGoodParams);
    % for i = 1:numel(cbProps)
    %     fprintf('cbProperty #%d: %s (%d/%d)\n', i, cbProps{i}, i, numel(cbProps));
    %     set(h.colorbar, cbProps{i}, cbGoodParams.(cbProps{i}))
    % end
    
    % Copy title
    h.colorbarTitle = title(h.colorbar, p.Results.colorbar.Title.String);
    [cbTtlGoodParams, cbTtlbadProps] = getGoodParams(p.Results.colorbar.Title, h.colorbarTitle, {'Parent'; 'Children';'Position'}, 'colorbarTitle'); 
    set(h.colorbarTitle, cbTtlGoodParams)
        
    % Copy ylabels
    h.colorbarYlabel = ylabel(h.colorbar, p.Results.colorbar.YLabel.String);
    [cbYLabGoodParams, cbYLabbadProps] = getGoodParams(p.Results.colorbar.YLabel, h.colorbarYlabel, {'Parent'; 'Children';'Position'},'colorbarYlabel'); 
    set(h.colorbarYlabel, cbYLabGoodParams)
    
    % Reset position of axes and legend and set position of colorbar if requirements are met. 
    posCopied = setPosition(p, h, destinationCreatedInternally, false, legHand, sourceFig); % see note [20]
        
else
    cbbadProps = table(); 
    cbTtlbadProps = table();
    cbYLabbadProps = table(); 
end
%% Toggle visibility 
% Destination figure visibilty is turned on under either of the following conditions. 
% 1) Figure was produced internally && source figure is visible.
% 2) Destination-figure handle (or axes handle) was provided in inputs and is visible (see note [22]).
if (destinationCreatedInternally(1) && strcmpi(sourceFig.Visible, 'on')) || ...
        (~destinationCreatedInternally(1) && strcmpi(h.figure.Visible','on'))
    h.figure.Visible = 'on'; 
end
drawnow(); pause(0.05); 
%% Show summary of properties that were not copied
if p.Results.listIgnoredProps
    % If position properties were copied, remove them here. 
    if posCopied(1) 
        uiaxbadProps(ismember(uiaxbadProps.axis,{'Position','InnerPosition','OuterPosition'}),:) = [];
    end
    if posCopied(2) 
        legbadProps(strcmpi(legbadProps.legend,'Position'),:) = [];
    end
    if posCopied(3) 
        cbbadProps(strcmpi(cbbadProps.colorbar,'Position'),:) = [];
    end
    % Put all badProps arrays into table
    badProps = {uiaxbadProps, axTtlbadProps, axSubtbadProps, axXLbadProps, axYLbadProps, ...
        axZLbadProps, legbadProps, legTtlbadProps, cbbadProps, cbTtlbadProps, cbYLabbadProps};
    numProps = cellfun(@numel,badProps);
    badProps(numProps==0) = []; 
    badParams = cellfun(@(c)[c;repmat({' '},max(numProps)-numel(c),1)],badProps,'UniformOutput',false);
    badParams = [badParams{:}];
    % Display table
    if ~isempty(badParams)
        fprintf(['\n%s\nThe following properties (rows) were not copied from the following objects (columns)\n'...
            'either because they are not editable or because they were intentionally ignored.\n'...
            'See %s for details.\n'], char(42*ones(1,70)), sprintf('<a href="matlab: help(''%s'') ">help(''%s'')</a>', which(mfilename), mfilename))
        disp(badParams)
        fprintf('%s\n\n',char(42*ones(1,70)))
    else
        fprintf('\n %s\n All properties were copied.\n See %s for details.\n%s\n\n ',char(42*ones(1,70)), ...
            sprintf('<a href="matlab: help(''%s'') ">help(''%s'')</a>', which(mfilename), mfilename), char(42*ones(1,70)))
    end
end
%% Check if a warning was thrown
pause(0.05) % set() needs time to complete prior to determining whether it caused an error.
if ~isempty(lastwarn())
    % A warning was thrown; tell user to inspect copied figure    
    msg = [newline, 'Based on the warning message(s) shown above, some properties may not have '...
        'been copied properly in copyUIAxes.  This can happen when a property change triggers changes '...
        'to other properties or due to the order in which the properties were set. Inspect ' ...
        'the copied axis closely and set neglected properties after running copyUIAxes().', newline];
    defaultMode = warning('query', 'backtrace');  %store default
    warning off backtrace                         %turn off backtrace
    fprintf('\n')
    warning(msg)
    warning(defaultMode.state, 'backtrace')       %turn back on default
    
else
    % A warning was not thrown; return the lastwarning.
    lastwarn(lastWarnMsg, lastWarnID)
   
end
%% Local functions
function [goodParams, badParams] = getGoodParams(hObjSource, hObjDest, badList, objName)
% The goal is to create a structure of parameter values to be changed in the destination object.  
% INPUTS
% hObjSource: handle to object being copied (ie, the old legend)
% hObjDest: handle to the new, existing copied object (ie, the new legend)
% badList: a nx1 cell array of strings that match property names (case sensitive).  These
%   properties will not be copied and will be removed from the structure.  
% objName: a char array identifying the object in hObhDest (for badParams output headers)
% OUTPUTS
% goodParams: a strcture of parameter names (fields) and their values that will be used
%   to update the parameters in hObjDest. To apply the parameter values, set(hObjDest, goodParams).
% badParams: a mx1 table of strings listing the properties that will not be copied either because
%   they were listed in badList or because they are not editable. Headers are defined by objName.
% List all parameters of the source object
params = get(hObjSource);
paramNames = fieldnames(hObjSource);
% Get list of editable params in destination object
editableParams = fieldnames(set(hObjDest));
% Remove the params that aren't editable or are unwanted in the destination obj
badParams = paramNames(~ismember(paramNames, editableParams));
badParams = unique([badParams; badList(ismember(badList,paramNames))]); % see note [4]
goodParams = rmfield(params,unique(badParams)); % see note [5]
badParams = table(sort(badParams),'VariableNames', {objName}); % see note [12]
% Change the order of properties (for troubleshooting purposes only)
% Properties in col1, if they exist in goodParams, should be set
% prior to their paired property in col2.  
%     AbeforeB = {
%         'ColorOrderIndex',           'ColorOrder'           % [15]
%         'LineStyleOrderIndex',       'LineStyleOrder'       % [15]
%         }; 
% 
%     params = fields(goodParams);
%     for i = 1:size(AbeforeB,1)
%         [~, paramIdx, ABidx] = intersect(params, AbeforeB(i,:));
%         if ABidx(1) > ABidx(2)
%             % Switch positions
%             paramOrder = (1:numel(params)).';
%             paramOrder(paramIdx) = flipud(paramIdx);
%             goodParams = orderfields(goodParams,paramOrder);
%         end
%     end
function posCopied = setPosition(p, h, destinationCreatedInternally, firstCall, legHand, sourceFig) % see note [20]
% By default, position properties are not copied unless copyPosition is set to true. This function
% sets the position of destination figure if the figure was produced internally and sets the position
% of axes if the axes was produced internally (both determined by destinationCreatedInternally).  If 
% there is a legend and a colorbar and those object 'Location' values are set to 'none' or 'manual',
% then their position properties will also be copied. Note that the position property has no effect
% on tiledlayout axes but since the axes this applies to will always be produced internally, we dont'
% need to worry about that. If copyPosition is true but the destination axes were provided as input, 
% there is no repositioning and and a warning is thrown when firstCall is true. If copyPosition is true
% and a destination figure was provided, the figure is resized to match source figure. 'p' is the input parser 
% object; 'h' is the structure of destination handles; 'legHand' is either empty of a legend handle.  
% 'sourceFig' is the handle to the source figure. 'posCopied' output is 1x3 logical vector indicating whether
% [axes, legend, colorbar] positions were copied.  
posCopied = false(1,3); 
if p.Results.copyPosition
    if firstCall && (destinationCreatedInternally(1) || isequal(destinationCreatedInternally,[false,true]))
        % Figure was created internally the figure handles was provided.
        % Set the figure to the same size as destination figure.
        h.figure.Units = sourceFig.Units;
        h.figure.Position(3:4) = sourceFig.Position(3:4);   % [21]
        movegui(h.figure) % Ensure it's on the screen
    end
    if ~destinationCreatedInternally(2) % If axes was not created internally (then the fig wasn't either)
        if firstCall
            % Axis position property is only copied when axes are produced internally.
            warning('The ''copyPosition'' property is ignored when the destination axes are created externally to %s.', mfilename())
        end
    else
        % Set destination axis position (but not units) [18, 20, 23].
        nestedParentFigPosition(p.Results.uiax, h.axes, firstCall);
        posCopied(1) = true; 
        
        % Set destination legend position (but not units) iff Location is 'none'
        if isfield(h,'legend') && ~isempty(legHand) && strcmpi(legHand.Location,'none')
            convertPosition(p, h, legHand, sourceFig)
            posCopied(2) = true; 
        end
        
        % Set destination colorbar position (but not units) iff Location is 'manual' 
        if isfield(h,'colorbar') && ~isempty(h.colorbar) && strcmpi(p.Results.colorbar.Location, 'manual')
            convertPosition(p, h, p.Results.colorbar, sourceFig)
            posCopied(3) = true; 
        end
    end
end
function convertPosition(p, h, obj, sourceFig) 
% The legend and colorbar position property for UIAxes is relative to the outerposition of the axes but for 
% regular axes the position properties are relative to the inner figure borders.  This function converts the 
% 1x4 position vector from uiaxis coordinates to destination-figure coordinates so that the object in the source
% figure is in the same position within the destination figure, assuming both figures are the same size.  
% Specifically, 'obj' is a legend or colorbar handle assigned to the dstination UIAxes.  obj.Position coordinates are relative 
% to p.Results.uiax.OuterPosition and are converted so that they are relative to the inner positoin of the figure 
% hosting h.axes.  The destination object is then moved.  It is assumed that the two figures are already the same size. 
% Get the handle to desination-object that corresponds to the obj which is in the source figure. 
drawnow(); pause(0.05)
assert(isprop(obj,'Type') && ismember(lower(obj.Type),{'legend','colorbar'}), 'Object repositioning only supports legends and colorbars.')
assert(isequal(ancestor(obj,'figure'),sourceFig), 'Input object must be an object within the source figure.')
if strcmpi(obj.Type, 'legend')
    obj(2) = h.legend; 
elseif strcmpi(obj.Type, 'colorbar')
    obj(2) = h.colorbar; 
else
    error('Object repositioning only supports legends and colorbars. Respositioning attempted with %s object.', obj.Type)
end
% Store original units and then change them temporarily
original.objUnits = obj(1).Units; 
original.sourceAxUnits = p.Results.uiax.Units; 
original.sourceFigUnits = sourceFig.Units; 
obj(1).Units = 'normalize'; 
obj(2).Units = 'normalize'; 
p.Results.uiax.Units = 'pixels'; 
sourceFig.Units = 'pixels';
% Coordinate conversion
normPos(1) = ((p.Results.uiax.OuterPosition(3) * obj(1).Position(1)) + p.Results.uiax.OuterPosition(1)) / sourceFig.InnerPosition(3); 
normPos(2) = ((p.Results.uiax.OuterPosition(4) * obj(1).Position(2)) + p.Results.uiax.OuterPosition(2)) / sourceFig.InnerPosition(4); 
normPos(3:4) = (p.Results.uiax.OuterPosition(3:4) .* obj(1).Position(3:4)) ./  sourceFig.InnerPosition(3:4); %assumes source figure size == desination fig size
% For testing only, confirm normalized coordinates by drawing a rectangle around the object in the source figure in norm coor.
% As of r2020a, annotation rectangles are buggy and do not appear so I'll use lines. Furthermore, applying annotation objects
% to uifigures in debug mode (only using live editor?) results in an unresponsive hang.  So we'll do this in an independent
% figure of equal size. 
%     testFig = figure('name', sprintf('%s_TEST_FIGURE',mfilename), 'Position', h.figure.Position);
%     annotation(testFig,'TextBox',normPos,'String','Legend');
% Move objects in destination figure
obj(2).Position = normPos;
% Return original units
obj(1).Units = original.objUnits; 
obj(2).Units = original.objUnits; % Use same units as destination
p.Results.uiax.Units = original.sourceAxUnits; 
sourceFig.Units = original.sourceFigUnits; 
function nestedParentFigPosition(source, destination, firstCall)
% The goal is to get the position of the source axes in hoObj relative to the source figure.
% That's not trivial since the axes could be tested in a Tab > TabGroup > Panel > GridLayout > Figure.  
% It computes the OuterPosition and InnerPosition of the source axes (source) relative to its parent figure
% in pixels units even if it's original units are not pixels.  In then moves the destination axes (destination)
% to the computed position within its host figure. The destination units are then returned to their original values. 
% Since this is potentially called >1 time per axes, to avoid redundant computation, the 3rd input (firstCall) is a 
% logical that is true only on the 1st time this is called.   When firstCall is false, there's no need to compute
% the positions again and the previous values are used.  
persistent outPos inPos
if firstCall
    drawnow(); pause(0.05)
    parentIsFig = false;
    innerPos = [];
    outerPos = [];
    h = source;
    while ~parentIsFig
        if isprop(h, 'Units')
            originalUnits = h.Units;
            h.Units = 'pixels';
            innerPos(end+1,:) = h.InnerPosition; %#ok<AGROW>
            outerPos(end+1,:) = h.OuterPosition; %#ok<AGROW>
            h.Units = originalUnits;
        end
        h = h.Parent;
        if isgraphics(h,'Figure')
            parentIsFig = true;
            continue
        end
    end
    % The first row of innerPos and outerPos is the input axes.  The last row is the child of the figure.
    % Add rows for col 1&2 to get horz & vert position in pixels.
    outPos = [sum(outerPos(:,1:2),1),outerPos(1,3:4)];
    inPos = [sum(innerPos(:,1:2),1),innerPos(1,3:4)];
end
% For testing only, confirm position by creating a figure and axes for comparison. 
%     testFig = figure('name', sprintf('%s_TEST_FIGURE',mfilename), 'Position', h.Position);
%     testAx = axes(testFig, 'Units', 'Pixels', 'OuterPosition', outPos, 'Position', inPos);  % see note [18]
originalDestinationUnits = destination.Units; 
destination.Units = 'pixels'; 
destination.OuterPosition = outPos; 
destination.Position = inPos;  % See note [18]
destination.Units = originalDestinationUnits; 
%% Notes
% [1] >580 views in past 30-days as of Oct 2019; >660 as of Feb 2019.
%   https://www.mathworks.com/matlabcentral/answers/281318-how-can-i-save-a-figure-within-app-designer
% [2] Some axes properties are read-only so we can't copy those over to the new axes.  Also, there are 
%   some axes properties that you don't want to copy over such as parent, children, X/Y/ZAxis, which will 
%   cause errors.  In addition to those, we will not copy the "Position" and "OuterPosition" properties
%   since the new axis position is pre-defined. 
% [3] Similar to [2], some legend and colorbar properties should not be copied.  I've chosen to not copy 
%   position due to the potential size differences between the new and old axes and we don't know what units 
%   the legend/cb contains so it could end up way off the plot.  User can adjust position afterwards. 
% [4] ZAxis (and maybe other properties) didn't exist in earlier releases (r2016a for example).
% [5] unique() needed because in earlier maltab release some of the properties in my additional bad field 
%   list are already included and duplicates cause error in rmfield().
% [6] Legend is currently a property of the UI axis. Previously it was not.  
%   Legends were not supported in app designer until r2016b:
%   https://www.mathworks.com/help/releases/R2016b/matlab/creating_guis/graphics-support-in-app-designer.html
% [7] In 2019a you can just call legend(h) with no strings and the legend obj will be created.  Testing in 
%   r2016b reveals that this method just returns an empty graphics obeject and a string is required.  This
%   is why we're copying the original legend strings over as soon as we create the legend. 
% [8] Unfortunately we can't use isprop(p.Results.uiax,'Legend') because it returns true in r2016a but then
%   when you try to access that property you get error: You cannot get the 'Legend' property of UIAxes. 
% [9] If I ever allow non-UIAxes, this test would have be a lot more flexible. The axis might be a polaraxis
%   or geoaxis etc but 'axes' will only confirm cartesian axes.  Here's a sloppy alternative:
%   @(h)~isempty(regexp(class(h),'^matlab.graphics.axis.','once'))
% [10] List of supported axes and plots for app designer: 
% 	https://www.mathworks.com/help/matlab/creating_guis/graphics-support-in-app-designer.html
% [11] Legends must be copied along with the uiaxes in copyobj() but since we can't use copyobj()
%   with UI axes, we cannot copy the legend. We created a new one but merely copying the 
%   legend properties does not preserve the 1:1 relationship between object handles and
%   thier legend strings.  For example, if the source plot has 2 objects but only the 2nd
%   one is represented in the legend, when the new legend is created it will detect that 
%   only 1 object is represented but it will show the displayname for the 1st object. 
%   Getting around this is really messy and the solution here may not always work.  If 
%   that's the case, it may be better to rebuild the legend. 
% [12] Several of the remaining fields contain sub-handles to graphics objects.  I'm not sure if copying
%   their values will results in any problems. To see which fields contain handles, 
%   goodParamIsHand = structfun(@ishghandle,goodParams,'UniformOutput',false);
% [13] The title, x/y/zLabel properties would copy fine from the source-axis but they have a 'parent' property 
%   that, when changed, would remove those objects from the source-axis rather than copying them.  There are also
%   some read-only properties of these objects that can't be copied.  Copying the position property can results in 
%   the label/title going off the figure border. Subtitle added in Matlab in r2020b. PositionConstraint added for 
%   2 reasons: 1-tiledLayout throws warning in 20b, 2-test 4.2 with results in a colapsed subplot.  
% [14] Below are links that show the rollout states of UIAxis support
%   * 2016a - 2017a : https://www.mathworks.com/help/releases/R2017a/matlab/creating_guis/graphics-support-in-app-designer.html
%   * 2017b - 2018a : https://www.mathworks.com/help/releases/R2018a/matlab/creating_guis/graphics-support-in-app-designer.html
%   * 2018b         : https://www.mathworks.com/help/releases/R2018b/matlab/creating_guis/graphics-support-in-app-designer.html
%   * 2019a         : https://www.mathworks.com/help/releases/R2019a/matlab/creating_guis/graphics-support-in-app-designer.html
%   * current       : https://www.mathworks.com/help/matlab/creating_guis/graphics-support-in-app-designer.html
% [15] See the link below regarding ColorOrderIndex & LineStyleOrderIndex.
%   https://www.mathworks.com/help/matlab/ref/matlab.graphics.axis.axes-properties.html#mw_23890c9d-6e04-4aab-92b5-f873dc229766
% [16] Prior to vs 1.2 user reported an error copying a categorical histogram (link A below). Categorical histograms use
%   a CategoricalRuler along the x or y axes [link B] while the default is a NumericRuler when an axis is created.  Categorical
%   data can be copied to a numericRuler axes but an error occurs when you copy the axis properties from a CategoricalRuler to
%   a NumericRuler.  Specifically, the XTick and XLim (or y-) values are problematic. Other categorical data plot (pie, scatter,
%   see link C) do not have this problem, only histogram() to my knowledge (as of vs 1.2.0).  To fix this, NumericRuler must be
%   converted to CategoricalRuler **prior to** copying the data to the axes.  The line that converts the axes was pulled from 
%   histogram (line 150) > categorical/categoricalHistogram (line 108) in r2019b, update 2. 
%   [A] https://www.mathworks.com/matlabcentral/answers/510081-how-to-export-xaxistick-labels-which-are-cell-arrays-to-figure-or-powerpoint-in-app-designer
%   [B] https://www.mathworks.com/help/matlab/ref/matlab.graphics.axis.decorator.categoricalruler-properties.html
%   [C] https://www.mathworks.com/help/matlab/matlab_prog/plot-categorical-data.html
% [17] Copying the XAxis and ZAxis come with all sorts of problems.  Many of the XAxis and YAxis properties are copied
%   when we copy the main axis handles properties, anyway.  For some reason, the YAxis is read-only and can't be copied
%   (see https://www.mathworks.com/matlabcentral/answers/438687-how-do-i-set-yaxis-of-axes-from-a-numeric-ruler-to-a-datetime-ruler).
%   The code below attempts to copy the X/ZAxes but I'm leaving it out for now. 
%     Try to copy the XAxis & ZAxis
%     YAxis is read-only for some reason [17].  Sometimes copying the X/ZAxis works but not always.
%     xzcopied = [false, false]; 
%     xzcopyProp = {'XAxis','ZAxis'};  % Must be exact prop names, case sensitive.
%     try
%         h.axes.XAxis = copy(p.Results.uiax.XAxis);
%         xzcopied(1) = true; 
%     catch
%         % No plan-B for now. Some XAxis properties will be copied later (see uiaxGoodParams).
%     end
%     try
%         h.axes.ZAxis = copy(p.Results.uiax.ZAxis);
%         xzcopied(2) = true; 
%     catch
%         % No plan-B for now. Some ZAxis properties will be copied later (see uiaxGoodParams). 
%     end
% 
%     % THIS LINE BELOW ALSO NEEDS TO BE ADDED LATER IN THE CODE
%     uiaxbadProps(ismember(uiaxbadProps.axis, xzcopyProp(xzcopied)),:) = []; %rm X|ZAxis if they were copied earlier
%
% [18] Inner/Outer/Position properties are defined differently for regular axes and UIAxes. According to the documentation, 
%   for regular axes, the innerposition==position which defines the axis box in 2D view; the outerPosition in regular axes
%   surrounds the axis labels and such.  For UIAxes, outerPosition==position and is the same as OuterPosition in regular
%   axes; the innerAxes property is read-only and is the same as the 'Position' property of regular axes. Furthermore, it 
%   seems like setting the OuterPos of the destination axis first is the right way to go.  
%   Also see: https://www.mathworks.com/help/matlab/creating_plots/automatic-axes-resize.html
% [19] https://www.mathworks.com/help/matlab/creating_plots/getting-and-setting-automatically-calculated-properties.html#buen4m6
% [20] When the legend location is set to 'none' or the colorbar location is set to 'manual', their positions are defined by 
%   the position property rather than the location property.  Here's a big problem:  For regular axes, the legend & colorbar 
%   position is relative to the figure borders.  For uiaxes, the legend & colorbar position is relative to the outerposition 
%   of the uiaxes.  I discovered this after trial and error, in r2020a.  Therefore, it doesn't make sense to copy the position
%   of colorbar and legend objects unless the axis position is also copied using copyPosition=true. So the legend and colorbar 
%   positions will follow these rules: 
%   A) If location is not none/manual, copy location as-is.  Position is set automatically.  
%   B) If copyPosition==false, replace none/manual with default locations. Position is set automatically.  
%   C) If copyPosition==true and the axes were produced internally, copy the none/manual location values and copy the units and position. 
%   D) If copyPosition==true and the axes were provided as input, throw a warning and do not copy any positions.
% [21] The innerposition and position properties of figures and UIFigures are always equal as of r2020a. 
% [22] In copyUIaxes_tester.mlx, tests that include a destination object as input (ie, test #2) produced figures with correct title
%   fontsize but when printed in live editor, the fontsize reduced to the default size (use "Contour & Catg Scat" for testing). 
%   After much troubleshooting and trial & error, for some reason, setting figure visibility to 'on', even if it was already on, solved
%   the problem most of the time. The problem was particularly reproducible when the input destination figure visibility was off.  
% [23] As of vs 3.0, when copyPositon = true, the axes position is copied relative to the figure even if the parent of the source
%   axes is not the figure (ie, a Tab or panel). 
% [24] Colormap wasn't a property of regular axes until r2018a (vs 9.4) which is when it first appeared in the documentation. Use colormap()
%   to set the colormap of an axes prior to r2018a.  
% [25] In Matlab r2020b when the colorbar is copied within a tiledlayout an error is thrown: "'Layout' value must be specified as a 
%   matlab.ui.layout.LayoutOptions object." Since the axes in tiledlayout are not uiaxes, we don't currently need to worry about copying
%   the layout option. As of r2020b, 'Layout' value must be specified as a matlab.ui.layout.LayoutOptions object.
%% Copyright (c) 2019, Adam Danz
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% * Redistributions of source code must retain the above copyright notice, this
% list of conditions and the following disclaimer.
% 
% * Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution
% * Neither the name of nor the names of itsoo
% contributors may be used to endorse or promote products derived from this
% software without specific prior written permission.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.