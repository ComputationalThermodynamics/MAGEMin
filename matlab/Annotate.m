classdef Annotate < handle %#ok<*INUSD,*INUSL,*CTCH>
%ANNOTATE Create an annotation that is pinned to the axes of a graph.
%
% Usage:
%
%   OBJ = ANNOTATE(AX, TYPE, X, Y)
%   OBJ = ANNOTATE(AX, TYPE, X, Y, 'PropertyName', Value, ...)
%
% Inputs:
%
%   AX <1x1 axes handle>
%     - Handle to axes where annotation is to be pinned
%
%   TYPE <character array>
%     - Type of annotation to create
%     - 'line', 'arrow', 'doublearrow', 'textarrow',
%       'rectangle', 'ellipse', or 'textbox'
%
%   X <1x2 or 2x1 numeric vector>
%     - Annotation X position in axes data coordinates
%     - Must lie within the visible axes limits
%     - For annotation types 'line', 'arrow', 'doublearrow', and 'textarrow', X(1)
%       specifies the tail end of the arrow and X(2) specifies the tip of the arrow head
%     - For annotation types 'rectangle', 'ellipse', and 'textbox', X(1)
%       specifies the left edge/extent and X(2) specifies the right edge/extent
%
%   Y <1x2 or 2x1 numeric vector>
%     - Annotation Y position in axes data coordinates
%     - Must lie within the visible axes limits
%     - For annotation types 'line', 'arrow', 'doublearrow', and 'textarrow', Y(1)
%       specifies the tail end of the arrow and Y(2) specifies the tip of the arrow head
%     - For annotation types 'rectangle', 'ellipse', and 'textbox', Y(1)
%       specifies the bottom edge/extent and Y(2) specifies the top edge/extent
%
%   Optional property name/value pairs can be specified to control the
%   appearance of the annotation object.  See 'annotation' documentation
%   and property pages for each annotation type for more details.
%
% Outputs:
%
%   OBJ <1x1 Annotate object handle>
%     - The Annotate object provides the following properties for
%       interaction with the annotation object.
%       * 'Position' allows for setting/getting the annotation object
%         position in axes data coordinates defined as [X(1),X(2),Y(1),Y(2)].
%         Note, rounding errors can occur.
%       * 'Primitive' is the handle to the annotation object, which allows
%         for further customization of its appearance.
%
% Description:
%
%   ANNOTATE creates an annotation object of type TYPE which is automatically
%   pinned to the axes with handle AX at the position specified by vectors
%   X and Y in axes data coordinates.
%
%   A UIContextMenu is attached to the annotation object which replicates the
%   standard annotation customizations (color, line width, line style, etc.),
%   normally accessed through edit plot mode, without having to enter edit plot
%   mode.  The UIContextMenu also provides move/resize and delete functionality.
%
%   ANNOTATE.ButtonDownFcn(SRC, EVNT, TYPE) is a static method available to
%   be executed as a button down callback function.  It interactively
%   creates an annotation object of type TYPE which is automatically pinned
%   to the axes ancestor of handle SRC.
%
%   ANNOTATE.ButtonDownFcn(___, 'PropertyName', Value, ...) interactively
%   creates the annotation object and applies the settings specified by the
%   property name/value pairs.
%
%   The property 'ParentAnnotation' is added to the annotation object so
%   its "parent annotation" Annotate object, constructed with the button
%   down callback function, can be recovered and used to set/get the
%   annotation position in axes data coordinates.
%
% Examples:
%
%   % programmatic placement through constructor
%   figure; x=0:0.1:6; y=sin(x); h=plot(x,y);
%   obj1 = Annotate(gca, 'ellipse', [1.5,2.5], [-0.8,0.4]);
%   obj2 = Annotate(gca, 'arrow', [x(1),x(10)]+0.2, [y(1),y(10)]);
%   obj3 = Annotate(gca, 'doublearrow', [4,3], [-0.6,0.4], 'linestyle', '--');
%   obj4 = Annotate(gca, 'line', [1,4], [0.6,0.6], 'color', 'r');
%   obj4.Position = [0.5,3.5,0.8,0.8]; % set new axes data coordinate position
%   obj5 = Annotate(gca, 'textbox', [2.5,4.5], [0.8,0.9], 'backgroundcolor', 'none', 'string', 'example textbox');
%   obj5.Primitive.EdgeColor = 'none'; obj5.Primitive.FontWeight = 'bold';
%   obj6 = Annotate(gca, 'textarrow', [4.5,x(55)], [0.2,y(55)], 'string', {'Programmatically','input text'});
%
%   % interactive placement through button down callback function
%   figure; x=0:0.1:6; y=sin(x); h=plot(x,y);
%   set(gca, 'ButtonDownFcn', {@Annotate.ButtonDownFcn, 'textarrow'}); % manually input text, hit escape when done
%   set(h, 'ButtonDownFcn', {@Annotate.ButtonDownFcn, 'rectangle'}); % can be set for axes child as well
%   % once created with the button down callback function, a primitive annotation could be set as the current object
%   % and a new axes data coordinate position could be specified with the following code
%   % obj = get(gco,'ParentAnnotation'); obj.Position = [1 3 0.0 0.5];
%
%   % support for uipanel hierarchy and multiple axes
%   figure; x=0:0.1:6; y=sin(x); p0=uipanel(gcf,'position',[0,0,1,1]);
%   p1=uipanel(p0,'position',[0.0,0.0,0.5,1.0]); ax1=axes('parent',p1); plot(ax1,x,y,'r');
%   p2=uipanel(p0,'position',[0.5,0.0,0.5,1.0]); ax2=axes('parent',p2); plot(ax2,x,y,'b');
%   obj1 = Annotate(ax1, 'ellipse', [1.5,2.5], [-0.8,0.4]);
%   obj2 = Annotate(ax2, 'arrow', [x(1),x(10)]+0.2, [y(1),y(10)]);
%   set(ax1, 'ButtonDownFcn', {@Annotate.ButtonDownFcn, 'arrow', 'color', 'b'});
%   set(ax2, 'ButtonDownFcn', {@Annotate.ButtonDownFcn, 'rectangle', 'color', 'r'});
%
% References:
%
%   http://undocumentedmatlab.com/blog/pinning-annotations-to-graphs
%   http://undocumentedmatlab.com/blog/inactive-control-tooltips-event-chaining
%   http://undocumentedmatlab.com/blog/adding-dynamic-properties-to-graphic-handles
%   http://undocumentedmatlab.com/blog/ishghandle-undocumented-input-parameter
%
% See also: annotation, getpixelposition, hgconvertunits, addlistener

% Author: Todd Baxter
% Date: 03-31-2020

    properties
        Primitive@handle scalar % primitive annotation object handle
    end

    properties (Dependent = true)
        Position@double vector % axes data coordinate position
    end

    properties (Dependent = true, Access = private)
        UIContextMenu
        WindowKeyPressFcn
    end

    properties (Access = private)
        AxesHandle
        FigureHandle
        ContainerHandle
        localUIContextMenuListener
        localCreateButtonMotionListener
        localCreateButtonUpListener
        localMoveResizeButtonMotionListener
        localMoveResizeButtonUpListener
        localKeyPressListener
        localMousePointerListener
        localDisableButtonDownListener
        localColorListener
        localLineWidthListener
        localPrimitivePositionListener
        localAxesPositionListener
        localEditPlotListener
        % callback execution helper properties
        prev_ptrloc % previous figure pointer location (in pixels)
        anchor_point % pointer location at start of creation (in pixels)
        callback_type % type of callback functionality current enabled
        prev_mousepointer % existing figure mouse pointer
        % private copy of dependent variables
        position_ % axes data coordinate position
        uicontextmenu_
    end

    properties (Constant = true, Access = private)
        % minimum annotation size after creation with button down callback
        MinimumWidth = 32; % pixels
        MinimumHeight = 18; % pixels
        % default text margin used for adjustments when line width changes
        DefaultTextMargin = 4; % pixels
    end

    methods
        % constructor
        function obj = Annotate(hAxes, annotation_type, xdata, ydata, varargin)
            %Annotate
            %   OBJ = ANNOTATE(AX,TYPE,X,Y) constructs a new Annotate
            %   object OBJ that creates a primitive annotation object of
            %   type TYPE which is automatically pinned to the axes with
            %   handle AX at the position specified by vectors X and Y in
            %   axes data coordinates.

            % input validation
            if ~ishghandle(hAxes,'axes')
                error('Annotate:invalidAxes', 'Invalid axes handle');
            end
            if ~isnumeric(xdata) || ~isvector(xdata) || numel(xdata) ~= 2
                error('Annotate:invalidX', 'X must be a 2-dimensional numeric vector');
            end
            if ~isnumeric(ydata) || ~isvector(ydata) || numel(ydata) ~= 2
                error('Annotate:invalidY', 'Y must be a 2-dimensional numeric vector');
            end
            xlims = get(hAxes, 'XLim');
            if any(xdata < xlims(1) | xdata > xlims(2))
                error('Annotate:outOfBoundsX', 'X must lie within its axes limits');
            end
            ylims = get(hAxes, 'YLim');
            if any(ydata < ylims(1) | ydata > ylims(2))
                error('Annotate:outOfBoundsY', 'Y must lie within its axes limits');
            end
            % make sure xdata and ydata are row vectors
            xdata = xdata(:)';
            ydata = ydata(:)';
            % get figure handle
            hFig = ancestor(hAxes, 'figure');
            % get primitive annotation container
            if ishg2()
                hContainer = get(hAxes,'Parent');
            else % HG1 backwards compatibility
                hContainer = hFig;
            end
            % convert from axes data space units to figure normalized units
            [figx, figy] = dsxy2figxy(hAxes, xdata, ydata);
            % create primitive annotation object
            switch lower(annotation_type)
                case {'line','arrow','doublearrow','textarrow'}
                    hThis = handle(annotation(hContainer, annotation_type, figx, figy));
                    if strcmpi(annotation_type,'textarrow')
                        set(hThis, 'HorizontalAlignment', 'left');
                        set(hThis, 'TextBackgroundColor', 'w');
                        set(hThis, 'TextEdgeColor', hThis.Color);
                        set(hThis, 'TextMargin', obj.DefaultTextMargin);
                    end
                case {'rectangle','ellipse','textbox'}
                    width = figx(2) - figx(1);
                    if width < 0
                        figx(1:2) = figx(2:-1:1);
                        width = -width;
                    end
                    height = figy(2) - figy(1);
                    if height < 0
                        figy(1:2) = figy(2:-1:1);
                        height = -height;
                    end
                    hThis = handle(annotation(hContainer, annotation_type, [figx(1),figy(1),width,height]));
                    if strcmpi(annotation_type,'textbox')
                        % convert normalized width/height to pixels
                        pixpos = hgconvertunits(hFig, [0,0,width,height], 'normalized', 'pixels', hContainer);
                        % allow textbox annotation to maintain its size if
                        % it is larger than the minimum
                        if pixpos(3) < obj.MinimumWidth && pixpos(4) < obj.MinimumHeight
                            % automatically size textbox annotation to fit
                            % text that is input (starting from a minimum
                            % size that is large enough to be visible)
                            normsize = hgconvertunits(hFig, [0,0,obj.MinimumWidth,obj.MinimumHeight], 'pixels', 'normalized', hContainer);
                            set(hThis, 'Position', [figx(1),figy(1),normsize(3),normsize(4)]);
                            set(hThis, 'FitBoxToText', 'on');
                        end
                        set(hThis, 'BackgroundColor', 'w');
                        set(hThis, 'EdgeColor', hThis.Color);
                        set(hThis, 'Margin', obj.DefaultTextMargin);
                    end
                otherwise
                    error('Annotate:unknownType', 'Unknown annotation type');
            end
            % apply optional primitive annotation property name/value pairs
            if numel(varargin) > 0
                set(hThis, varargin{:});
            end
            % pin the annotation to the axes
            for affNum = hThis.PinAff(:)'
                if ishg2()
                    pinAtAffordance(hThis, affNum);
                else % HG1 backwards compatibility
                    Annotate.pinAtAffordance(hThis, affNum);
                end
            end
            % enable textbox editing (if applicable) when the string
            % property has not been specified as an optional argument and
            % the constructor has not been called by button down callback
            s = dbstack;
            if isprop(hThis,'Editing') && ~any(strcmpi('String',varargin)) && ~any(strcmpi('Annotate.ButtonDownFcn',{s(:).name}))
                % click mouse outside annotation or hit escape key to exit
                hThis.Editing = 'on';
            end
            % set object handles
            obj.Primitive = hThis;
            obj.AxesHandle = hAxes;
            obj.FigureHandle = hFig;
            obj.ContainerHandle = hContainer;
            % set current figure mouse pointer
            obj.prev_mousepointer = get(hFig, 'Pointer');
            % set private copy of axes data coordinate position
            obj.position_ = [xdata,ydata];
            % create window button listeners
            obj.localUIContextMenuListener = addlistener(hFig, 'WindowMousePress', @obj.localUIContextMenuFcn);
            obj.localMoveResizeButtonMotionListener = addlistener(hFig, 'WindowMouseMotion', @obj.localMoveResizeButtonMotionFcn);
            set_listener_enabled(obj.localMoveResizeButtonMotionListener, 'off');
            obj.localMoveResizeButtonUpListener = addlistener(hFig, 'WindowMouseRelease', @obj.localMoveResizeButtonUpFcn);
            set_listener_enabled(obj.localMoveResizeButtonUpListener, 'off');
            obj.localMousePointerListener = addlistener(hFig, 'WindowMouseMotion', @obj.localMousePointerFcn);
            set_listener_enabled(obj.localMousePointerListener, 'off');
            obj.localDisableButtonDownListener = addlistener(hFig, 'WindowMousePress', @obj.localDisableButtonDownFcn);
            set_listener_enabled(obj.localDisableButtonDownListener, 'off');
            obj.localKeyPressListener = addlistener(hFig, 'WindowKeyPress', @obj.localKeyPressFcn);
            set_listener_enabled(obj.localKeyPressListener, 'off');
            % for consistency, be sure a figure key press function is set,
            % so the figure retains focus when keys are pressed
            if isempty(obj.WindowKeyPressFcn)
                set(obj.FigureHandle, 'WindowKeyPressFcn', @obj.localHelperKeyPressFcn);
            end
            % create property listeners
            if strcmpi(annotation_type,'textarrow')
                obj.localColorListener = addlistener(hThis, 'Color', 'PostSet', @obj.localColorFcn);
            end
            if ismember(lower(annotation_type), {'textbox','textarrow'})
                obj.localLineWidthListener = addlistener(hThis, 'LineWidth', 'PostSet', @obj.localTextMarginFcn);
            end
            if ismember(lower(annotation_type), {'rectangle','ellipse','textbox'})
                if ishg2()
                    % 'MarkedClean' event occurs when the axes changes
                    % location, size, or its x/y limits (even during a
                    % zoom/pan reset view)
                    obj.localAxesPositionListener = addlistener(hAxes, 'MarkedClean', @obj.localPositionFcn);
                else % HG1 backwards compatibility
                    obj.localPrimitivePositionListener = addlistener(hThis, 'Position', 'PostSet', @obj.localPositionFcn);
                end
                % listener to allow edit plot mode to move annotation
                obj.localEditPlotListener = addlistener(hFig, 'WindowMouseMotion', @obj.localEditPlotFcn);
            end
            % set primitive annotation object context menu
            if ishg2()
                obj.UIContextMenu = obj.Primitive.getScribeMenus();
            else % HG1 backwards compatibility
                obj.UIContextMenu = obj.Primitive.createScribeContextMenu(hFig);
            end
            % add property to primitive annotation object that references
            % this "parent annotation" object being constructed
            if ishg2()
                hProp = addprop(hThis, 'ParentAnnotation');
                hThis.ParentAnnotation = obj;
                hProp.SetAccess = 'private';
            else % HG1 backwards compatibility
                hProp = schema.prop(hThis, 'ParentAnnotation', 'mxArray');
                hThis.ParentAnnotation = obj;
                hProp.AccessFlags.PublicSet = 'off';
            end
            % set primitive annotation delete callback function to call the
            % class destructor so object does not persist when it is
            % created with a button down callback function
            set(hThis, 'DeleteFcn', @(src,evnt)delete(obj));
        end

        % destructor
        function delete(obj)
            try
                delete(obj.Primitive);
                set(obj.FigureHandle, 'Pointer', obj.prev_mousepointer);
                set(obj.FigureHandle, 'WindowKeyPressFcn', obj.WindowKeyPressFcn);
                % delete listeners
                delete(obj.localUIContextMenuListener);
                if ~isempty(obj.localCreateButtonMotionListener)
                    delete(obj.localCreateButtonMotionListener);
                end
                if ~isempty(obj.localCreateButtonUpListener)
                    delete(obj.localCreateButtonUpListener);
                end
                delete(obj.localMoveResizeButtonMotionListener);
                delete(obj.localMoveResizeButtonUpListener);
                delete(obj.localKeyPressListener);
                delete(obj.localMousePointerListener);
                delete(obj.localDisableButtonDownListener);
                if ~isempty(obj.localColorListener)
                    delete(obj.localColorListener);
                end
                if ~isempty(obj.localLineWidthListener)
                    delete(obj.localLineWidthListener);
                end
                if ~isempty(obj.localPrimitivePositionListener)
                    delete(obj.localPrimitivePositionListener);
                end
                if ~isempty(obj.localAxesPositionListener)
                    delete(obj.localAxesPositionListener);
                end
                if ~isempty(obj.localEditPlotListener)
                    delete(obj.localEditPlotListener);
                end
            catch
                % nothing to do
            end
        end

        function set.Position(obj, xy)
            % set primitive annotation object position according to input
            % specified as [X(1),X(2),Y(1),Y(2)] in axes data coordinates

            % input validation
            if ~isvector(xy) || numel(xy) ~= 4
                error('Annotate:invalidPosition', 'Position property must be specified as 4-dimensional vector [X(1),X(2),Y(1),Y(2)]');
            end
            % convert from axes data space units to figure normalized units
            [figx, figy] = dsxy2figxy(obj.AxesHandle, xy(1:2), xy(3:4));
            % disable annotation position listener
            if ~isempty(obj.localPrimitivePositionListener)
                set_listener_enabled(obj.localPrimitivePositionListener, 'off');
            end
            % update primitive annotation object position
            if isprop(obj.Primitive,'X')
                obj.Primitive.X = figx;
                obj.Primitive.Y = figy;
            else
                left = figx(1);
                width = figx(2) - figx(1);
                if width < 0
                    % reverse order of X position arguments
                    left = figx(2);
                    width = abs(width);
                    xy(1:2) = xy(2:-1:1);
                end
                bottom = figy(1);
                height = figy(2) - figy(1);
                if height < 0
                    % reverse order of Y position arguments
                    bottom = figy(2);
                    height = abs(height);
                    xy(3:4) = xy(4:-1:3);
                end
                obj.Primitive.Position = [left,bottom,width,height];
            end
            % re-enable annotation position listener
            if ~isempty(obj.localPrimitivePositionListener)
                set_listener_enabled(obj.localPrimitivePositionListener, 'on');
            end
            % set private copy of axes data coordinate position
            obj.position_ = xy;
        end

        function xy = get.Position(obj)
            % return the primitive annotation object position specified as
            % [X(1),X(2),Y(1),Y(2)] in axes data coordinates
            % note: rounding errors can occur

            % get figure position of annotation in pixels
            hFig = obj.FigureHandle;
            figpos = hgconvertunits(hFig, obj.Primitive.Position, obj.Primitive.Units, 'pixels', obj.ContainerHandle);
            % get pixel position of axes
            hAxes = obj.AxesHandle;
            if ishg2()
                % w.r.t. axes parent
                axpos = getpixelposition(hAxes);
            else % HG1 backwards compatibility
                % w.r.t. figure
                axpos = getpixelposition(hAxes, true);
            end
            % get the axis limits [xlim ylim (zlim)]
            axlim = axis(hAxes);
            % reverse order of axes limits (if necessary)
            if strcmp(get(hAxes,'xdir'),'reverse')
                axlim(1:2) = axlim(2:-1:1);
            end
            if strcmp(get(hAxes,'ydir'),'reverse')
                axlim(3:4) = axlim(4:-1:3);
            end
            % convert axes limits to linear scale (if necessary)
            if strcmp(get(hAxes,'xscale'),'log')
                axlim(1:2) = log10(axlim(1:2));
            end
            if strcmp(get(hAxes,'yscale'),'log')
                axlim(3:4) = log10(axlim(3:4));
            end
            % compute axes width and height
            axwidth = diff(axlim(1:2));
            axheight = diff(axlim(3:4));
            % compute the annotation position in axes data coordinatess
            pos(1) = (figpos(1) - axpos(1)) / axpos(3) * axwidth + axlim(1);
            pos(2) = (figpos(2) - axpos(2)) / axpos(4) * axheight + axlim(3);
            pos(3) = figpos(3) * axwidth / axpos(3);
            pos(4) = figpos(4) * axheight / axpos(4);
            % return X,Y position pair
            xy = [pos(1), pos(1)+pos(3), pos(2), pos(2)+pos(4)];
            % convert X,Y position pair to logarithmic scale (if necessary)
            if strcmp(get(hAxes,'xscale'),'log')
                xy(1:2) = 10.^xy(1:2);
            end
            if strcmp(get(hAxes,'yscale'),'log')
                xy(3:4) = 10.^xy(3:4);
            end
        end

        function set.UIContextMenu(obj, hMenus)
            % build a UIContextMenu for the primitive annotation object so
            % it can be used without having to enable edit plot mode

            uicm_obj = uicontextmenu;
            uimenu(uicm_obj, 'Label', 'Move/Resize', 'Callback', @obj.enableLocalCallbacks);
            % reparent the built-in annotation context menu items to the
            % local UIContextMenu and turn on the visibility so they
            % actually get displayed
            for i = 1:length(hMenus)
                hMenu = hMenus(i);
                set(hMenu, 'Parent', uicm_obj, 'Visible', 'on');
                % built-in arrow and textarrow reverse direction menu does
                % not work without edit plot mode enabled
                if strcmp(get(hMenu,'Label'), 'Reverse Direction')
                    set(hMenu, 'Callback', @obj.localReverseDirectionFcn);
                end
                if i == 1
                    set(hMenu, 'Separator', 'on');
                end
                subMenus = allchild(hMenu); % find all hidden menu children
                set(subMenus(end:-1:1), 'Parent', hMenu, 'Visible', 'on');
            end
            uimenu(uicm_obj, 'Label', 'Delete', 'Callback', @(src,evnt)delete(obj), 'Separator', 'on');
            % set the UIContextMenu for the primitive annotation object
            obj.Primitive.UIContextMenu = uicm_obj;
            % save this UIContextMenu object so it can be enabled and
            % disabled based on whether edit plot mode is active or not
            obj.uicontextmenu_ = uicm_obj;
        end

        function uicm_obj = get.UIContextMenu(obj)
            % return the local UIContextMenu unless edit plot mode is
            % active, because those context menu items would be duplicated

            if isactiveuimode(obj.FigureHandle,'Standard.EditPlot')
                uicm_obj = [];
                return;
            end
            uicm_obj = obj.uicontextmenu_;
        end

        function set.WindowKeyPressFcn(obj, keypressfcn)
            % set the figure key press function

            set(obj.FigureHandle, 'WindowKeyPressFcn', keypressfcn);
        end

        function keypressfcn = get.WindowKeyPressFcn(obj)
            % return the current figure key press function, ignoring the
            % local helper key press function that could be set

            keypressfcn = get(obj.FigureHandle, 'WindowKeyPressFcn');
            if isa(keypressfcn,'function_handle') && ~isempty(regexp(func2str(keypressfcn),'obj.localHelperKeyPressFcn','once'))
                keypressfcn = '';
            end
        end
    end

    methods (Static = true)
        function obj = ButtonDownFcn(src, evnt, annotation_type, varargin)
            %Annotate.ButtonDownFcn
            %   OBJ = ANNOTATE.ButtonDownFcn(SRC,EVNT,TYPE) interactively
            %   creates a primitive annotation object of type TYPE which is
            %   automatically pinned to the axes ancestor of handle SRC.
            %   This static method can be directly set as a button down
            %   callback function, or it can be wrapped in another button
            %   down callback function to get access to the returned
            %   Annotate object OBJ that is created.

            % get axes handle
            hAxes = ancestor(src, 'axes');
            % get clicked point in axes coordinates
            currPt = get(hAxes, 'CurrentPoint');
            xdata = currPt(1,1);
            ydata = currPt(1,2);
            % replicate this current point for annotation creation
            xdata = repmat(xdata,1,2);
            ydata = repmat(ydata,1,2);
            % construct Annotate object
            obj = Annotate(hAxes, annotation_type, xdata, ydata, varargin{:});
            % set creation anchor point as pointer location w.r.t. figure
            hFig = ancestor(hAxes, 'figure');
            ptrloc = hgconvertunits(hFig, [get(hFig,'CurrentPoint'),0,0], get(hFig,'Units'), 'Pixels', hFig);
            obj.anchor_point = floor(ptrloc(1:2));
            % creation listener callbacks
            if isprop(obj.Primitive,'MoveStyle')
                obj.Primitive.MoveStyle = 'topright';
            else % HG1 backwards compatibility
                obj.Primitive.MoveMode = 'topright';
            end
            % disable annotation position listener
            if ~isempty(obj.localPrimitivePositionListener)
                set_listener_enabled(obj.localPrimitivePositionListener, 'off');
            end
            % add creation listeners
            obj.localCreateButtonMotionListener = addlistener(hFig, 'WindowMouseMotion', @obj.localCreateButtonMotionFcn);
            % enable textbox editing (if applicable) when the string
            % property has not been specified as an optional argument
            enable_edit = false;
            if isprop(obj.Primitive,'Editing') && ~any(strcmpi('String',varargin))
                enable_edit = true;
            end
            obj.localCreateButtonUpListener = addlistener(hFig, 'WindowMouseRelease', @(src,evnt)obj.localCreateButtonUpFcn(src,evnt,enable_edit));
        end
    end

    methods (Access = private)
        function localUIContextMenuFcn(obj, src, evnt)
            % detect if mouse has been clicked on the annotation that
            % enabled this listener callback function
            if isprop(evnt,'HitObject')
                HitObject = evnt.HitObject;
            else % HG1 backwards compatibility
                HitObject = evnt.CurrentObject;
            end
            if isequal(obj.Primitive,HitObject)
                obj.Primitive.UIContextMenu = obj.UIContextMenu;
            end
        end

        function localEditPlotFcn(obj, src, evnt)
            % allow annotation position to change in edit plot mode
            hFig = obj.FigureHandle;
            if isactiveuimode(hFig,'Standard.EditPlot') && isequal(obj.Primitive,get(hFig,'CurrentObject'))
                % set private copy of axes data coordinate position
                obj.position_ = obj.Position;
            end
        end

        function localReverseDirectionFcn(obj, src, evnt)
            obj.Primitive.X = fliplr(obj.Primitive.X);
            obj.Primitive.Y = fliplr(obj.Primitive.Y);
        end

        function enableLocalCallbacks(obj, src, evnt)
            set(obj.Primitive, 'Selected', 'on');
            set_listener_enabled(obj.localMousePointerListener, 'on');
            set_listener_enabled(obj.localDisableButtonDownListener, 'on');
            set_listener_enabled(obj.localKeyPressListener, 'on');
            if isempty(obj.WindowKeyPressFcn)
                % be sure a figure key press function is set, so the figure
                % retains focus when keys are pressed
                set(obj.FigureHandle, 'WindowKeyPressFcn', @obj.localHelperKeyPressFcn);
            end
            set(obj.Primitive, 'ButtonDownFcn', @obj.localMoveResizeButtonDownFcn);
        end

        function localMousePointerFcn(obj, src, evnt)
            if ~isprop(evnt,'CurrentPoint')
                mode_type = obj.Primitive.findMoveMode(evnt);
            else % HG1 backwards compatibility
                mode_type = obj.Primitive.findMoveMode(evnt.CurrentPoint);
            end
            % set the appropriate mouse pointer
            switch lower(mode_type)
                case 'none'
                    set(obj.FigureHandle, 'Pointer', obj.prev_mousepointer);
                case 'mouseover'
                    set(obj.FigureHandle, 'Pointer', 'fleur');
                case {'topleft','bottomright'}
                    set(obj.FigureHandle, 'Pointer', 'topl');
                case {'topright','bottomleft'}
                    set(obj.FigureHandle, 'Pointer', 'topr');
                case {'left','right'}
                    set(obj.FigureHandle, 'Pointer', 'left');
                case {'top','bottom'}
                    set(obj.FigureHandle, 'Pointer', 'bottom');
                otherwise
                    set(obj.FigureHandle, 'Pointer', obj.prev_mousepointer);
            end
        end

        function localDisableButtonDownFcn(obj, src, evnt)
            % detect if mouse has been clicked somewhere other than the
            % annotation that enabled this listener callback function
            if isprop(evnt,'HitObject')
                HitObject = evnt.HitObject;
            else % HG1 backwards compatibility
                HitObject = evnt.CurrentObject;
            end
            if ~isequal(obj.Primitive,HitObject)
                disableLocalCallbacks(obj);
            end
        end

        function localHelperKeyPressFcn(obj, src, evnt)
            % a figure window key press listener has been defined to handle
            % the key press functionality, this callback is only needed to
            % ensure the figure retains focus when keys are pressed
            %
            % note: this, or a similar, approach is needed because in HG2
            % disabling a figure key press function from within that
            % function itself will cause the figure to lose focus
        end

        function localKeyPressFcn(obj, src, evnt)
            % move the primitive annotation object in the direction of the
            % arrow key that was pressed, or disable callbacks when the
            % escape (ESC) key is pressed
            if ishg2()
                key = lower(evnt.Key);
            else % HG1 backwards compatibility
                key = get(obj.FigureHandle, 'CurrentCharacter');
            end
            switch key
                case {27,'escape'}
                    disableLocalCallbacks(obj);
                    return;
                case {28,'leftarrow'}
                    delta = [-1 0];
                case {29,'rightarrow'}
                    delta = [+1 0];
                case {30,'uparrow'}
                    delta = [0 +1];
                case {31,'downarrow'}
                    delta = [0 -1];
                otherwise
                    return;
            end
            % disable annotation position listener
            if ~isempty(obj.localPrimitivePositionListener)
                set_listener_enabled(obj.localPrimitivePositionListener, 'off');
            end
            % move primitive annotation object in direction of arrow key
            obj.Primitive.move(delta);
            % set private copy of axes data coordinate position
            obj.position_ = obj.Position;
            % re-enable annotation position listener
            if ~isempty(obj.localPrimitivePositionListener)
                set_listener_enabled(obj.localPrimitivePositionListener, 'on');
            end
        end

        function disableLocalCallbacks(obj)
            % disable callbacks when mouse is clicked somewhere other than
            % the annotation that enabled the listener callback function,
            % or when the escape (ESC) key is pressed
            set(obj.Primitive, 'Selected', 'off');
            set(obj.Primitive, 'ButtonDownFcn', '');
            set_listener_enabled(obj.localMousePointerListener, 'off');
            set_listener_enabled(obj.localDisableButtonDownListener, 'off');
            set_listener_enabled(obj.localKeyPressListener, 'off');
            set(obj.FigureHandle, 'Pointer', obj.prev_mousepointer);
        end

        function localCreateButtonMotionFcn(obj, src, evnt)
            if isprop(evnt,'Point')
                ptrloc = evnt.Point;
            else % HG1 backwards compatibility
                ptrloc = evnt.CurrentPoint;
            end
            hFig = obj.FigureHandle;
            if ~strcmpi(get(hFig,'units'),'pixels')
                ptrloc = hgconvertunits(hFig, [ptrloc,0,0], get(hFig,'units'), 'pixels', hFig);
                ptrloc = ptrloc(1:2);
            end
            hContainer = obj.ContainerHandle;
            if ~isequal(hFig,hContainer)
                containerpos = getpixelposition(hContainer, true);
                ptrloc = ptrloc - containerpos(1:2);
            end
            obj.Primitive.resize(ptrloc);
        end

        function localCreateButtonUpFcn(obj, src, evnt, enable_edit)
            set_listener_enabled(obj.localCreateButtonMotionListener, 'off');
            set_listener_enabled(obj.localCreateButtonUpListener, 'off');
            % calculate how far mouse has moved from the anchor point,
            % which indicates the size of the annotation
            if isprop(evnt,'Point')
                ptrloc = evnt.Point;
            else % HG1 backwards compatibility
                ptrloc = evnt.CurrentPoint;
            end
            hFig = obj.FigureHandle;
            if ~strcmpi(get(hFig,'units'),'pixels')
                ptrloc = hgconvertunits(hFig, [ptrloc,0,0], get(hFig,'units'), 'pixels', hFig);
                ptrloc = ptrloc(1:2);
            end
            deltaPt = ptrloc - obj.anchor_point;
            % protect against the user accidentally creating an annotation
            % that is too small to be easily visible
            normsize = hgconvertunits(hFig, [0,0,obj.MinimumWidth,obj.MinimumHeight], 'pixels', 'normalized', obj.ContainerHandle);
            if isprop(obj.Primitive,'X')
                if norm(deltaPt) < norm([obj.MinimumWidth,obj.MinimumHeight])
                    % set minimum sizes in same direction as mouse moved
                    direction = sign(deltaPt);
                    direction(direction == 0) = 1;
                    obj.Primitive.Position = [obj.Primitive.Position(1:2),normsize(3:4).*direction];
                end
            else
                if abs(deltaPt(1)) < obj.MinimumWidth && abs(deltaPt(2)) < obj.MinimumHeight
                    % set minimum sizes in same direction as mouse moved
                    if obj.anchor_point(1) < ptrloc(1)
                        left = obj.Primitive.Position(1);
                    else
                        left = obj.Primitive.Position(1) + obj.Primitive.Position(3) - normsize(3);
                    end
                    if obj.anchor_point(2) < ptrloc(2)
                        bottom = obj.Primitive.Position(2);
                    else
                        bottom = obj.Primitive.Position(2) + obj.Primitive.Position(4) - normsize(4);
                    end
                    obj.Primitive.Position = [left,bottom,normsize(3:4)];
                end
            end
            % set private copy of axes data coordinate position
            obj.position_ = obj.Position;
            % re-enable annotation position listener
            if ~isempty(obj.localPrimitivePositionListener)
                set_listener_enabled(obj.localPrimitivePositionListener, 'on');
            end
            % enable text editing only when annotation creation is done
            if enable_edit
                % click mouse outside annotation or hit escape key to exit
                obj.Primitive.Editing = 'on';
            end
        end

        function localMoveResizeButtonDownFcn(obj, src, evnt)
            % initialize the previous pointer location
            hFig = obj.FigureHandle;
            currPt = hgconvertunits(hFig, [get(hFig,'CurrentPoint'),0,0], get(hFig,'Units'), 'Pixels', hFig);
            obj.prev_ptrloc = currPt(1:2);
            % enable window button motion/up listeners
            set_listener_enabled(obj.localMoveResizeButtonMotionListener, 'on');
            set_listener_enabled(obj.localMoveResizeButtonUpListener, 'on');
            % disable mouse pointer listener
            set_listener_enabled(obj.localMousePointerListener, 'off');
            % disable annotation position listener
            if ~isempty(obj.localPrimitivePositionListener)
                set_listener_enabled(obj.localPrimitivePositionListener, 'off');
            end
            % set callback type
            if isprop(obj.Primitive,'MoveStyle')
                mode_type = obj.Primitive.MoveStyle;
            else % HG1 backwards compatibility
                mode_type = obj.Primitive.MoveMode;
            end
            if strcmpi(mode_type,'mouseover')
                obj.callback_type = 'move';
            else
                obj.callback_type = 'resize';
            end
        end

        function localMoveResizeButtonMotionFcn(obj, src, evnt)
            % get pointer location w.r.t. figure
            if isprop(evnt,'Point')
                ptrloc = evnt.Point;
            else % HG1 backwards compatibility
                ptrloc = evnt.CurrentPoint;
            end
            hFig = obj.FigureHandle;
            if ~strcmpi(get(hFig,'units'),'pixels')
                ptrloc = hgconvertunits(hFig, [ptrloc,0,0], get(hFig,'units'), 'pixels', hFig);
                ptrloc = ptrloc(1:2);
            end
            % move or resize primitive annotation object
            switch obj.callback_type
                case 'move'
                    obj.Primitive.move(ptrloc-obj.prev_ptrloc);
                case 'resize'
                    hContainer = obj.ContainerHandle;
                    if ~isequal(hFig,hContainer)
                        containerpos = getpixelposition(hContainer, true);
                        ptrloc = ptrloc - containerpos(1:2);
                    end
                    obj.Primitive.resize(ptrloc);
            end
            % set previous pointer location
            obj.prev_ptrloc = ptrloc;
            % set private copy of axes data coordinate position
            obj.position_ = obj.Position;
        end

        function localMoveResizeButtonUpFcn(obj, src, evnt)
            % disable window button motion/up listeners
            set_listener_enabled(obj.localMoveResizeButtonMotionListener, 'off');
            set_listener_enabled(obj.localMoveResizeButtonUpListener, 'off');
            % re-enable mouse pointer listener
            set_listener_enabled(obj.localMousePointerListener, 'on');
            % re-enable annotation position listener
            if ~isempty(obj.localPrimitivePositionListener)
                set_listener_enabled(obj.localPrimitivePositionListener, 'on');
            end
        end

        function localColorFcn(obj, src, evnt)
            % this function is needed for textarrow annotations to keep the
            % textedgecolor property in sync with the color property
            obj.Primitive.TextEdgeColor = obj.Primitive.Color;
        end

        function localTextMarginFcn(obj, src, evnt)
            % adjust text margin as line width changes
            if isprop(obj.Primitive,'TextMargin')
                obj.Primitive.TextMargin = obj.DefaultTextMargin + fix((obj.Primitive.LineWidth - 1)/2);
            end
            if isprop(obj.Primitive,'Margin')
                obj.Primitive.Margin = obj.DefaultTextMargin + fix((obj.Primitive.LineWidth - 1)/2);
            end
        end

        function localPositionFcn(obj, src, evnt)
            % this function is needed for annotations with only one
            % affordance ('rectangle' and 'ellipse') in order to maintain
            % the same axes data coordinate position
            if ~isempty(obj.localPrimitivePositionListener)
                set_listener_enabled(obj.localPrimitivePositionListener, 'off');
            end
            obj.Position = obj.position_;
            if ~isempty(obj.localPrimitivePositionListener)
                set_listener_enabled(obj.localPrimitivePositionListener, 'on');
            end
        end
    end

    % The following is the HG1 version of the 'pinAtAffordance' function code,
    % copied from 'toolbox\matlab\scribe\@scribe\@scribeobject\pinAtAffordance.m',
    % with a minor update to its 'pointinaxes' local function that fixes a bug it
    % had with axes embedded in a uipanel hierarchy.
    methods (Access = private, Static = true)
        function pinAtAffordance(hThis, affNum)
            % Pin a scribe object at the given affordance

            %   Copyright 2006 The MathWorks, Inc.

            % First, unpin any pins at a given affordance:
            hThis.unpinAtAffordance(affNum);

            % First, convert the affordance position into pixels representing a the
            % point in the figure
            hAff = hThis.Srect(affNum);
            point = [get(hAff,'XData') get(hAff,'YData')];
            hFig = ancestor(hThis,'Figure');
            point = hgconvertunits(hFig,[point 0 0],'Normalized','Pixels',hFig);
            point = point(1:2);

            % Before we create a pin, make sure we are over a pinnable object
            % First check for a valid axes:
            axlist = findobj(hFig,'type','axes');
            pinax = [];
            for i=1:length(axlist)
                if Annotate.pointinaxes(axlist(i),point) && ~isappdata(axlist(i),'NonDataObject')
                    pinax = axlist(i);
                end
            end

            % If there is no axes that the annotation can be pinned to, return early
            % and don't even create a pin.
            if isempty(pinax)
                return;
            end

            pinobj = [];
            % Turn off the hittest property of this group to find out if we are over an
            % object
            hitState = get(hThis,'HitTest');
            set(hThis,'HitTest','off');
            obj = handle(hittest(hFig,point));
            set(hThis,'HitTest',hitState);

            if ~isempty(obj) && ~isa(obj,'scribe.scribeobject') && ...
                    ~strcmpi(get(obj,'tag'),'DataTipMarker')
                type = get(obj,'type');
                if strcmpi(type,'surface')||strcmpi(type,'patch')||strcmpi(type,'line')
                    pinobj = obj;
                end
            end

            % Create the pin
            hPin = scribe.scribepin('Parent',hThis.Parent,'Target',hThis,'DataAxes',pinax,'Affordance',affNum);
            % Create the correspondence between the axes and the pin
            repin(hPin(end),point,pinax,pinobj);
            set(hPin,'Enable','on');

            % Link the pin to the scribe object
            hThis.Pin(end+1) = hPin;
        end

        function isin = pointinaxes(ax,p)
            pos = getpixelposition(ax,true);
            isin = false;
            if p(1)>=pos(1) && p(2)>=pos(2) && ...
                    p(1)<=(pos(1)+pos(3)) && p(2)<=(pos(2)+pos(4))
                isin = true;
            end
        end
    end
end
