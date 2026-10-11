function varargout = paper_figure_utils(action, varargin)
%PAPER_FIGURE_UTILS Shared paper-figure layout, legend and export utilities.
%
% This helper centralizes the figure behavior that was previously implemented
% only in generate_paper_figures.m.  It is intentionally presentation-only:
% it does not alter simulation data, controller outputs, statistics or seeds.
%
% Usage
%   paper_figure_utils('defaults');
%   [fig,layout] = paper_figure_utils('create', name, position, nrow, ncol, titleText);
%   ax = paper_figure_utils('tile', layout, tileIndex);
%   paper_figure_utils('panel_title', ax, '(a) Panel title');
%   loc = paper_figure_utils('align_legends', axesVector, legendVector);
%   paper_figure_utils('export', fig, outputDir, baseName); % PDF + PNG with outer white margin
%   paper_figure_utils('export', fig, outputDir, baseName, ...
%       'SaveJPEG',true,'JPEGResolution',1000); % add 1000-dpi JPEG
%
% PDF export uses print -dpdf -painters with an explicitly enlarged page.
% Raster export requests an explicit outer padding when supported (R2025a+)
% and falls back to print on older MATLAB releases.  This prevents titles,
% labels and legends from visually touching surrounding manuscript text.

    if nargin < 1 || isempty(action)
        error('paper_figure_utils:MissingAction','An action is required.');
    end

    switch lower(char(action))
        case 'defaults'
            apply_defaults();
        case 'create'
            [varargout{1:nargout}] = create_layout(varargin{:});
        case 'tile'
            varargout{1} = create_tile(varargin{:});
        case 'panel_title'
            set_panel_title(varargin{:});
        case 'align_legends'
            [varargout{1:nargout}] = align_legends(varargin{:});
        case 'export'
            export_pair(varargin{:});
        otherwise
            error('paper_figure_utils:UnknownAction','Unknown action: %s',char(action));
    end
end

function apply_defaults()
    set(groot, 'DefaultAxesFontSize', 11);
    set(groot, 'DefaultTextFontSize', 12);
    set(groot, 'DefaultLineLineWidth', 1.5);
    set(groot, 'DefaultAxesLabelFontSizeMultiplier', 1.05);
    set(groot, 'DefaultAxesTitleFontWeight', 'normal');
    set(groot, 'DefaultAxesXGrid', 'on');
    set(groot, 'DefaultAxesYGrid', 'on');
end

function [fig,layout] = create_layout(fig_name, fig_pos, nrow, ncol, main_title, varargin)
    p = inputParser;
    addParameter(p,'Visible','on',@(x)ischar(x)||isstring(x));
    addParameter(p,'ShowMainTitle',true,@(x)islogical(x)||isnumeric(x));
    parse(p,varargin{:});
    o = p.Results;

    apply_defaults();
    fig = figure('Name',char(fig_name),'Position',fig_pos,'Color','w', ...
                 'Visible',char(o.Visible),'InvertHardcopy','off');

    % Explicit axes geometry prevents subplot/tiledlayout from squeezing text
    % into adjacent panels or changing inner plot widths when tick labels differ.
    if nrow == 3 && ncol == 1
        left = 0.115; right = 0.045; bottom = 0.075; top = 0.118;
        row_gap = 0.105; col_gap = 0.0;
    elseif nrow == 4 && ncol == 1
        left = 0.115; right = 0.045; bottom = 0.070; top = 0.115;
        row_gap = 0.080; col_gap = 0.0;
    elseif nrow == 2 && ncol == 1
        left = 0.115; right = 0.045; bottom = 0.090; top = 0.125;
        row_gap = 0.135; col_gap = 0.0;
    elseif nrow == 1 && ncol == 3
        left = 0.070; right = 0.035; bottom = 0.155; top = 0.118;
        row_gap = 0.0; col_gap = 0.080;
    elseif nrow == 1 && ncol == 2
        left = 0.085; right = 0.040; bottom = 0.145; top = 0.125;
        row_gap = 0.0; col_gap = 0.090;
    elseif nrow == 2 && ncol == 2
        left = 0.085; right = 0.040; bottom = 0.095; top = 0.125;
        row_gap = 0.135; col_gap = 0.090;
    elseif nrow == 3 && ncol == 3
        left = 0.070; right = 0.035; bottom = 0.075; top = 0.118;
        row_gap = 0.105; col_gap = 0.065;
    else
        left = 0.095; right = 0.045; bottom = 0.090; top = 0.125;
        row_gap = 0.100; col_gap = 0.070;
    end

    axes_w = (1-left-right-(ncol-1)*col_gap)/ncol;
    axes_h = (1-top-bottom-(nrow-1)*row_gap)/nrow;
    if axes_w <= 0 || axes_h <= 0
        error('paper_figure_utils:InvalidLayout','Margins leave no room for axes.');
    end

    layout = struct('fig',fig,'nrow',nrow,'ncol',ncol, ...
                    'left',left,'right',right,'bottom',bottom,'top',top, ...
                    'row_gap',row_gap,'col_gap',col_gap, ...
                    'axes_w',axes_w,'axes_h',axes_h);

    if logical(o.ShowMainTitle) && nargin >= 5 && ~isempty(main_title)
        title_gap = 0.008;
        title_bottom = 1-top+title_gap;
        title_height = max(0.028,top-title_gap-0.010);
        annotation(fig,'textbox',[left,title_bottom,1-left-right,title_height], ...
            'String',char(main_title),'FontSize',13,'FontWeight','bold', ...
            'HorizontalAlignment','center','VerticalAlignment','middle', ...
            'EdgeColor','none','Interpreter','none','FitBoxToText','off');
    end
end

function ax = create_tile(layout, tile_idx)
    if tile_idx < 1 || tile_idx > layout.nrow*layout.ncol
        error('paper_figure_utils:BadTile','tile_idx is outside the layout.');
    end
    row = ceil(tile_idx/layout.ncol);
    col = tile_idx-(row-1)*layout.ncol;
    left = layout.left+(col-1)*(layout.axes_w+layout.col_gap);
    top_edge = 1-layout.top-(row-1)*(layout.axes_h+layout.row_gap);
    bottom = top_edge-layout.axes_h;

    ax = axes('Parent',layout.fig,'Units','normalized', ...
              'Position',[left,bottom,layout.axes_w,layout.axes_h]);
    set(ax,'FontSize',11,'Box','on','Layer','top','TickDir','out', ...
           'XGrid','on','YGrid','on');
    try
        ax.PositionConstraint = 'innerposition';
    catch
        try
            ax.ActivePositionProperty = 'position'; %#ok<SETPROP>
        catch
        end
    end
end

function set_panel_title(ax, txt)
    grid(ax,'on');
    set(ax,'XGrid','on','YGrid','on');
    h = title(ax,char(txt),'FontSize',11,'FontWeight','bold','Interpreter','none');
    try
        h.Units = 'normalized';
        p = h.Position;
        p(2) = 1.018;
        h.Position = p;
    catch
    end
end

function [location, legend_handles] = align_legends(ax_list, legend_handles, varargin)
%ALIGN_LEGENDS Choose one low-occupancy corner and use it in all subplots.
% The scoring is performed jointly across the supplied axes, so legends in a
% multi-panel figure stay in the same corner whenever possible.
    p = inputParser;
    addParameter(p,'NumColumns',[],@(x)isempty(x)||(isnumeric(x)&&isscalar(x)&&x>=1));
    addParameter(p,'FontSize',9.0,@(x)isnumeric(x)&&isscalar(x)&&x>0);
    parse(p,varargin{:});
    o = p.Results;

    ax_list = ax_list(isgraphics(ax_list));
    legend_handles = legend_handles(isgraphics(legend_handles));
    if isempty(legend_handles)
        location = 'northeast';
        return;
    end

    candidates = {'northeast','northwest','southeast','southwest'};
    scores = zeros(1,numel(candidates));
    for a = 1:numel(ax_list)
        scores = scores + corner_occupancy(ax_list(a));
    end
    [~,idx] = min(scores);
    location = candidates{idx};

    for k = 1:numel(legend_handles)
        lgd = legend_handles(k);
        try, lgd.Location = location; catch, set(lgd,'Location',location); end
        try, lgd.AutoUpdate = 'off'; catch, end
        try, lgd.Box = 'on'; catch, end
        try, lgd.FontSize = o.FontSize; catch, end
        try, lgd.Interpreter = 'none'; catch, end

        ncol = o.NumColumns;
        if isempty(ncol)
            try
                s = lgd.String;
                if ischar(s), nitems = 1; else, nitems = numel(s); end
                if nitems >= 5
                    ncol = 2;
                elseif nitems == 4
                    ncol = 2;
                else
                    ncol = 1;
                end
            catch
                ncol = 1;
            end
        end
        try, lgd.NumColumns = ncol; catch, end
    end
    drawnow;
end

function score = corner_occupancy(ax)
    % Return occupancy scores [NE NW SE SW]; lower is better.
    score = zeros(1,4);
    try
        xl = xlim(ax); yl = ylim(ax);
        xscale = get(ax,'XScale'); yscale = get(ax,'YScale');
        xdir = get(ax,'XDir'); ydir = get(ax,'YDir');
    catch
        score(:) = 1;
        return;
    end

    children = get(ax,'Children');
    series_count = 0;
    for k = 1:numel(children)
        h = children(k);
        if ~(isprop(h,'XData') && isprop(h,'YData'))
            continue;
        end
        try
            x = double(h.XData(:)); y = double(h.YData(:));
        catch
            continue;
        end
        good = isfinite(x) & isfinite(y);
        x = x(good); y = y(good);
        if isempty(x), continue; end

        % Bar objects occupy the full vertical interval from BaseValue to
        % YData, not just their top endpoints. Sample that interval so a
        % lower-corner legend is not incorrectly judged to be empty.
        if isprop(h,'BaseValue') && numel(x)==numel(y)
            try
                bv = double(h.BaseValue);
                ns = 9;
                xx = repelem(x,ns);
                yy = zeros(numel(y)*ns,1);
                for ib = 1:numel(y)
                    ii = (ib-1)*ns + (1:ns);
                    yy(ii) = linspace(bv,y(ib),ns);
                end
                x = xx(:); y = yy(:);
            catch
            end
        end

        if numel(x) > 2500
            ii = unique(round(linspace(1,numel(x),2500)));
            x = x(ii); y = y(ii);
        end

        [xn,xok] = normalize_axis_values(x,xl,xscale);
        [yn,yok] = normalize_axis_values(y,yl,yscale);
        ok = xok & yok & isfinite(xn) & isfinite(yn);
        xn = xn(ok); yn = yn(ok);
        if isempty(xn), continue; end
        if strcmpi(xdir,'reverse'), xn = 1-xn; end
        if strcmpi(ydir,'reverse'), yn = 1-yn; end

        % Legend footprint proxy: roughly 42% of width x 34% of height.
        rect = [xn>=0.58 & yn>=0.66, ... % NE
                xn<=0.42 & yn>=0.66, ... % NW
                xn>=0.58 & yn<=0.34, ... % SE
                xn<=0.42 & yn<=0.34];    % SW
        score = score + mean(rect,1);

        % Also penalize points close to the four corners, which helps when a
        % curve passes through a legend box but contributes few sampled points.
        corners = [1 1; 0 1; 1 0; 0 0];
        for c = 1:4
            d2 = (xn-corners(c,1)).^2 + (yn-corners(c,2)).^2;
            score(c) = score(c) + 0.10*mean(exp(-d2/0.075));
        end
        series_count = series_count + 1;
    end

    if series_count == 0
        score(:) = 0;
    end
end

function [vnorm,ok] = normalize_axis_values(v,lim,scale)
    ok = true(size(v));
    if strcmpi(scale,'log')
        ok = v>0 & lim(1)>0 & lim(2)>0;
        vv = nan(size(v));
        vv(ok) = log10(v(ok));
        ll = log10(lim);
    else
        vv = v;
        ll = lim;
    end
    den = ll(2)-ll(1);
    if ~isfinite(den) || abs(den) < eps
        vnorm = 0.5*ones(size(v));
    else
        vnorm = (vv-ll(1))/den;
    end
end

function export_pair(fig, out_dir, base_name, varargin)
    p = inputParser;
    addParameter(p,'PdfPadPt',18,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
    addParameter(p,'RasterPadPt',18,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
    addParameter(p,'SavePNG',true,@(x)islogical(x)||isnumeric(x));
    addParameter(p,'SaveJPEG',false,@(x)islogical(x)||isnumeric(x));
    addParameter(p,'Resolution',300,@(x)isnumeric(x)&&isscalar(x)&&x>0);
    addParameter(p,'JPEGResolution',1000,@(x)isnumeric(x)&&isscalar(x)&&x>0);
    parse(p,varargin{:});
    o = p.Results;

    if ~exist(out_dir,'dir'), mkdir(out_dir); end
    drawnow;
    set(fig,'Renderer','painters','Color','w','InvertHardcopy','off');

    % Preserve the original figure units.  print-based fallbacks use an
    % explicit physical page, while exportgraphics R2025a+ receives an
    % explicit Padding value in points.
    old_units = get(fig,'Units');
    cleanup = onCleanup(@() set(fig,'Units',old_units));
    set(fig,'Units','points');
    pos = get(fig,'Position');
    fw = max(pos(3),1); fh = max(pos(4),1);

    if logical(o.SavePNG)
        png_path = fullfile(out_dir,[char(base_name) '.png']);
        try
            exportgraphics(fig,png_path, ...
                'Resolution',o.Resolution, ...
                'BackgroundColor','white', ...
                'Units','points', ...
                'Padding',o.RasterPadPt);
        catch
            % Compatibility path for MATLAB releases without Padding.
            rpad = o.RasterPadPt;
            set(fig,'PaperUnits','points','PaperPositionMode','manual', ...
                'PaperPosition',[rpad rpad fw fh], ...
                'PaperSize',[fw+2*rpad fh+2*rpad]);
            print(fig,png_path,'-dpng',sprintf('-r%d',round(o.Resolution)));
        end
    end

    if logical(o.SaveJPEG)
        jpg_path = fullfile(out_dir,[char(base_name) '.jpg']);
        jpg_dpi = round(o.JPEGResolution);
        try
            exportgraphics(fig,jpg_path, ...
                'Resolution',jpg_dpi, ...
                'BackgroundColor','white', ...
                'Units','points', ...
                'Padding',o.RasterPadPt);
        catch
            % Compatibility path for MATLAB releases without Padding.
            rpad = o.RasterPadPt;
            set(fig,'PaperUnits','points','PaperPositionMode','manual', ...
                'PaperPosition',[rpad rpad fw fh], ...
                'PaperSize',[fw+2*rpad fh+2*rpad]);
            try
                print(fig,jpg_path,'-djpeg100',sprintf('-r%d',jpg_dpi));
            catch
                print(fig,jpg_path,'-djpeg',sprintf('-r%d',jpg_dpi));
            end
        end
    end

    % Vector PDF with a real physical page border.
    ppad = o.PdfPadPt;
    set(fig,'PaperUnits','points','PaperPositionMode','manual', ...
        'PaperPosition',[ppad ppad fw fh], ...
        'PaperSize',[fw+2*ppad fh+2*ppad]);

    pdf_path = fullfile(out_dir,[char(base_name) '.pdf']);
    try
        print(fig,pdf_path,'-dpdf','-painters',sprintf('-r%d',round(o.Resolution)));
    catch ME
        error('paper_figure_utils:PDFExportFailed', ...
              'PDF export failed for %s: %s',char(base_name),ME.message);
    end
    clear cleanup;
    set(fig,'Units',old_units);
end
