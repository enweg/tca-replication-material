function fig = plot_irf_comparison(irf_data, varargin)
% PLOT_IRF_COMPARISON  Publication-quality 2×2 panel plot comparing three
% IRF types (total, consumption, investment) across four estimation methods.
%
% The function creates a clean 2×2 tiled panel figure with:
%   - optional confidence bands,
%   - a single shared legend placed below the panels,
%   - a global title,
%   - an optional global subtitle,
%   - and configurable y-axis grouping across panels.
%
% The global title and subtitle are left-aligned and drawn at figure level
% rather than through tiledlayout-managed title objects. This gives more
% precise control over typography, alignment, and spacing, and ensures that
% there is always sufficient room between the header block and the first
% row of panels.
%
% Y-axis grouping lets you decide which panels should share the same y-axis
% scale. The grouping is specified by a vector of length 4, ordered in the
% same way as the panels are drawn:
%
%       [panel 1, panel 2, panel 3, panel 4]
%
% corresponding to:
%
%       [top-left, top-right, bottom-left, bottom-right]
%
% Panels with the same integer are assigned to the same y-axis group and
% therefore use the same y-limits. Panels with different integers use
% independently determined y-limits. In all cases, every panel shows its
% own y-axis ticks and tick labels for readability.
%
% -------------------------------------------------------------------------
% INPUT FORMAT
% -------------------------------------------------------------------------
% irf_data : required
%   Either
%     (i) a numeric array of size [T × 3 × 4], where
%            T     = number of horizons,
%            dim 2 = IRF type:
%                      1 = Total
%                      2 = Consumption channel
%                      3 = Investment channel
%            dim 3 = estimation method / panel:
%                      1 = top-left
%                      2 = top-right
%                      3 = bottom-left
%                      4 = bottom-right
%
%   or
%
%     (ii) a struct with fields
%            .data   : [T × 3 × 4] numeric array (required)
%            .lower  : [T × 3 × 4] lower confidence band (optional)
%            .upper  : [T × 3 × 4] upper confidence band (optional)
%
% -------------------------------------------------------------------------
% OPTIONAL NAME-VALUE ARGUMENTS
% -------------------------------------------------------------------------
% 'MethodNames' : cell(1,4) or string(1,4)
%   Titles shown inside each panel.
%   Default:
%       {'Method 1','Method 2','Method 3','Method 4'}
%
% 'IRFNames' : cell(1,3) or string(1,3)
%   Legend entries for the three IRF types.
%   Default:
%       {'Total','Consumption','Investment'}
%
% 'Horizon' : numeric vector of length T
%   Horizontal axis values.
%   Example:
%       0:T-1
%   Default:
%       (1:T)'
%
% 'YLabel' : char or string scalar
%   Common y-axis label shown on the left column.
%   Default:
%       'Response (%)'
%
% 'XLabel' : char or string scalar
%   Common x-axis label shown on the bottom row.
%   Default:
%       'Horizon (quarters)'
%
% 'Title' : char or string scalar
%   Global figure title, left-aligned above the tiled layout.
%   Default:
%       ''
%
% 'Subtitle' : char or string scalar
%   Optional global subtitle, left-aligned below the main title in a smaller
%   muted font. Intended for sample notes, identification notes, scaling
%   notes, or brief methodological context.
%
%   Examples:
%       'Sample: 1998M2–2019M12'
%       'Shock normalised to a 25 bps increase in the 1-year rate'
%       'Bands denote 68% posterior credible intervals'
%
%   Default:
%       ''
%
% 'YAxisGroups' : numeric vector of length 4 containing integers
%   Defines which panels share the same y-axis scale.
%
%   Panels with the same integer are grouped together and receive the same
%   automatically computed y-limits. Panels with different integers receive
%   separate automatically computed y-limits.
%
%   Examples:
%       [1 1 1 1]   -> all four panels share one y-axis scale
%       [1 2 3 4]   -> all four panels have separate y-axis scales
%       [1 1 2 2]   -> top row shares one scale, bottom row another
%       [1 2 1 2]   -> left column shares one scale, right column another
%
%   The order is:
%       [top-left, top-right, bottom-left, bottom-right]
%
%   Default:
%       [1 1 1 1]
%
% 'YLim' : numeric vector [lo hi]
%   Manual y-axis limits applied to all panels.
%   If empty, limits are determined automatically. In that automatic case,
%   limits are computed separately within each y-axis group defined by
%   'YAxisGroups'.
%
%   Default:
%       []
%
% 'Colors' : numeric [3 × 3]
%   RGB colour matrix for the three IRF types.
%   Each row corresponds to one IRF series.
%   Default:
%       colourblind-safe palette
%
% 'LineWidth' : positive scalar
%   Line width of the IRF lines.
%   Default:
%       2.2
%
% 'BandAlpha' : scalar in [0,1]
%   Transparency level of confidence bands.
%   Only used when lower/upper bands are provided.
%   Default:
%       0.15
%
% 'ZeroLine' : logical scalar
%   If true, a horizontal zero line is drawn in each panel.
%   Default:
%       true
%
% 'MaxHorizon' : positive integer scalar
%   Truncates the plot to the first MaxHorizon periods.
%   Must satisfy MaxHorizon <= T.
%   Default:
%       T
%
% 'FigureSize' : numeric [w h]
%   Figure size in pixels.
%   Default:
%       [1100 700]
%
% -------------------------------------------------------------------------
% OUTPUT
% -------------------------------------------------------------------------
% fig : figure handle
%   Handle to the created MATLAB figure.
%
% -------------------------------------------------------------------------
% EXAMPLES
% -------------------------------------------------------------------------
% % --- Minimal call ------------------------------------------------------
% T = 20;
% data = randn(T, 3, 4) * 0.3;
% plot_irf_comparison(data);
%
% % --- Add title and subtitle -------------------------------------------
% plot_irf_comparison(data, ...
%     'Title',    'IRF Comparison across Identification Schemes', ...
%     'Subtitle', 'Sample: 1998M2–2019M12');
%
% % --- Separate y-axis scale for each panel ------------------------------
% plot_irf_comparison(data, ...
%     'YAxisGroups', [1 2 3 4], ...
%     'Title',       'IRFs by Method', ...
%     'Subtitle',    'Each panel uses its own y-axis scale');
%
% % --- Share y-axis within rows -----------------------------------------
% plot_irf_comparison(data, ...
%     'YAxisGroups', [1 1 2 2], ...
%     'Title',       'IRFs by Method', ...
%     'Subtitle',    'Top row and bottom row use different shared scales');
%
% % --- With confidence bands and custom labels --------------------------
% s.data  = data;
% s.lower = data - 0.15;
% s.upper = data + 0.15;
%
% plot_irf_comparison(s, ...
%     'MethodNames', {'OLS','2SLS','SVAR','Local Proj.'}, ...
%     'IRFNames',    {'Total','Consumption','Investment'}, ...
%     'Horizon',     0:T-1, ...
%     'Title',       'IRF Comparison across Identification Schemes', ...
%     'Subtitle',    'Bands denote 68% confidence intervals', ...
%     'YAxisGroups', [1 1 2 2], ...
%     'YLabel',      'Response (pp)');
%
% -------------------------------------------------------------------------
% NOTES
% -------------------------------------------------------------------------
% 1. The function assumes exactly three IRF series and exactly four methods.
% 2. Confidence bands, if supplied, must have the same size as the data.
% 3. Title and subtitle spacing is handled explicitly by reserving header
%    room above the tiled layout. This avoids crowding near the first row
%    of panels.
% 4. When 'YLim' is provided manually, it overrides automatic group-specific
%    limit computation and is applied to all panels.
% 5. Every panel always shows its own y-axis ticks and tick labels, even if
%    it belongs to a shared y-axis group.
%
% -------------------------------------------------------------------------
% AUTHOR   : auto-generated, revised
% REQUIRES : MATLAB R2019b or later (uses tiledlayout)
% -------------------------------------------------------------------------

    %% 1. Unpack input data

    if isstruct(irf_data)
        if ~isfield(irf_data, 'data')
            error('plot_irf_comparison:missingDataField', ...
                'When irf_data is a struct, it must contain a .data field.');
        end

        lower_band = [];
        upper_band = [];

        if isfield(irf_data, 'lower')
            lower_band = irf_data.lower;
        end
        if isfield(irf_data, 'upper')
            upper_band = irf_data.upper;
        end

        irf_data = irf_data.data;
    else
        lower_band = [];
        upper_band = [];
    end

    validateattributes(irf_data, {'numeric'}, {'real','finite','nonnan', ...
        'ndims',3,'size',[NaN,3,4]}, 'plot_irf_comparison', 'irf_data');

    T = size(irf_data, 1);

    if ~isempty(lower_band)
        validateattributes(lower_band, {'numeric'}, ...
            {'real','finite','nonnan','size',size(irf_data)}, ...
            'plot_irf_comparison', 'lower');
    end

    if ~isempty(upper_band)
        validateattributes(upper_band, {'numeric'}, ...
            {'real','finite','nonnan','size',size(irf_data)}, ...
            'plot_irf_comparison', 'upper');
    end

    %% 2. Defaults and name-value parsing

    % Colourblind-safe palette (Wong-style palette with slight tuning)
    default_colors = [
        0.00  0.45  0.70;   % blue   - Total
        0.84  0.37  0.00;   % orange - Consumption
        0.00  0.62  0.45;   % green  - Investment
    ];

    p = inputParser();
    p.FunctionName = 'plot_irf_comparison';

    p.addParameter('MethodNames', ...
        {'Method 1','Method 2','Method 3','Method 4'}, ...
        @(x) validateTextList(x, 4));

    p.addParameter('IRFNames', ...
        {'Total','Consumption','Investment'}, ...
        @(x) validateTextList(x, 3));

    p.addParameter('Horizon', (1:T)', ...
        @(x) isnumeric(x) && isvector(x) && numel(x) == T);

    p.addParameter('YLabel', 'Response (%)', @validateTextScalar);
    p.addParameter('XLabel', 'Horizon (quarters)', @validateTextScalar);
    p.addParameter('Title', '', @validateTextScalar);
    p.addParameter('Subtitle', '', @validateTextScalar);

    p.addParameter('YAxisGroups', [1 1 1 1], @validateYAxisGroups);

    p.addParameter('YLim', [], ...
        @(x) isempty(x) || (isnumeric(x) && numel(x) == 2 && x(1) < x(2)));

    p.addParameter('Colors', default_colors, ...
        @(x) isnumeric(x) && isequal(size(x), [3 3]));

    p.addParameter('LineWidth', 2.2, ...
        @(x) isnumeric(x) && isscalar(x) && x > 0);

    p.addParameter('BandAlpha', 0.15, ...
        @(x) isnumeric(x) && isscalar(x) && x >= 0 && x <= 1);

    p.addParameter('ZeroLine', true, ...
        @(x) islogical(x) && isscalar(x));

    p.addParameter('MaxHorizon', T, ...
        @(x) isnumeric(x) && isscalar(x) && x >= 1 && x <= T && ...
              x == floor(x));

    p.addParameter('FigureSize', [1100 700], ...
        @(x) isnumeric(x) && numel(x) == 2 && all(x > 0));

    p.parse(varargin{:});
    opt = p.Results;

    % Convert text inputs to character vectors for consistent downstream use
    opt.YLabel      = convertToChar(opt.YLabel);
    opt.XLabel      = convertToChar(opt.XLabel);
    opt.Title       = convertToChar(opt.Title);
    opt.Subtitle    = convertToChar(opt.Subtitle);
    opt.MethodNames = cellstr(string(opt.MethodNames));
    opt.IRFNames    = cellstr(string(opt.IRFNames));

    horizon      = opt.Horizon(:);
    clrs         = opt.Colors;
    lw           = opt.LineWidth;
    H            = opt.MaxHorizon;
    yaxis_groups = opt.YAxisGroups(:).';
    n_methods    = 4;
    n_types      = 3;

    %% 3. Apply horizon truncation

    irf_data = irf_data(1:H, :, :);
    horizon  = horizon(1:H);

    if ~isempty(lower_band)
        lower_band = lower_band(1:H, :, :);
    end
    if ~isempty(upper_band)
        upper_band = upper_band(1:H, :, :);
    end

    %% 4. Determine y-limits panel-by-panel or group-by-group

    panel_ylims = zeros(n_methods, 2);

    if isempty(opt.YLim)
        unique_groups = unique(yaxis_groups, 'stable');

        for g = unique_groups
            members = find(yaxis_groups == g);

            % Gather all values in the current y-axis group
            all_vals = irf_data(:, :, members);
            all_vals = all_vals(:);

            if ~isempty(lower_band)
                tmp = lower_band(:, :, members);
                all_vals = [all_vals; tmp(:)];
            end

            if ~isempty(upper_band)
                tmp = upper_band(:, :, members);
                all_vals = [all_vals; tmp(:)];
            end

            ylims_g = computeNiceSymmetricLimits(all_vals);
            panel_ylims(members, :) = repmat(ylims_g, numel(members), 1);
        end
    else
        panel_ylims = repmat(opt.YLim(:).', n_methods, 1);
    end

    %% 5. Create figure, reserve header space, and add global title block

    fig = figure( ...
        'Color', 'white', ...
        'Position', [100 100 opt.FigureSize(1) opt.FigureSize(2)], ...
        'PaperPositionMode', 'auto');

    hasTitle    = ~isempty(strtrim(opt.Title));
    hasSubtitle = ~isempty(strtrim(opt.Subtitle));

    % Explicitly reserve top space so the header never feels squeezed
    if hasTitle && hasSubtitle
        tlo_pos = [0.075 0.095 0.89 0.76];
    elseif hasTitle || hasSubtitle
        tlo_pos = [0.075 0.095 0.89 0.80];
    else
        tlo_pos = [0.075 0.095 0.89 0.84];
    end

    tlo = tiledlayout(fig, 2, 2, ...
        'TileSpacing', 'compact', ...
        'Padding', 'compact');

    tlo.Units    = 'normalized';
    tlo.Position = tlo_pos;

    addGlobalTitleBlock(fig, opt.Title, opt.Subtitle);

    %% 6. Draw panels

    ax = gobjects(n_methods, 1);

    for m = 1:n_methods
        ax(m) = nexttile(tlo);
        hold(ax(m), 'on');

        % Subtle background
        set(ax(m), 'Color', [0.975 0.977 0.982]);

        % Light horizontal grid
        set(ax(m), ...
            'XGrid',         'off', ...
            'YGrid',         'on', ...
            'GridColor',     [1 1 1], ...
            'GridAlpha',     0.9, ...
            'GridLineStyle', '-', ...
            'LineWidth',     0.5);

        % Zero line
        if opt.ZeroLine
            yline(ax(m), 0, ...
                'Color',     [0.5 0.5 0.5], ...
                'LineWidth', 0.9, ...
                'LineStyle', '--', ...
                'HandleVisibility', 'off');
        end

        % Confidence bands
        if ~isempty(lower_band) && ~isempty(upper_band)
            for k = 1:n_types
                lo = lower_band(:, k, m);
                hi = upper_band(:, k, m);

                fill(ax(m), ...
                    [horizon; flipud(horizon)], ...
                    [lo;      flipud(hi)], ...
                    clrs(k,:), ...
                    'FaceAlpha',        opt.BandAlpha, ...
                    'EdgeColor',        'none', ...
                    'HandleVisibility', 'off');
            end
        end

        % IRF lines
        line_styles = {'-','--',':'};
        markers     = {'none','none','none'};

        for k = 1:n_types
            plot(ax(m), ...
                horizon, irf_data(:, k, m), ...
                'Color',        clrs(k,:), ...
                'LineWidth',    lw, ...
                'LineStyle',    line_styles{k}, ...
                'Marker',       markers{k}, ...
                'DisplayName',  opt.IRFNames{k});
        end

        % Panel title inside axes
        ylims_m = panel_ylims(m, :);
        text(ax(m), ...
            horizon(1) + 0.03 * (horizon(end) - horizon(1)), ...
            ylims_m(2) - 0.06 * (ylims_m(2) - ylims_m(1)), ...
            opt.MethodNames{m}, ...
            'FontSize',   13, ...
            'FontWeight', 'bold', ...
            'FontName',   'Helvetica', ...
            'Color',      [0.20 0.20 0.20], ...
            'VerticalAlignment', 'top', ...
            'Interpreter', 'none');

        % Axes formatting
        set(ax(m), ...
            'XLim',       [horizon(1) horizon(end)], ...
            'YLim',       ylims_m, ...
            'FontSize',   12, ...
            'FontName',   'Helvetica', ...
            'TickDir',    'out', ...
            'TickLength', [0.012 0.012], ...
            'XColor',     [0.35 0.35 0.35], ...
            'YColor',     [0.35 0.35 0.35], ...
            'Box',        'off', ...
            'Layer',      'bottom');

        % Remove redundant x-tick labels only on top row
        if m == 1 || m == 2
            set(ax(m), 'XTickLabel', {});
        end

        % Keep y-axis tick labels on every panel for readability

        % Axis labels
        if m == 3 || m == 4
            xlabel(ax(m), opt.XLabel, ...
                'FontSize', 13, ...
                'FontName', 'Helvetica', ...
                'Color', [0.35 0.35 0.35], ...
                'Interpreter', 'none');
        end

        if m == 1 || m == 3
            ylabel(ax(m), opt.YLabel, ...
                'FontSize', 13, ...
                'FontName', 'Helvetica', ...
                'Color', [0.35 0.35 0.35], ...
                'Interpreter', 'none');
        end

        hold(ax(m), 'off');
    end

    %% 7. Link y-axes within groups

    unique_groups = unique(yaxis_groups, 'stable');
    for g = unique_groups
        members = find(yaxis_groups == g);
        if numel(members) > 1
            linkaxes(ax(members), 'y');
        end
    end

    %% 8. Shared legend

    line_handles = gobjects(n_types, 1);
    for k = 1:n_types
        line_handles(k) = findobj(ax(1), 'Type', 'Line', ...
            'DisplayName', opt.IRFNames{k});
    end

    leg = legend(line_handles, opt.IRFNames, ...
        'Orientation', 'horizontal', ...
        'FontSize',    12, ...
        'FontName',    'Helvetica', ...
        'TextColor',   [0.2 0.2 0.2], ...
        'Box',         'off', ...
        'Interpreter', 'none');

    leg.Layout.Tile = 'south';

    %% 9. Final polish

    for m = 1:n_methods
        ax(m).XAxis.LineWidth = 0.8;
        ax(m).YAxis.LineWidth = 0.8;
    end

end


function tf = validateTextScalar(x)
    tf = ischar(x) || (isstring(x) && isscalar(x));
end


function tf = validateTextList(x, n)
    tf = (iscellstr(x) && numel(x) == n) || ...
         (isstring(x) && numel(x) == n);
end


function tf = validateYAxisGroups(x)
    tf = isnumeric(x) && isvector(x) && numel(x) == 4 && ...
         all(isfinite(x)) && all(x == floor(x));
end


function c = convertToChar(x)
    if isstring(x)
        c = char(x);
    else
        c = x;
    end
end


function ylims = computeNiceSymmetricLimits(all_vals)
% COMPUTENICESYMMETRICLIMITS  Computes symmetric rounded y-limits.

    raw_max = max(abs(all_vals));

    if isempty(raw_max) || raw_max == 0
        ylims = [-1 1];
        return
    end

    raw_max = 1.15 * raw_max;   % breathing room
    mag     = 10^floor(log10(raw_max));
    y_abs   = ceil(raw_max / mag) * mag;
    ylims   = [-y_abs, y_abs];
end


function addGlobalTitleBlock(fig, titleText, subtitleText)
% ADDGLOBALTITLEBLOCK  Draws a left-aligned global title/subtitle block.
%
% This helper places title and subtitle at figure level rather than relying
% on tiledlayout-managed title objects. This gives better control over:
%   - left alignment,
%   - vertical spacing,
%   - font sizing,
%   - and separation from the first row of panels.

    titleText    = convertToChar(titleText);
    subtitleText = convertToChar(subtitleText);

    hasTitle    = ~isempty(strtrim(titleText));
    hasSubtitle = ~isempty(strtrim(subtitleText));

    x0 = 0.075;   % align with tiled layout left edge
    w  = 0.89;

    if hasTitle
        annotation(fig, 'textbox', [x0 0.948 w 0.036], ...
            'String', titleText, ...
            'EdgeColor', 'none', ...
            'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'middle', ...
            'FontSize', 17, ...
            'FontWeight', 'bold', ...
            'FontName', 'Helvetica', ...
            'Color', [0.14 0.14 0.14], ...
            'Interpreter', 'none');
    end

    if hasSubtitle
        if hasTitle
            y_sub = 0.912;
        else
            y_sub = 0.942;
        end

        annotation(fig, 'textbox', [x0 y_sub w 0.032], ...
            'String', subtitleText, ...
            'EdgeColor', 'none', ...
            'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'middle', ...
            'FontSize', 13, ...
            'FontWeight', 'normal', ...
            'FontName', 'Helvetica', ...
            'Color', [0.38 0.40 0.43], ...
            'Interpreter', 'none');
    end
end
