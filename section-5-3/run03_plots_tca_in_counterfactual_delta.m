clear; clc;
% load outputs/tca-in-counterfactual/resultsDelta-from0p01-by0p01-to0p99.mat;
load outputs/tca-in-counterfactual/resultsDelta-from0p001-by0p001-to0p05.mat;

deltaValues = results.deltaValues; 
inflationDecompositions = results.inflationDecompostions; 

% FILTERING: comment if all should be plotted
minDelta = 0.0125; 
maxDelta = 0.05; 
idxInside = (results.deltaValues >= minDelta) & (results.deltaValues <= maxDelta); 
deltaValues = deltaValues(idxInside); 
inflationDecompositions = inflationDecompositions(:, :, idxInside);

%% --- PLOTTING


figure('Color', 'w', 'Units', 'inches', 'Position', [1, 1, 18, 8.4]);

tl = tiledlayout(2, 3, ...
    'TileSpacing', 'loose', ...
    'Padding', 'loose');

colIdxs = [1, 2, 3, 4, 5];
plotTitles = {'Consumption Response', 'Investment Response', 'Inflation Response', 'Consumption Channel', 'Investment Channel'};
ax = gobjects(1, numel(plotTitles));

for j = 1:numel(plotTitles)
    ax(j) = nexttile;
    hold(ax(j), 'on');
    plotPanel(ax(j), inflationDecompositions, deltaValues, colIdxs(j), plotTitles{j});
    title(ax(j), plotTitles{j}, ...
        'FontSize', 16, ...
        'FontWeight', 'bold');
    hold(ax(j), 'off');
end

xlabel(tl, 'Horizon', ...
    'Interpreter', 'latex', ...
    'FontSize', 15, ...
    'FontWeight', 'bold');
% ylabel(tl, 'Response', ...
%     'Interpreter', 'latex', ...
%     'FontSize', 15, ...
%     'FontWeight', 'bold');

% colour bar
wongBlue   = [0, 114, 178] / 255;
wongYellow = [240, 228, 66] / 255;

nColors = 256;
alpha = linspace(0, 1, nColors)';

% low detla = blue, high delta = yellow
colorbarMap = (1 - alpha) .* wongBlue + alpha .* wongYellow;
colormap(gcf, colorbarMap);

% Tell MATLAB what numerical range the colour scale represents
for j = 1:numel(ax)
    clim(ax(j), [min(deltaValues), max(deltaValues)]);
end

cb = colorbar(ax(end));
cb.Layout.Tile = 'east';
cb.Label.Interpreter = 'latex';
cb.Label.String = '$\delta$';
cb.Label.FontWeight = 'bold';
cb.Ticks = linspace(min(deltaValues), max(deltaValues), 5);
cb.TickLabels = compose('%.1f%%', cb.Ticks * 100);
cb.TickDirection = 'out';
cb.Box = 'off';
cb.FontName = 'Helvetica';
cb.FontSize = 14;
cb.Label.FontSize = 20;

%% --- SAVING

fig = gcf;
drawnow;
set(fig, ...
    'PaperUnits', 'inches', ...
    'PaperPosition', [0 0 18 8.4], ...
    'PaperSize', [18.6 8.4], ...
    'PaperPositionMode', 'manual', ...
    'InvertHardcopy', 'off');
print(fig, fullfile('outputs/tca-in-counterfactual', 'tca-in-counterfactual-delta.pdf'), ...
    '-dpdf', '-painters');

%% --- FUNCTIONS

function plotPanel(ax, inflationDecompositions, deltaValues, colIdx, panelTitle)
    nHorizons = size(inflationDecompositions, 1);
    horizons = 0:(nHorizons-1);

    [sortedDeltaValues, idx] = sort(deltaValues, 'ascend');

    wongBlue   = [0, 114, 178] / 255;
    wongYellow = [240, 228, 66] / 255;

    nLines = numel(sortedDeltaValues);
    alpha = linspace(0, 1, nLines)';
    lineColors = (1 - alpha) .* wongBlue + alpha .* wongYellow;

    for k = 1:nLines
        i = idx(k);

        plot(ax, horizons, inflationDecompositions(:, colIdx, i), ...
            'Color', lineColors(k, :), ...
            'LineWidth', 1.15);
    end

    xlim(ax, [horizons(1), horizons(end)]);
    yline(ax, 0, '-', ...
        'Color', [0.65, 0.65, 0.65], ...
        'LineWidth', 0.6, ...
        'HandleVisibility', 'off');
    grid(ax, 'on');
    box(ax, 'off');
    set(ax, ...
        'Color', 'w', ...
        'FontName', 'Helvetica', ...
        'FontSize', 14, ...
        'FontWeight', 'bold', ...
        'GridAlpha', 0.12, ...
        'GridColor', [0, 0, 0], ...
        'Layer', 'top', ...
        'LineWidth', 0.8, ...
        'TickDir', 'out', ...
        'TickLabelInterpreter', 'latex', ...
        'XMinorTick', 'off', ...
        'YMinorTick', 'off');
end
