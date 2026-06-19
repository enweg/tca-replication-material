%% SCRIPT OPTIONS
clear; clc;

MAX_HORIZON = 9;

%% --- BASELINE RESULTS

load outputs/tca-vs-counterfactual/resultsBaseline;
irfs = resultsBaseline.irfsVarma; 

effectConsumptionChannel = resultsBaseline.TCAThroughConsumptionInvestmentInvariant.orderCBeforeITrue.effectsConsumption; 
effectInvestmentChannel = resultsBaseline.TCAThroughConsumptionInvestmentInvariant.orderCBeforeITrue.effectsInvestment; 

% computing the effects by subtracting the not channel from the total effect
effectConsClosestCounterfactual = resultsBaseline.TCAClosestCounterfactual.orderCBeforeITrue.effectNotConsumption; 
effectConsClosestCounterfactual = irfs(:, 3, :) - effectConsClosestCounterfactual;
effectInvClosestCounterfactual = resultsBaseline.TCAClosestCounterfactual.orderCBeforeITrue.effectNotInvestment; 
effectInvClosestCounterfactual = irfs(:, 3, :) - effectInvClosestCounterfactual; 

%% --- COUNTERFACTUAL: NO CONSUMPTION RESPONSE

load outputs/tca-vs-counterfactual/resultsCounterfactualNoConsumption; 

irfsNoConsumption = resultsCounterfactualNoConsumption.irfsVarma; 
consCounterfactual = irfs - irfsNoConsumption;

%% --- COUNTERFACTUAL: NO INVESTMENT RESPONSE

load outputs/tca-vs-counterfactual/resultsCounterfactualNoInvestment; 

irfsNoInvestment = resultsCounterfactualNoInvestment.irfsVarma;
invCounterfactual = irfs - irfsNoInvestment; 

%% --- PLOTTING

total = irfs(3, 3, :);
total = total(:);

consChannel = effectConsumptionChannel(3, 1, :);
consChannel = consChannel(:);

invChannel = effectInvestmentChannel(3, 1, :);
invChannel = invChannel(:); 

consClosestCounterfactual = effectConsClosestCounterfactual(3, 1, :); 
consClosestCounterfactual = consClosestCounterfactual(:); 

invClosestCounterfactual = effectInvClosestCounterfactual(3, 1, :);
invClosestCounterfactual = invClosestCounterfactual(:);

consCounterfactual = consCounterfactual(3, 3, :); 
consCounterfactual = consCounterfactual(:);

invCounterfactual = invCounterfactual(3, 3, :); 
invCounterfactual = invCounterfactual(:);

horizons = 0:(length(total(:))-1);
horizons = horizons';


blue = [0.0000 0.4470 0.6980];    % Wong blue
yellow = [0.9020 0.6240 0.0000];  % Wong orange/yellow
black = [0 0 0];

% fig = figure('Color', 'white', 'Position', [100 100 1400 620]);
fig = figure('Color', 'white', 'Position', [100 100 1400 620], 'Renderer', 'painters');
tlo = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'loose');
tlo.Units = 'normalized';
tlo.Position = [0.06 0.13 0.89 0.80];

[hTotal, hConsumption, hInvestment] = plotPanel(tlo, horizons, total, consCounterfactual, invCounterfactual, ...
    "Counterfactual", MAX_HORIZON, blue, yellow, black);
plotPanel(tlo, horizons, total, consChannel, invChannel, ...
    "TCA", MAX_HORIZON, blue, yellow, black);
plotPanel(tlo, horizons, total, consClosestCounterfactual, invClosestCounterfactual, ...
    "TCA Closest to Counterfactual", MAX_HORIZON, blue, yellow, black);

% xlabel(tlo, "Quarters", 'FontSize', 16, 'FontWeight', 'bold');
% ylabel(tlo, "% of real potential GDP", 'FontSize', 16, 'FontWeight', 'bold');

leg = legend([hTotal hConsumption hInvestment], {'Total', 'Consumption channel', 'Investment channel'}, ...
    'Orientation', 'horizontal', ...
    'Box', 'off', ...
    'FontSize', 15, ...
    'Location', 'southoutside');
leg.Layout.Tile = 'south';

drawnow;
set(fig, ...
    'PaperUnits', 'inches', ...
    'PaperPosition', [0 0 14 6.2], ...
    'PaperSize', [14 6.2], ...
    'PaperPositionMode', 'manual', ...
    'InvertHardcopy', 'off');
print(fig, fullfile('outputs/tca-vs-counterfactual', 'sec2-tca-vs-counterfactual.pdf'), ...
    '-dpdf', '-painters');

function [hTotal, hConsumption, hInvestment] = plotPanel(tlo, horizons, total, consumptionChannel, investmentChannel, ...
    panelTitle, maxHorizon, blue, yellow, black)

    ax = nexttile(tlo);
    plotRange = 1:maxHorizon;

    hold(ax, 'on');
    plot(ax, horizons(plotRange), zeros(size(horizons(plotRange))), ...
        'LineStyle', ':', ...
        'Color', [0.25 0.25 0.25], ...
        'LineWidth', 2.0, ...
        'MarkerSize', 1.5, ...
        'HandleVisibility', 'off');
    hTotal = plot(ax, horizons(plotRange), total(plotRange), ...
        'Color', black, ...
        'LineWidth', 1.8, ...
        'Marker', 'o', ...
        'MarkerFaceColor', black, ...
        'MarkerSize', 4.5, ...
        'DisplayName', 'Total');  % total
    hConsumption = plot(ax, horizons(plotRange), consumptionChannel(plotRange), ...
        'Color', blue, ...
        'LineWidth', 2.4, ...
        'DisplayName', 'Consumption channel');  % consumption channel
    hInvestment = plot(ax, horizons(plotRange), investmentChannel(plotRange), ...
        'Color', yellow, ...
        'LineWidth', 2.4, ...
        'DisplayName', 'Investment channel');  % investment channel
    hold(ax, 'off');

    title(ax, panelTitle, 'FontSize', 16, 'FontWeight', 'bold');
    set(ax, ...
        'Box', 'on', ...
        'FontSize', 14, ...
        'GridAlpha', 0.25, ...
        'Layer', 'top', ...
        'LineWidth', 0.8, ...
        'XGrid', 'on', ...
        'YGrid', 'on');
    xlim(ax, [horizons(plotRange(1)) horizons(plotRange(end))]);
    ylim(ax, [-0.6 0.2]);
    yticks(ax, -0.6:0.2:0.2);
end
