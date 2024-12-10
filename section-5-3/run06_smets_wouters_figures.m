%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figures for alternative-2 channels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
clc;

horizons = [16, 32, 40, 52, 100];
for i=1:length(horizons)
    horizon = horizons(i);
    df = readtable(sprintf("output/alternative2-effects-horizons-%d", horizon));

    stackdata = [df.interest_rate df.expectations df.output_wage, df.expectations_output_wage];
    fig = figure;
    hold on
    bar(0:horizon, stackdata, 'stacked', 'EdgeColor', 'none')
    plot(0:horizon, df.total, "-*", "LineWidth", 3, color="black");
    hold off

    % title('Decomposition of Monetary Policy Shock in Smets & Wouters 2007', 'FontSize', 15)
    xlabel('Horizon', 'FontSize', 15)
    ylabel('Response of Inflation', 'FontSize', 15)
    lgd = legend( ...
        "Pure Interest Rate", ...
        "Expectations (no Output-Wage)", ...
        "Output-Wage (no Expectations)", ...
        "Expectations and Output-Wage", ...
        "Total", ...
        'FontSize', 15, ...
        'Location', 'southoutside', ...
        'Orientation', 'horizontal', ...
        'Box', 'off' ...
    );
    ax = gca;
    ax.FontSize = 15;
    ax.Box = 'on';
    ax.LineWidth = 2;
    grid on;

    width = 30;
    height = 15;
    set(fig, 'PaperSize', [width height]);
    set(fig, 'PaperPosition', [0 0 width height]);
    filename = sprintf("plots/alternative2-smets-wouters-decomposition-horizon-%d.pdf", horizon);
    saveas(fig, filename);
end
