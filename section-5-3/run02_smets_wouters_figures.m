clear;
clc;

horizons = [16, 32, 40, 52];
for i=1:length(horizons)
    horizon = horizons(i);
    df = readtable(sprintf("output/effects-horizons-%d", horizon));

    stackdata = [df.interest_rate df.expectations df.output_wage];
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
        "Expectations (no Output or Wage)", ...
        "Output or Wage", ...
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
    filename = sprintf("plots/smets-wouters-decomposition-horizon-%d.pdf", horizon);
    saveas(fig, filename);
end
