projectdir(parts...) = normpath(joinpath(@__DIR__, parts...))

using CSV, DataFrames
using CairoMakie
using Colors
using Makie

my_theme = Theme(
    Axis = (
        xlabelsize = 25,
        ylabelsize = 25,
        xticklabelsize = 25,
        yticklabelsize = 25,
        titlesize = 25,
    ), 
    Label = (
        fontsize = 25,
    ), 
    Lines = (
        linewidth = 2,
    ), 
    Legend = (
        labelsize = 25,
        titleposition = :left, 
        titlesize = 25,
        framevisible = false
    )
);

################################################################################
# Decomposition Figure (original ordering)
################################################################################

horizon = 40
effects = CSV.read(projectdir("output", "effects-horizons-$horizon.csv"), DataFrame)
first(effects, 5)

total = effects[!, :total]
interest_channel = effects[!, :interest_rate]
expectations_channel = effects[!, :expectations]
output_wage_channel = effects[!, :output_wage]

tbl = (
    x = 0:horizon,
    total = total,
    x_barplot = vcat(collect(0:horizon), collect(0:horizon), collect(0:horizon)),
    composition = vcat(interest_channel, expectations_channel, output_wage_channel),
    grp = vcat(fill(1, horizon+1), fill(2, horizon+1), fill(3, horizon+1))
)

fig = with_theme(my_theme) do
    colors = Makie.wong_colors()
    colors = [colors[1], RGBA{Float32}(217/255, 83/255, 26/255, 1.0f0), colors[2]]
    fig = Figure(;size=(1200, 500));
    ax = Axis(fig[1,1]; title="", xlabel="Horizon", ylabel="Response of Inflation", ytickformat="{:.2f}");
    bar = barplot!(ax, tbl.x_barplot, tbl.composition; stack = tbl.grp, color = colors[tbl.grp]);
    line = lines!(ax, tbl.x, tbl.total; color = :black, linewidth=3);
    scat = scatter!(ax, tbl.x, tbl.total; color = :black, markersize=15);
    elements = [PolyElement(polycolor = colors[i]) for i = 1:3];
    elements = vcat([[LineElement(color = :black, linewidth=3), MarkerElement(marker = :circle, markercolor = :black, markersize=15)]], elements);
    labels = ["Total", "Pure Interest Rate", "Expectations (no Output or Wage)", "Output or Wage"];
    Legend(fig[2, :], elements, labels, "Effect"; orientation = :horizontal, framevisible = false, titleposition = :left);
    fig
end;

save(projectdir("plots", "smets-wouters-decomposition-horizon-$horizon-julia.pdf"), fig)

################################################################################
# Decomposition Figure (alternative ordering)
################################################################################

horizon = 40
effects = CSV.read(projectdir("output", "alternative-effects-horizons-$horizon.csv"), DataFrame)
first(effects, 5)

total = effects[!, :total]
interest_channel = effects[!, :interest_rate]
expectations_channel = effects[!, :expectations]
output_wage_channel = effects[!, :output_wage]

tbl = (
    x = 0:horizon,
    total = total,
    x_barplot = vcat(collect(0:horizon), collect(0:horizon), collect(0:horizon)),
    composition = vcat(interest_channel, expectations_channel, output_wage_channel),
    grp = vcat(fill(1, horizon+1), fill(2, horizon+1), fill(3, horizon+1))
)

fig = with_theme(my_theme) do
    colors = Makie.wong_colors()
    colors = [colors[1], RGBA{Float32}(217/255, 83/255, 26/255, 1.0f0), colors[2]]
    fig = Figure(;size=(1200, 500));
    ax = Axis(fig[1,1]; title="", xlabel="Horizon", ylabel="Response of Inflation", ytickformat="{:.2f}");
    bar = barplot!(ax, tbl.x_barplot, tbl.composition; stack = tbl.grp, color = colors[tbl.grp]);
    line = lines!(ax, tbl.x, tbl.total; color = :black, linewidth=3);
    scat = scatter!(ax, tbl.x, tbl.total; color = :black, markersize=15);
    elements = [PolyElement(polycolor = colors[i]) for i = 1:3];
    elements = vcat([[LineElement(color = :black, linewidth=3), MarkerElement(marker = :circle, markercolor = :black, markersize=15)]], elements);
    labels = ["Total", "Pure Interest Rate", "Expectations", "Output or Wage (No Expectations)"];
    Legend(fig[2, :], elements, labels, "Effect"; orientation = :horizontal, framevisible = false, titleposition = :left);
    fig
end;

save(projectdir("plots", "smets-wouters-alternative-decomposition-horizon-$horizon-julia.pdf"), fig)

################################################################################
# All total IRFs
################################################################################

horizon = 40
irfs = CSV.read(projectdir("output", "irfs-horizons-$horizon.csv"), DataFrame)
first(irfs, 5)

fig = with_theme(my_theme) do
    horizons = 0:horizon
    fig = Figure(;size=(1500, 800));
    g1 = GridLayout(fig[1, 1], alignmode = Outside(15))
    g2 = GridLayout(fig[2, 1], alignmode = Outside(15))
    box = Box(fig[2, 1], color = (:black, 0.1), strokecolor = :transparent)
    Makie.translate!(box.blockscene, 0, 0, -100)
    Label(fig[2, 1, Top()], text="Output-Wage System (y)", font=:bold, padding = (0, 0, 5, 0), fontsize=25)

    ax = Axis(g1[1,1]; title="Interest Rate (r)")
    lines!(ax, horizons, irfs[!, :robs]; color=:black)
    hlines!(ax, [0]; color=:black, linestyle=:dash)

    ax = Axis(g1[1,2]; title="Inflation (p)")
    lines!(ax, horizons, irfs[!, :pinfobs]; color=:black)
    hlines!(ax, [0]; color=:black, linestyle=:dash)

    ax = Axis(g1[1,3]; title="Expected Inflation (e)")
    lines!(ax, horizons, irfs[!, :piexp]; color=:black)
    hlines!(ax, [0]; color=:black, linestyle=:dash)

    ax = Axis(g2[1,1]; title="Consumption")
    lines!(ax, horizons, irfs[!, :dc]; color=:black)
    hlines!(ax, [0]; color=:black, linestyle=:dash)

    ax = Axis(g2[1,2]; title="Investment")
    lines!(ax, horizons, irfs[!, :dinve]; color=:black)
    hlines!(ax, [0]; color=:black, linestyle=:dash)

    ax = Axis(g2[1,3]; title="Output")
    lines!(ax, horizons, irfs[!, :dy]; color=:black)
    hlines!(ax, [0]; color=:black, linestyle=:dash)

    ax = Axis(g2[1,4]; title="Wages")
    lines!(ax, horizons, irfs[!, :dw]; color=:black)
    hlines!(ax, [0]; color=:black, linestyle=:dash)

    fig
end

fig
save(projectdir("plots", "smets-wouters-irfs-horizon-$horizon-julia.pdf"), fig)
