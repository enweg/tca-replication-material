################################################################################
# REPLICATION INFORMATION
# 1. Check the current path by running line XXX. 
# 2. If the current path is different from the folder `section-5-3`, then 
#    uncomment lines 13-14, adjust the path in quotation marks on line 13 and 
#    run the line. 
# 3. Adjust the `horizon` on line 15 to the horizon for which the 
#    decomposition was computed. 
# 4. No adjustments are needed for the rest of the code. 
################################################################################

pwd()
# cd("./section-5-3/")
# pwd()
horizon = 20

################################################################################
# NO ADJUSTMENTS ARE NEEDED BELOW THIS LINE
################################################################################

projectdir(parts...) = normpath(joinpath(@__DIR__, parts...))
using Pkg;
Pkg.activate(".");  # assumes working directory is section-5-3
Pkg.instantiate();
using CSV, DataFrames
using CairoMakie
using Colors
using Makie

my_theme = Theme(
    Axis=(
        xlabelsize=25,
        ylabelsize=25,
        xticklabelsize=25,
        yticklabelsize=25,
        titlesize=25,
    ),
    Label=(
        fontsize=25,
    ),
    Lines=(
        linewidth=2,
    ),
    Legend=(
        labelsize=25,
        titleposition=:left,
        titlesize=25,
        framevisible=false
    )
);

################################################################################
# Decomposition Figure (original ordering)
################################################################################

effects = CSV.read(projectdir("output", "effects-horizons-$horizon.csv"), DataFrame)
first(effects, 5)

total = effects[!, :total]
through_wages = effects[!, :w]
not_through_wages = effects[!, :not_w]

tbl = (
    x=0:horizon,
    total=total,
    x_barplot=vcat(collect(0:horizon), collect(0:horizon)),
    composition=vcat(not_through_wages, through_wages),
    grp=vcat(fill(1, horizon + 1), fill(2, horizon + 1))
)

fig = with_theme(my_theme) do
    colors = Makie.wong_colors()
    colors = [colors[1], RGBA{Float32}(217 / 255, 83 / 255, 26 / 255, 1.0f0), colors[2]]
    fig = Figure(; size=(1200, 500))
    ax = Axis(fig[1, 1]; title="", xlabel="Horizon", ylabel="Response of Inflation", ytickformat="{:.2f}")
    ylims!(ax, -0.055, 0.02)
    bar = barplot!(ax, tbl.x_barplot, tbl.composition; stack=tbl.grp, color=colors[tbl.grp])
    line = lines!(ax, tbl.x, tbl.total; color=:black, linewidth=3)
    scat = scatter!(ax, tbl.x, tbl.total; color=:black, markersize=15)
    elements = [PolyElement(polycolor=colors[i]) for i = 1:2]
    elements = vcat([[LineElement(color=:black, linewidth=3), MarkerElement(marker=:circle, markercolor=:black, markersize=15)]], elements)
    labels = ["Total", "Demand Channel", "Wage Channel"]
    Legend(fig[2, :], elements, labels, "Effect"; orientation=:horizontal, framevisible=false, titleposition=:left)
    fig
end;

save(projectdir("plots", "smets-wouters-decomposition-horizon-$horizon-julia.pdf"), fig)


################################################################################
# Decomposition Figure (alternative ordering)
################################################################################

effects = CSV.read(projectdir("output", "effects-horizons-$horizon-alternative.csv"), DataFrame)

total = effects[!, :total]
through_wages = effects[!, :w]
not_through_wages = effects[!, :not_w]

tbl = (
    x=0:horizon,
    total=total,
    x_barplot=vcat(collect(0:horizon), collect(0:horizon)),
    composition=vcat(not_through_wages, through_wages),
    grp=vcat(fill(1, horizon + 1), fill(2, horizon + 1))
)

fig = with_theme(my_theme) do
    colors = Makie.wong_colors()
    colors = [colors[1], RGBA{Float32}(217 / 255, 83 / 255, 26 / 255, 1.0f0), colors[2]]
    fig = Figure(; size=(1200, 500))
    ax = Axis(fig[1, 1]; title="", xlabel="Horizon", ylabel="Response of Inflation", ytickformat="{:.2f}")
    ylims!(ax, -0.055, 0.02)
    bar = barplot!(ax, tbl.x_barplot, tbl.composition; stack=tbl.grp, color=colors[tbl.grp])
    line = lines!(ax, tbl.x, tbl.total; color=:black, linewidth=3)
    scat = scatter!(ax, tbl.x, tbl.total; color=:black, markersize=15)
    elements = [PolyElement(polycolor=colors[i]) for i = 1:2]
    elements = vcat([[LineElement(color=:black, linewidth=3), MarkerElement(marker=:circle, markercolor=:black, markersize=15)]], elements)
    labels = ["Total", "Demand Channel", "Wage Channel"]
    Legend(fig[2, :], elements, labels, "Effect"; orientation=:horizontal, framevisible=false, titleposition=:left)
    fig
end;

save(projectdir("plots", "smets-wouters-decomposition-horizon-$horizon-julia-alternative.pdf"), fig)

################################################################################
# The differences between the actual decompositions are very small. 
# The next graph makes these differences more clear
################################################################################

effects_orig = CSV.read(projectdir("output", "effects-horizons-$horizon.csv"), DataFrame)
effects_alt = CSV.read(projectdir("output", "effects-horizons-$horizon-alternative.csv"), DataFrame)

total = effects[!, :total]
through_wages_orig = effects_orig[!, :w]
not_through_wages_orig = effects_orig[!, :not_w]
through_wages_alt = effects_alt[!, :w]
not_through_wages_alt = effects_alt[!, :not_w]

tbl = (
    x=0:horizon,
    total=total,
    x_barplot=vcat(collect(0:horizon), collect(0:horizon), collect(0:horizon), collect(0:horizon)),
    y_barplot=vcat(not_through_wages_orig, through_wages_orig, not_through_wages_alt, through_wages_alt),
    grp_stack=vcat(fill(1, horizon + 1), fill(2, horizon + 1), fill(1, horizon + 1), fill(2, horizon + 1)),
    grp_dodge=vcat(fill(1, horizon + 1), fill(1, horizon + 1), fill(2, horizon + 1), fill(2, horizon + 1)),
    grp=vcat(fill(1, horizon + 1), fill(2, horizon + 1), fill(3, horizon + 1), fill(4, horizon + 1))
)

mix_colours(col1, col2, gamma) = gamma * col1 + (1 - gamma) * col2
add_black_tint(col1; black=Colors.RGBA(0.0f0, 0.0f0, 0.0f0, 1.0f0), gamma=0.5) = mix_colours(black, col1, gamma)

colors_orig = [
    RGBA{Float32}(0.0f0, 0.44705883f0, 0.69803923f0, 1.0f0),
    RGBA{Float32}(217 / 255, 83 / 255, 26 / 255, 1.0f0)
]
colors_alt = mix_colours.(colors_orig, Colors.RGBA(1.0f0, 1.0f0, 1.0f0, 1.0f0), 0.5)
colors = vcat(colors_orig, colors_alt)

fig = with_theme(my_theme) do
    fig = Figure(; size=(1200, 500));
    ax = Axis(fig[1, 1]; title="", xlabel="Horizon", ylabel="Response of Inflation", ytickformat="{:.2f}");
    ylims!(ax, -0.055, 0.02);
    bar = barplot!(ax, tbl.x_barplot, tbl.y_barplot; stack=tbl.grp_stack, dodge=tbl.grp_dodge, color=colors[tbl.grp], width=(1 ./ tbl.grp_dodge),);
    line = lines!(ax, tbl.x, tbl.total; color=:black, linewidth=3);
    scat = scatter!(ax, tbl.x, tbl.total; color=:black, markersize=15);
    elements = [PolyElement(polycolor=col) for col in colors];
    elements = vcat([[LineElement(color=:black, linewidth=3), MarkerElement(marker=:circle, markercolor=:black, markersize=15)]], elements);
    labels = ["Total", "Demand Channel", "Wage Channel", "Demand Channel", "Wage Channel"];
    Legend(fig[2, :], [[elements[1]], elements[2:3], elements[4:end]], [[labels[1]], labels[2:3], labels[4:end]], ["", "first-round\nordering", "second-round\nordering"]; orientation=:horizontal, framevisible=false, titleposition=:left, nbanks=2)
    fig
end;

save(projectdir("plots", "smets-wouters-decomposition-horizon-$horizon-julia-both.pdf"), fig)

################################################################################
# All total IRFs
################################################################################

irfs = CSV.read(projectdir("output", "irfs-horizons-$horizon.csv"), DataFrame)
first(irfs, 5)

fig = with_theme(my_theme) do
    horizons = 0:horizon
    fig = Figure(; size=(1500, 800))
    g1 = GridLayout(fig[1, 1], alignmode=Outside(15))
    g2 = GridLayout(fig[2, 1], alignmode=Outside(15))
    box = Box(fig[2, 1], color=(:black, 0.1), strokecolor=:transparent)
    Makie.translate!(box.blockscene, 0, 0, -100)
    Label(fig[2, 1, Top()], text="Output System (y)", font=:bold, padding=(0, 0, 5, 0), fontsize=25)

    ax = Axis(g1[1, 1]; title="Interest Rate (r)")
    lines!(ax, horizons, irfs[!, :robs]; color=:black)
    hlines!(ax, [0]; color=:black, linestyle=:dash)

    ax = Axis(g1[1, 2]; title="Inflation (p)")
    lines!(ax, horizons, irfs[!, :pinfobs]; color=:black)
    hlines!(ax, [0]; color=:black, linestyle=:dash)

    ax = Axis(g1[1, 3]; title="Wages (w)")
    lines!(ax, horizons, irfs[!, :dw]; color=:black)
    hlines!(ax, [0]; color=:black, linestyle=:dash)

    ax = Axis(g2[1, 1]; title="Consumption")
    lines!(ax, horizons, irfs[!, :dc]; color=:black)
    hlines!(ax, [0]; color=:black, linestyle=:dash)

    ax = Axis(g2[1, 2]; title="Investment")
    lines!(ax, horizons, irfs[!, :dinve]; color=:black)
    hlines!(ax, [0]; color=:black, linestyle=:dash)

    ax = Axis(g2[1, 3]; title="Output")
    lines!(ax, horizons, irfs[!, :dy]; color=:black)
    hlines!(ax, [0]; color=:black, linestyle=:dash)

    ax = Axis(g2[1, 4]; title="Hours Worked")
    lines!(ax, horizons, irfs[!, :labobs]; color=:black)
    hlines!(ax, [0]; color=:black, linestyle=:dash)

    fig
end;

save(projectdir("plots", "smets-wouters-irfs-horizon-$horizon-julia.pdf"), fig)
