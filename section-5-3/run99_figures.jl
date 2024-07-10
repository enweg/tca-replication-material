projectdir(parts...) = normpath(joinpath(@__DIR__, parts...))

using CSV, DataFrames
using CairoMakie
using Colors

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
