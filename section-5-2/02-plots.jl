using DrWatson
@quickactivate "gov-cons-anticipation"
using Pkg; Pkg.instantiate()
using DataFrames, CSV, DataFramesMeta
using CairoMakie
CairoMakie.activate!()
using JLD2

output = JLD2.load(projectdir("output", "estimation-results.jld2"))
irfs = output["irfs"]
effects = output["effects"]

include(scriptsdir("decomposition-plot.jl"))

colors = Makie.wong_colors();
shock_size = 1 / irfs[1, 1, 1];
fig = Figure(;size = (1000, 400));
_ = plot_decomposition_ramey_noinstrument!(fig[1, 1], irfs .* shock_size, effects .* shock_size, 4, title="GDP", ylabel="% of real potential GDP");
_ = plot_decomposition_ramey_noinstrument!(fig[1, 2], irfs .* shock_size, effects .* shock_size, 3, title="Total Government Spending", xlabel="Quarters");
_ = plot_decomposition_ramey_noinstrument!(fig[1, 3], irfs .* shock_size, effects .* shock_size, 2, title="Defense Spending");
elements = [PolyElement(polycolor = colors[i]) for i = 1:2];
elements = vcat([[LineElement(color = :black), MarkerElement(marker = :circle, markercolor = :black)]], elements);
labels = ["Total", "Anticipation", "Implementation"];
Legend(fig[2, :], elements, labels, "Effect"; orientation = :horizontal, framevisible = false, titleposition = :left);
save(plotsdir("ramey-anticipation-direct-measure.pdf"), fig)
