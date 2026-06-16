
using Pkg;
Pkg.activate(".")  # assumes current working directory is 'section-5-1'
Pkg.instantiate()

using JLD2
using CairoMakie
CairoMakie.activate!()

include("./scripts/plot.jl");

################################################################################
# Gertler and Karadi
################################################################################

gk = load("output/gk.jld2")
gk_total = gk["gk_total"]
gk_non_contemporaneous = gk["gk_non_contemporaneous"]
gk_contemporaneous = gk["gk_contemporaneous"]

fig = with_theme(my_theme) do
    plot_decomposition(
        gk_total,
        gk_non_contemporaneous,
        gk_contemporaneous;
        names = ["FFR", "Output Gap", "Inflation"],
        select_vars = [1, 3],
    )
end;
save("./plots/instrument-comparison-GK.pdf", fig)

################################################################################
# Romer and Romer
################################################################################

rr = load("output/rr.jld2")
rr_total = rr["rr_total"]
rr_non_contemporaneous = rr["rr_non_contemporaneous"]
rr_contemporaneous = rr["rr_contemporaneous"]

fig = with_theme(my_theme) do
    plot_decomposition(
        rr_total,
        rr_non_contemporaneous,
        rr_contemporaneous;
        names = ["FFR", "Output Gap", "Inflation"],
        select_vars = [1, 3],
    )
end;
save("./plots/instrument-comparison-RR.pdf", fig)
