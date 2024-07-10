using DrWatson
@quickactivate "gov-cons-anticipation"
using DataFrames, CSV, DataFramesMeta
using StatsBase, LinearAlgebra
using JLD2

using Pkg
Pkg.add(url="https://github.com/enweg/TransmissionChannelAnalysis.jl.git#b03801a");
using TransmissionChannelAnalysis

include(scriptsdir("utilities.jl"))

################################################################################
# Estimation by Local Projections
################################################################################
include(scriptsdir("lp-irfs.jl"))

data_lp = CSV.read(datadir("data_lp.csv"), DataFrame)

k = size(data_lp, 2)
order = [:newsy, :gdef, :g, :y]
data_lp = select(data_lp, order...)
horizon = 20
horizons = 0:horizon
irfs = zeros(k, k, length(horizons))
for (i, shock) in enumerate(order)
    irfs_tmp, _ = _lp_irf(data_lp, shock, order[1:(i-1)], 4; include_constant=true)
    irfs[:, i, :] = irfs_tmp[:, 1, :]
end

irfs_stacked = TransmissionChannelAnalysis.to_transmission_irfs(irfs)
# same since we already use a Cholesky identification
irfs_ortho_stacked = TransmissionChannelAnalysis.to_transmission_irfs(irfs)

################################################################################
# Computing transmission effects
################################################################################

effects = fill(NaN, size(irfs))
k = size(irfs, 2)
not_until = 20
for h = 0:min((size(irfs, 3)-1), not_until)
    s = join(["!x$(2+(i*k))" for i=0:h], " & ")
    # @show h, s
    cond = make_condition(s)
    effect = transmission(1, irfs_stacked, irfs_ortho_stacked, cond; method=:irfs)
    for i in 1:size(effects, 1)
        effects[i, 1, h+1] = effect[i+(h*size(irfs, 2))]
        if h == not_until
            effects[i, 1, (h+1):end] = effect[(i+(h*k)):k:end]
        end
    end
end

################################################################################
# Saving results
################################################################################

JLD2.save(projectdir("output", "estimation-results.jld2"), Dict(
    "irfs" => irfs, 
    "effects" => effects
))


