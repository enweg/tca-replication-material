################################################################################
# Preliminaries
################################################################################

using Pkg; 
Pkg.add("DrWatson")  # replication manager
using DrWatson
@quickactivate "instrument-comparison"  # activating the replication environment
Pkg.instantiate()  # installing all required packages

using JLD2
using DataFrames, CSV, DataFramesMeta, TSFrames
using Random, Distributions
using LinearAlgebra
using Dates 
using CairoMakie
CairoMakie.activate!()

# loading utility functions
include(scriptsdir("utils-data-wrangling.jl"))
include(scriptsdir("svar-utils.jl"))
include(scriptsdir("svar-internal-instrument.jl"))
include(scriptsdir("svar-external-instrument.jl"));
include(scriptsdir("plot.jl"))

################################################################################
# Loading data
################################################################################

data_mckay_wolf = CSV.read(datadir("mckay_wolf_2023", "_data.csv"), DataFrame)
data_mckay_wolf = @chain data_mckay_wolf begin
    @rename begin
        :ygap = :ygap_hp
        :ad_shock = :ad
    end
    @transform :lpgdp = log.(:pgdp)
    @transform :dlpgdp = diff(:lpgdp)
    @transform :infl = 4*:dlpgdp*100
end

vardata = @chain data_mckay_wolf begin
    @rsubset 1969 <= :date <= 2007.75
    # GK instrument, RR instrument, Fed Funds Rate, Output Gap, Inflation, Commodity Prices
    @select :mp1_tc :rr_3 :ffr :ygap :infl :lpcom
end
vardata = mapcols(_replace_missing_zero, vardata)

################################################################################
# Gertler and Karadi
################################################################################

data_gk = select(vardata, :mp1_tc, :ffr, :ygap, :infl, :lpcom);
data_gk = Matrix(Float64.(data_gk));

# computing structural IRFs via internal instruments
relative_irfs_gk = internal_instrument_SVAR(data_gk, 4, 0:40; include_constant = true, include_linear_trend = true);

# getting the orthogonal IRFs / the Cholesky IRFs 
phi_tilde = orthogonal_irfs(data_gk[:, 2:end], 4, 0:40; include_constant = true, include_linear_trend = true);

# We are interested in the effect not going through the fed funds rate 
# contemporaneously. So we are interested in the transmission query 
# Q(!ffr) = Q(TRUE) - Q(ffr) according to the theory. 
# Q(TRUE) are the total impulse responses from the monetary policy shock 
# to the outcome variables
# Q(ffr) can be computed by multiplying the structural irf of the monetary policy
# shock on the ffr with the cholesky irf of a shock to the ffr on other variables, 
# and dividing this (normalising this) by the response of a cholesky ffr shock 
# on ffr itself.
irfs_transmission_gk = relative_irfs_gk[:, 1:1, :] .- relative_irfs_gk[1, 1, 1]*phi_tilde[:, 1:1, :]/phi_tilde[1, 1, 1];

# setting the shock size to 25bps
gk_total = relative_irfs_gk[:, 1:1, :] * 0.25
gk_non_contemporaneous = irfs_transmission_gk[:, 1:1, :] * 0.25
gk_contemporaneous = gk_total .- gk_non_contemporaneous;

save(projectdir("output", "gk.jld2"), Dict(
    "gk_total" => gk_total, 
    "gk_non_contemporaneous" => gk_non_contemporaneous, 
    "gk_contemporaneous" => gk_contemporaneous
))

################################################################################
# Romer and Romer
################################################################################

data_rr = select(vardata, :rr_3, :ffr, :ygap, :infl, :lpcom);
names_rr = names(data_rr)
data_rr = Matrix(Float64.(data_rr));

# structural irfs
relative_irfs_rr = internal_instrument_SVAR(data_rr, 4, 0:40; include_constant = true, include_linear_trend = true);

# cholesky irfs
phi_tilde = orthogonal_irfs(data_rr[:, 2:end], 4, 0:40; include_constant = true, include_linear_trend = true);
# computation is the same as above 
irfs_transmission_rr = relative_irfs_rr[:, 1:1, :] .- relative_irfs_rr[1, 1, 1]*phi_tilde[:, 1:1, :]/phi_tilde[1, 1, 1];

rr_total = relative_irfs_rr[:, 1:1, :] * 0.25
rr_non_contemporaneous = irfs_transmission_rr[:, 1:1, :] * 0.25
rr_contemporaneous = rr_total .- rr_non_contemporaneous;

save(projectdir("output", "rr.jld2"), Dict(
    "rr_total" => rr_total, 
    "rr_non_contemporaneous" => rr_non_contemporaneous, 
    "rr_contemporaneous" => rr_contemporaneous
))

