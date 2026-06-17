using Pkg
Pkg.activate(".")  # assumes current working directory is 'section-5-2'
Pkg.resolve()
Pkg.instantiate()

using TransmissionChannelAnalysis
using DataFrames, CSV
using CairoMakie  # for plotting

data = CSV.read("./data/data_lp.csv", DataFrame)
data = select(data, :newsy, :gdef, :y)

# 1. Defining the model
model = LP(data, :newsy, 4, 0:20; include_constant = true)
fit!(model)

# 2. Obtaining total effects
method = Recursive()
irfs = IRF(model, method, 20)
irfs = irfs.irfs[:, :, :]

# 3. Defining the transmission matrix
transmission_order = [:newsy, :gdef, :y]

# 4. Defining transmission channels
anticipation_channel = not_through(model, :gdef, 0:20, transmission_order)
implementation_channel = !anticipation_channel

# 5. Computing transmission effects
anticipation_effects =
    transmission(model, method, 1, anticipation_channel, transmission_order, 20)
implementation_effects =
    transmission(model, method, 1, implementation_channel, transmission_order, 20)

# 6. Visualising the effects
teffects = [anticipation_effects, implementation_effects]
channel_names = ["Anticipation", "Implementation"]

fig = Figure(; size = (800, 400));
ax1 =
    Axis(fig[1, 1]; title = "GDP", ylabel = "% of real potential GDP", xlabel = "Quarters")
plot_decomposition!(ax1, 3, irfs, teffects)
ax2 = Axis(fig[1, 2]; title = "Defense Spending", xlabel = "Quarters")
plot_decomposition!(ax2, 2, irfs, teffects)
add_decomposition_legend!(fig[2, :], channel_names)

# uncomment following line to show plot
# fig

save("./plots/ramey-anticipation-direct-measure.pdf", fig)
