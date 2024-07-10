function plot_decomposition_ramey_noinstrument!(fig, irfs, effects, i; title="", xlabel="", ylabel="", colors=Makie.wong_colors())

    total = irfs[i, 1, :]
    anticipation = effects[i, 1, :]
    implementation = total - anticipation
    horizon = length(total) - 1
    tbl = (
        x = 0:horizon,
        total = total,
        x_barplot = vcat(collect(0:horizon), collect(0:horizon)),
        composition = vcat(anticipation, implementation),
        grp = vcat(fill(1, size(anticipation)), fill(2, size(implementation)))
    )

    ax = Axis(fig; title=title, xlabel=xlabel, ylabel=ylabel, ytickformat="{:.2f}")
    bar = barplot!(ax, tbl.x_barplot, tbl.composition; stack = tbl.grp, color = colors[tbl.grp])
    line = lines!(ax, tbl.x, tbl.total; color = :black)
    scat = scatter!(ax, tbl.x, tbl.total; color = :black)

    return fig
end
