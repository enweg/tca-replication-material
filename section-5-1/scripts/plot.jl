function plot_decomposition(
        total, 
        non_contemporanous, 
        contemporanous; 
        plot_width = 600, 
        plot_height = 400, 
        names = nothing, 
        select_vars = 1:size(total, 1)
)
    
    colors = Makie.wong_colors()
    k = size(total, 1)
    horizons = 0:(size(total, 3) - 1)
    
    fig = Figure(; resolution = (length(select_vars)*plot_width, plot_height))
    for (j, i) in enumerate(select_vars)
        tbl = (;
            x = horizons, 
            total = total[i, 1, :], 
            x_bar = vcat(horizons, horizons),
            decomposition = vcat(contemporanous[i, 1, :], non_contemporanous[i, 1, :]),
            grp = vcat(fill(1, length(horizons)), fill(2, length(horizons)))
        )        
        title = !isnothing(names) ? names[i] : ""
        ax = Axis(fig[1, j]; title = title, ytickformat="{:.2f}")
        barplot!(ax, tbl.x_bar, tbl.decomposition; stack = tbl.grp, color = colors[tbl.grp])
        lines!(ax, tbl.x, tbl.total; color = :black)
        scatter!(ax, tbl.x, tbl.total; color = :black)
        xlims!(ax, -1, 41)
    end
    elements = [PolyElement(polycolor = colors[i]) for i in 1:2]
    elements = vcat([[LineElement(color = :black), MarkerElement(marker = :circle, color = :black)]], elements)
    Legend(
        fig[2, :], 
        elements, 
        ["Total", "Contemporaneous", "Non-Contemporaneous"], 
        "Effect"; 
        orientation = :horizontal
    )
    return fig
end

const my_theme = Theme(
    Axis = (
        xlabelsize = 25,
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
