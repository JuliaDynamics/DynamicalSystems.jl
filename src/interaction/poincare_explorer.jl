using AbstractPlotting, GLMakie
using StatsBase
using StatsMakie
using Colors: color

"""
    plot_sim(sim; colors=axes(sim, 1), idxs=[1,2])

This function considers the `sim` variable as a vector of arrays and creates
a 2D scatter plot where each element of `sim` is considered as a separate series.
In order to change the transparency(α) of a given plot a vector of `Observable`s is used.
Each series has an individual plot and an entry in the `series_alpha` vector.
The function also adds a slider to the scatter plot, which controls the size of
the points. The function returns the `scene` and `series_alpha`.
"""

function plot_sim(sim; colors=axes(sim, 1), idxs=[1,2])
    ui, ms = AbstractPlotting.textslider(range(0.001, stop=1., length=1000), "scale", start=0.05)
    data = Scene(resolution=(1000, 1000))
    colormap = to_colormap(:viridis, size(sim, 1))
    get_color(i) = AbstractPlotting.interpolated_getindex(colormap, colors[i], extrema(colors))
    series_alpha = map(eachindex(sim)) do i
        simᵢ = sim[i]
        alpha = Node(1.0)
        if length(sim[i]) ≠ 0
            cmap = lift(α-> RGBAf0.(color.(fill(get_color(i), size(simᵢ, 1))), α), alpha)
            scatter!(data, [Point2f0(simᵢ[i, idxs[1]], simᵢ[i, idxs[2]]) for i in axes(simᵢ, 1)],
            colormap=cmap, color=fill(colors[i], size(simᵢ, 1)), markersize=ms)
        end
        alpha
    end

    scene = Scene()

    hbox(ui, data, parent=scene)
    return scene, series_alpha
end

"""
    plot_hist(hist)

Plot a histogram where the transparency(α) of each bin can be changed
and return the scene together with the αs.
"""
function plot_hist(hist)
    cmap = to_colormap(:viridis, length(hist.weights))
    hist_α = [Node(1.) for i in cmap]
    bincolor(αs...) = RGBAf0.(color.(cmap), αs)
    colors = lift(bincolor, hist_α...)
    hist_sc = plot(hist, color=colors)

    return hist_sc, hist_α
end

"""
    change_α(series_alpha, idxs, α=0.005)

Given a vector of `Observable`s that represent the αs for some series, change
the elements with indices given by `idxs` to the value `α`.
This can be used to hide some series by using a low α value (default).
To restore the initial color, use `α = 1`.
"""
function change_α(series_alpha, idxs, α=0.005)
    foreach(i-> series_alpha[i][] = α, idxs)
end

"""
    get_series_idx(selected_plot, scene)

Get the index of the `selected_plot` in `scene`.
"""
function get_series_idx(selected_plot, scene)
    # TODO: There is probably a better or more efficient way of doing this.
    plot_idx = findfirst(map(p->selected_plot === p, scene.plots))
    # println("scatter ", plot_idx)

    plot_idx
end

"""
    setup_click(scene, idx=1)

Give a `scene` return a `Node` that listens to left clicks inside the scene.
The `idx` argument is used to index the tuple `(plt, click_idx)` which gives
the selected plot and the index of the selected element in the plot.
"""
function setup_click(scene, idx=1)
    selection = Node{Any}(0)
    on(scene.events.mousebuttons) do buttons
        if ispressed(scene, Mouse.left) && AbstractPlotting.is_mouseinside(scene)
            plt, click_idx = AbstractPlotting.mouse_selection(scene)
            selection[] = (plt, click_idx)[idx]
        end
    end
    return selection
end

"""
    bin_with_val(val, hist)

Get the index of the bin in the histogram(`hist`) that contains the given value (`val`).
"""
bin_with_val(val, hist) = searchsortedfirst(hist.edges[1], val) - 1

"""
    idxs_in_bin(i, hist, val; closed=:left)

Given the values(`val`) which are histogramed in `hist`, find all the indices
which correspond to the values in the `i`-th bin.
"""
function idxs_in_bin(i, hist, val; closed=:left)
    h = hist.edges[1]
    inbin(x) = (closed == :left) ? h[i] ≤ x < h[i+1] : h[i] < x ≤ h[i+1]
    idx = findall(inbin, val)
    return idx
end


"""
    select_series(scene, selected_plot, scatter_α, hist_α, data, hist)

Setup selection of a series in a scatter plot and the corresponding histogram.
When a point of the scatter plot is cliked, the corresponding series is
highlighted (or selected) by changing the transparency of all the other series
(and corresponding histogram bins) to a very low value.
When the click is outside, the series is deselected, that is all the αs are
set back to 1.
"""
function select_series(scene, selected_plot, scatter_α, hist_α, data, hist)
    series_idx = map(get_series_idx, selected_plot, scene)
    on(series_idx) do i
        if !isa(i, Nothing)
            scatter_α[i - 1][] = 1.
            change_α(scatter_α, setdiff(axes(scatter_α, 1), i - 1))
            selected_bin = bin_with_val(data[i-1], hist)
            hist_α[selected_bin][] = 1.
            change_α(hist_α, setdiff(axes(hist_α, 1), selected_bin))
        else
            change_α(scatter_α, axes(scatter_α, 1), 1.)
            change_α(hist_α, axes(hist_α, 1), 1.)
        end
        return nothing
    end
end

"""
    select_bin(hist_idx, hist, hist_α, scatter_α, data; closed=:left)

Setup a selection of a histogram bin and the corresponding series in the
scatter plot. See also [`select_series`](@ref).
"""
function select_bin(hist_idx, hist, hist_α, scatter_α, data; closed=:left)
    on(hist_idx) do i
        if i ≠ 0
            hist_α[i][] = 1.
            change_α(hist_α, setdiff(axes(hist.weights, 1), i))
            change_α(scatter_α, idxs_in_bin(i, hist, data, closed=closed), 1.)
            change_α(scatter_α, setdiff(axes(scatter_α, 1), idxs_in_bin(i, hist, data, closed=closed)))
        else
            change_α(scatter_α, axes(scatter_α, 1), 1.)
            change_α(hist_α, axes(hist_α, 1), 1.)
        end
        return nothing
    end
end

function poincare_explorer(sim, vals; nbins=50, closed=:left)
    hist = fit(StatsBase.Histogram, vals, nbins=nbins, closed=closed)

    colors = Float32.(axes(hist.weights, 1))

    scatter_sc_with_ui, scatter_α = plot_sim(sim, colors=vals)
    scatter_sc = scatter_sc_with_ui.children[2]

    hist_sc, hist_α = plot_hist(hist)

    # set axis names
    # if axis == 3
    #     scatter_sc[Axis][:names, :axisnames][] = ("q₂","p₂")
    # else
    #     scatter_sc[Axis][:names, :axisnames][] = ("q₁","p₁")
    # end
    sc = AbstractPlotting.vbox(scatter_sc_with_ui, hist_sc)

    selected_plot = setup_click(scatter_sc, 1)
    hist_idx = setup_click(hist_sc, 2)

    select_series(scatter_sc, selected_plot, scatter_α, hist_α, vals, hist)
    select_bin(hist_idx, hist, hist_α, scatter_α, vals, closed=closed)

    return sc
end

sim = [rand(50,50) for i=1:10]
vals = rand(10)

poincare_explorer(sim, vals)
