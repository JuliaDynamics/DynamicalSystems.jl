using Makie, Observables
using StatsBase
using StatsMakie
export trajectory_highlighter

const DEFAULT_α = 0.01

"""
    trajectory_highlighter(datasets, vals; kwargs...)
Open an interactive application for highlighting specific datasets
and properties of these datasets. `datasets` is a vector of `Dataset` from
**DynamicalSystems.jl**. Each dataset corresponds to a specific value from `vals`
(a `Vector{<:Real}`). The value of `vals` gives each dataset
a specific color based on a colormap.

The application is composed of two scenes: the left scene plots the datasets,
while the right scene plots the histogram of the `vals`. The function
returns the two scenes `data_scene, hist_scene`.

Two dimensional datasets are plotted as scatter plots (and are assumed to be discrete
in nature), while three dimensional are plotted as lines (and are assummed continuous).

## Interaction
Clicking on a bin of the histogram plot will "highlight" all data
whose value belongs in that bin. Here highlighting actually means "hidding"
(i.e. reducing their alpha value) all other data besides the ones you want
to highlight. Clicking on empty space on the histogram plot will reset
highlighting.

Clicking on a plotted series in the left window will highlight this series
as well as the histogram bin that contains its value. Clicking on empty
space will reset the highlighting.

## Keyword Arguments
* `nbins = 10, closed = :left` : used in producing the histogram.
* `α = 0.05` : the alpha value of the hidden data.
* `hα = 0.2` : the alpha value of the hidden histogram bins.
* `cmap = :viridis` : the colormap used.
* `hname = "value"` : name for the histogram axis.
* `kwargs...` : Anything else is propagated to `plot_dataset`.
"""
function trajectory_highlighter(datasets, vals;
    nbins = 10, closed=:left, α = 0.05,
    cmap = :viridis, hα = 0.2, hname = "value", kwargs...)

    N = length(datasets)
    N == length(vals) || error("data and value must have equal length")

    # First prepare the colors of the datasets:
    colormap = Makie.categorical_colors(cmap, length(datasets))
    get_color(i) = Makie.color(Makie.interpolated_getindex(
        colormap, vals[i], extrema(vals)
    ))
    # The colors are observables; the transparency can be changed
    data_αs = [Observable(1.0) for i in 1:N]
    colors = [lift(α -> RGBAf(get_color(i), α), data_αs[i]) for i ∈ 1:N]
    data_scene = plot_datasets(datasets, colors; kwargs...)

    # now time for the histogram:
    hist = fit(StatsBase.Histogram, vals, nbins=nbins, closed=closed)
    hist_scene, hist_αs = plot_histogram(hist, cmap)

    hist_scene[Axis][:names, :axisnames] = (hname, "count")

    sc = Makie.vbox(data_scene, hist_scene)

    selected_plot = setup_click(data_scene, 1)
    hist_idx = setup_click(hist_scene, 2)

    select_series(data_scene, selected_plot, data_αs, hist_αs, vals, hist, α, hα)
    select_bin(hist_idx, hist, hist_αs, data_αs, vals, closed=closed, α = α, hα = hα)
    display(sc)
    return data_scene, hist_scene
end


"""
    plot_histogram(hist, cmap) -> hist_scene, hist_αs
Plot a histogram where the transparency (α) of each bin can be changed
and return the scene together with the αs. The bins are colored
according to a colormap.
"""
function plot_histogram(hist::StatsBase.Histogram, cmap)
    c = Makie.categorical_colors(cmap, length(hist.weights))
    hist_αs = [Observable(1.) for i in c]
    bincolor(αs...) = RGBAf.(c, αs)
    colors = lift(bincolor, hist_αs...)
    hist_scene = plot(hist, color=colors)
    return hist_scene, hist_αs
end

"""
    change_α(series_alpha, idxs, α = DEFAULT_α)
Given a vector of `Observable`s that represent the αs for some series, change
the elements with indices given by `idxs` to the value `α`.
This can be used to hide some series by using a low α value (default).
To restore the initial color, use `α = 1`.
"""
function change_α(series_alpha, idxs, α = DEFAULT_α)
    foreach(i -> series_alpha[i][] = α, idxs)
end

"""
    get_series_idx(selected_plot, scene)
Get the index of the `selected_plot` in `scene`.
"""
function get_series_idx(selected_plot, scene)
    # TODO: Find a less kacky way of doing this
    # Note that the first index will be for the axis
    plot_idx = findfirst(p->selected_plot === p, scene.plots)
    return isnothing(plot_idx) ? nothing : plot_idx - 1
end

"""
    setup_click(scene, idx=1)
Given a `scene` return a `Observable` that listens to left clicks inside the scene.
The `idx` argument is used to index the tuple `(plt, click_idx)` which gives
the selected plot and the index of the selected element in the plot.
"""
function setup_click(scene, idx=1)
    selection = Observable{Any}(0)
    on(scene.events.mousebuttons) do event
        if event.button == Mouse.left && event.action == Mouse.press && Makie.is_mouseinside(scene)
            plt, click_idx = mouse_selection(scene)
            if toplevelparent(plt) isa BarPlot
                click_idx = (click_idx + 1) ÷ 4
            end
            selection[] = (plt, click_idx)[idx]
        end
        return false
    end
    return selection
end

toplevelparent(plt) = isa(plt.parent, Scene) ? plt : toplevelparent(plt.parent)
toplevelparent(::Nothing) = nothing

"""
    bin_with_val(val, hist)
Get the index of the bin in the histogram(`hist`) that contains the given value (`val`).
"""
bin_with_val(val, hist) = searchsortedfirst(hist.edges[1], val) - 1

"""
    idxs_in_bin(i, hist, val; closed=:left)
Given the values (`val`) which are histogramed in `hist`, find all the indices
which correspond to the values in the `i`-th bin.
"""
function idxs_in_bin(i, hist, val; closed=:left)
    h = hist.edges[1]
    inbin(x) = (closed == :left) ? h[i] ≤ x < h[i+1] : h[i] < x ≤ h[i+1]
    idx = findall(inbin, val)
    return idx
end


"""
    select_series(scene, selected_plot, data_αs, hist_αs, data, hist)
Setup selection of a series in a scatter plot and the corresponding histogram.
When a point of the scatter plot is clicked, the corresponding series is
highlighted (or selected) by changing the transparency of all the other series
(and corresponding histogram bins) to a very low value.
When the click is outside, the series is deselected, that is all the αs are
set back to 1.
"""
function select_series(scene, selected_plot, data_αs,
                       hist_αs, data, hist,
                       α = DEFAULT_α, hα = 0.2
                       )
    series_idx = map(get_series_idx, selected_plot, scene)
    on(series_idx) do i
        if !isa(i, Nothing)
            data_αs[i][] = 1.0
            change_α(data_αs, setdiff(axes(data_αs, 1), i), α)
            selected_bin = bin_with_val(data[i], hist)
            hist_αs[selected_bin][] = 1.0
            change_α(hist_αs, setdiff(axes(hist_αs, 1), selected_bin), hα)
        else
            change_α(data_αs, axes(data_αs, 1), 1.0)
            change_α(hist_αs, axes(hist_αs, 1), 1.0)
        end
        return nothing
    end
end

"""
    select_bin(hist_idx, hist, hist_αs, data_αs, data; closed=:left, α = DEFAULT_α)
Setup a selection of a histogram bin and the corresponding series in the
scatter plot. See also [`select_series`](@ref).
"""
function select_bin(hist_idx, hist, hist_αs, data_αs, data;
    closed=:left, α = DEFAULT_α, hα = 0.2)

    on(hist_idx) do i
        if i ≠ 0
            hist_αs[i][] = 1.0
            change_α(hist_αs, setdiff(axes(hist.weights, 1), i), hα)
            change_α(data_αs, idxs_in_bin(i, hist, data, closed=closed), 1.0)
            change_α(data_αs, setdiff(
                axes(data_αs, 1), idxs_in_bin(i, hist, data, closed=closed)
            ), α)
        else
            change_α(data_αs, axes(data_αs, 1), 1.0)
            change_α(hist_αs, axes(hist_αs, 1), 1.0)
        end
        return nothing
    end
end
