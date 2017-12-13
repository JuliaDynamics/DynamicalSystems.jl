export phasespace, plot_linear_regions
using PyPlot

"""
    phasespace(ds, limits, density::Int, t::Int; kwargs...)
    phasespace(ds, ics, t::Int; kwargs...)
Plot a phasespace of a 2D discrete dynamical system, by evolving it for `t` time.

Initial conditions can be given in the form of `limits, density` with
`limits = [(x0, xf), (y0, yf)]` or directly
as a vector of initial conditions using the second method.

## Keyword Arguments:
* `ax=PyPlot.gca()` : the axes to plot on.
* `kwargs...` : Keyword arguments passed into `ax[:plot]()`.
"""
function phasespace(ds::DiscreteDS{2, T, F, J}, limits,
    density, t;
    ax = PyPlot.gca(), kwargs...) where {T, F, J}

    data = Dataset{2, T}()
    for x in linspace(limits[1][1], limits[1][2], density)
        for y in linspace(limits[2][1], limits[2][2], density)
            ds.state = SVector{2,T}(x,y)
            append!(data, trajectory(ds, t))
        end
    end
    ax[:plot](data[:,1], data[:,2];
    marker = "s", ms = 0.5, color="black", kwargs..., lw=0)
    return nothing
end

function phasespace(ds::DiscreteDS{2, T, F, J}, ics::Vector{<:SVector}, t;
    ax = PyPlot.gca(), kwargs...) where {T, F, J}

    data = Dataset{2, T}()
    for ic in ics
        ds.state = ic
        append!(data, trajectory(ds, t))
    end
    ax[:plot](data[:,1], data[:,2];
    marker = "s", ms = 0.5, color="black", kwargs..., lw=0)
    return nothing
end



# This function exists ONLY FOR TESTING! Do not use it elsewhere!
function _plot_lrs(x, y, lrs, tangents)
  for i ∈ 1:length(lrs)-1
    PyPlot.plot(x[lrs[i]:lrs[i+1]], y[lrs[i]:lrs[i+1]])
  end
end

"""
    plot_linear_regions(x, y; dxi = 1, tol = 0.2)
Visualize the outcome of calling `linear_regions` by plotting each
linear segment with a different color (on the current axes).
"""
function plot_linear_regions(x, y; dxi = 1, tol = 0.2)
    PyPlot.plot(x,y, lw=0, ms= 5, marker="o", color = "black")
  _plot_lrs(x, y, linear_regions(x, y;  dxi = dxi, tol = tol)...)
end
