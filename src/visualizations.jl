export phasespace, plot_linear_regions
export plot_dataset, plot_trajectory
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
function phasespace(ds::DiscreteDynamicalSystem{IIP, S, 2}, limits,
    density, t;
    ax = PyPlot.gca(), kwargs...) where {IIP, S}

    data = Dataset{2, eltype(S)}()
    for x in range(limits[1][1], stop=limits[1][2], length=density)
        for y in range(limits[2][1], stop=limits[2][2], length=density)
            append!(data, trajectory(ds, t, SVector{2}(x,y)))
        end
    end
    ax[:plot](data[:,1], data[:,2];
    marker = "s", ms = 0.5, color="black", kwargs..., lw=0)
    return nothing
end
#
function phasespace(ds::DiscreteDynamicalSystem{IIP, S, 2},
    ics::Vector{<:AbstractVector}, t;
    ax = PyPlot.gca(), kwargs...) where {IIP, S}

    data = Dataset{2, T}()
    for ic in ics
        append!(data, trajectory(ds, t, SVector{2}(ic)))
    end
    ax[:plot](data[:,1], data[:,2];
    marker = "s", ms = 0.5, color="black", kwargs..., lw=0)
    return nothing
end



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
    PyPlot.plot(x, y; ms= 5, marker="o", color = "black", lw=0)
    _plot_lrs(x, y, linear_regions(x, y;  dxi = dxi, tol = tol)...)
end

"""
    plot_dataset(data::Dataset, tvec)
Plot each column of a dataset as a timeseries vs `tvec`.
"""
function plot_dataset(data::Dataset, tvec)
    for (i, col) in enumerate(columns(data))
        PyPlot.plot(tvec, col, label = "var. $i")
    end
    PyPlot.legend()
end

"""
    plot_trajectory(ds::DynamicalSystem, T [, dt])
Plot each column of the `trajectory` as a timeseries vs `tvec`.
"""
function plot_trajectory(ds::DDS, T, dt = 1)
    data = trajectory(ds, T, dt=dt)
    tvec = inittime(ds):dt:(T+inittime(ds))
    plot_dataset(data, tvec)
end
function plot_trajectory(ds::CDS, T, dt = 0.01)
    data = trajectory(ds, T, dt=dt)
    tvec = inittime(ds):dt:(T+inittime(ds))
    plot_dataset(data, tvec)
end
