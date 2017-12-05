export phasespace

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
            dss = set_state(ds, SVector{2,T}(x,y))
            append!(data, trajectory(dss, t))
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
        dss = set_state(ds, ic)
        append!(data, trajectory(dss, t))
    end
    ax[:plot](data[:,1], data[:,2];
    marker = "s", ms = 0.5, color="black", kwargs..., lw=0)
    return nothing
end
