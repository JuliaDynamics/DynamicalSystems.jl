# Attractor Basins, Tipping Points

In this page we list several functions related with basins of attraction and tipping points. In the example [Basins in Higher Dimensions](@ref) we try to apply every single function listed below, so check this for an example application of everything listed here!

## Computing Basins of Attraction

In DynamicalSystems.jl we provide performant methods for estimating basins of attraction of various attractors.
The performance of these methods comes from a _constraint on a 2D plane_, as you will see below.

```@docs
basins_of_attraction
match_attractors!
```

## Final state sensitivity
```@docs
uncertainty_exponent
```

## Tipping points
```@docs
basin_fractions
tipping_probabilities
```

## Discrete system example
```@example MAIN
function newton_map(dz, z, p, n)
    z1 = z[1] + im*z[2]
    dz1 = newton_f(z1, p[1])/newton_df(z1, p[1])
    z1 = z1 - dz1
    dz[1]=real(z1)
    dz[2]=imag(z1)
    return
end
newton_f(x, p) = x^p - 1
newton_df(x, p)= p*x^(p-1)

# dummy Jacobian function due to https://github.com/JuliaDiff/ForwardDiff.jl/issues/520
function newton_map_J(J,z0, p, n) end

ds = DiscreteDynamicalSystem(newton_map, [0.1, 0.2], [3.0], newton_map_J)
xg = yg = range(-1.5,1.5,length=400)
basins, attractors = basins_of_attraction((xg, yg), ds)
basins
```
```@example MAIN
attractors
```

Now let's plot this as a heatmap
```@example MAIN
# Set up some code for plotting attractors
function scatter_attractors(attractors)
    for k ∈ keys(attractors)
        x, y = columns(attractors[k])
        scatter(x, y; color = "C$(k-1)", alpha = 0.75, edgecolor = "white")
    end
end
LC =  matplotlib.colors.ListedColormap
cmap = LC([matplotlib.colors.to_rgb("C$k") for k in 0:3])
vmin = 1; vmax = 4

fig = figure()
pcolormesh(xg, yg, basins'; cmap, vmin, vmax)
scatter_attractors(attractors)
fig.tight_layout(pad=0.3); fig
```


## Stroboscopic map example
This example targets periodically driven 2D continuous dynamical systems, like the Duffing oscillator:
```@example MAIN
using DynamicalSystems, PyPlot
ω=1.0; f = 0.2
ds = Systems.duffing([0.1, 0.25]; ω, f, d = 0.15, β = -1)
```
For stroboscopic maps, we strongly recommend using a higher precision integrator from OrdinaryDiffEq.jl.

```@example MAIN
using OrdinaryDiffEq
diffeq = (alg = Tsit5(), reltol = 1e-9)
xg = yg = range(-2.2,2.2,length=200)
basins, attractors = basins_of_attraction((xg, yg), ds; T=2π/ω, diffeq)
basins
```

And visualize the result as a heatmap, scattering the found attractors via scatter.

```@example MAIN
fig = figure()
pcolormesh(xg, yg, basins'; cmap, vmin, vmax)
scatter_attractors(attractors)
fig.tight_layout(pad=0.3); fig
```

## Basins in Higher Dimensions
In this section we will calculate the basins of attraction of as system projected on a 2D plane. See the [Three dimensional basins](@ref) section for something even more complex.

For this example we will compute the basins of attraction of the magnetic pendulum, projected on the x-y plane, and then analyze their properties using all functions stated initially in this page.
### Computing the basins

```@example MAIN
ds = Systems.magnetic_pendulum(d=0.2, α=0.2, ω=0.8, N=3)
xg = yg = range(-4, 4, length=150)
@time basins, attractors = basins_of_attraction((xg, yg), ds; diffeq = (reltol = 1e-9,))
attractors
```
Alright, so far so good, we found 3 attractors (the 3 magnets).
Let's visualize this beauty now

```@example MAIN
fig = figure()
pcolormesh(xg, yg, basins'; cmap, vmin, vmax)
scatter_attractors(attractors)
fig.tight_layout(pad=0.3); fig
```

### Computing the uncertainty exponent
Let's now calculate the [`uncertainty_exponent`](@ref) for this system as well.
The calculation is straightforward:
```@example MAIN
ε, f_ε, α = uncertainty_exponent(xg, yg, basins)
fig = figure()
plot(log.(ε), log.(f_ε))
plot(log.(ε), log.(ε) .* α)
fig.tight_layout(pad=0.3); fig
```

### Computing the tipping probabilities
We will compute the tipping probabilities using the magnetic pendulum's example
as the "before" state. For the "after" state we will change the `γ` parameter of the
third magnet to be so small, it's basin of attraction will virtually disappear.
```@example MAIN
ds = Systems.magnetic_pendulum(d=0.2, α=0.2, ω=0.8, N=3, γs = [1.0, 1.0, 0.1])
basins_after, attractors_after = basins_of_attraction((xg, yg), ds; diffeq = (reltol = 1e-9,))
match_attractors!(basins, attractors, basins_after, attractors_after)
fig = figure()
pcolormesh(xg, yg, basins_after'; vmin, vmax, cmap)
scatter_attractors(attractors_after)
fig.tight_layout(pad=0.3); fig
```

```@example MAIN
P = tipping_probabilities(basins, basins_after)
```
As you can see `P` has size 3×2, as after the change only 2 attractors have been identified in the system (3 still exist but our state space discretization isn't accurate enough to find the 3rd because it has such a small basin).
Also, the first row of `P` is 50% probability to each other magnet, as it should be due to the system's symmetry.

## Higher-dimensional basins
To showcase the power of [`basins_of_attraction`](@ref) we need to use a system whose attractors span higher-dimensional space. An example is 
```@example MAIN
ds = Systems.thomas_cyclical(b = 0.1665)
```
which, for this parameter, contains 5 coexisting attractors. 3 of them are entangled periodic orbits that span across all three dimensions, and the remaining 2 are fixed points.

To compute the basins we define a three-dimensional grid and call on it [`basins_of_attraction`](@ref).

```julia
# This computation takes about an hour
xg = yg = zg = range(-6.0, 6.0; length = 250)
basins, attractors = basins_of_attraction((xg, yg, zg), ds)
attractors
```
```
Dict{Int16, Dataset{3, Float64}} with 5 entries:
  5 => 3-dimensional Dataset{Float64} with 1 points
  4 => 3-dimensional Dataset{Float64} with 362 points
  2 => 3-dimensional Dataset{Float64} with 364 points
  3 => 3-dimensional Dataset{Float64} with 362 points
  1 => 3-dimensional Dataset{Float64} with 1 points
```

The basins of attraction are very complicated. We can try to visualize them by animating the 2D slices at each z value, to obtain:

```@raw html
<video width="75%" height="auto" controls loop>
<source src="https://raw.githubusercontent.com/JuliaDynamics/JuliaDynamics/master/videos/chaos/cyclical_basins.mp4?raw=true" type="video/mp4">
</video>
```

Then, we visualize the attractors to obtain:

```@raw html
<video width="75%" height="auto" controls loop>
<source src="https://raw.githubusercontent.com/JuliaDynamics/JuliaDynamics/master/videos/chaos/cyclical_attractors.mp4?raw=true" type="video/mp4">
</video>
```

In the animation above, the scattered points are the attractor values the function [`basins_of_attraction`](@ref) found by itself. Of course, for the periodic orbits these points are incomplete. Once the function's logic understood we are on an attractor, it stops computing. However, we also simulated lines, by evolving initial conditions colored appropriately with the basins output.

The animation was produced with the code:
```julia
using GLMakie
fig = Figure()
display(fig)
ax = fig[1,1] = Axis3(fig; title = "found attractors")
cmap = cgrad(COLORSCHEME[1:5], 5; categorical = true)

for i in keys(attractors)
    tr = attractors[i]
    markersize = length(attractors[i]) > 10 ? 2000 : 6000
    marker = length(attractors[i]) > 10 ? :circle : :rect
    scatter!(ax, columns(tr)...; markersize, marker, transparency = true, color = cmap[i])
    j = findfirst(isequal(i), bsn)
    x = xg[j[1]]
    y = yg[j[2]]
    z = zg[j[3]]
    tr = trajectory(ds, 100, SVector(x,y,z); Ttr = 100)
    lines!(ax, columns(tr)...; linewidth = 1.0, color = cmap[i])
end

a = range(0, 2π; length = 200) .+ π/4

record(fig, "cyclical_attractors.mp4", 1:length(a)) do i
    ax.azimuth = a[i]
end
```