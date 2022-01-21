# Attractor Basins, Tipping Points

In this page we list several functions related with basins of attraction and tipping points. In the example [2D basins of higher dimensional system](@ref) we try to apply every single function listed below, so check this for an example application of everything listed here!

## Computing basins of attraction

```@docs
basins_of_attraction
automatic_Δt_basins
match_attractors!
```

## Final state sensitivity / fractal boundaries
```@docs
basins_fractal_dimension
basin_entropy
basins_fractal_test
uncertainty_exponent
```

## Tipping points
```@docs
basin_fractions
tipping_probabilities
```

## Discrete system example
```@example MAIN
using DynamicalSystems, PyPlot
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
basins, attractors = basins_of_attraction((xg, yg), ds; show_progress = false)
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
        scatter(x, y; color = "C$(k-1)", edgecolor = "white")
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
basins, attractors = basins_of_attraction((xg, yg), ds; T=2π/ω, diffeq, show_progress = false)
basins
```

And visualize the result as a heatmap, scattering the found attractors via scatter.

```@example MAIN
fig = figure()
pcolormesh(xg, yg, basins'; cmap, vmin, vmax)
scatter_attractors(attractors)
fig.tight_layout(pad=0.3); fig
```

## 2D basins of higher dimensional system
In this section we will calculate the basins of attraction of the four-dimensional magnetic pendulum. We know that the attractors of this system are all individual fixed points on the (x, y) plane so we will only compute the basins there. See the [3D basins](@ref) section for something even more complex.

### Computing the basins

```@example MAIN
ds = Systems.magnetic_pendulum(d=0.2, α=0.2, ω=0.8, N=3)
xg = yg = range(-4, 4, length=150)
@time basins, attractors = basins_of_attraction((xg, yg), ds; diffeq = (reltol = 1e-9,), show_progress = false)
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
ε, f_ε, α = uncertainty_exponent(basins)
fig = figure()
plot(log.(ε), log.(f_ε))
plot(log.(ε), log.(ε) .* α; label = "α = $(round(α; digits=3))")
legend()
fig.tight_layout(pad=0.3); fig
```
The actual uncertainty exponent is the slope of the curve.

### Computing the tipping probabilities
We will compute the tipping probabilities using the magnetic pendulum's example
as the "before" state. For the "after" state we will change the `γ` parameter of the
third magnet to be so small, its basin of attraction will virtually disappear.
```@example MAIN
ds = Systems.magnetic_pendulum(d=0.2, α=0.2, ω=0.8, N=3, γs = [1.0, 1.0, 0.1])
basins_after, attractors_after = basins_of_attraction(
    (xg, yg), ds; diffeq = (reltol = 1e-9,), show_progress = false
)
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

### Computing the fractions of initial conditions going to each attractor
We can also compute the fraction of initial conditions that go to each of the three
attractors. First, we need to define the initial conditions to be analyzed. Using the
`sampler` function, we can simply provide the bounds of the region we're interested in. With
the uniform sampling method, the fractions we'll compute can also be interpreted as the
volume of the basins of attraction.
```@example MAIN
s = ChaosTools.sampler(min_bounds=[-4, -4, -1, -1], max_bounds=[4, 4, 1, 1], method="uniform")
ics = Dataset([s() for i=1:500])
```

Next, we need to provide the function that extracts features from each trajectory. These
features should be similar for trajectories belonging to the same attractor because then the
algorithm will identify clusters of features, and consider each cluster to be a unique
attractor. This is the key step for this method, and unfortunately requires some knowledge
of the system. For the magnetic pendulum, we know the attractors are in different regions in
space, so the final position is already enough to distinguish between them.
```@example MAIN
feature_extraction(t, y) = y[end][1:2] #3 and 4 are the velocities, which go to 0
fs, class_labels = basin_fractions_clustering(
    ds, feature_extraction, ics; T=2000, Ttr=1998, show_progress=true)
fs
```

Due to the rotational symmetry in this system, the attractors are also symmetric, and should
attract the same amount of initial conditions. We can see this in the results, in which the
fractions of initial conditions are roughly the same. As the number of initial conditions is
increased, the numbers get even closer to the theoretical `1/3`.

In this first example, we assumed we didn't know where the attractors were, so we identified
them (and the size of their basins of attraction) simply from the features via the
unsupervised classification mode of `basin_fractions_clustering`. But we happen to know the
positions of the attractors for this specific system, so we can also use the supervised
classification mode. For it, we define where points where the attractors are, and the
algorithm simply classifies the trajectory according to its closest attractor. In this
system, the attractors are fixed points uniformly distributed in the unit circle (as also
shown in the results found via `basins_of_attraction`):
```@example MAIN
attractors_ic = Dataset([[cos(2π*i/3) sin(2π*i/3) 0 0] for i=1:3])
fs, class_labels = basin_fractions_clustering(
    ds, feature_extraction, ics, attractors_ic; T=2000, Ttr=1998, show_progress=true)
fs
```

The results are similar to before, with the fractions close to `1/3`. This supervised mode
is much quicker than the unsupervised one, but of course requires additional information.
Let's use this improved speed to increase the number of initial conditions and then plot the
basins of attraction. We can do this because the method also returns the labels for each
initial condition, identifying to which attractor they were classified.
```@example MAIN
ics = Dataset([s() for i=1:30000])
fs, class_labels = basin_fractions_clustering(
    ds, feature_extraction, ics, attractors_ic; T=2000, Ttr=1998, show_progress=true)
cmap2 = LC([matplotlib.colors.to_rgb("C$k") for k in 0:2])
fig = figure()
scatter(ics[:,1], ics[:,2], c=class_labels, cmap=cmap2, s=80)
fig.tight_layout(pad=0.3); fig
```

The integration of the 30000 trajectories can be performed much faster if it is parallelized.
You can read more about this in `basin_fractions_clustering`'s documentation.


## 3D basins
To showcase the true power of [`basins_of_attraction`](@ref) we need to use a system whose attractors span higher-dimensional space. An example is 
```@example MAIN
ds = Systems.thomas_cyclical(b = 0.1665)
```
which, for this parameter, contains 5 coexisting attractors. 3 of them are entangled periodic orbits that span across all three dimensions, and the remaining 2 are fixed points.

To compute the basins we define a three-dimensional grid and call on it [`basins_of_attraction`](@ref).

```julia
# This computation takes about an hour
xg = yg = zg = range(-6.0, 6.0; length = 251)
basins, attractors = basins_of_attraction((xg, yg, zg), ds)
attractors
```
```
Dict{Int16, Dataset{3, Float64}} with 5 entries:
  5 => 3-dimensional Dataset{Float64} with 1 points
  4 => 3-dimensional Dataset{Float64} with 379 points
  6 => 3-dimensional Dataset{Float64} with 1 points
  2 => 3-dimensional Dataset{Float64} with 538 points
  3 => 3-dimensional Dataset{Float64} with 537 points
  1 => 3-dimensional Dataset{Float64} with 1 points
```

The basins of attraction are very complicated. We can try to visualize them by animating the 2D slices at each z value, to obtain:

```@raw html
<video width="75%" height="auto" controls autoplay loop>
<source src="https://raw.githubusercontent.com/JuliaDynamics/JuliaDynamics/master/videos/chaos/cyclical_basins.mp4?raw=true" type="video/mp4">
</video>
```

Then, we visualize the attractors to obtain:

```@raw html
<video width="75%" height="auto" controls autoplay loop>
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
cmap = cgrad(:dense, 6; categorical = true)

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

## Poincaré map example
In the previous example we saw that this system has periodic attractors. In the Poincaré map these periodic attractors become points. We can use the functionality of [`basins_of_attraction`](@ref) and [`poincaremap`](@ref) to find basins of attraction on the Poincaré surface of section.

```@example MAIN
ds = Systems.thomas_cyclical(b = 0.1665)
xg = yg = range(-6.0, 6.0; length = 100)
pmap = poincaremap(ds, (3, 0.0), 1e6; 
    rootkw = (xrtol = 1e-8, atol = 1e-8), reltol=1e-9
)
basins, attractors = basins_of_attraction((xg,yg), pmap; show_progress = false)
attractors
```
```@example MAIN
attractors[1]
```
Looks good so far, but let's plot it as well:
```@example MAIN
fig = figure()
pcolormesh(xg, yg, basins'; cmap, vmin, vmax)
scatter_attractors(attractors)
fig.tight_layout(pad=0.3); fig
```
This aligns perfectly with the video we produced above.