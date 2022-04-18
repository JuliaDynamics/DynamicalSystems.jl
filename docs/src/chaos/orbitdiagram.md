# Orbit Diagrams & PSOS
## Orbit Diagrams of Maps
An orbit diagram (also called bifurcation diagram) is a way to visualize the asymptotic
behavior of a map, when a parameter of the system is changed
```@docs
orbitdiagram
```

For example, let's compute the famous orbit diagram of the logistic map:
```@example MAIN
using DynamicalSystems, CairoMakie

ds = Systems.logistic()
i = 1
pvalues = 3:0.005:4
ics = [rand() for m in 1:10]
n = 2000
Ttr = 2000
p_index = 1
output = orbitdiagram(ds, i, p_index, pvalues; n = n, Ttr = Ttr)

L = length(pvalues)
x = Vector{Float64}(undef, n*L)
y = copy(x)
for j in 1:L
    x[(1 + (j-1)*n):j*n] .= pvalues[j]
    y[(1 + (j-1)*n):j*n] .= output[j]
end

fig, ax = scatter(x, y; axis = (xlabel = L"r", ylabel = L"x"),
    markersize = 0.8, color = ("black", 0.05),
)
ax.title = "total points: $(L*n)"
xlims!(ax, pvalues[1], pvalues[end]); ylims!(ax,0,1)
fig
```

Notice that if you are using `PyPlot`, the plotting process will be slow, since it is slow at plotting big numbers of points.

The function is not limited to 1D maps, and can be applied just as well to any
discrete system.

## Poincaré Surface of Section
Also called [Poincaré map](https://en.wikipedia.org/wiki/Poincar%C3%A9_map) is a
technique to reduce a continuous system into a discrete map with 1 less dimension.
We are doing this using the function:
```@docs
poincaresos
```

Here is an example of the [Henon-Heiles](https://en.wikipedia.org/wiki/H%C3%A9non%E2%80%93Heiles_system) system showing the mixed nature of the phase space
```@example MAIN
using DynamicalSystems, CairoMakie

hh = Systems.henonheiles()

plane = (1, 0.0)
u0s = [[0.0, -0.25, 0.42081, 0.0],
[0.0, -0.31596, 0.354461, 0.0591255],
[0.0, 0.1, 0.5, 0.0],
[0.0, -0.0910355, 0.459522, -0.173339],
[0.0, -0.205144, 0.449328, -0.0162098]]

fig = Figure(resolution = (500,500))
ax = Axis(fig[1,1]; xlabel = L"q_2", ylabel = L"p_2")
for u0 in u0s
    psos = poincaresos(hh, plane, 20000.0; u0 = u0)
    scatter!(ax, psos[:, 2], psos[:, 4]; markersize = 2.0)
end
fig
```

Here the surface of section was the (hyper-) plane that $q_1 = 0$. Some chaotic and regular orbits can be seen in the plot. You can tell the regular orbits apart because they look like a single connected curve. This is the result of cutting a 2-torus by a plane!

### Advanced hyperplane
Finally here is one more example with a more complex hyperplane:
```@example MAIN
gis = Systems.gissinger([2.32865, 2.02514, 1.98312])

# Define appropriate hyperplane for gissinger system
const ν = 0.1
const Γ = 0.9 # default parameters of the system

# I want hyperperplane defined by these two points:
Np(μ) = SVector{3}(sqrt(ν + Γ*sqrt(ν/μ)), -sqrt(μ + Γ*sqrt(μ/ν)), -sqrt(μ*ν))
Nm(μ) = SVector{3}(-sqrt(ν + Γ*sqrt(ν/μ)), sqrt(μ + Γ*sqrt(μ/ν)), -sqrt(μ*ν))

# Create hyperplane passing through Np, Nm and 0:
using LinearAlgebra
gis_plane(μ) = [cross(Np(μ), Nm(μ))..., 0]

μ = 0.119
set_parameter!(gis, 1, μ)
fig = Figure()
ax = Axis3(fig[1,1]; xlabel="Q", ylabel="D", zlabel="V")
psos = poincaresos(gis, gis_plane(μ), 10000.0, Ttr = 200.0)
scatter!(ax, columns(psos)...)
fig
```

### Stroboscopic Map
A special case of a PSOS is a stroboscopic map, which is defined for non-autonomous
systems with periodic time dependence, like e.g. the [`Systems.duffing`](@ref) oscillator.

A "cut" through the phase-space can be produced at every period $T = 2\pi/\omega$. There is no
reason to use `poincaresos` for this though, because you can simply use
[`trajectory`](@ref) and get the solution with a certain time sampling rate.
For example, this piece of code:
```julia
using DynamicalSystems, Plots

ds = Systems.duffing(β = -1, ω = 1, f = 0.3) # non-autonomous chaotic system

frames=120
a = trajectory(ds, 100000.0, Δt = 2π/frames, Ttr=20π) # every period T = 2π/ω

orbit_length = div(size(a)[1], frames)
a = Matrix(a)

@gif for i in 1:frames
    orbit_points = i:frames:(orbit_length*frames)
    scatter(a[orbit_points, 1], a[orbit_points, 2], markersize=1, html_output_format=:png,
        leg=false, framestyle=:none, xlims=extrema(a[:,1]), ylims=extrema(a[:,2]))
end
```

Produces this nice animation:

![](https://raw.githubusercontent.com/JuliaDynamics/JuliaDynamics/master/videos/chaos/Duffing_stroboscopic_plot.gif?raw=true)


## Producing Orbit Diagrams for Flows
The [`orbitdiagram`](@ref) does not make much sense for continuous systems, besides the
trivial case where the system is at a fixed point. In order for [`orbitdiagram`](@ref) to have meaning one must have a map.

If only there was a way to turn a continuous system into a map... **OH WAIT!** That is
what [`poincaresos`](@ref) does! By performing successive surfaces of section at different parameter values, one can indeed "produce" an orbit diagram for a flow.

We have bundled this process in the following function:
```@docs
produce_orbitdiagram
```

For example, we will calculate the orbit diagram of the Shinriki oscillator, a continuous system that undergoes a period doubling route to chaos, much like the logistic map!

```@example MAIN
ds = Systems.shinriki([-2, 0, 0.2])

pvalues = range(19, stop = 22, length = 201)
i = 1
plane = (2, 0.0)
tf = 200.0
p_index = 1

output = produce_orbitdiagram(ds, plane, i, p_index, pvalues;
                              tfinal = tf, Ttr = 200.0)

fig = Figure()
ax = Axis(fig[1,1]; xlabel = L"R_1", ylabel = L"V_1")
for (j, p) in enumerate(pvalues)
    scatter!(ax, fill(p, length(output[j])), output[j]; 
        color = ("black", 0.5), markersize = 1
    )
end
fig
```
