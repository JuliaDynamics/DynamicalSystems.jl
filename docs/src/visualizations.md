# [Interactive GUIs, animations, visualizations](@id visualizations)

Using the functionality of package extensions in Julia v1.9+, DynamicalSystems.jl provides various visualization tools as soon as the [Makie](https://makie.juliaplots.org/stable/) package comes into scope (i.e., when `using Makie` or any of its backends like `GLMakie`).

The main functionality is [`interactive_trajectory`](@ref) that allows building custom GUI apps for visualizing the time evolution of dynamical systems. The remaining GUI applications in this page are dedicated to more specialized scenarios.

## Interactive- or animated trajectory evolution

The following GUI is obtained with the function [`interactive_trajectory_timeseries`](@ref) and the code snippet below it!

```@raw html
<video width="100%" height="auto" controls autoplay loop>
<source src="https://raw.githubusercontent.com/JuliaDynamics/JuliaDynamics/master/videos/interact/interactive_trajectory.mp4?raw=true" type="video/mp4">
</video>
```

```@setup MAIN
using DynamicalSystems
struct ShowFile
    file::String
end
function Base.show(io::IO, ::MIME"text/plain", f::ShowFile)
    write(io, read(f.file))
end
```
```@example MAIN
ShowFile(joinpath(dirname(pathof(DynamicalSystems)), "../test/mtk_gui.jl")) # hide
```

```@docs
interactive_trajectory_timeseries
interactive_trajectory
```

### Example 1: interactive trajectory animation

```@example MAIN
using DynamicalSystems, CairoMakie
F, G, a, b = 6.886, 1.347, 0.255, 4.0
ds = PredefinedDynamicalSystems.lorenz84(; F, G, a, b)

u1 = [0.1, 0.1, 0.1] # periodic
u2 = u1 .+ 1e-3     # fixed point
u3 = [-1.5, 1.2, 1.3] .+ 1e-9 # chaotic
u4 = [-1.5, 1.2, 1.3] .+ 21e-9 # chaotic 2
u0s = [u1, u2, u3, u4]

fig, dsobs = interactive_trajectory(
    ds, u0s; tail = 1000, fade = true,
    idxs = [1,3],
)

fig
```

We could interact with this plot live, like in the example video above. We can also progress the visuals via code as instructed by [`interactive_trajectory`](@ref) utilizing the second returned argument `dsobs`:

```@example MAIN
step!(dsobs, 2000)
fig
```

(if you progress the visuals via code you probably want to give `add_controls = false` as a keyword to [`interactive_trajectory`](@ref))

### Example 2: Adding parameter-dependent elements to a plot

In this advanced example we add plot elements to the provided figure, and also utilize the parameter observable in `dsobs` to add animated plot elements that update whenever a parameter updates. The final product of this snippet is in fact the animation at the top of the docstring of [`interactive_trajectory_panel`](@ref).

We start with an interactive trajectory panel of the Lorenz63 system, in which we also add sliders for interactively changing parameter values

```@example MAIN
using DynamicalSystems, CairoMakie

ps = Dict(
    1 => 1:0.1:30,
    2 => 10:0.1:50,
    3 => 1:0.01:10.0,
)
pnames = Dict(1 => "σ", 2 => "ρ", 3 => "β")

lims = ((-30, 30), (-30, 30), (0, 100))

ds = PredefinedDynamicalSystems.lorenz()

u1 = [10,20,40.0]
u3 = [20,10,40.0]
u0s = [u1, u3]

fig, dsobs = interactive_trajectory(
    ds, u0s; parameter_sliders = ps, pnames, lims
)

fig
```

If now one interactively clicked (if using GLMakie) the parameter sliders and then update, the system parameters would be updated accordingly. We can also add new plot elements that depend on the parameter values using the `dsobs`:


```@example MAIN
# Fixed points of the lorenz system (without the origin)
lorenz_fixedpoints(ρ,β) = [
    Point3f(sqrt(β*(ρ-1)), sqrt(β*(ρ-1)), ρ-1),
    Point3f(-sqrt(β*(ρ-1)), -sqrt(β*(ρ-1)), ρ-1),
]

# add an observable trigger to the system parameters
fpobs = map(dsobs.param_observable) do params
    σ, ρ, β = params
    return lorenz_fixedpoints(ρ, β)
end

# If we want to plot directly on the trajectory axis, we need to
# extract it from the figure. The first entry of the figure is a grid layout
# containing the axis and the GUI controls. The [1,1] entry of the layout
# is the axis containing the trajectory plot

ax = content(fig[1,1][1,1])
scatter!(ax, fpobs; markersize = 10, marker = :diamond, color = :red)

fig
```

Now, after the live animation "run" button is pressed, we can interactively change the parameter ρ and click update, in which case both the dynamical system's ρ parameter will change, but also the location of the red diamonds.

We can also change the parameters non-interactively using `set_parameter!`

```@example MAIN
set_parameter!(dsobs, 2, 50.0)

fig
```

```@example MAIN
set_parameter!(dsobs, 2, 10.0)

fig
```

Note that the sliders themselves did not change, as this functionality is for "offline" creation of animations where one doesn't interact with sliders. The keyword `add_controls` should be given as `false` in such scenarios.

### Example 3: Observed timeseries of the system

```@example MAIN
using DynamicalSystems, CairoMakie
using LinearAlgebra: norm, dot

# Dynamical system and initial conditions
ds = Systems.thomas_cyclical(b = 0.2)
u0s = [[3, 1, 1.], [1, 3, 1.], [1, 1, 3.]] # must be a vector of states!

# Observables we get timeseries of:
function distance_from_symmetry(u)
    v = SVector{3}(1/√3, 1/√3, 1/√3)
    t = dot(v, u)
    return norm(u - t*v)
end
fs = [3, distance_from_symmetry]

fig, dsobs = interactive_trajectory_timeseries(ds, fs, u0s;
    idxs = [1, 2], Δt = 0.05, tail = 500,
    lims = ((-2, 4), (-2, 4)),
    timeseries_ylims = [(-2, 4), (0, 5)],
    add_controls = false,
    figure = (size = (800, 400),)
)

fig
```

we can progress the simulation:

```@example MAIN
step!(dsobs, 200)
fig
```

or we can even make a nice video out of it:

```@example MAIN

record(fig, "thomas_cycl.mp4", 1:100) do i
    step!(dsobs, 10)
end
```

```@raw html
<video width="auto" controls autoplay loop>
<source src="../thomas_cycl.mp4" type="video/mp4">
</video>
```

## Cobweb Diagrams
```@raw html
<video width="100%" height="auto" controls autoplay loop>
<source src="https://raw.githubusercontent.com/JuliaDynamics/JuliaDynamics/master/videos/interact/cobweb.mp4?raw=true" type="video/mp4">
</video>
```

```@docs
interactive_cobweb
```

The animation at the top of this section was done with

```julia
using DynamicalSystems, GLMakie

# the second range is a convenience for intermittency example of logistic
rrange = 1:0.001:4.0
# rrange = (rc = 1 + sqrt(8); [rc, rc - 1e-5, rc - 1e-3])

lo = Systems.logistic(0.4; r = rrange[1])
interactive_cobweb(lo, rrange, 5)
```

## Orbit Diagrams
*Notice that orbit diagrams and bifurcation diagrams are different things in DynamicalSystems.jl*

```@raw html
<video width="100%" height="auto" controls autoplay loop>
<source src="https://raw.githubusercontent.com/JuliaDynamics/JuliaDynamics/master/videos/interact/odhenon.mp4?raw=true" type="video/mp4">
</video>
```

```@docs
interactive_orbitdiagram
scaleod
```

The animation at the top of this section was done with

```julia
i = p_index = 1
ds, p_min, p_max, parname = Systems.henon(), 0.8, 1.4, "a"
t = "orbit diagram for the Hénon map"

oddata = interactive_orbitdiagram(ds, p_index, p_min, p_max, i;
                                  parname = parname, title = t)

ps, us = scaleod(oddata)
```

## Interactive Poincaré Surface of Section

```@raw html
<video width="100%" height="auto" controls autoplay loop>
<source src="https://raw.githubusercontent.com/JuliaDynamics/JuliaDynamics/master/videos/interact/interactive_psos.mp4?raw=true" type="video/mp4">
</video>
```

```@docs
interactive_poincaresos
```

To generate the animation at the start of this section you can run
```julia
using InteractiveDynamics, GLMakie, OrdinaryDiffEq, DynamicalSystems
diffeq = (alg = Vern9(), abstol = 1e-9, reltol = 1e-9)

hh = Systems.henonheiles()

potential(x, y) = 0.5(x^2 + y^2) + (x^2*y - (y^3)/3)
energy(x,y,px,py) = 0.5(px^2 + py^2) + potential(x,y)
const E = energy(get_state(hh)...)

function complete(y, py, x)
    V = potential(x, y)
    Ky = 0.5*(py^2)
    Ky + V ≥ E && error("Point has more energy!")
    px = sqrt(2(E - V - Ky))
    ic = [x, y, px, py]
    return ic
end

plane = (1, 0.0) # first variable crossing 0

# Coloring points using the Lyapunov exponent
function λcolor(u)
    λ = lyapunovs(hh, 4000; u0 = u)[1]
    λmax = 0.1
    return RGBf(0, 0, clamp(λ/λmax, 0, 1))
end

state, scene = interactive_poincaresos(hh, plane, (2, 4), complete;
labels = ("q₂" , "p₂"),  color = λcolor, diffeq...)
```

## Scanning a Poincaré Surface of Section

```@raw html
<video width="100%" height="auto" controls autoplay loop>
<source src="https://raw.githubusercontent.com/JuliaDynamics/JuliaDynamics/master/videos/interact/psos_brainscan.mp4?raw=true" type="video/mp4">
</video>
```

```@docs
interactive_poincaresos_scan
```

The animation at the top of this page was done with

```julia
using GLMakie, DynamicalSystems
using OrdinaryDiffEq: Vern9

diffeq = (alg = Vern9(), abstol = 1e-9, reltol = 1e-9)
ds = PredefinedDynamicalSystems.henonheiles()
ds = CoupledODEs(ds, diffeq)

u0s = [
    [0.0, -0.25, 0.42081, 0.0],
    [0.0, 0.1, 0.5, 0.0],
    [0.0, -0.31596, 0.354461, 0.0591255]
]
# inputs
trs = [trajectory(ds, 10000, u0)[1][:, SVector(1,2,3)] for u0 ∈ u0s]
j = 2 # the dimension of the plane

interactive_poincaresos_scan(trs, j; linekw = (transparency = true,))
```

## Interactive 2D dynamical system

```@docs
interactive_clicker
```

The `interactive_clicker` function can be used to spin up a GUI
for interactively exploring the state space of a 2D dynamical system.

For example, the following code show how to interactively explore a
[`ProjectedDynamicalSystem`](@ref):

```julia
using GLMakie, DynamicalSystems

# This is the 3D Lorenz model
lorenz = Systems.lorenz()

projection = [1, 2]
complete_state = [0.0]
projected_ds = ProjectedDynamicalSystem(lorenz, projection, complete_state)

interactive_clicker(projected_ds; tfinal = (10.0, 150.0))
```
