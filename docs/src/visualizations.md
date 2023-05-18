# Interactive GUIs, animations, visualizations

Using the functionality of package extensions in Julia v1.9+, DynamicalSystems.jl provides various visualization tools as soon as the [Makie](https://makie.juliaplots.org/stable/) package comes into scope (i.e., when `using Makie` or any of its backends like `GLMakie`).

The main functionality is [`interactive_trajectory`](@ref) that allows building custom GUI apps for visualizing the time evolution of dynamical systems. The remaining GUI applications in this page are dedicated to more specialized scenarios.


## Interactive trajectory evolution

```@raw html
<video width="100%" height="auto" controls autoplay loop>
<source src="https://raw.githubusercontent.com/JuliaDynamics/JuliaDynamics/master/videos/interact/interactive_trajectory.mp4?raw=true" type="video/mp4">
</video>
```

```@docs
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

### Example 2: Adding parameter-dependent elements to the plot
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

If now one interactively clicked (if using GLMakie) the parameter sliders and then update, the system parameters would be updated accordingly. We do it here manually via code

lalala.


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
