# Interactive GUIs, animations, visualizations

Using the functionality of package extensions in Julia v1.9+, DynamicalSystems.jl provides various visualization tools as soon as the [Makie](https://makie.juliaplots.org/stable/) package comes into scope (i.e., when `using Makie` or any of its backends like `GLMakie`).

The main functionality is [`interactive_trajectory_panel`](@ref) that allows building custom GUI apps for visualizing the time evolution of dynamical systems. The remaining GUI applications in this page are dedicated to more specialized scenarios.


## Interactive trajectory evolution

```@raw html
<video width="100%" height="auto" controls autoplay loop>
<source src="https://raw.githubusercontent.com/JuliaDynamics/JuliaDynamics/master/videos/interact/interactive_trajectory.mp4?raw=true" type="video/mp4">
</video>
```

```@docs
interactive_trajectory_panel
```


### Example 1: interactive trajectory animation

```@example MAIN
F, G, a, b = 6.886, 1.347, 0.255, 4.0
ds = PredefinedDynamicalSystems.lorenz84(; F, G, a, b)

u1 = [0.1, 0.1, 0.1] # periodic
u2 = u1 .+ 1e-3     # fixed point
u3 = [-1.5, 1.2, 1.3] .+ 1e-9 # chaotic
u4 = [-1.5, 1.2, 1.3] .+ 21e-9 # chaotic 2
u0s = [u1, u2, u3, u4]

fig, dso = interactive_trajectory_panel(
    ds, u0s; tail = 1000, fade = true,
    idxs = [1,3],
)

fig
```

We could interact with this plot live, like in the example video above. We can also progress the visuals via code as instructed by [`intractive_trajectory_panel`](@ref) utilizing the second returned argument `dso`:

```@example MAIN
step!(dso, 2000)
fig
```

(if you progress the visuals via code you probably want to give `add_controls = false` as a keyword to [`interactive_trajectory_panel`](@ref))