# [Overarching tutorial for DynamicalSystems.jl](@id tutorial)

This page serves as a short but to-the-point introduction to the **DynamicalSystems.jl** library. It outlines the core components, and how they establish an interface that is used by the rest of the library. It also provides a couple of usage examples to connect the various packages of the library together.

Going through this tutorial should take you about 20 minutes.

## Installation

To install **DynamicalSystems.jl**, simply do:
```julia
using Pkg; Pkg.add("DynamicalSystems")
```

As discussed in the [contents](@ref contents) page, this installs several packages for the Julia language, that are all exported under a common name. To use them, simply do:
```julia
using DynamicalSystems
```
in your Julia session.

## Core components

The individual packages that compose `DynamicalSystems` interact flawlessly with each other because of the following two components:

1. The [`StateSpaceSet`](@ref), which represents numerical data. They can be observed or measured from experiments, sampled trajectories of dynamical systems, or just unordered sets in a state space. A `StateSpaceSet` is a container of equally-sized points, representing multivariate timeseries or multivariate datasets. Timeseries, which are univariate sets, are represented by the `AbstractVector{<:Real}` Julia base type.
2. The [`DynamicalSystem`](@ref), which is the abstract representation of a dynamical system with a known dynamic evolution rule. `DynamicalSystem` defines an extendable interface, but typically one uses existing implementations such as [`DeterministicIteratedMap`](@ref) or [`CoupledODEs`](@ref).

## Making dynamical systems

In the majority of cases, to make a dynamical system one needs three things:

1. The dynamic rule `f`: A Julia function that provides the instructions of how to evolve the dynamical system in time.
2. The state `u`: An array-like container that contains the variables of the dynamical system and also defines the starting state of the system.
3. The parameters `p`: An arbitrary container that parameterizes `f`.

For most concrete implementations of `DynamicalSystem` there are two ways of defining `f, u`.
The distinction is done on whether `f` is defined as an in-place (iip) function or out-of-place (oop) function.

* **oop** : `f` **must** be in the form `f(u, p, t) -> out`
    which means that given a state `u::SVector{<:Real}` and some parameter container
    `p` it returns the output of `f` as an `SVector{<:Real}` (static vector).
* **iip** : `f` **must** be in the form `f!(out, u, p, t)`
    which means that given a state `u::AbstractArray{<:Real}` and some parameter container `p`,
    it writes in-place the output of `f` in `out::AbstractArray{<:Real}`.
    The function **must** return `nothing` as a final statement.

`t` stands for current time in both cases.
**iip** is suggested for systems with high dimension and **oop** for small.
The break-even point is between 10 to 100 dimensions but should be benchmarked
on a case-by-case basis as it depends on the complexity of `f`.

!!! note "Autonomous vs non-autonomous systems"
    Whether the dynamical system is autonomous (`f` doesn't depend on time) or not, it is still necessary to include `t` as an argument to `f`. Some algorithms utilize this information, some do not, but we prefer to keep a consistent interface either way.


### Example: Henon map

Let's make the Henon map, defined as
```math
\begin{aligned}
x_{n+1} &= 1 - ax^2_n+y_n \\
y_{n+1} & = bx_n
\end{aligned}
```
with parameters $a = 1.4, b = 0.3$.

First, we define the dynamic rule as a standard Julia function. Since the dynamical system is only two-dimensional, we should use the _out-of-place_ form that returns an `SVector` with the next state:

```@example MAIN
using DynamicalSystems

function henon_rule(u, p, n) # here `n` is "time", but we don't use it.
    x, y = u # system state
    a, b = p # system parameters
    xn = 1.0 - a*x^2 + y
    yn = b*x
    return SVector(xn, yn)
end
```

Then, we define initial state and parameters

```@example MAIN
u0 = [0.2, 0.3]
p0 = [1.4, 0.3]
```

Lastly, we give these three to the `DeterministicIteratedMap`:

```@example MAIN
henon = DeterministicIteratedMap(henon_rule, u0, p0)
```

`henon` is a `DynamicalSystem`, one of the two core structures of the library.
They can evolved interactively, and queried, using the interface defined by [`DynamicalSystem`](@ref). The simplest thing you can do with a `DynamicalSystem` is to get its trajectory:

```@example MAIN
total_time = 10_000
X, t = trajectory(henon, total_time)
```

```@example MAIN
X
```

`X` is a `StateSpaceSet`, the second of the core structures of the library. We'll see below how, and where, to use a `StateSpaceset`, but for now let's just do a scatter plot

```@example MAIN
using CairoMakie
scatter(X[:, 1], X[:, 2])
```

### Example: Lorenz96

Let's also make another dynamical system, the Lorenz96 model:
```math
\frac{dx_i}{dt} = (x_{i+1}-x_{i-2})x_{i-1} - x_i + F
```
for $i \in \{1, \ldots, N\}$ and $N+j=j$.

Here, instead of a discrete time map we have $N$ coupled ordinary differential equations. However, creating the dynamical system works out just like above, but using `CoupledODEs` instead of `DeterministicIteratedMap`.

First, we make the dynamic rule function. Since this dynamical system can be arbitrarily high-dimensional, we prefer to use the _in-place_ form for `f`, overwriting in place the rate of change in a pre-allocated container. It is [customary](https://docs.julialang.org/en/v1/manual/style-guide/#bang-convention) to append the name of functions that modify their arguments in-place with a bang (`!`).

```@example MAIN
function lorenz96_rule!(du, u, p, t)
    F = p[1]; N = length(u)
    # 3 edge cases
    du[1] = (u[2] - u[N - 1]) * u[N] - u[1] + F
    du[2] = (u[3] - u[N]) * u[1] - u[2] + F
    du[N] = (u[1] - u[N - 2]) * u[N - 1] - u[N] + F
    # then the general case
    for n in 3:(N - 1)
        du[n] = (u[n + 1] - u[n - 2]) * u[n - 1] - u[n] + F
    end
    return nothing # always `return nothing` for in-place form!
end
```

then, like before, we define an initial state and parameters, and initialize the system

```@example MAIN
N = 6
u0 = range(0.1, 1; length = N)
p0 = [8.0]
lorenz96 = CoupledODEs(lorenz96_rule!, u0, p0)
```

and, again like before, we may obtain a trajectory the same way

```@example MAIN
total_time = 12.5
sampling_time = 0.02
Y, t = trajectory(lorenz96, total_time; Ttr = 2.2, Δt = sampling_time)
Y
```

We can't scatterplot something 6-dimensional but we can visualize all timeseries

```@example MAIN
fig = Figure()
ax = Axis(fig[1, 1]; xlabel = "time", ylabel = "variable")
for var in columns(Y)
    lines!(ax, t, var)
end
fig
```

### ODE solving

Continuous time dynamical systems are evolved through DifferentialEquations.jl.
When initializing a `CoupledODEs` you can tune the solver properties to your heart's content using any of the [ODE solvers](https://diffeq.sciml.ai/latest/solvers/ode_solve/) and any of the [common solver options](https://diffeq.sciml.ai/latest/basics/common_solver_opts/). For example:

```@example MAIN
using OrdinaryDiffEq # accessing the ODE solvers
diffeq = (alg = Vern9(), abstol = 1e-9, reltol = 1e-9)
lorenz96_vern = ContinuousDynamicalSystem(lorenz96_rule!, u0, p0; diffeq)
```

```@example MAIN
Y, t = trajectory(lorenz96_vern, total_time; Ttr = 2.2, Δt = sampling_time)
Y[end]
```

## Using dynamical systems

You may use the [`DynamicalSystem`](@ref) interface to develop algorithms that utilize dynamical systems with a known evolution rule. The two main packages of the library that do this are [`ChaosTools`](@ref) and [`Attractors`](@ref). For example, you may want to compute the Lyapunov spectrum of the Lorenz96 system from above. This is as easy as calling the `lyapunovspectrum` function with `lorenz96`

```@example MAIN
steps = 10_000
lyapunovspectrum(lorenz96, steps)
```
As expected, there is at least one positive Lyapunov exponent, because the system is chaotic, and at least one zero Lyapunov exponent, because the system is continuous time.

Alternatively, you may want to estimate the basins of attraction of a multistable dynamical system. The Henon map is "multistable" in the sense that some initial conditions diverge to infinity, and some others converge to a chaotic attractor. Computing these basins of attraction is simple with [`Attractors`](@ref), and would work as follows:

```@example MAIN
# define a state space grid to compute the basins on:
xg = yg = range(-2, 2; length = 201)
# find attractors using recurrences in state space:
mapper = AttractorsViaRecurrences(henon, (xg, yg); sparse = false)
# compute the full basins of attraction:
basins, attractors = basins_of_attraction(mapper; show_progress = false)
```

```@example MAIN
fig, ax = heatmap(xg, yg, basins)
x, y = columns(X) # attractor of Henon map
scatter!(ax, x, y; color = "black")
fig
```

You could also be using a `DynamicalSystem` instance directly to build your own algorithm if it isn't already implemented (and then later contribute it so it _is_ implemented ;) ). A dynamical system can be evolved forwards in time using `step!`:

```@example MAIN
henon
```

Notice how the time is not 0, because `henon` has already been stepped when we called the function `basins_of_attraction` with it. We can step it more:

```@example MAIN
step!(henon)
```

```@example MAIN
step!(henon, 2)
```

For more information on how to directly use `DynamicalSystem` instances, see the documentation of [`DynamicalSystemsBase`](@ref).

## State space sets

Let's recall that the output of the `trajectory` function is a `StateSpaceSet`:
```@example MAIN
X
```

It is printed like a matrix where each column is the timeseries of each dynamic variable. In reality, it is a vector of statically sized vectors (for performance reasons). When indexed with 1 index, it behaves like a vector of vectors
```@example MAIN
X[1]
```

```@example MAIN
X[2:5]
```

When indexed with two indices, it behaves like a matrix

```@example MAIN
X[2:5, 2]
```

When iterated, it iterates over the contained points
```@example MAIN
for (i, point) in enumerate(X)
    @show point
    i > 5 && break
end
```

```@example MAIN
map(point -> point[1] + point[2], X)
```

The columns of the set are obtained with the convenience `columns` function

```@example MAIN
x, y = columns(X)
summary.((x, y))
```

## Using state space sets

Several packages of the library deal with `StateSpaceSets`.

You could use [`ComplexityMeasures`](@ref) to obtain the entropy, or other complexity measures, of a given set. Below, we obtain the entropy of the natural density of the chaotic attractor by partitioning into a histogram of approximately `50` bins per dimension:
```@example MAIN
prob_est = ValueHistogram(50)
entropy(prob_est, X)
```

Alternatively, you could use [`FractalDimensions`](@ref) to get the fractal dimensions of the chaotic attractor of the henon map using the Grassberger-Procaccia algorithm:
```@example MAIN
grassberger_proccacia_dim(X)
```

Or, you could obtain a recurrence matrix of a state space set with [`RecurrenceAnalysis`](@ref)
```@example MAIN
R = RecurrenceMatrix(Y, 8.0)
Rg = grayscale(R)
rr = recurrencerate(R)
heatmap(Rg; colormap = :grays,
    axis = (title = "recurrence rate = $(rr)", aspect = 1,)
)
```


## More nonlinear timeseries analysis

A `trajectory` of a known dynamical system is one way to obtain a `StateSpaceSet`. However, another common way is via a delay coordinates embedding of a measured/observed timeseries. For example, we could use `optimal_separated_de` from [`DelayEmbeddings`](@ref) to create an optimized delay coordinates embedding of a timeseries

```@example MAIN
w = Y[:, 1] # first variable of Lorenz96
𝒟, τ, e = optimal_separated_de(w)
𝒟
```

and compare

```@example MAIN
fig = Figure()
axs = [Axis3(fig[1, i]) for i in 1:2]
for (S, ax) in zip((Y, 𝒟), axs)
    lines!(ax, S[:, 1], S[:, 2], S[:, 3])
end
fig
```

Since `𝒟` is just another state space set, we could be using any of the above analysis pipelines on it just as easily.

The last package to mention here is [`TimeseriesSurrogates`](@ref), which ties with all other observed/measured data analysis by providing a framework for confidence/hypothesis testing. For example, if we had a measured timeseries but we were not sure whether it represents a deterministic system with structure in the state space, or mostly noise, we could do a surrogate test. For this, we use `surrogenerator` and `RandomFourier` from [`TimeseriesSurrogates`](@ref), and the `generalized_dim` from [`FractalDimensions`](@ref) (because it performs better in noisy sets)

```@example MAIN
x = X[:, 1] # Henon map timeseries
# contaminate with noise
using Random: Xoshiro
rng = Xoshiro(1234)
x .+= randn(rng, length(x))/100
# compute noise-contaminated fractal dim.
Δ_orig = generalized_dim(embed(x, 2, 1))
```

And we do the surrogate test
```@example MAIN
surrogate_method = RandomFourier()
sgen = surrogenerator(x, surrogate_method, rng)

Δ_surr = map(1:1000) do i
    s = sgen()
    generalized_dim(embed(s, 2, 1))
end
```

and visualize the test result

```@example MAIN
fig, ax = hist(Δ_surr)
vlines!(ax, Δ_orig)
fig
```
since the real value is outside the distribution we have confidence the data are not pure noise.


## Integration with ModelingToolkit.jl

DynamicalSystems.jl understands when a model has been generated via [ModelingToolkit.jl](https://docs.sciml.ai/ModelingToolkit/stable/). The symbolic variables used in ModelingToolkit.jl can be used to access the state or parameters of the dynamical system.

To access this functionality, the `DynamicalSystem` must be created from a `DEProblem` of the SciML ecosystem, and the `DEProblem` itself must be created from a ModelingToolkit.jl model.

Let's create a the Roessler system as an MTK model:
```@example MAIN
using ModelingToolkit

@mtkmodel Roessler begin
    @variables t
    D = Differential(t)
    @parameters begin
        a = 0.2
        b = 0.2
        c = 5.7
    end
    @variables begin
        x(t) = 1.0
        y(t) = 0.0
        z(t) = 0.0
        nlt(t) # nonlinear term
    end
    @equations begin
        D(x) ~ -y -z
        D(y) ~ x + a*y
        D(z) ~ b + nlt
        nlt ~ z*(x - c)
    end
end

@mtkbuild model = Roessler()
```
this model can then be made into an `ODEProblem`:
```@example MAIN
prob = ODEProblem(model)
```
(notice that because we specified initial values for all parameters and variables during the model creation  we do need to provide additional initial values)

Now, this problem can be made into a [`CoupledODEs`](@ref):

```@example MAIN
roessler = CoupledODEs(prob)
```

This dynamical system instance can be used in the rest of the library like anything else. Additionally, you can "observe" referenced symbolic variables:
```@example MAIN
observe_state(roessler, model.x)
```

```@example MAIN
observe_state(roessler, model.nlt)
```

These observables can also be used in the GUI visualization [`interactive_trajectory_timeseries`](@ref).

You can also symbolically alter parameters

```@example MAIN
current_parameter(roessler, model.c)
```

```@example MAIN
set_parameter!(roessler, 5.0, model.c)
```

and this symbolic indexing can be given anywhere in the ecosystem where you would be altering the parameters.

## Core components reference

```@docs
StateSpaceSet
DynamicalSystem
```

## Dynamical system implementations

```@docs
DeterministicIteratedMap
CoupledODEs
StroboscopicMap
PoincareMap
ProjectedDynamicalSystem
ArbitrarySteppable
```

## Learn more

To learn more, you need to visit the documentation pages of the modules that compose DynamicalSystems.jl. See the [contents](@ref contents) page for more!
