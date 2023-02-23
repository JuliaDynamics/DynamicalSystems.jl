# [Overarching tutorial for DynamicalSystems.jl](@id tutorial)

This page serves as a short, but to-the-point, introduction in the **DynamicalSystems.jl** library. It outlines the core components, and how they establish an interface that is used by the rest of the library. It also provides a couple of usage examples to connect the various packages of the library together.

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
2. The [`DynamicalSystem`](@ref), which is the abstract representation of a dynamical system with a known dynamic evolution rule. `DynamicalSystem` defines an extendable interface, but typically one uses concrete implementations such as [`DeterministicIteratedMap`](@ref) or [`CoupledODEs`](@ref).

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

### Example: Henon map
Let's make the Henon map, defined as
$$
\begin{aligned}
x_{n+1} &= 1 - ax^2_n+y_n \\
y_{n+1} & = bx_n
\end{aligned}
$$
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

`henon` is a `DynamicalSystem`, of the two core structures of the library.
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
$$
\frac{dx_i}{dt} = (x_{i+1}-x_{i-2})x_{i-1} - x_i + F
$$
for $i \in \{1, \ldots, N\}$ and $N+j=j$.

Here, instead of a discrete time map we have $N$ coupled ordinary differential equations. However, creating the dynamical system works out just like above, but using `CoupledODEs` instead of `DeterministicIteratedMap`.

First, we make the dynamic rule function. Since this dynamical system can be arbitrarily high-dimensional, we prefer to use the _in-place_ form for `f`, overwriting in place the rate of change in a pre-allocated container.

```@example MAIN
function lorenz96_rule(du, u, p, t)
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
lorenz96 = CoupledODEs(lorenz96_rule, u0, p0)
```

and, again like before, we may obtain a trajectory the same way

```@example MAIN
total_time = 7.5
sampling_time = 0.02
Y, t = trajectory(lorenz96, total_time; Ttr = 2.2, Δt = sampling_time)
Y
```

We can't scatterplot something 6-dimensional but we can visualize all timeseries

```@example MAIN
fig = Figure()
ax = Axis(fig[1,1]; xlabel = "time", ylabel = "variable")
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
lorenz96_vern = ContinuousDynamicalSystem(lorenz96_rule, u0, p0; diffeq)
```

```@example MAIN
X, t = trajectory(lorenz96_vern, total_time; Ttr = 2.2, Δt = sampling_time)
X[end]
```


## Using dynamical systems

You may use the [`DynamicalSystem`](@ref) interface to develop algorithms that utilize dynamical systems with a known evolution rule. The two main packages of the library that do this are `ChaosTools` and `Attractors`. For example, you may want to compute the Lyapunov spectrum of the Lorenz96 system from above. This is as easy as calling the `lyapunovspectrum` function with `lorenz96`

```@example MAIN
steps = 10_000
lyapunovspectrum(lorenz96, steps)
```
As expected, there is at least one positive Lyapunov exponent (before the system is chaotic) and at least one zero Lyapunov exponent, because the system is continuous time.

Alternatively, you may want to estimate the basins of attraction of a multistable dynamical system. The Henon map is "multistable" in the sense that some initial conditions go to infinity, and some others converge to a chaotic attractor. Computing these basins of attraction is simple with `Attractors`, and would work as follows:

```@example MAIN
# define a state space grid to compute the basins on:
xg = yg = range(-3, 3; length = 400)
# find attractors using recurrences in state space:
mapper = AttractorsViaRecurrences(henon, (xg, yg); sparse = false)
# compute the full basins of attraction:
basins, attractors = basins_of_attraction(mapper; show_progress = false)
```

```@example MAIN
heatmap(xg, yg, basins)
```

## State space sets

Let's recall that the output of the `trajectory` function is a `StateSpaceSet`:
```@example MAIN
X
```

It is printed like a matrix where each column is the timeseries of each dynamic variable. In reality, it is a vector of statically sized vectors (for performance reasons). When iterated or indexed with 1 index, it behaves like a vector of vectors.
```@example MAIN
X[1]
```

```@example MAIN
X[2:5]
```

When indexed with two indices, it behaves like a matrix.

```@example MAIN
A[2:5, 2]
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
```

## Using state space sets



## More nonlinear timeseries analysis

- DelayEmbeddings.jl


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

To learn more, you need to visit the documentation pages of the individual packages. See the [contents](@ref contents) page for more!