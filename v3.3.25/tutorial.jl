# # [Overarching tutorial for **DynamicalSystems.jl**](@id tutorial)

# This page serves as a short but to-the-point introduction to the **DynamicalSystems.jl**
# library. It outlines the core components, and how they establish an interface that
# is used by the rest of the library. It also provides a couple of usage examples
# to connect the various packages of the library together.

# Going through this tutorial should take you about 20 minutes.

# !!! note "Also available as a Jupyter notebook"
#     This tutorial is also available online as a [Jupyter notebook](https://github.com/JuliaDynamics/DynamicalSystems.jl/blob/gh-pages/dev/tutorial.ipynb).

# ## Installation

# To install **DynamicalSystems.jl**, simply do:
# ```julia
# using Pkg; Pkg.add("DynamicalSystems")
# ```

# This installs several packages for the Julia language. These are the sub-modules/packages
# that comprise **DynamicalSystems.jl**, see [contents](@ref contents) for more.
# All of the functionality is brought into scope when doing:

using DynamicalSystems

# in your Julia session.

# ### Package versions used

import Pkg

#nb # Activate an environment in the folder containing the notebook
#nb Pkg.activate(dirname(@__DIR__))
#nb Pkg.add(["DynamicalSystems", "CairoMakie", "GLMakie", "OrdinaryDiffEq", "BenchmarkTools"])
Pkg.status(["DynamicalSystems", "CairoMakie", "GLMakie", "OrdinaryDiffEq", "BenchmarkTools"]; mode = Pkg.PKGMODE_MANIFEST)

#nb # ## **DynamicalSystems.jl** summary

#nb @doc DynamicalSystems

# ## Core components

# The individual packages that compose `DynamicalSystems` interact flawlessly with each other because of the following two components:

# 1. The [`StateSpaceSet`](@ref), which represents numerical data. They can be observed or measured from experiments, sampled trajectories of dynamical systems, or just unordered sets in a state space. A `StateSpaceSet` is a container of equally-sized points, representing multivariate timeseries or multivariate datasets. Timeseries, which are univariate sets, are represented by the `AbstractVector{<:Real}` Julia base type.
# 2. The [`DynamicalSystem`](@ref), which is the abstract representation of a dynamical system with a known dynamic evolution rule. `DynamicalSystem` defines an extendable interface, but typically one uses existing implementations such as [`DeterministicIteratedMap`](@ref) or [`CoupledODEs`](@ref).

# ## Making dynamical systems

# In the majority of cases, to make a dynamical system one needs three things:

# 1. The dynamic rule `f`: A Julia function that provides the instructions of how to evolve the dynamical system in time.
# 2. The state `u`: An array-like container that contains the variables of the dynamical system and also defines the starting state of the system.
# 3. The parameters `p`: An arbitrary container that parameterizes `f`.

# For most concrete implementations of `DynamicalSystem` there are two ways of defining `f, u`.
# The distinction is done on whether `f` is defined as an in-place (iip) function or out-of-place (oop) function.

# * **oop** : `f` **must** be in the form `f(u, p, t) -> out`
#     which means that given a state `u::SVector{<:Real}` and some parameter container
#     `p` it returns the output of `f` as an `SVector{<:Real}` (static vector).
# * **iip** : `f` **must** be in the form `f!(out, u, p, t)`
#     which means that given a state `u::AbstractArray{<:Real}` and some parameter container `p`,
#     it writes in-place the output of `f` in `out::AbstractArray{<:Real}`.
#     The function **must** return `nothing` as a final statement.

# `t` stands for current time in both cases.
# **iip** is suggested for systems with high dimension and **oop** for small.
# The break-even point is between 10 to 100 dimensions but should be benchmarked
# on a case-by-case basis as it depends on the complexity of `f`.

# !!! note "Autonomous vs non-autonomous systems"
#     Whether the dynamical system is autonomous (`f` doesn't depend on time) or not, it is still necessary to include `t` as an argument to `f`. Some algorithms utilize this information, some do not, but we prefer to keep a consistent interface either way.


# ### Example: Henon map

# Let's make the Henon map, defined as
# ```math
# \begin{aligned}
# x_{n+1} &= 1 - ax^2_n+y_n \\
# y_{n+1} & = bx_n
# \end{aligned}
# ```
# with parameters $a = 1.4, b = 0.3$.

# First, we define the dynamic rule as a standard Julia function. Since the dynamical system is only two-dimensional, we should use the _out-of-place_ form that returns an `SVector` with the next state:

using DynamicalSystems

function henon_rule(u, p, n) # here `n` is "time", but we don't use it.
    x, y = u # system state
    a, b = p # system parameters
    xn = 1.0 - a*x^2 + y
    yn = b*x
    return SVector(xn, yn)
end

# Then, we define initial state and parameters

u0 = [0.2, 0.3]
p0 = [1.4, 0.3]

# Lastly, we give these three to the `DeterministicIteratedMap`:

henon = DeterministicIteratedMap(henon_rule, u0, p0)

# `henon` is a `DynamicalSystem`, one of the two core structures of the library.
# They can evolved interactively, and queried, using the interface defined by [`DynamicalSystem`](@ref).
# The simplest thing you can do with a `DynamicalSystem` is to get its trajectory:

total_time = 10_000
X, t = trajectory(henon, total_time)
X

# `X` is a `StateSpaceSet`, the second of the core structures of the library.
# We'll see below how, and where, to use a `StateSpaceset`, but for now let's just do a scatter plot

using CairoMakie
scatter(X)

# ### Example: Lorenz96

# Let's also make another dynamical system, the Lorenz96 model:
# ```math
# \frac{dx_i}{dt} = (x_{i+1}-x_{i-2})x_{i-1} - x_i + F
# ```
# for $i \in \{1, \ldots, N\}$ and $N+j=j$.

# Here, instead of a discrete time map we have $N$ coupled ordinary differential equations. However, creating the dynamical system works out just like above, but using `CoupledODEs` instead of `DeterministicIteratedMap`.

# First, we make the dynamic rule function. Since this dynamical system can be arbitrarily high-dimensional, we prefer to use the _in-place_ form for `f`, overwriting in place the rate of change in a pre-allocated container. It is [customary](https://docs.julialang.org/en/v1/manual/style-guide/#bang-convention) to append the name of functions that modify their arguments in-place with a bang (`!`).

function lorenz96_rule!(du, u, p, t)
    F = p[1]; N = length(u)
    ## 3 edge cases
    du[1] = (u[2] - u[N - 1]) * u[N] - u[1] + F
    du[2] = (u[3] - u[N]) * u[1] - u[2] + F
    du[N] = (u[1] - u[N - 2]) * u[N - 1] - u[N] + F
    ## then the general case
    for n in 3:(N - 1)
        du[n] = (u[n + 1] - u[n - 2]) * u[n - 1] - u[n] + F
    end
    return nothing # always `return nothing` for in-place form!
end

# then, like before, we define an initial state and parameters, and initialize the system

N = 6
u0 = range(0.1, 1; length = N)
p0 = [8.0]
lorenz96 = CoupledODEs(lorenz96_rule!, u0, p0)

# and, again like before, we may obtain a trajectory the same way

total_time = 12.5
sampling_time = 0.02
Y, t = trajectory(lorenz96, total_time; Ttr = 2.2, Δt = sampling_time)
Y

# We can't scatterplot something 6-dimensional but we can visualize all timeseries

fig = Figure()
ax = Axis(fig[1, 1]; xlabel = "time", ylabel = "variable")
for var in columns(Y)
    lines!(ax, t, var)
end
fig

# ### ODE solving and choosing solver

# Continuous time dynamical systems are evolved through DifferentialEquations.jl.
# In this sense, the above `trajectory` function is a simplified version of `DifferentialEquations.solve`.
# If you only care about evolving a dynamical system forwards in time, you are probably better off using
# DifferentialEquations.jl directly. **DynamicalSystems.jl** can be used to do many other things that either occur during
# the time evolution or after it, see the section below on [using dynamical systems](@ref using).

# When initializing a `CoupledODEs` you can tune the solver properties to your heart's
# content using any of the [ODE solvers](https://diffeq.sciml.ai/latest/solvers/ode_solve/)
# and any of the [common solver options](https://diffeq.sciml.ai/latest/basics/common_solver_opts/).
# For example:

using OrdinaryDiffEq: Vern9 # accessing the ODE solvers
diffeq = (alg = Vern9(), abstol = 1e-9, reltol = 1e-9)
lorenz96_vern = ContinuousDynamicalSystem(lorenz96_rule!, u0, p0; diffeq)

#

Y, t = trajectory(lorenz96_vern, total_time; Ttr = 2.2, Δt = sampling_time)
Y[end]

# The choice of the solver algorithm can have **huge impact on the performance and stability of the ODE integration!**
# We will showcase this with two simple examples

# #### Higher accuracy, higher order

# The solver `Tsit5` (the default solver) is most performant when medium-low error
# tolerances are requested. When we require very small error tolerances, choosing a different solver
# can be more accurate. This can be especially impactful for chaotic dynamical systems.
# Let's first expliclty ask for a given accuracy when solving the ODE by passing the
# keywords `abstol, reltol` (for absolute and relative tolerance respectively),
# and compare performance to a naive solver one would use:

using BenchmarkTools: @btime
using OrdinaryDiffEq: BS3 # 3rd order solver

for alg in (BS3(), Vern9())
    diffeq = (; alg, abstol = 1e-12, reltol = 1e-12)
    lorenz96 = CoupledODEs(lorenz96_rule!, u0, p0; diffeq)
    @btime step!($lorenz96, 100.0) # evolve for 100 time units
end

# The performance difference is dramatic!

# #### Stiff problems

# A "stiff" ODE problem is one that can be numerically unstable unless the step size (or equivalently, the step error tolerances) are extremely small. There are several situations where a problem may be come "stiff":

# - The derivative values can get very large for some state values.
# - There is a large _timescale separation_ between the dynamics of the variables
# - There is a large _speed separation_ between different state space regions

# One must be aware whether this is possible for their system and choose a solver that is better suited to tackle stiff problems. If not, a solution may diverge and the ODE integrator will throw an error or a warning.

# Many of the problems in DifferentialEquations.jl are suitable for dealing with stiff problems. We can create a stiff problem by using the well known Van der Pol  oscillator _with a timescale separation_:

# ```math
# \begin{aligned}
# \dot{x} & = y \\
# \dot{y} /  \mu &= (1-x^2)y - x
# \end{aligned}
# ```

# with $\mu$ being the timescale of the $y$ variable in units of the timescale of the $x$ variable. For very large values of $\mu$ this problem becomes stiff.

# Let's compare

using OrdinaryDiffEq: Tsit5, Rodas5P

function vanderpol_rule(u, μ, t)
    x, y = u
    dx = y
    dy = μ*((1-x^2)*y - x)
    return SVector(dx, dy)
end

μ = 1e6

for alg in (Tsit5(), Rodas5P()) # default vs specialized solver
    diffeq = (; alg, abstol = 1e-12, reltol = 1e-12, maxiters = typemax(Int))
    vdp = CoupledODEs(vanderpol_rule, SVector(1.0, 1.0), μ; diffeq)
    @btime step!($vdp, 100.0)
end

# We see that the stiff solver `Rodas5P` is much faster than the default `Tsit5` when there is a large timescale separation. This happened because `Rodas5P` required much less steps to integrated the same total amount of time. In fact, there are cases where regular solvers will _fail_ to integrate the ODE if the problem is very stiff, e.g. in the [ROBER example](https://docs.sciml.ai/SciMLBenchmarksOutput/stable/StiffODE/ROBER/).

# So using an appropriate solver really does matter!
# For more information on choosing solvers consult the DifferentialEquations.jl documentation.

# ## Interacting with dynamical systems

# The [`DynamicalSystem`](@ref) type defines an extensive interface for what it means to be a
# "dynamical system". This interface can be used to (1) define fundamentally new types
# of dynamical systems, (2) to develop algorithms that utilize
# dynamical systems with a known evolution rule. It can also be used to simply
# query and alter properties of a given dynamical system. For example, we have

lorenz96

# which we can evolve forwards in time using `step!`

step!(lorenz96, 100.0) # progress for `100.0` units of time

# and we can then query what is the current state that the dynamical system was brought into

current_state(lorenz96)

# we can also restart the system at a different state using `set_state!`

set_state!(lorenz96, rand(6))

# or we can alter system parameters given the index of the parameter and the value to set it to

set_parameter!(lorenz96, 1, 9.6) # change first parameter of the parameter container
current_parameters(lorenz96)

# For more functions that query or alter a dynamical system see its docstring: [`DynamicalSystem`](@ref).

# ## [Using dynamical systems](@id using)

# Now, as an end-user, you are most likely to be giving a `DynamicalSystem` instance to a library function.
# For example, you may want to compute the Lyapunov spectrum of the Lorenz96 system from above,
# which is a functionality offered by [`ChaosTools`](@ref).
# This is as easy as calling the `lyapunovspectrum` function with `lorenz96`

steps = 10_000
lyapunovspectrum(lorenz96, steps)

# As expected, there is at least one positive Lyapunov exponent, because the system is chaotic,
# and at least one zero Lyapunov exponent, because the system is continuous time.

# A fantastic feature of **DynamicalSystems.jl** is that all library functions work for any
# applicable dynamical system. The exact same `lyapunovspectrum` function would also work
# for the Henon map.

lyapunovspectrum(henon, steps)

# Something else that uses a dynamical system is estimating the basins of attraction
# of a multistable dynamical system.
# The Henon map is "multistable" in the sense that some initial conditions diverge to
# infinity, and some others converge to a chaotic attractor.
# Computing these basins of attraction is simple with [`Attractors`](@ref), and would work as follows:

## define a state space grid to compute the basins on:
xg = yg = range(-2, 2; length = 201)
## find attractors using recurrences in state space:
mapper = AttractorsViaRecurrences(henon, (xg, yg); sparse = false)
## compute the full basins of attraction:
basins, attractors = basins_of_attraction(mapper; show_progress = false)

# Let's visualize the result

heatmap_basins_attractors((xg, yg), basins, attractors)

# ## Stochastic systems

# DynamicalSystems.jl has some support for stochastic systems
# in the form of Stochastic Differential Equations (SDEs).
# Just like `CoupledODEs`, one can make `CoupledSDEs`!
# For example here is a stochastic version of a FitzHugh-Nagumo model

using StochasticDiffEq # load extention for `CoupledSDEs`

function fitzhugh_nagumo(u, p, t)
    x, y = u
    ϵ, β, α, γ, κ, I = p
    dx = (-α * x^3 + γ * x - κ * y + I) / ϵ
    dy = -β * y + x
    return SVector(dx, dy)
end
p = [1.,3.,1.,1.,1.,0.]
sde = CoupledSDEs(fitzhugh_nagumo, zeros(2), p; noise_strength = 0.05)

# In this particular example the SDE noise is white noise (Wiener process)
# with strength (σ) of 0.05. See the documentation of `CoupledSDEs` for alternatives.

# In any case, in DynamicalSystems.jl all dynamical systems are part of the same
# interace, stochastic or not. As long as the algorithm is not influenced by stochasticity,
# we can apply it to `CoupledSDEs` just as well. For example, we can study multistability
# in a stochastic system. In contrast to the previous example of the Henon map,
# we have to use an alternative algorithm, because `AttractorsViaRecurrences`
# only works for deterministic systems. So instead we'll use `AttractorsViaFeaturizing`:

featurizer(X, t) = X[end]

mapper = AttractorsViaFeaturizing(sde, featurizer; Ttr = 200, T = 10)

xg = yg = range(-1, 1; length = 101)

sampler, _ = statespace_sampler((xg, yg))

fs = basins_fractions(mapper, sampler)

# and we can see the stored "attractors"

attractors = extract_attractors(mapper)
fig, ax = scatter(attractors[1])
scatter!(attractors[2])
fig

# The mathematical concept of attractors
# doesn't translate trivially to stochastic systems but thankfully
# this system has two fixed point attractors that are only mildly perturbed
# by the noise.

# ## Parallelization

# !!! warning "Dynamical systems are modified!"
#     It is not immediatelly obvious, but all library functions that obtain as an input a
#     `DynamicalSystem` instance will modify it, in-place. For example the `current_state`
#     of the system before and after giving it to a function such as `basins_of_attraction`
#     will not be the same! This also affects parallelization, see below.


# Since `DynamicalSystem`s are mutable, one needs to copy them before parallelizing,
# to avoid having to deal with complicated race conditions etc. The simplest way is with
# `deepcopy`. Here is an example block that shows how to parallelize calling some expensive
# function (e.g., calculating the Lyapunov exponent) over a parameter range
# (or alternatively, over different initial conditions) using `Threads`:


# ```julia
# ds = DynamicalSystem(f, u, p) # some concrete implementation
# parameters = 0:0.01:1
# outputs = zeros(length(parameters))
#
# # Since `DynamicalSystem`s are mutable, we need to copy to parallelize
# systems = [deepcopy(ds) for _ in 1:Threads.nthreads()-1]
# pushfirst!(systems, ds) # we can save 1 copy
#
# Threads.@threads for (i, p) in enumerate(parameters)
#     system = systems[Threads.threadid()]
#     set_parameter!(system, index, parameters[i])
#     outputs[i] = expensive_function(system, args...)
# end
# ```

# ## Interactive GUIs

# A particularly useful feature are interactive GUI apps one
# can launch to examine a `DynamicalSystem`. The simplest is [`interactive_trajectory_timeseries`](@ref).
# To actually make it interactive one needs to enable GLMakie.jl as a backend:

# ```julia
# import GLMakie
# GLMakie.activate!()
# ```

# and then launch the app:

u0s = [10rand(5) for _ in 1:3]
parameter_sliders = Dict(1 => 0:0.01:32)

fig, dsobs = interactive_trajectory_timeseries(
    lorenz96, [1, 2, 3, 4, 5], u0s;
    Δt = 0.02, parameter_sliders
)
step!(dsobs, 500) # hide
fig

# ### Developing new algorithms

# You could also be using a `DynamicalSystem` instance directly to build your own algorithm if it isn't already implemented (and then later contribute it so it _is_ implemented ;) ). A dynamical system can be evolved forwards in time using `step!`:

henon

# Notice how the time is not 0, because `henon` has already been stepped when we called the function `basins_of_attraction` with it. We can step it more:

step!(henon)

#

step!(henon, 2)

# For more information on how to directly use `DynamicalSystem` instances, see the documentation of [`DynamicalSystemsBase`](@ref).

# ## State space sets

# Let's recall that the output of the `trajectory` function is a `StateSpaceSet`:

X

# This is the main data structure used in **DynamicalSystems.jl** to handle numerical data.
# It is printed like a matrix where each column is the timeseries of each dynamic variable.
# In reality, it is a vector equally-sized vectors representing state space points.
# _(For advanced users: `StateSpaceSet` directly subtypes `AbstractVector{<:AbstractVector}`)_

# When indexed with 1 index, it behaves like a vector of vectors

X[1]

#

X[2:5]

# When indexed with two indices, it behaves like a matrix

X[7:13, 2] # 2nd column

# When iterated, it iterates over the contained points

for (i, point) in enumerate(X)
    @show point
    i > 5 && break
end

#

map(point -> point[1] + 1/(point[2]+0.1), X)

# The columns of the set are obtained with the convenience `columns` function

x, y = columns(X)
summary.((x, y))

# Because `StateSpaceSet` really is a vector of vectors, it can be given
# to any Julia function that accepts such an input. For example,
# the Makie plotting ecosystem knows how to plot vectors of vectors.
# That's why this works:

scatter(X)

# even though Makie has no knowledge of the specifics of `StateSpaceSet`.

# ## Using state space sets

# Several packages of the library deal with `StateSpaceSets`.

# You could use [`ComplexityMeasures`](@ref) to obtain the entropy, or other complexity measures, of a given set. Below, we obtain the entropy of the natural density of the chaotic attractor by partitioning into a histogram of approximately `50` bins per dimension:

prob_est = ValueHistogram(50)
entropy(prob_est, X)

# Or, obtain the permutation and sample entropies of the two columns of `X`:

pex = entropy_permutation(x; m = 4)
sey = entropy_sample(y; m = 2)
pex, sey

# Alternatively, you could use [`FractalDimensions`](@ref) to get the fractal dimensions of the chaotic attractor of the henon map using the Grassberger-Procaccia algorithm:

grassberger_proccacia_dim(X)

# Or, you could obtain a recurrence matrix of a state space set with [`RecurrenceAnalysis`](@ref)

R = RecurrenceMatrix(Y, 8.0)
Rg = grayscale(R)
rr = recurrencerate(R)
heatmap(Rg; colormap = :grays,
    axis = (title = "recurrence rate = $(round(rr; digits = 3))", aspect = 1)
)


# ## More nonlinear timeseries analysis

# A `trajectory` of a known dynamical system is one way to obtain a `StateSpaceSet`. However, another common way is via a delay coordinates embedding of a measured/observed timeseries. For example, we could use `optimal_separated_de` from [`DelayEmbeddings`](@ref) to create an optimized delay coordinates embedding of a timeseries

w = Y[:, 1] # first variable of Lorenz96
𝒟, τ, e = optimal_separated_de(w)
𝒟

# and compare

fig = Figure()
axs = [Axis3(fig[1, i]) for i in 1:2]
for (S, ax) in zip((Y, 𝒟), axs)
    lines!(ax, S[:, 1], S[:, 2], S[:, 3])
end
fig

# Since `𝒟` is just another state space set, we could be using any of the above analysis pipelines on it just as easily.

# The last package to mention here is [`TimeseriesSurrogates`](@ref), which ties with all other observed/measured data analysis by providing a framework for confidence/hypothesis testing. For example, if we had a measured timeseries but we were not sure whether it represents a deterministic system with structure in the state space, or mostly noise, we could do a surrogate test. For this, we use `surrogenerator` and `RandomFourier` from [`TimeseriesSurrogates`](@ref), and the `generalized_dim` from [`FractalDimensions`](@ref) (because it performs better in noisy sets)

x # Henon map timeseries
## contaminate with noise
using Random: Xoshiro
rng = Xoshiro(1234)
x .+= randn(rng, length(x))/100
## compute noise-contaminated fractal dim.
Δ_orig = generalized_dim(embed(x, 2, 1))

# And we do the surrogate test

surrogate_method = RandomFourier()
sgen = surrogenerator(x, surrogate_method, rng)
Δ_surr = map(1:1000) do i
    s = sgen()
    generalized_dim(embed(s, 2, 1))
end

# and visualize the test result

fig, ax = hist(Δ_surr)
vlines!(ax, Δ_orig)
fig

# since the real value is outside the distribution we have confidence the data are not pure noise.

# ## Integration with ModelingToolkit.jl

# **DynamicalSystems.jl** understands when a model has been generated via [ModelingToolkit.jl](https://docs.sciml.ai/ModelingToolkit/stable/). The symbolic variables used in ModelingToolkit.jl can be used to access the state or parameters of the dynamical system.

# To access this functionality, the `DynamicalSystem` must be created from a `DEProblem` of the SciML ecosystem, and the `DEProblem` itself must be created from a ModelingToolkit.jl model.

# !!! note "ProcessBasedModelling.jl"
#     ProcessBasedModelling.jl is an extension to ModelingToolkit.jl for creating
#     models from a set of equations. It has been designed to be useful for scenarios
#     applicable to a typical nonlinear dynamics analysis workflow,
#     and provides better error messages during system construction than MTK.
#     Have a look [at its docs](https://juliadynamics.github.io/ProcessBasedModelling.jl/stable/)!


# Let's create a the Roessler system as an MTK model:

using ModelingToolkit

@variables t # use unitless time
D = Differential(t)
@mtkmodel Roessler begin
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

# this model can then be made into an `ODEProblem`:

prob = ODEProblem(model)

# (notice that because we specified initial values for all parameters and variables during the model creation  we do need to provide additional initial values)

# Now, this problem can be made into a [`CoupledODEs`](@ref):

roessler = CoupledODEs(prob)

# This dynamical system instance can be used in the rest of the library like anything else.
# Additionally, you can "observe" referenced symbolic variables:

observe_state(roessler, model.x)

# or, more commonly, you can observe a `Symbol` that has the name of the symbolic variable:

observe_state(roessler, :nlt)

# These observables can also be used in the GUI visualization [`interactive_trajectory_timeseries`](@ref).

# You can also symbolically alter parameters

current_parameter(roessler, :c)

#

set_parameter!(roessler, :c, 5.0)

#

current_parameter(roessler, :c)

# This symbolic indexing can be given anywhere in the ecosystem where you would be
# altering the parameters, such as the function `global_continuation` from [`Attractors`](@ref).

# ## Core components reference

# ```@docs
# StateSpaceSet
# DynamicalSystem
# ```

# ## Dynamical system implementations

# ```@docs
# DeterministicIteratedMap
# CoupledODEs
# StroboscopicMap
# PoincareMap
# ProjectedDynamicalSystem
# ArbitrarySteppable
# CoupledSDEs
# ```

# ## Dynamical system interface
# ```@docs
# current_state
# initial_state
# observe_state
# state_name
# current_parameters
# current_parameter
# parameter_name
# initial_parameters
# isdeterministic
# isdiscretetime
# dynamic_rule
# current_time
# initial_time
# isinplace(::DynamicalSystem)
# successful_step
# referrenced_sciml_model
# ```

# ## Learn more

# To learn more, you need to visit the documentation pages of the modules that compose **DynamicalSystems.jl**. See the [contents](@ref contents) page for more!
