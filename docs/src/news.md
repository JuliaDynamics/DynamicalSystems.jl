# News
## Timeseries Prediction
A new module has been added to **DynamicalSystems.jl**: `TimeseriesPrediction` (version `v0.2.0`), which
tries to predict timeseries using methods from nonlinear dynamics and chaos!

The first available method is `localmodel_tsp` that uses local averages! See
the new documentation page for more!

## Cao's Method
With `ChaosTools v0.8.0`, the well-known method for estimating dimension for a
[`Reconstruction`](@ref) is now implemented and exported!

## Multi-time, Multi-Diensional Reconstructions
With the latest version of `DynamicalSystemsBase v0.8.0` we now have the possibility
for both multi-time and multi-dimensional delay reconstructions! The new documentation
string for [`Reconstruction`](@ref) has all the relevant information.

## Release v0.11.0
With new version v0.11.0 for **DynamicalSytems.jl** (`DynamicalSystemsBase` and `ChaosTools` at version 0.6) we have some major improvements of the library all around. Here I list the syntactic changes, the internal changes, the prospect for other developers and the gains we have from making all these changes!

### Syntax changes
There are no syntax changes (or *any* changes) to functions that handle numerical
data (`Dataset` and the likes).

The syntax for both discrete and continuous systems has changed to
```julia
DiscreteDynamicalSystem(eom, state, p [, jacobian [, J0]]; t0::Int = 0)
ContinuousDynamicalSystem(eom, state, p [, jacobian [, J0]]; t0 = 0.0)
```
The equations of motion function `eom` can be one of two forms:

* **iip** : The `eom` **must** be in the form `eom(x, p, t) -> SVector`
  which means that given a state `x::SVector` and some parameter container
  `p` it returns an [`SVector`](http://juliaarrays.github.io/StaticArrays.jl/stable/pages/api.html#SVector-1)
  containing the next state.
* **oop** : The `eom` **must** be in the form `eom!(xnew, x, p, t)`
  which means that given a state `x::Vector` and some parameter container `p`,
  it writes in-place the new state in `xnew`.

There is no constructor that takes an `ODEProblem` anymore.

In addition, `DynamicalSystem` and subtypes are now **immutable**. One cannot set
their state in place, or anything like that. Instead, all high-level functions
allow you to choose an initial state.

In summary:

* All discrete systems are now simply `DiscreteDynamicalSystem`.
* Continuous systems have been renamed to `ContinuousDynamicalSystem`.
* Don't use `set_state!`, etc. Instead use the keyword argument `u0` of
  methods like e.g. `gali`.

There are some syntax changes to high-level functions from `ChaosTools` as well.
For example, `lyapunovs` now has call signature
```julia
lyapunovs(ds::DynamicalSystem, N, k::Int | Q0; kwargs...) -> λs
```
It is advised to first look the documentation string of a function you want to use
before usage!

### Internal changes & Prospects
The internals of `DynamicalSystemsBase` have been completely re-worked from the ground up.
Here are the highlights:

* All [`DynamicalSystem`](@ref) objects are immutable, and contain a problem `prob`
  the jacobian function and an initialized Jacobian matrix.
* All functions that use a `DynamicalSystem` have changed behavior.
  The way the functions work now is that
  when given a `DynamicalSystem` the create an appropriate *integrator* from it.
  Then we use `step!(integrator, dt)` and use the integrator state to perform
  calculations.
* Eight possible system types are now available:
    * Continuous or Discrete.
    * In-place or out-of-place (large versus small systems).
    * Auto-differentiated or not (for the Jacobian function).
    * This is only possible due to the strictness of defining the `eom` function.
    * Robust multiple dispatch on all 8 types (again, only possible due to the strictness of the `eom` function).


* Three low-lever integrator constructing functions are available, that only need
  a `DynamicalSystem` and (optionally) an initial state:
    1. [`integrator`](@ref) for "normal" integration of a system.
    2. [`tangent_integrator`](@ref) for integration of a state and deviation vectors (that live on the tangent space).
    3. [`parallel_integrator`](@ref) for integrating in "parallel" a number of states. Notice that the states are integrated at *exact* same times, even for continuous systems.


* All three of the above integrators work perfectly fine for all eight combinations
  of systems and also have performant implementations.
* Simple internal implementation for Discrete system integrator that is tailored
  to the needs of **DynamicalSystems.jl**. It follows the high-level syntax of DifferentialEquations.jl: there is
  an implementation of a minimal discrete problem, as well as a minimal discrete
  integrator that steps via `step!(integrator, dt)`.
* `DynamicalSystem` is type-parameterized in such a way that allows for easy multiple
  dispatch and clear source code.

Because the resulting behavior is very robust and efficient, it allows
`DynamicalSystemsBase` to also be a library used by other developers that want to
develop techniques/functions for dynamical systems.

In addition, there is absolutely no drawback in having a `ContinuousDynamicalSystem`
instead of an `ODEProblem`, since the field `.prob` of the system is exactly this
`ODEProblem`, and can be passed directly to things like `solve` from DifferentialEquations.jl.

### Gains

* Out-of-place continuous systems are now possible!
* Auto-differentiated methods compute the vector field only once.
* Safe, robust implementations due to the immutability of the central structure `DynamicalSystem`.
* No problems with parallelization/threading/etc.
* Even clearer source code! Most `ChaosTools` functions are now composed of a
  set-up part and an implementation part, both of which are clear to read and understand.
      * Also clarity on discrete systems, since they are all fused into one structure.
      * Low-level functions can be used easily by users that want performance for loops.

* Lyapunov exponent calculating functions now have full flexibility in all aspects
  (initial deviation vectors/transient times/pretty much anything).
* Big performance gains all around, and especially in methods that propagate tangent space.
  For example, the function that calculates the Lyapunov spectrum of the towel map, needs:
  ```julia
  using DynamicalSystems
  ds = Systems.towel()
  l = lyapunovs(ds, 1000)
  ```
  ```
  Float64[3]
  0.42883526635723973
  0.36501911701374234
  -3.2835393321781092
  ```
  ```julia
  using BenchmarkTools
  @btime lyapunovs($ds, 1000);
  228.265 μs (176 allocations: 11.28 KiB)
  ```

### Still need to be done
We still need to create higher level functions like `set_state!` or `set_deviations!`
that set the state or deviation vectors on the integrator. These functions
will use multiple dispatch and thus work for all 24 combinations.

At the moment a "sloppy" implementation for everything is present in the source code
of `ChaosTools`, but this can be massively reduced into well-thought functions
and multiple dispatch usage.
