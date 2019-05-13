# News
## Predictability of a dynamical system in v1.3
See the function [`predictability`](@ref)!

## New default solver in v1.2; `OrdinaryDiffEq` dependency dropped!
In the newest version of **DynamicalSystems.jl**, which includes `DynamicalSystemsBase v1.2`, the default solver for all continuous systems has been changed to **`SimpleATsit5`** from the package [SimpleDiffEq.jl](https://github.com/JuliaDiffEq/SimpleDiffEq.jl). The previous default solver was `Vern9` from `OrdinaryDiffEq` This has reduced the "first run time" massively for functions that use a `ContinuousDynamicalSystem`!

In addition the `OrdinaryDiffEq` dependency was dropped. This is a big benefit for precompilation times! This does *not* mean that you can't use all the solvers from `OrdinaryDiffEq`! You can still use any solver you want, provided that you do `using OrdinaryDiffEq` to access the solvers!

!!! warning "Numeric values will change slightly"
    Although there was no API change the default solver did change. This means that if you were using the default solver for some code, the numeric values you will now obtain will slightly change as a result. This does not mean that any functionality broke, but be aware that if you were depending on the *exact value* of e.g. `lyapunovs`, you will now have a different value instead.

    Provided that you have already been using long enough integration times, the convergence of your results has not been affected though.

## `RecurrenceAnalysis` joins **DynamicalSystems.jl** in v1.1!
The excellent Julia package `RecurrenceAnalysis` (authored by Helios de Rosario, `@heliosdrm`) is now part of **DynamicalSystems.jl** and reexported by it. Besides adding all the amazing functionality of `RecurrenceAnalysis` to **DynamicalSystems.jl**, other benefits also sprung up:
* New package `DelayEmbeddings` that defines `Dataset` and provides all delay embedding functionality.
* Method to estimate delay embedding dimension (now in `DelayEmbeddings`) was enriched with two more algorithms!
* Mutual information (currently wip) is also being improved.
* `RecurrenceAnalysis` now also has a documentation page, soon to be updated with more real world examples.

## DynamicalSystems v1.0 - Julia 1.0 release
All support for any Julia version before 0.7 is dropped.

Please be sure to check out the `CHANGELOG.md` files of the individual repositories. There all changes are listed in detail. Here we note only the most important ones.

* `TimeseriesPrediction` is *not* installed with `DynamicalSystems` for version 1.0, because it is undergoing major changes. It is not even ready for Julia 1.0 actually.

* `DynamicalSystem` has been totally reworked for better clarity: it does not store a "problem" anymore, only the absolutely necessary ingredients to create one. The API did not change though!
* `Reconstruction` has been renamed to `reconstruct`, and now always returns a `Dataset`. In addition, now the parameter `D` (now renamed to `γ`) stands for the number of temporal neighbors. **This is a breaking change!**. The change allows more intuition across the different versions of `reconstruct`.
* The various offered integrators became more robust, and now allow passing callbacks etc. for the DifferentialEquations.jl event handling.
* Brand new algorithm for computing Poincare surfaces of section. It is not also more clear and understandable from the old one, but also much faster as well!!!
* Mutual information computation method. Also new method for optimal delay time using the Mutual information!

As always: be sure to read the documentation string before using a function: the docstrings are always updated and will show latest changes even if we (mistakenly) missed them when writing the documentation pages!


## Timeseries Prediction
A new module has been added to **DynamicalSystems.jl**: `TimeseriesPrediction` (version `v0.2.0`), which
tries to predict timeseries using methods from nonlinear dynamics and chaos!

The first available method is `localmodel_tsp` that uses local averages! See
the new documentation page for more!

## Cao's Method
With `ChaosTools v0.8.0`, the well-known method for estimating dimension for a
`Reconstruction` is now implemented and exported! See `estimate_dimension`.

## Multi-time, Multi-Diensional Reconstructions
With the latest version of `DynamicalSystemsBase v0.8.0` we now have the possibility
for both multi-time and multi-dimensional delay reconstructions! The new documentation
string for `Reconstruction` has all the relevant information.

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

* All `DynamicalSystem` objects are immutable, and contain a problem `prob`
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
    1. `integrator` for "normal" integration of a system.
    2. `tangent_integrator` for integration of a state and deviation vectors (that live on the tangent space).
    3. `parallel_integrator` for integrating in "parallel" a number of states. Notice that the states are integrated at *exact* same times, even for continuous systems.


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
