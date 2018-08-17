# Advanced documentation
This section overviews the various integrators available from `DynamicalSystemsBase`, as well as gives some insight into the internals, so that other developers that want to use this library can build upon it.

# ADD CALLBACK EXAMPLE
# ADD SPARSE JACOBIAN IN CSM

## Integrators
```@docs
integrator
parallel_integrator
tangent_integrator
```
---
Notice that the state type `integrator.u` of each integrator is quite different and *does change* between the possible versions of a [`DynamicalSystem`](@ref)!

## Integrator state functions
There are four functions associated with the integrators that we export:
```@docs
get_state
set_state!
get_deviations
set_deviations!
```
These functions work with *any* possible integrator and it is best to use the
to change states robustly!

## Re-initializing an integrator
It is more efficient to re-initialize an integrator using `reinit!`
than to create a new one.
This can be very helpful when looping over initial conditions and/or parameter values.

All high-level functions from `ChaosTools` have a set-up part that creates
an integrator, and a low-level part that does the computation.

The low level part is your friend! Use it! See the [Using `gali`](chaos/chaos_detection/#using-gali) page for an example.

The `reinit!` call signature is the same for continuous and discrete systems.
In the following, `state` is supposed to be a `D` dimensional vector (state of the dynamical system).

1. `reinit!(integ, state)` : to be used with standard [`integrator`](@ref).
3. `reinit!(integ, Vector_of_states)` : to be used with the [`parallel_integrator`](@ref).
2. `reinit!(integ, state, Q0::AbstractMatrix)` : to be used with the [`tangent_integrator`](@ref). This three argument version of `reinit!` is exported from `DynamicalSystemsBase`.

See below for a couple of examples.

### Re-init of continuous tangent integrator
Here we compute the [`lyapunovs`](@ref) for many different initial conditions.
```julia
ds = Systems.lorenz()
tinteg = tangent_integrator(ds, 2)
ics = [rand(3) for i in 1:100]
for ic in ics
  reinit!(tinteg, ic, orthonormal(3, 2))
  λ = lyapunovs(tinteg, 1000, 0.1, 10.0)
  # reminder: lyapunovs(tinteg, N, dt::Real, Ttr::Real = 0.0)
end
```


### Re-init of discrete parallel integrator
Here we compute the [`lyapunov`](@ref) for many different parameters.
```julia
ds = Systems.henon()
u0 = rand(SVector{2})
ps = 1.2:0.01:1.4
pinteg = parallel_integrator(ds, [u0, u0 + 1e-9rand(SVector{2})])
for p in ps
  set_parameter!(ds, 1, p)
  reinit!(pinteg, [u0, u0 + 1e-9rand(SVector{2})])
  λ = lyapunov(pinteg, 1000, 10, 1, 1e-9, 1e-6, 1e-12)
  # reminder: lyapunov(pinteg, T, Ttr, dt, d0, ut, lt)
end
```

## `DynamicalSystem` implementation
```julia
abstract type DynamicalSystem{
        IIP,     # is in place , for dispatch purposes and clarity
        S,       # state type
        D,       # dimension
        F,       # equations of motion
        P,       # parameters
        JAC,     # jacobian
        JM,      # jacobian matrix
        IAD}     # is auto-differentiated
    # one-liner: {IIP, S, D, F, P, JAC, JM, IAD}
    # Subtypes of DynamicalSystem have fields:
    # 1. f
    # 2. u0
    # 3. p
    # 4. t0
    # 5. jacobian (function)
    # 6. J (matrix)
end
```
The `DynamicalSystem` stores only the absolutely necessary information. Every other functionality of **DynamicalSystems.jl** initializes an integrator.

The final type-parameter `IAD` is useful when creating the `tangent_integrator`, so that the vector field is not computed twice!
