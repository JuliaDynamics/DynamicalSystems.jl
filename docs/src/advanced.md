# Advanced documentation
This section overviews the various integrators available from `DynamicalSystemsBase`, as well as gives some insight into the internals, so that other developers that want to use this library can build upon it.

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
to change states robustly

## Re-initializing an integrator
It is more efficient to re-initialize an integrator than to create a new one.
This can be very helpful when looping over initial conditions and/or parameter values.

All high-level functions from `ChaosTools` have a set-up part that creates
an integrator, and a low-level part that does the computation.

The low level part is your friend! Use it! See the [Using `gali`](chaos/chaos_detection/#using-gali) page for an example.

### Discrete
For discrete systems, there is no "re-initialization". Just change the
state or deviation vectors on the spot using `set_state!` or `set_deviations!`.
E.g.
```julia
ds = Systems.towel()
integ = integrator(ds)
set_state!(integ, rand(3))
tinteg = tangent_integrator(ds, 3)
set_state!(tinteg, rand(3))
```
etc, etc.

To change parameters simply change the field `p` of an integrator.
```julia
ds = Systems.henon()
pinteg = parallel_integrator(ds, [rand(2), rand(2)])
pinteg.p[2] = 0.45
```

### Continuous
For continuous systems one needs to properly re-initialize the integrator instance
so that derivatives are re-computed. We still advise to use `set_state!` and
`set_deviations!` though, so you do not have to deal with the fact that
different integrators have different types of `u`.

You should do something like
```julia
ds = Systems.lorenz()
integ = integrator(ds)
set_state!(integ, rand(3))
reinit!(integ, integ.u)
```
*(it is important to use `integ.u` as the second argument to `reinit!`)*

The above code would work if the integrator was `parallel_integrator` or `tangent_integrator`!

The full documentation for `reinit!(::ODEIntegrator)` is [here](http://docs.juliadiffeq.org/latest/basics/integrator.html#Reinit-1). Although,
for usage within **DynamicalSystems.jl** the other arguments do not matter, because
steps are never saved.

In the continuous case you **must** `reinit!` even after changing a parameter value,
because the derivatives need to be re-computed.

## Implementation of `DynamicalSystem`
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
    # 1. prob
    # 2. jacobian (function)
    # 3. J (matrix)  <- will allow Sparse implementation in the future
end
```
This allows easily using multiple dispatch on the first three type parameters,
which are the most important for dispatching purposes.

The final type-parameter `IAD` is useful when creating the `tangent_integrator`,
so that the vector field is not computed twice!
