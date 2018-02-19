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

If you want to use `reinit!` (see below), first check the state type of the `integrator.u`.

For `tangent_integrator` the `u` is a matrix with the first column being the
actual system state, and all the rest being deviation vectors.

For `parallel_integrator` the `u` is a vector of vectors, where each vector is a
a system state. Only in the case of in-place continuous systems is the state actually a matrix with each column being a state, because DifferentialEquations at the moment does not support `Vector[Vector]` as integrator state.

For the "standard" integrator the `u` is always a vector of course.

## Re-initializing an integrator
It is much more efficient to re-inialize an integrator than to create a new one.
This can be very helpful when looping over initial conditions and/or parameter values.

All high-level functions from `ChaosTools` have a set-up part that creates
an integrator, and a low-level part that does the computation.

The low level part is your friend! Use it! See the [Using `gali`](chaos/chaos_detection/#using-gali) page for an example.

### Discrete
The `reinit!` signature is:
```julia
reinit!(integ::MinimalDiscreteIntegrator, u; t0 = integ.t0)
```
`reinit!` is not necessary when one changes parameters `p`. Just change
them in place, using `integ.p[index] = value`.

### Continuous
The call signature is identical, with a bunch of extra keyword arguments.
```julia
reinit!(integ::ODEIntegrator, u = integ.sol.prob.u0;  t0 = integ.sol.prob.tspan[1])
```
The full documentation for `reinit!(::ODEIntegrator)` is [here](http://docs.juliadiffeq.org/latest/basics/integrator.html#Reinit-1). Although,
for usage within **DynamicalSystems.jl** the other arguments do not matter, because
steps are never saved.

In the continuous case you should `reinit!` even after changing a parameter value,
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
