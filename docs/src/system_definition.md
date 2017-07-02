# System Definition
For `DynamicalSystems.jl` a system is simple a structure that contains the system's state, the equations of motion and the Jacobian. The last two are *functions* that take as an input a state.

By taking advantage of the package [`ForwardDiff.jl`](https://github.com/JuliaDiff/ForwardDiff.jl) an automated Jacobian
function can always be supplemented by the package. More details are enclosed in the indivivdual sections, however the documentation strings of all the constructors are
also self-contained.

!!! warning "Non-autonomous systems"
    This package does **not** accept non-autonomous systems. To use such systems with this package increase
    the dimensionality of your system by 1, by introducing an additional variable
    ``\tau`` such that ``d\tau/dt = 1`` (or ``\tau_{n+1} = \tau_n + 1``). This additional variable will serve as
    the "time" in your equations of motion.

## Discrete Systems
Discrete systems are of the form: $\vec{x}_{n+1} = \vec{f}(\vec{x}_n)$.
The Type representing such a system is called `DiscreteDS` and it is immutable. The reason for the choice of immutable type is simply speed: it is faster than the mutable
when evolving it.

The constructor is:
```julia
DiscreteDS(state, eom [, jacob])
```
Here `state` is simply the state the system starts (a.k.a. initial conditions) and
`eom` is a *function* that takes a state ``\vec{x}_{n}`` as an input and returns the next state ``\vec{x}_{n+1}`` as an output.

The `jacob` is also a *function* that takes a state as an input and returns the
Jacobian matrix of the system. This however is optional and if not provided by the user, will be calculated automatically using the package `ForwardDiff.jl`.

!!! tip "Return form of the `eom` function"
    It is **heavilty** advised that the equations of motion `eom` function returns an `SVector` from
    the julia package [`StaticArrays.jl`](https://www.google.de/search?q=julia+staticarrays&ie=utf-8&oe=utf-8&client=firefox-b-ab&gfe_rd=cr&ei=J0dSWdXTObLPXqvhj9AE) and similarly the `jacob` function returns an `SMatrix`. [Numerous benchmarks](https://github.com/Datseris/DynamicalSystems.jl/tree/master/test/benchmarks) have been made in order to deduce the most efficient possible way to define
    a system, and this way was proved to be the best.

### 1-dimensional Discrete Systems
asdf
## Continuous Systems
asf
## Pre-Defined Systems
asf
## Numerical Data
asf
## Convenience Functions
asf
