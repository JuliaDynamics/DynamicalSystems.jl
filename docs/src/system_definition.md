# System Definition
For `DynamicalSystems.jl` a system is simple a structure that contains the system's state, the equations of motion and the Jacobian. The last two are *functions* that take as an input a state. This page treats systems where one already **knows the equations of motion**.

By taking advantage of the package [`ForwardDiff.jl`](https://github.com/JuliaDiff/ForwardDiff.jl) an automated Jacobian
function can always be supplemented by the package. More details are enclosed in the indivivdual sections, however the documentation strings of all the constructors are
also self-contained.

!!! warning "Non-autonomous systems"
    This package does **not** accept non-autonomous systems. To use such systems with this package increase
    the dimensionality of your system by 1, by introducing an additional variable
    ``\tau`` such that ``d\tau/dt = 1`` (or ``\tau_{n+1} = \tau_n + 1``). This additional variable will serve as
    the "time" in your equations of motion.

## Discrete Systems
Discrete systems are of the form:
```math
\vec{x}_{n+1} = \vec{f}(\vec{x}_n).
```
The Type representing such systems for dimensionality `D ≤ 10` is called `DiscreteDS`.

The constructor is:
```julia
DiscreteDS(state, eom [, jacob])
```
Here `state` is simply the state the system starts (a.k.a. initial conditions) and
`eom` is a *function* that takes a `state` as an input and returns the next state
as an output.

The `jacob` is also a *function* that takes a `state` as an input and returns the
Jacobian matrix of the system (at this state). This however is optional and if not provided by the user, will be calculated automatically using the package `ForwardDiff.jl`.

!!! tip "Return form of the `eom` function for small `D`"
    It is **heavily** advised that the equations of motion `eom` function returns an `SVector` from
    the julia package [`StaticArrays.jl`](https://www.google.de/search?q=julia+staticarrays&ie=utf-8&oe=utf-8&client=firefox-b-ab&gfe_rd=cr&ei=J0dSWdXTObLPXqvhj9AE) and similarly the `jacob` function returns an `SMatrix`. [Numerous benchmarks](https://github.com/Datseris/DynamicalSystems.jl/tree/master/test/benchmarks) have been made in order to deduce the most efficient possible way to define
    a system, and this way was proved to be the best.

For example, let's create one of the [Predefined Systems](#Predefined-systems) offered by this package, the Hénon map:
```julia
using DynamicalSystems
using StaticArrays #only necessary when defining a system

eom_henon(x) = SVector{2}(1.0 - a*x[1]^2 + x[2], b*x[1])
jacob_henon(x) = @SMatrix [-2*a*x[1] 1.0; b 0.0]

ds = DiscreteDS(rand(2), eom_henon, jacob_henon)
```
If we did not want to write a Jacobian (due to e.g. unending laziness), we could
```julia
ds_nojac = DiscreteDS(rand(2), eom_henon)
```
and the Jacobian function would be created automatically.


### Large Discrete Systems
TBA

### 1-dimensional Discrete Systems
In the case of maps, there a special structure for one-dimensional systems, since
they are commonly used in scientific research. The syntax is `DiscreteDS1D(state, eom [, deriv])`. In this one-dimensional case, you don't need to worry about `StaticArrays.jl`
because everything is in plain numbers. For example:
```julia
using DynamicalSystems

@inline eom_logistic(r) = (x) -> r*x*(1-x)  # this is a closure
@inline deriv_logistic(r) = (x) -> r*(1-2x) # this is a closure
r = 3.7
logistic = DiscreteDS1D(rand(), eom_logistic(r), deriv_logistic(r))
```
Once again, if you skip the derivative functions it will be calculated automatically
using `ForwardDiff.jl`.

## Continuous Systems
TBA
## Predefined Systems
TBA
## Numerical Data
TBA
## Convenience Functions
TBA
