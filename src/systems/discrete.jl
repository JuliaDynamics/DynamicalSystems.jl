using StaticArrays, ForwardDiff, Requires

export DiscreteDS, DiscreteDS1D, evolve, evolve!, trajectory, dimension
export BigDiscreteDS

#####################################################################################
#                                   Constructors                                    #
#####################################################################################
"Abstract type representing discrete systems."
abstract type DiscreteDynamicalSystem <: DynamicalSystem end
"""
    DiscreteDS(state, eom [, jacob]; name="") <: DynamicalSystem
`D`-dimensional discrete dynamical system (used for `D ≤ 10`).
## Fields:
* `state::SVector{D}` : Current state-vector of the system, stored in the data format
  of `StaticArray`'s `SVector`.
* `eom` (function) : The function that represents the system's equations of motion
  (also called vector field). The function is of the format: `eom(u) -> SVector`
  which means that given a state-vector `u` it returns an `SVector` containing the
  next state.
* `jacob` (function) : A function that calculates the system's jacobian matrix,
  based on the format: `jacob(u) -> SMatrix` which means that given a state-vector
  `u` it returns an `SMatrix` containing the Jacobian at that state.
* `name::String` : A name for the dynamical system (possibly including parameter
  values), solely for pretty-printing purposes. Always passed to the constructors
  as a keyword.

If the `jacob` is not provided by the user, it is created automatically
using the module [`ForwardDiff`](http://www.juliadiff.org/ForwardDiff.jl/stable/).
"""
mutable struct DiscreteDS{D, T<:Number, F, J} <: DiscreteDynamicalSystem
    state::SVector{D,T}
    eom::F
    jacob::J
    name::String
end
# constructor without jacobian (uses ForwardDiff)
function DiscreteDS(u0::AbstractVector, eom;name="")
    su0 = SVector{length(u0)}(u0)
    cfg = ForwardDiff.JacobianConfig(eom, u0)
    @inline ForwardDiff_jac(x) = ForwardDiff.jacobian(eom, x, cfg)
    return DiscreteDS(su0, eom, ForwardDiff_jac, name)
end
function DiscreteDS(u0::AbstractVector, eom, jac;name="")
    D = length(u0)
    su0 = SVector{D}(u0)
    T = eltype(su0); F = typeof(eom); J = typeof(jac)
    return DiscreteDS{D, T, F, J}(su0, eom, jac,name)
end

"""
    DiscreteDS1D(state, eom [, deriv]; name="") <: DynamicalSystem
One-dimensional discrete dynamical system.
## Fields:
* `state::Real` : Current state of the system.
* `eom` (function) : The function that represents the system's equation of motion:
  `eom(x) -> Real`.
* `deriv` (function) : A function that calculates the system's derivative given
  a state: `deriv(x) -> Real`. If it is not provided by the user
  it is created automatically using the module
  [`ForwardDiff`](http://www.juliadiff.org/ForwardDiff.jl/stable/).
* `name::String` : A name for the dynamical system (possibly including parameter
  values), solely for pretty-printing purposes. Always passed to the constructors
  as a keyword.
"""
mutable struct DiscreteDS1D{S<:Real, F, D} <: DiscreteDynamicalSystem
    state::S
    eom::F
    deriv::D
    name::String
end
function DiscreteDS1D(x0, eom;name="")
    ForwardDiff_der(x) = ForwardDiff.derivative(eom, x)
    DiscreteDS1D(x0, eom, ForwardDiff_der,name)
end
DiscreteDS1D(a,b,c;name="")=DiscreteDS1D(a,b,c,name)

"""
    BigDiscreteDS(state, eom! [, jacob! [, J]]; name="") <: DynamicalSystem
`D`-dimensional discrete dynamical system (used for `D > 10`). The equations
for this system
perform all operations `in-place`.
## Fields:
* `state::Vector{T}` : Current state-vector of the system.
* `eom!` (function) : The function that represents the system's equations of motion
  (also called vector field). The function is of the format: `eom!(xnew, x)`
  which means that given a state-vector `x` and another similar one `xnew`,
  it writes in-place the new state in `xnew`.
* `jacob!` (function) : A function that calculates the system's jacobian matrix,
  based on the format: `jacob!(J, x)` which means that given a state-vector
  `x` it writes in-place the Jacobian in `J`.
* `J::Matrix{T}` : Initialized Jacobian matrix (optional). This field is not
  displayed.
* `dummystate::Vector{T}` : Dummy vector, which most of the time fills the
  role of the previous state in e.g. [`evolve`](@ref). This field is not
  displayed.
* `name::String` : A name for the dynamical system (possibly including parameter
  values), solely for pretty-printing purposes. Always passed to the constructors
  as a keyword.

If the `jacob!` is not provided by the user, it is created automatically
using the module [`ForwardDiff`](http://www.juliadiff.org/ForwardDiff.jl/stable/).
"""
mutable struct BigDiscreteDS{T<:Number, F, J} <: DiscreteDynamicalSystem
    state::Vector{T}
    eom!::F
    jacob!::J
    J::Matrix{T}
    dummystate::Vector{T}
    name::String
end
function BigDiscreteDS(u0, f!, j!,
    J::Matrix = zeros(eltype(u0), length(u0), length(u0)); name="")

    dum = copy(u0)
    BigDiscreteDS(u0, f!, j!, J, dum, name)
end
function BigDiscreteDS(u0, f!,
    J::Matrix = zeros(eltype(u0), length(u0), length(u0)); name="")
    dum = copy(u0)

    cfg = ForwardDiff.JacobianConfig(f!, dum, u0)
    FD_jacob!(J, x) = ForwardDiff.jacobian!(J, f!, dum, x, cfg)
    return BigDiscreteDS(u0, f!, FD_jacob!, J, dum, name)
end

dimension(::DiscreteDS{D, T, F, J}) where {D, T, F, J} = D
dimension(::DiscreteDS1D) = 1
dimension(ds::BigDiscreteDS) = length(ds.state)

#####################################################################################
#                               System Evolution                                    #
#####################################################################################
"""
    evolve(ds::DynamicalSystem, T=1; diff_eq_kwargs = Dict()) -> final_state
Evolve a `ds` for total "time" `T` and return the `final_state` (does not change
`ds.state`).
For discrete systems `T` corresponds to steps and
thus it must be integer. See [`trajectory`](@ref) for using `diff_eq_kwargs`.

This function *does not store* any information about intermediate steps.
Use [`trajectory`](@ref) if you want to produce a trajectory of the system.
If you want to
perform step-by-step evolution of a continuous system, use
`ODEIntegrator(ds, args...)` and
the `step!(integrator)` function provided by `DifferentialEquations`.

See also [`evolve!`](@ref).
"""
function evolve(ds::DiscreteDynamicalSystem, N::Int = 1)
    st = ds.state
    for i in 1:N
        st = ds.eom(st)
    end
    return st
end

function evolve(ds::BigDiscreteDS, N::Int = 1)
    st = copy(ds.state)
    for i in 1:N
        ds.dummystate .= st
        ds.eom!(st, ds.dummystate)
    end
    return st
end

"""
    evolve!(ds::DynamicalSystem, T; diff_eq_kwargs = Dict())
Same as [`evolve`](@ref), but also updates the system's `state` field with the final
state after evolution.
"""
function evolve!(ds::DiscreteDynamicalSystem, N::Int = 1)
    ds.state = evolve(ds, N)
    return ds
end

function evolve!(ds::BigDiscreteDS, N::Int = 1)
    for i in 1:N
        ds.dummystate .= ds.state
        ds.eom!(ds.state, ds.dummystate)
    end
    return ds
end


"""
```julia
trajectory(ds::DynamicalSystem, T; kwargs...) -> dataset
```
Return a dataset what will contain the trajectory of the sytem,
after evolving it for time `T`. See [`Dataset`](@ref) for info on how to
manipulate this object.

For the discrete case, `T` is an integer and a `T×D` dataset is returned
(`D` is the system dimensionality). For the
continuous case, a `W×D` dataset is returned, with `W = length(0:dt:T)` with
`0:dt:T` representing the time vector (*not* returned).
## Keywords:
* `dt = 0.05` : (only for continuous) Time step of value output during the solving
  of the continuous system.
* `diff_eq_kwargs = Dict()` : (only for continuous) A dictionary `Dict{Symbol, ANY}`
  of keyword arguments
  passed into the `solve` of the `DifferentialEquations.jl` package,
  for example `Dict(:abstol => 1e-9)`. If you want to specify a solver,
  do so by using the symbol `:solver`, e.g.:
  `Dict(:solver => DP5(), :maxiters => 1e9)`. This requires you to have been first
  `using OrdinaryDiffEq` to access the solvers.
"""
function trajectory(ds::DiscreteDS, N::Real)
    st = ds.state
    ts = [st]
    f = ds.eom
    for i in 2:N
        st = f(st)
        push!(ts, st)
    end
    return Dataset(ts)
end

function trajectory(ds::DiscreteDS1D, N::Int)
    x = deepcopy(ds.state)
    f = ds.eom
    ts = Vector{typeof(x)}(N)
    ts[1] = x
    for i in 2:N
        x = f(x)
        ts[i] = x
    end
    return ts
end

function trajectory(ds::BigDiscreteDS, N::Int)
    x = copy(ds.state)
    f! = ds.eom!
    ts = [zeros(eltype(x), dimension(ds)) for i in 1:N]
    ts[1] = x
    for i in 2:N
        ds.dummystate .= ts[i-1]
        f!(ts[i], ds.dummystate)
    end
    return Dataset(ts)
end





#####################################################################################
#                                Pretty-Printing                                    #
#####################################################################################
import Base.show
function Base.show(io::IO, ds::DiscreteDS{N, S, F, J}) where
    {N<:ANY, S<:ANY, F<:ANY, J<:ANY}
    if ds.name == ""
        text = "$(dimension(ds))-dimensional discrete system"
    else
        text = ds.name
    end
    print(io, text*"\n",
    " state: $(ds.state)\n", " e.o.m.: $F\n", " jacobian: $J")
end

@require Juno begin
    function Juno.render(i::Juno.Inline, ds::DiscreteDS{N, S, F, J}) where
        {N<:ANY, S<:ANY, F<:ANY, J<:ANY}
        t = Juno.render(i, Juno.defaultrepr(ds))
        if ds.name == ""
            text = "$(dimension(ds))-dimensional discrete system"
        else
            text = ds.name
        end
        t[:head] = Juno.render(i, Text(text))
        t[:children] = t[:children][1:3] # remove showing field dummystate
        t
    end
end


### Big Discrete
function Base.show(io::IO, ds::BigDiscreteDS{T, F, J}) where
    {T, F<:ANY, J<:ANY}
    if ds.name == ""
        text = "$(dimension(ds))-dimensional Big discrete system"
    else
        text = ds.name
    end
    print(io, text*"\n",
    " state: $(ds.state)\n", " e.o.m.: $F\n", " jacobian: $J")
end

@require Juno begin
    function Juno.render(i::Juno.Inline, ds::BigDiscreteDS{T, F, J}) where
        {T<:ANY, F<:ANY, J<:ANY}
        if ds.name == ""
            text = "$(dimension(ds))-dimensional Big discrete system"
        else
            text = ds.name
        end
        t = Juno.render(i, Juno.defaultrepr(ds))
        t[:head] = Juno.render(i, Text(text))
        t[:children] = t[:children][1:3] # remove showing field dummystate
        t
    end
end

### 1-D
function Base.show(io::IO, s::DiscreteDS1D{S, F, J}) where {S<:ANY, F<:ANY, J<:ANY}
    if ds.name == ""
        text = "1-dimensional discrete system"
    else
        text = s.name
    end
    print(io, "1-dimensional discrete dynamical system:\n",
    "state: $(s.state)\n", "e.o.m.: $F\n", "jacobian: $J")
end
@require Juno begin
    function Juno.render(i::Juno.Inline, s::DiscreteDS1D{S, F, J}) where
        {S<:ANY, F<:ANY, J<:ANY}
        t = Juno.render(i, Juno.defaultrepr(s))
        if ds.name == ""
            text = "1-dimensional discrete system"
        else
            text = s.name
        end
        t[:head] = Juno.render(i, Text(text))
        t[:children] = t[:children][1:3] # remove showing field dummystate
        t
    end
end
