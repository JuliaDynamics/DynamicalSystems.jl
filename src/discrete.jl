using StaticArrays, ForwardDiff, Requires

export DiscreteDS, DiscreteDS1D, evolve, evolve!, timeseries, dimension, jacobian

abstract type DiscreteDynamicalSystem <: DynamicalSystem end
#######################################################################################
#                                     Constructors                                    #
#######################################################################################
function test_functions(u0, eom, jac)
  length(size(u0)) == 1 || throw(ArgumentError("Initial condition must an AbstractVector"))
  D = length(u0)
  su0 = SVector{D}(u0); sun = eom(u0);
  length(sun) == length(u0) ||
  throw(DimensionMismatch("E.o.m. does not give same sized vector as initial condition"))
  if !issubtype((typeof(sun)), SVector)
    throw(ArgumentError("E.o.m. should create an SVector (from StaticArrays)"))
  end
  J1 = jac(u0); J2 = jac(SVector{length(u0)}(u0))
  if !issubtype((typeof(J1)), SMatrix) || !issubtype((typeof(J2)), SMatrix)
    throw(ArgumentError("Jacobian function should create an SMatrix (from StaticArrays)!"))
  end
  return true
end
function test_functions(u0, eom)
  jac = (x) -> ForwardDiff.jacobian(eom, x)
  test_discrete(u0, eom, fd_jac)
end

"""
    DiscreteDS <: DynamicalSystem
`D`-dimensional discrete dynamical system (used for `D ≤ 10`).
# Fields:
* `state::SVector{D}` : Current state-vector of the system, stored in the data format
  of `StaticArray`'s `SVector`.
* `eom::F` (function) : The function that represents the system's equations of motion
  (also called vector field). The function is of the format: `eom(u) -> SVector`
  which means that given a state-vector `u` it returns an `SVector` containing the
  next state.
* `jacob::J` (function) : A function that calculates the system's jacobian matrix,
  based on the format: `jacob(u) -> SMatrix` which means that given a state-vector
  `u` it returns an `SMatrix` containing the Jacobian at that state.

# Constructors:
* `DiscreteDS(u0, eom, jac)` : The default constructor.
* `DiscreteDS(u0, eom)` : The Jacobian function is created with *tremendous* efficiency
  using the module `ForwardDiff`. Most of the time, for low dimensional systems, this
  Jacobian is within a few % of speed of a user-defined one.
"""
mutable struct DiscreteDS{D, T<:Real, F, J} <: DiscreteDynamicalSystem
  state::SVector{D,T}
  eom::F
  jacob::J
end
# constsructor without jacobian (uses ForwardDiff)
function DiscreteDS(u0::AbstractVector, eom)
  su0 = SVector{length(u0)}(u0)
  @inline ForwardDiff_jac(x) = ForwardDiff.jacobian(eom, x)
  # test_functions(su0, eom, ForwardDiff_jac)
  return DiscreteDS(su0, eom, ForwardDiff_jac)
end
function DiscreteDS(u0::AbstractVector, eom, jac)
  su0 = SVector{length(u0)}(u0)
  # test_functions(su0, eom, jac)
  return DiscreteDS(su0, eom, jac)
end

"""
    DiscreteDS1D <: DynamicalSystem
One-dimensional discrete dynamical system.
# Fields:
* `state::Real` : Current state of the system.
* `eom::F` (function) : The function that represents the system's equations of motion:
  `eom(x) -> Real`.
* `deriv::D` (function) : A function that calculates the system's derivative given
  a state: `deriv(x) -> Real`.
# Constructors:
* `DiscreteDS1D(x0, eom, deriv)` : The default constructor with user-provided
  derivative function (most efficient)
* `DiscreteDS1d(x0, eom)` : The derivative function is created
  automatically using the module `ForwardDiff`.
"""
mutable struct DiscreteDS1D{S<:Real, F, D} <: DiscreteDynamicalSystem
  state::S
  eom::F
  deriv::D
end
function DiscreteDS1D(x0, eom)
  fd_deriv(x) = ForwardDiff.derivative(eom, x)
  DiscreteDS1D(x0, eom, fd_deriv)
end



dimension(::DiscreteDS{D, T, F, J})  where {D<:ANY, T<:ANY, F<:ANY, J<:ANY} = D
dimension(::DiscreteDS1D) = 1
is1D(::DiscreteDS1D) = true
jacobian(ds::DynamicalSystem) = ds.jacob(ds.state)
#######################################################################################
#                                 System Evolution                                    #
#######################################################################################
"""
```julia
evolve([state, ] ds::DynamicalSystem, T=1; diff_eq_kwargs = Dict()) -> new_state
```
Evolve a `state` (or the system's state) under the dynamics
of `ds` for total "time" `T`. For discrete systems `T` corresponds to steps and
thus it must be integer. Returns the final state after evolution.

The **keyword** argument `diff_eq_kwargs` is a `Dict{Symbol, Any}` and is only
applicable for continuous systems. It contains keyword arguments passed into the
`solve` of the `DifferentialEquations` package, like for
example `:abstol => 1e-9`. If you want to specify the solving algorithm,
do so by using `:solver` as one of your keywords, like `:solver => DP5()`.

This function *does not store* any information about intermediate steps.
Use `timeseries` if you want to produce timeseries of the system.
"""
function evolve(ds::DiscreteDynamicalSystem, N::Int = 1)
  st = ds.state
  st = evolve(st, ds, N)
end
function evolve(state, ds::DiscreteDynamicalSystem, N::Int = 1)
  f = ds.eom
  for i in 1:N
    state = f(state)
  end
  return state
end

"""
```julia
evolve!(ds::DynamicalSystem, T=1; diff_eq_kwargs = Dict()) -> ds
```
Evolve (in-place) a dynamical system for total "time" `T`.
For discrete systems `T` corresponds to steps and
thus it must be integer.

The last **keyword** argument `diff_eq_kwargs` is a `Dict{Symbol, Any}` and is only
applicable for continuous systems. It contains keyword arguments passed into the
`solve` of the `DifferentialEquations` package, like for
example `:abstol => 1e-9`. If you want to specify the solving algorithm,
do so by using `:solver` as one of your keywords, like `:solver => DP5()`.

This function *does not store* any information about intermediate steps.
Use `timeseries` if you want to produce timeseries of the system.
"""
function evolve!(ds::DiscreteDynamicalSystem, N::Int = 1)
  st = ds.state
  ds.state = evolve(st, ds, N)
  return ds
end


"""
```julia
timeseries(ds::DiscreteDS, N::Int; mutate = true)
```
Create a `N×D` matrix that will contain the timeseries of the sytem, after evolving it
for `N` steps (`D` is the system dimensionality).

*Each column corresponds to one dynamic variable.*

The keyword argument `mutate` controls whether the system's state is updated to
the final state of the timeseries.
"""
function timeseries(ds::DiscreteDS, N::Int; mutate = true)
  st = ds.state
  T = eltype(st)
  D = length(st)
  ts = Array{T}(N, D)
  f = ds.eom
  ts[1,:] .= ds.state
  for i in 2:N
    st = f(st)
    ts[i, :] .= st
  end
  if mutate
    ds.state = ts[end, :]
  end
  return ts
end

function timeseries(ds::DiscreteDS1D, N::Int; mutate = true)
  x = deepcopy(ds.state)
  f = ds.eom
  ts = Vector{eltype(x)}(N)
  ts[1] = x
  for i in 2:N
    x = f(x)
    ts[i] = x
  end
  if mutate
    ds.state = x
  end
  return ts
end

#######################################################################################
#                                 Pretty-Printing                                     #
#######################################################################################
import Base.show
function Base.show(io::IO, s::DiscreteDS{N, S, F, J}) where {N<:ANY, S<:ANY, F<:ANY, J<:ANY}
  print(io, "$N-dimensional discrete dynamical system:\n",
  "state: $(s.state)\n", "e.o.m.: $F\n", "jacobian: $J")
end

@require Juno begin
  function Juno.render(i::Juno.Inline, s::DiscreteDS{N, S, F, J}) where {N<:ANY, S<:ANY, F<:ANY, J<:ANY}
    t = Juno.render(i, Juno.defaultrepr(s))
    t[:head] = Juno.render(i, Text("$N-dimensional discrete dynamical system"))
    t
  end
end

# 1-D
function Base.show(io::IO, s::DiscreteDS1D{S, F, J}) where {S<:ANY, F<:ANY, J<:ANY}
  print(io, "1-dimensional discrete dynamical system:\n",
  "state: $(s.state)\n", "e.o.m.: $F\n", "jacobian: $J")
end
@require Juno begin
  function Juno.render(i::Juno.Inline, s::DiscreteDS1D{S, F, J}) where {S<:ANY, F<:ANY, J<:ANY}
    t = Juno.render(i, Juno.defaultrepr(s))
    t[:head] = Juno.render(i, Text("1-dimensional discrete dynamical system"))
    t
  end
end
