using StaticArrays, ForwardDiff, Requires

export DiscreteDS, DiscreteDS1D, evolve, jacobian, timeseries, setu

abstract type DiscreteDynamicalSystem <: DynamicalSystem end
#######################################################################################
#                                     Constructors                                    #
#######################################################################################
function test_discrete(u0, eom, jac)
  is1D(u0) || throw(ArgumentError("Initial condition must a Vector"))
  D = length(u0)
  su0 = SVector{D}(u0); sun = eom(u0);
  length(sun) == length(s) ||
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
function test_discrete(u0, eom)
  jac = (x) -> ForwardDiff.jacobian(eom, x)
  test_discrete(u0, eom, fd_jac)
end

"""
    DiscreteDS <: DynamicalSystem
Immutable (for efficiency reasons) structure representing a `D`-dimensional
Discrete dynamical system.
# Fields:
* `state::SVector{D}` : Current state-vector of the system, stored in the data format
  of `StaticArray`'s `SVector`.
* `eom::F` (function) : The function that represents the system's equations of motion
  (also called vector field). The function is of the format: `eom(u) -> SVector`
  which means that   given a state-vector `u` it returns an `SVector` containing the
  next state.
* `jacob::J` (function) : A function that calculates the system's jacobian matrix,
  based on the format: `jacob(u) -> SMatrix` which means that given a state-vector
  `u` it returns an `SMatrix` containing the Jacobian at that state.

The function `DynamicalBilliards.test_discrete(u0, eom[, jac])` is provided to help
you ensure that your setup is correct.
# Constructors:
* `DiscreteDS(u0, eom, jac)` : The default constructor.
* `DiscreteDS(u0, eom)` : The Jacobian function is created with *tremendous* efficiency
  using the module `ForwardDiff`. Most of the time, for low dimensional systems, this
  Jacobian is within a few % of speed of a user-defined one.
"""
struct DiscreteDS{D, T<:Real, F, J} <: DiscreteDynamicalSystem
  state::SVector{D,T}
  eom::F
  jacob::J
end
# constsructor without jacobian (uses ForwardDiff)
function DiscreteDS(u0::AbstractVector, eom)
  su0 = SVector{length(u0)}(u0)
  @inline ForwardDiff_jac(x) = ForwardDiff.jacobian(eom, x)
  return DiscreteDS(su0, eom, fd_jac)
end
function DiscreteDS(u0::AbstractVector, eom, jac)
  su0 = SVector{length(u0)}(u0)
  return DiscreteDS(su0, eom, jac)
end

"""
    DiscreteDS1D <: DynamicalSystem
Immutable structure representing an one-dimensional Discrete dynamical system.
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
struct DiscreteDS1D{S<:Real, F, D} <: DiscreteDynamicalSystem
  state::S
  eom::F
  deriv::D
end
function DiscreteDS1D(x0, eom)
  fd_deriv(x) = ForwardDiff.derivative(eom, x)
  DiscreteDS1D(x0, eom, fd_deriv)
end

"""
    setu(u, ds::DynamicalSystem) -> ds_new
Create a new system, identical to `ds` but with state `u`.
"""
setu(u0, ds::DiscreteDS) = DiscreteDS(u0, ds.eom, ds.jacob)
setu(x0, ds::DiscreteDS1D) = DiscreteDS1D(x0, ds.eom, ds.deriv)

"""
    jacobian(ds::DynamicalSystem)
Return the Jacobian matrix of the system at the current state.
"""
jacobian(s::DynamicalSystem) = s.jacob(s.state)

is1D(::DiscreteDS1D) = true
#######################################################################################
#                                 System Evolution                                    #
#######################################################################################
"""
```julia
evolve(st, ds::DiscreteDynamicalSystem, N = 1)
evolve(ds::DiscreteDynamicalSystem, N = 1)
```
Evolve a state `st` (or system `ds`) under the dynamics
of `ds` for `N` steps.
Because both `st` and `ds` are immutable, call as: `st = evolve(st, ds, N)` or
`ds = evolve(ds, N)`.

This function does not store any information about intermediate steps.
Use `timeseries` if you want to produce timeseries of the system.
"""
function evolve(ds::DiscreteDynamicalSystem, N::Int = 1)
  st = deepcopy(ds.state)
  st = evolve(st, ds, N)
  return setu(st, ds)
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
timeseries(ds::DiscreteDynamicalSystem, N::Int)
```
Create an `N×D` matrix that will contain the timeseries of the sytem, after evolving it
for `N` steps. (`D` is the system dimensionality)

Returns a vector for 1-dimensional systems.
"""
function timeseries(s::DiscreteDS, N::Int)
  d = s
  T = eltype(d.state)
  D = length(d.state)
  ts = Array{T}(N, D)
  ts[1,:] .= d.state
  for i in 2:N
    d = evolve(d)
    ts[i, :] .= d.state
  end
  return ts
end

function timeseries(s::DiscreteDS1D, N::Int)
  x = deepcopy(s.state)
  f = s.eom
  ts = Vector{eltype(x)}(N)
  ts[1] = x
  for i in 2:N
    x = f(x)
    ts[i] = x
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
  print(io, "1-dimensional Discrete dynamical system:\n",
  "state: $(s.state)\n", "e.o.m.: $F\n", "jacobian: $J")
end
@require Juno begin
  function Juno.render(i::Juno.Inline, s::DiscreteDS1D{S, F, J}) where {S<:ANY, F<:ANY, J<:ANY}
    t = Juno.render(i, Juno.defaultrepr(s))
    t[:head] = Juno.render(i, Text("1-dimensional Discrete dynamical system"))
    t
  end
end
