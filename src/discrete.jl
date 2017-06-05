using StaticArrays, ForwardDiff, Requires

export DiscreteDS, evolve, jacobian, timeseries
#######################################################################################
#                                     Constructors                                    #
#######################################################################################
"""
    DiscreteDS <: DynamicalSystem
Immutable structure representing a D-dimensional Discrete dynamical system.
# Fields:
* `state::SVector{D}` : Current state-vector of the system, stored in the data format
  of StaticArray's `SVector`.
* `eom::F` (function) : The function that represents the system's equations of motion (also called
  vector field). The function is of the format: `eom(u) -> SVector` which means that
  given a state-vector `u` it returns an `SVector` containing the next state.
* `jacob::J` (function) : A function that calculates the system's jacobian matrix, based on the
  format: `jacob(u) -> SMatrx` which means that given a state-vector `u` it returns
  an SMatrix containing the Jacobian at that state.
# Constructors:
* `DiscreteDS(u0, eom, jac)` : The default constructor. **Ensures that the functions
  given are at the form described here**, or results in error otherwise.
* `DiscreteDS(u0, eom)` : The Jacobian function is created with tremendous efficiency
  using the module `ForwardDiff`. Most of the time, for low dimensional systems, this
  Jacobian is within a few % of speed with a user-defined one.
"""
struct DiscreteDS{D, S<:Number, F, J} <: DynamicalSystem
  state::SVector{D,S}
  eom::F
  jacob::J
end
# constsructor without jacobian (uses ForwardDiff)
function DiscreteDS(u0, eom)
  @assert length(size(u0)) == 1 "Initial condition must an one-dimensional vector"
  D = length(u0)
  su0 = SVector{D}(u0); sun = eom(u0);
  ssu0 = eom(sun) #test that eom works with SVector as well

  if !issubtype((typeof(sun)), SVector)
    error("E.o.m. should create an SVector (from StaticArrays)")
  end

  @inline jac(x) = ForwardDiff.jacobian(eom, x)
  J = jac(su0)
  if !issubtype((typeof(J)), SMatrix)
    error("Fatal error in DiscreteDS: Jacobian is not SMatrix")
  end
  return DiscreteDS(su0, eom, jac)
end

function DiscreteDS(u0, eom, jac)
  @assert length(size(u0)) == 1 "Initial condition must an one-dimensional vector"
  D = length(u0)
  su0 = SVector{D}(u0); sun = eom(u0);
  ssu0 = eom(sun) #test that eom works with SVector as well

  if !issubtype((typeof(sun)), SVector)
    error("E.o.m. should create an SVector (from StaticArrays)!")
  end

  J1 = jac(u0); J2 = jac(su0)
  if !issubtype((typeof(J1)), SMatrix) || !issubtype((typeof(J2)), SMatrix)
    error("jacobian function should create an SMatrix (from StaticArrays)!")
  end
  return DiscreteDS(su0, eom, jac)
end

#######################################################################################
#                                 System Evolution                                    #
#######################################################################################
"""
```julia
evolve(s::DiscreteDS, N::Int = 1)
```
Evolve a discrete system for `N` steps. Because `DiscreteDS` is immutable, `evolve`
has to be called as `s = evolve(s, N)`.

This function does not store any information about intermediate steps. Use `timeseries`
if you want to keep intermediate information.
"""
@inline function evolve(s::DiscreteDS, N::Int)
  d = s
  for i in 1:N
    d = DiscreteDS(d.eom(d.state), d.eom, d.jacob)
  end
  return d
end
@inline function evolve(s::DiscreteDS)
  DiscreteDS(s.eom(s.state), s.eom, s.jacob)
end

"""
```julia
timeseries(s::DiscreteDS, N::Int)
```
Create an `N×D` matrix that will contain the timeseries of the sytem, after evolving it
for `N` steps. (`D` is the system dimensionality)
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

jacobian(s::DiscreteDS) = s.jacob(s.state)

#######################################################################################
#                                 Pretty-Printing                                     #
#######################################################################################
import Base.show

function Base.show(io::IO, s::DiscreteDS{N, S, F, J}) where {N<:ANY, S<:ANY, F<:ANY, J<:ANY}
  print(io, "$N-dimensional Discrete dynamical system:\n",
  "state: $(s.state)\n", "e.o.m.: $F\n", "jacobian: $J")
end

 @require Juno begin
  import Juno.render
  function Juno.render(i::Juno.Inline, s::DiscreteDS{N, S, F, J}) where {N<:ANY, S<:ANY, F<:ANY, J<:ANY}
    t = Juno.render(i, Juno.defaultrepr(s))
    t[:head] = Juno.render(i, Text("$N-dimensional Discrete dynamical system"))
    t
  end
end
