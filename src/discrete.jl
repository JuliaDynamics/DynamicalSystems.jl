using StaticArrays, ForwardDiff
is1D(u) = length(size(u)) == 1
#######################################################
#                     Constructors                    #
#######################################################
struct DiscreteDS{N, S<:Real, F, J} #<: DynamicalSystem
  state::SVector{N,S}
  eom::F
  jacob::J
end
# constsructor without jacobian (uses ForwardDiff)
function DiscreteDS(u0, eom)
  @assert is1D(u) "Initial condition must an one-dimensional vector"
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
  @assert is1D(u) "Initial condition must an one-dimensional vector"
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

#######################################################
#                      Evolutions                     #
#######################################################
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
