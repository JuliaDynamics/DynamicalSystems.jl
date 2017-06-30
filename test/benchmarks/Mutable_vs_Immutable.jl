using BenchmarkTools, StaticArrays
abstract type System end
@inline function eom_towel(x)
  x1, x2, x3 = x[1], x[2], x[3]
  SVector( 3.8*x1*(1-x1) - 0.05*(x2+0.35)*(1-2*x3),
  0.1*( (x2+0.35)*(1-2*x3) - 1 )*(1 - 1.9*x1),
  3.78*x3*(1-x3)+0.2*x2 )
end
@inline function jacob_towel(x)
  @SMatrix [3.8*(1 - 2x[1]) -0.05*(1-2x[3]) 0.1*(x[2] + 0.35);
  -0.19((x[2] + 0.35)*(1-2x[3]) - 1)  0.1*(1-2x[3])*(1-1.9x[1])  -0.2*(x[2] + 0.35)*(1-1.9x[1]);
  0.0  0.2  3.78(1-2x[3]) ]
end
u0 = @SVector [0.085, -0.121, 0.075]
struct DiscreteDS{D, T<:Real, F, J} <: System
  state::SVector{D,T}
  eom::F
  jacob::J
end

mutable struct DiscreteDSmut{D, T<:Real, F, J} <: System
  state::SVector{D,T}
  eom::F
  jacob::J
end

function evolve(state, ds::System, N::Int = 1)
  f = ds.eom
  for i in 1:N
    state = f(state)
  end
  return state
end

function evolve!(ds::DiscreteDS, N::Int = 1)
  st = ds.state
  st = evolve(st, ds, N)
  return DiscreteDS(st, ds.eom, ds.jacob)
end

function evolve!(ds::DiscreteDSmut, N::Int = 1)
  st = ds.state
  st = evolve(st, ds, N)
  ds.state = st
  return ds
end

ds = DiscreteDS(u0, eom_towel, jacob_towel)
dsmut = DiscreteDSmut(u0, eom_towel, jacob_towel)

ds = evolve!(ds, 100)
evolve!(dsmut, 100)

a = @btime evolve!($ds, 10000)
b = @btime evolve!($dsmut, 10000)
