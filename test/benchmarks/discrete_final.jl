using StaticArrays, BenchmarkTools, ForwardDiff
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 10
show_ind_bench = true
show_jac = false
N = 1000000
println("Benchmarking evolution of Full Discrete systems (with Jac and eom and state)")
println("The eom will always return an SVector, since this is the fastest for")
println("the Jacobian. Evolution for total of $N steps")
u = rand(3); un = rand(3); J = rand(3,3)
su = SVector{3}(u); sun = SVector{3}(un);
mu = MVector{3}(u); mun = MVector{3}(un)
@inline function eom_towel(x)
  x1, x2, x3 = x[1], x[2], x[3]
  SVector(3.8*x1*(1-x1)-0.05*(x2+0.35)*(1-2*x3),
  0.1*((x2+0.35)*(1-2*x3)-1 )*(1-1.9*x1),
  3.78*x3*(1-x3)+0.2*x2)
end

println("V1: immutable struct with state SVector (new system at every evolve)")
struct dds1{S<:Real, F, J} #<: DynamicalSystem
  state::SVector{3,S}
  eom::F
  jacob!::J
end
function dds1(u0, eom)
  @inline jac!(J, x) = ForwardDiff.jacobian!(J, eom, x)
  return dds1(u0, eom, jac!)
end
@inline function evolve1(s::dds1, N::Int)
  for i in 1:N
    s = (dds1(s.eom(s.state), s.eom))
  end
  return s
end
s1 = dds1(su, eom_towel)
s1 = evolve1(s1, 1000)
b1 = @benchmark (s = evolve1($s1, $N))
j1 = @benchmark ()
display(b1)

println("V2: immutable struct with state Vector (changed in-place)")
struct dds2{S<:Real, F, J} #<: DynamicalSystem
  state::Vector{S}
  eom::F
  jacob!::J
end
function dds2(u0, eom)
  @inline jac!(J, x) = ForwardDiff.jacobian!(J, eom, x)
  return dds2(u0, eom, jac!)
end
@inline function evolve2!(s::dds2, N::Int)
  eom = s.eom
  for i in 1:N
    s.state .= eom(s.state)
  end
end
s2 = dds2(u, eom_towel)
evolve2!(s2, 1000)
b2 = @benchmark evolve2!($s2, $N)
display(b2)


println("V3: immutable struct with state MVector (changed in-place)")
struct dds3{S<:Real, F, J} #<: DynamicalSystem
  state::MVector{3, S}
  eom::F
  jacob!::J
end
function dds3(u0, eom)
  @inline jac!(J, x) = ForwardDiff.jacobian!(J, eom, x)
  return dds3(u0, eom, jac!)
end
@inline function evolve3!(s::dds3, N::Int)
  eom = s.eom
  for i in 1:N
    s.state .= eom(s.state)
  end
end
s3 = dds3(mu, eom_towel)
evolve3!(s3, 1000)
b3 = @benchmark evolve3!($s3, $N)
display(b3)

beches = [minimum(b) for b in [b1,b2,b3]]
i = argmin(beches)
println("Conclusions:")
println("Minimum time of evolve-call for N=$(N): v$i")
println("Judged vs v$i:")
for (j,b) in enumerate(beches)
  j == i && continue
  println("v$i vs v$j")
  sleep(0.01)
  display(judge(beches[i],beches[j]))
  sleep(0.01)
end




#### This is broken in Base (julia;s fault.) The immutable should be faster normally.
using BenchmarkTools
println("One dimensional benchmarks")
struct ds1d1{S<:Real, F}
  state::S
  eom::F
end
function evolve1d1(s::ds1d1, N::Int)
  x = deepcopy(s.state)
  f = s.eom
  for i in 1:N
    x = f(x)
  end
  return ds1d1(x, s.eom)
end

mutable struct ds2d2{S<:Real, F}
  state::S
  eom::F
end
function evolve2d2(s::ds2d2, N::Int)
  x = deepcopy(s.state)
  f = s.eom
  for i in 1:N
    x = f(x)
  end
  s.state = x
end

function bench(N)
  x0 = rand()

  @inline eom_logistic(x) = 4*x*(1-x)
  d1 = ds1d1(x0, eom_logistic)
  println("V1: Immutable struct")
  b1 = @benchmark evolve1d1($d1, $N)
  display(b1)

  println("V2: mutable struct")
  d2 = ds2d2(x0, eom_logistic)
  b2 = @benchmark evolve2d2($d2, $N)
  display(b2)
end

println("For N = 1000")
bench(1000)
println("For N = 1000000")
bench(1000000)

## Profiling:
Profile.clear()
N = 1000000000
x0 = rand()
@inline eom_logistic(x) = 4*x*(1-x)
d1 = ds1d1(x0, eom_logistic)
evolve1d1(d1, N)
@profile evolve1d1(d1, N)
Profile.print()

Profile.clear()
N = 1000000000
x0 = rand()
@inline eom_logistic(x) = 4*x*(1-x)
d2 = ds2d2(x0, eom_logistic)
evolve2d2(d2, N)
@profile evolve2d2(d2, N)
Profile.print()

## Profiling Juno
# Profile.clear()
# N = 10000000
# x0 = rand()
# eom_logistic(x) = 4*x*(1-x)
# d1 = ds1d1(x0, eom_logistic)
# @profile evolve1d1(d1, N)
# Atom.Profiler.tree()
#
# N = 10000000
# x0 = rand()
# eom_logistic(x) = 4*x*(1-x)
# d1 = ds1d1(x0, eom_logistic)
# @profile evolve1d1(d1, N)
# Juno.profiletree()
