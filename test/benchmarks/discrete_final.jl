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
i = indmin(beches)
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




####
@inline fj(u) = ForwardDiff.jacobian(eom_towel, u)
@inline fj!(J, u) = ForwardDiff.jacobian!(J, eom_towel, u)

j1 = @benchmark fj($su)
display(j1)
j2 = @benchmark fj!($J, $su)
display(j2)
