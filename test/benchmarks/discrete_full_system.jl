using StaticArrays, BenchmarkTools, ForwardDiff
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 10
show_ind_bench = true
println("Benchmarking evolution of Full Discrete systems (with Jac)")
println("The jacobian method will always be the in-place one")
u = rand(3); un = rand(3); J = rand(3,3)
su = SVector{3}(u); sun = SVector{3}(un);
mu = MVector{3}(u); mun = MVector{3}(un)


println("\nVersion 1: e.o.m. return static array: eom_1(u) -> un::SVector")
println("and system.u is also an SVector")
mutable struct system_1f{T<:Real, F<:Function}
  u::SVector{3, T}
  eom::F
end
@inline @inbounds function eom_1(x)
  x1 = x[1]; x2=x[2]; x3 = x[3]; T = eltype(x)
  SVector{3, T}(
  3.8*x1*(1-x1)-0.05*(x2+0.35)*(1-2*x3),
  0.1*( (x2+0.35)*(1-2*x3)-1 )*(1-1.9*x1),
  3.78*x3*(1-x3)+0.2*x2
  )
end
@inline jac1!(J, x) = ForwardDiff.jacobian!(J, eom_1, x)
s1 = system_1f(SVector{3}(u), eom_1)

evolve1(s::system_1f) = (un = s.eom(s.u); s.u = un)
show_ind_bench && println("system evolution:")
b1 = @benchmark evolve1($s1)
show_ind_bench && display(b1)
show_ind_bench && println("in-place jacobian")
bj1 = @benchmark (jac1!($J, $(s1.u)))
show_ind_bench && display(bj1)
show_ind_bench && sleep(0.1)



println("\nVersion 2: e.o.m. are in-place with Base Arrays")
println("and system.u is also a Base Array")
mutable struct system_2f{T<:Real, F<:Function}
  u::Vector{T}
  eom::F
end
@inline @inbounds function eom_2(xn, x)
  x1 = x[1]; x2=x[2]; x3 = x[3]
  xn[1] = 3.8*x1*(1-x1)-0.05*(x2+0.35)*(1-2*x3)
  xn[2] = 0.1*( (x2+0.35)*(1-2*x3)-1 )*(1-1.9*x1)
  xn[3] = 3.78*x3*(1-x3)+0.2*x2
end
s2 = system_2f(u, eom_2)

evolve2(s::system_2f) = (un = copy(s.u); s.eom(s.u, un))
show_ind_bench && println("system evolution:")
b2 = @benchmark evolve2($s2)
show_ind_bench && display(b2)

println("in-place jacobian")
cfg = ForwardDiff.JacobianConfig(eom_2, un, u)
@inline jac2!(J, y, x, cfg) = ForwardDiff.jacobian!(J, eom_2, y, x, cfg)
bj2 = @benchmark jac2!($J, $un, $u, $cfg)
show_ind_bench && display(bj2)
show_ind_bench && sleep(0.1)


println("\nVersion 3: system is an immutable struct with SVector")
struct system_3f{T<:Real, F<:Function}
    u::SVector{3,T}
    eom::F
end
@inline function eom_3(x)
    x1, x2, x3 = x[1], x[2], x[3]
    return SVector(3.8*x1*(1-x1)-0.05*(x2+0.35)*(1-2*x3),
                   0.1*((x2+0.35)*(1-2*x3)-1 )*(1-1.9*x1),
                   3.78*x3*(1-x3)+0.2*x2)
end
@inline evolve3(s::system_3f) = system_3f(s.eom(s.u), s.eom)
s3 = system_3f(su, eom_3)
b3 = @benchmark evolve3(s) setup=(s = system_3f(SVector{3}(rand(3)), eom_3))
show_ind_bench && println("system evolution:")
show_ind_bench && display(b3)
show_ind_bench && println("in-place jacobian (same as Version 1)")
@inline jac3!(J, x) = ForwardDiff.jacobian!(J, eom_3, x)
bj3 = @benchmark (jac3!($J, $(s3.u)))
show_ind_bench && display(bj3)
show_ind_bench && sleep(0.1)

println("\nVersion 4: e.o.m. are in-place with Mutable StaticArrays")
println("and system.u is also a MutableStaticArray")
mutable struct system_4f{T<:Real, F<:Function}
  u::MVector{3, T}
  eom::F
end
@inline @inbounds function eom_4(xn, x)
  x1 = x[1]; x2=x[2]; x3 = x[3]
  xn[1] = 3.8*x1*(1-x1)-0.05*(x2+0.35)*(1-2*x3)
  xn[2] = 0.1*( (x2+0.35)*(1-2*x3)-1 )*(1-1.9*x1)
  xn[3] = 3.78*x3*(1-x3)+0.2*x2
end
s4 = system_4f(mu, eom_4)

evolve4(s::system_4f) = (un = copy(s.u); s.eom(s.u, un))
show_ind_bench && println("system evolution:")
b4 = @benchmark evolve4($s4)
show_ind_bench && display(b4)

println("in-place jacobian")
cfg = ForwardDiff.JacobianConfig(eom_4, mun, mu)
@inline jac4!(J, y, x, cfg) = ForwardDiff.jacobian!(J, eom_4, y, x, cfg)
bj4 = @benchmark jac4!($J, $un, $u, $cfg)
show_ind_bench && display(bj4)
show_ind_bench && sleep(0.1)



println("-----Comparison of evolve call:-----")
for (i, b) in enumerate([b2,b3,b4])
  println("v$(i+1) versus v1")
  sleep(0.1)
  display(judge(median(b), median(b1)))
  sleep(0.1)
end
println("-----Comparison of in-place jacobian call:-----")
for (i, b) in enumerate([bj2,bj3,bj4])
  println("v$(i+1) versus v1")
  sleep(0.1)
  display(judge(median(b), median(bj1)))
  sleep(0.1)
end

println("Conclusions:")
evcall = [median(b) for b in [b1,b2,b3,b4]]
jaccall = [median(b) for b in [bj1,bj2,bj3,bj4]]
println("Minimum time of evolve call for method v$(indmin(evcall))")
println("Minimum time of jacob call for method v$(indmin(jaccall))")
