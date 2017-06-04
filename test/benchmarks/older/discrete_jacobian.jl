using StaticArrays, BenchmarkTools, ForwardDiff
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 10
show_ind_bench = true
println("Benchmarking Jacobian calculation of 3D Discrete systems for different types of e.o.m.")

u = rand(3); un = rand(3); su = SVector{3}(u); mu = MVector{3}(u)
mun = MVector{3}(un)
J = rand(3,3); MJ = MMatrix{3,3}(J); MI = MMatrix{3,3}(J);
mutable struct system1
  J::StaticArrays.MArray{Tuple{3,3},Float64,2,9}
end
sys = system1(MI)



println("\nVersion 1: StaticArray returned: eom_1(u) -> un::SVector")
@inline @inbounds function eom_1(x)
  x1 = x[1]; x2=x[2]; x3 = x[3]
  SVector{3}(
  3.8*x1*(1-x1)-0.05*(x2+0.35)*(1-2*x3),
  0.1*( (x2+0.35)*(1-2*x3)-1 )*(1-1.9*x1),
  3.78*x3*(1-x3)+0.2*x2
  )
end
println("--In-place jacobian")
@inline jac1!(J, x) = ForwardDiff.jacobian!(J, eom_1, x)
bj1 = @benchmark jac1!($sys.J, $su)
show_ind_bench && display(bj1)
println("--Return-value jacobian")
@inline jac1(x) = (j = ForwardDiff.jacobian(eom_1, x); sys.J .= j)
bjr1 = @benchmark jac1($su)
show_ind_bench && display(bjr1)




println("\nVersion 2: Base Array returned eom_2(x) -> xn::Vector ")
@inline @inbounds function eom_2(x)
  x1 = x[1]; x2=x[2]; x3 = x[3];
  xn = zeros(x)
  xn[1] = 3.8*x1*(1-x1)-0.05*(x2+0.35)*(1-2*x3)
  xn[2] = 0.1*( (x2+0.35)*(1-2*x3)-1 )*(1-1.9*x1)
  xn[3] = 3.78*x3*(1-x3)+0.2*x2
  return xn
end
cfg = ForwardDiff.JacobianConfig(eom_2, u)
println("--In-place jacobian")
@inline jac2!(J, x, cfg) = ForwardDiff.jacobian!(J, eom_2, x, cfg)
bj2 = @benchmark jac2!($sys.J, $su, $cfg)
show_ind_bench && display(bj2)
println("--Return-value jacobian")
@inline jac2(x, cfg) = (j = ForwardDiff.jacobian(eom_2, x, cfg); sys.J .= j)
bjr2 = @benchmark jac2($su, $cfg)
show_ind_bench && display(bjr2)




println("\nVersion 3: MVector returned eom_3(x) -> xn::MVector ")
@inline @inbounds function eom_3(x)
  x1 = x[1]; x2=x[2]; x3 = x[3];
  xn = MVector{3}(x)
  xn[1] = 3.8*x1*(1-x1)-0.05*(x2+0.35)*(1-2*x3)
  xn[2] = 0.1*( (x2+0.35)*(1-2*x3)-1 )*(1-1.9*x1)
  xn[3] = 3.78*x3*(1-x3)+0.2*x2
  return xn
end
cfg = ForwardDiff.JacobianConfig(eom_3, mu)
println("--In-place jacobian")
@inline jac3!(J, x, cfg) = ForwardDiff.jacobian!(J, eom_3, x, cfg)
bj3 = @benchmark jac3!($MJ, $mu, $cfg)
show_ind_bench && display(bj3)
println("--Return-value jacobian")
@inline jac3(x, cfg) = (j = ForwardDiff.jacobian(eom_3, x, cfg); sys.J .= j)
bjr3 = @benchmark jac3($mu, $cfg)
show_ind_bench && display(bjr3)



println("\nVersion 4: In place with 2 arguments: eom_towel!(xn, x) (changes xn)")
@inline @inbounds function eom_4!(xn, x)
  xn[1] = 3.8*x[1]*(1-x[1])-0.05*(x[2]+0.35)*(1-2*x[3])
  xn[2] = 0.1*( (x[2]+0.35)*(1-2*x[3])-1 )*(1-1.9*x[1])
  xn[3] = 3.78*x[3]*(1-x[3])+0.2*x[2]
end
cfg = ForwardDiff.JacobianConfig(eom_4!, un, u)
println("--In-place jacobian")
@inline jac4!(J, y, x, cfg) = ForwardDiff.jacobian!(J, eom_4!, y, x, cfg)
bj4 = @benchmark jac4!($J, $un, $u, $cfg)
show_ind_bench && display(bj4)
println("--Return-value jacobian")
@inline jac4(y, x, cfg) = (j = ForwardDiff.jacobian(eom_4!, y, x, cfg); sys.J .= j)
bjr4 = @benchmark jac4($un, $u, $cfg)
show_ind_bench && display(bjr4)



println("\nVersion 5: In place with 2 arguments & using MVector")
@inline @inbounds function eom_5!(xn::MVector, x::MVector)
  xn[1] = 3.8*x[1]*(1-x[1])-0.05*(x[2]+0.35)*(1-2*x[3])
  xn[2] = 0.1*( (x[2]+0.35)*(1-2*x[3])-1 )*(1-1.9*x[1])
  xn[3] = 3.78*x[3]*(1-x[3])+0.2*x[2]
end
cfg5 = ForwardDiff.JacobianConfig(eom_5!, mun, mu)
println("--In-place jacobian")
@inline jac5!(J, y, x, cfg) = ForwardDiff.jacobian!(J, eom_5!, y, x, cfg)
bj5 = @benchmark jac5!($MJ, $mun, $mu, $cfg5)
show_ind_bench && display(bj5)
println("--Return-value jacobian")
@inline jac5(y, x, cfg) = (j = ForwardDiff.jacobian(eom_5!, y, x, cfg); sys.J .= j)
bjr5 = @benchmark jac5($mun, $mu, $cfg5)
show_ind_bench && display(bjr5)


println("comparison of in-place jacobian call:")
for (i, b) in enumerate([bj2,bj3,bj4,bj5])
  println("v$(i+1) versus v1")
  sleep(0.1)
  display(judge(median(b), median(bj1)))
  sleep(0.1)
end

println("Conclusion about ForwardDiff jacobians:")
println("All methods for Jacobian are slower than the first (which uses SVector)")
println("However method 4 comes extremely close to it!")
println("It is actually even faster than the one with MVector")
