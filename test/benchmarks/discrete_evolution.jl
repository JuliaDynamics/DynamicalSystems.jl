using StaticArrays, BenchmarkTools
println("Benchmarking evolution of Discrete systems for different types of e.o.m.")
u = rand(3); un = rand(3)


println("\nVersion 1: StaticArray returned: eom_1(u) -> un::SVector")

@inline @inbounds function eom_1(x)::SVector{3, <:Real}
  x1 = x[1]; x2=x[2]; x3 = x[3]
  SVector{3}(
  3.8*x1*(1-x1)-0.05*(x2+0.35)*(1-2*x3),
  0.1*( (x2+0.35)*(1-2*x3)-1 )*(1-1.9*x1),
  3.78*x3*(1-x3)+0.2*x2
  )
end
println("eom call:")
e1 = @benchmark eom_1($u)
display(e1)

mutable struct system_1
  u::SVector{3, <:Real}
end
evolve1(s::system_1) = (un = eom_1(s.u); s.u = un)
s1 = system_1(SVector{3}(u))
println("system evolution:")
b1 = @benchmark evolve1($s1)
display(b1)
sleep(0.1)
println("-------------------------------------------")

println("\nVersion 2: Base Array returned eom_2(x) -> xn::Vector ")
@inline @inbounds function eom_2(x)
  x1 = x[1]; x2=x[2]; x3 = x[3];
  xn = zeros(x)
  xn[1] = 3.8*x1*(1-x1)-0.05*(x2+0.35)*(1-2*x3)
  xn[2] = 0.1*( (x2+0.35)*(1-2*x3)-1 )*(1-1.9*x1)
  xn[3] = 3.78*x3*(1-x3)+0.2*x2
  return xn
end
println("eom call:")
e2 = @benchmark eom_2($u)
display(e2)

mutable struct system_2
  u::Vector{<:Real}
end
evolve2(s::system_2) = (un = eom_2(s.u); s.u .= un)
s2 = system_2(u)
println("system evolution:")
b2 = @benchmark evolve2($s2)
display(b2)
sleep(0.1)
println("-------------------------------------------")


println("\nVersion 3: MVector returned eom_3(x) -> xn::MVector ")
@inline @inbounds function eom_3(x)::MVector{3, <:Real}
  x1 = x[1]; x2=x[2]; x3 = x[3];
  xn = MVector{3}(x)
  xn[1] = 3.8*x1*(1-x1)-0.05*(x2+0.35)*(1-2*x3)
  xn[2] = 0.1*( (x2+0.35)*(1-2*x3)-1 )*(1-1.9*x1)
  xn[3] = 3.78*x3*(1-x3)+0.2*x2
  return xn
end
println("call:")
e3 = @benchmark eom_3($u)
display(e3)

mutable struct system_3
  u::MVector{3, <:Real}
end
evolve3(s::system_3) = (un = eom_3(s.u); s.u .= un)
s3 = system_3(MVector{3}(u))
println("system evolution:")
b3 = @benchmark evolve3($s3)
display(b3)
sleep(0.1)
println("-------------------------------------------")

println("\nVersion 4: In place with 2 arguments: eom_towel!(xn, x) (changes xn)")
@inline @inbounds function eom_4!(xn, x)
  xn[1] = 3.8*x[1]*(1-x[1])-0.05*(x[2]+0.35)*(1-2*x[3])
  xn[2] = 0.1*( (x[2]+0.35)*(1-2*x[3])-1 )*(1-1.9*x[1])
  xn[3] = 3.78*x[3]*(1-x[3])+0.2*x[2]
end
println("eom call:")
e4 = @benchmark eom_4!($un, $u)
display(e4)

mutable struct system_4
  u::Vector{<:Real}
end
evolve4(s::system_4) = (u0 = copy(s.u); eom_4!(s.u, u0))
s4 = system_4((u))
println("system evolution:")
b4 = @benchmark evolve4($s4)
display(b4)
sleep(0.1)
println("-------------------------------------------")


println("\nVersion 5: In place with 2 arguments & using MVector")
@inline @inbounds function eom_5!(xn::MVector, x)
  xn[1] = 3.8*x[1]*(1-x[1])-0.05*(x[2]+0.35)*(1-2*x[3])
  xn[2] = 0.1*( (x[2]+0.35)*(1-2*x[3])-1 )*(1-1.9*x[1])
  xn[3] = 3.78*x[3]*(1-x[3])+0.2*x[2]
end
println("eom call:")
mun = MVector{3}(un)
e5 = @benchmark eom_5!($mun, $u)
display(e5)

mutable struct system_5
  u::MVector{3, <:Real}
end
evolve5(s::system_5) = (u0 = copy(s.u); eom_5!(s.u, u0))
s5 = system_5(MVector{3}(u))
println("system evolution:")
b5 = @benchmark evolve5($s5)
display(b5)
sleep(0.1)
println("-------------------------------------------")

println("Judging version 1 with the rest... (using median)")
println("comparison of e.o.m. call:")
for (i, b) in enumerate([e2,e3,e4,e5])
  println("v$(i+1) versus v1")
  sleep(0.1)
  display(judge(median(b), median(e1)))
  sleep(0.1)
end

println("comparison of evolve call:")
for (i, b) in enumerate([b2,b3,b4,b5])
  println("v$(i+1) versus v1")
  sleep(0.1)
  display(judge(median(b), median(b1)))
  sleep(0.1)
end
println("Conclusion: Versions 3,4 and 5 are by far the best for `evolve!`.")
println("Why isn't evolve performing well for version 1? Because the system field `u`")
println("has to be updated. In all other versions you can do s.u .= u, but not for v1!")
