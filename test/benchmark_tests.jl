using BenchmarkTools, ForwardDiff

#######################################################################################
println("Benchmark: having a Dummy field/array or giving zeros(u)")
println("to ForwardDiff.jacobian!")
d = 20
eom!(xn, x) = (xn .= 2x)
function f1(d)
  x = rand(d); y = copy(x)
  benchdum = rand(d,d)
  jacob1! = (J, x) -> ForwardDiff.jacobian!(J, eom!, y, x)
  for i in 1:1000
    jacob1!(benchdum, x)
  end
end
function f2(d)
  x = rand(d); y = copy(x)
  benchdum = rand(d,d)
  jacob2! = (J, x) -> ForwardDiff.jacobian!(J, eom!, zeros(x), x)
  for i in 1:1000
    jacob2!(benchdum, x)
  end
end

println("For d = $d")
b1 = @benchmark f1($d)
b2 = @benchmark f2($d)
println("Mean b1 = $(mean(b1))")
println("Mean b2 = $(mean(b2))")
#=
Result: Creating and giving a zeros(u) in the Jacobian of ForwardDiff
makes no difference with giving something pre-allocated
=#

#######################################################################################
println("\n\nTest difference of discrete evolve using two vectors or only one")
@inline @inbounds function eom_towel!(x)
  x1 = x[1]; x2=x[2]; x3 = x[3]
  x[1] = 3.8*x1*(1-x1)-0.05*(x2+0.35)*(1-2*x3)
  x[2] = 0.1*( (x2+0.35)*(1-2*x3)-1 )*(1-1.9*x1)
  x[3] = 3.78*x3*(1-x3)+0.2*x2
  return x
end
@inline @inbounds function eom_towel!(xn, x)
  xn[1] = 3.8*x[1]*(1-x[1])-0.05*(x[2]+0.35)*(1-2*x[3])
  xn[2] = 0.1*( (x[2]+0.35)*(1-2*x[3])-1 )*(1-1.9*x[1])
  xn[3] = 3.78*x[3]*(1-x[3])+0.2*x[2]
end

function inte1(x0, n)
  for i in 1:n
    eom_towel!(x0)
  end
  return x0
end
function inte2(x0, n)
  x1 = copy(x0)
  for i in 1:n
    eom_towel!(x1, x0)
    copy!(x0, x1) #This line is necessary because for all real calculations the
    # "old" state has to be updated somehow.
  end
  return x1
end

for i in [100, 10000, 1000000]
  x0 = rand(3)
  a = @benchmark inte1($x0, $i)
  b = @benchmark inte2($x0, $i)
  println("\nn = $i")
  println("inte1 mean: $(mean(a))")
  println("inte2 mean: $(mean(b))")
end
#=
Result: The version with only one vector is more than twice as fast
=#

#######################################################################################
println("\n\nTest difference of jacobian with having the 2 argument form, versus")
println("creating an anonymous function that does one-argument form (by creating vec)")
@inbounds function eom_towel!(x)
  x1 = x[1]; x2=x[2]; x3 = x[3]
  x[1] = 3.8*x1*(1-x1)-0.05*(x2+0.35)*(1-2*x3)
  x[2] = 0.1*( (x2+0.35)*(1-2*x3)-1 )*(1-1.9*x1)
  x[3] = 3.78*x3*(1-x3)+0.2*x2
  return x
end
@inbounds function eom_towel!(xn, x)
  xn[1] = 3.8*x[1]*(1-x[1])-0.05*(x[2]+0.35)*(1-2*x[3])
  xn[2] = 0.1*( (x[2]+0.35)*(1-2*x[3])-1 )*(1-1.9*x[1])
  xn[3] = 3.78*x[3]*(1-x[3])+0.2*x[2]
end

function jac1(x0, n)
  fakef = (u) -> (un = copy(u); eom_towel!(un); un)
  jacob! = (J, u) -> ForwardDiff.jacobian!(J, fakef, u)
  XX = rand(3,3)
  for i in 1:n
    jacob!(XX, x0)
  end
end
function jac2(x0, n)
  x = rand(3)
  jacob! = (J, u) -> ForwardDiff.jacobian!(J, eom_towel!, x, u)
  XX = rand(3,3)
  for i in 1:n
    jacob!(XX, x0)
  end
end

for i in [100, 10000, 100000]
  x0 = rand(3)
  a = @benchmark jac1($x0, $i)
  b = @benchmark jac2($x0, $i)
  println("\nn = $i")
  println("jac1 mean: $(mean(a))")
  println("jac2 mean: $(mean(b))")
end
#=
Result: both versions have identical speed
=#
