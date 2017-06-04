using ForwardDiff, BenchmarkTools, StaticArrays

# Using StaticArrays #
#--------------------#

@inline eom(u) = @SVector [3.8*u[1]*(1-u[1])-0.05*(u[2]+0.35)*(1-2*u[3]),
0.1*((u[2]+0.35)*(1-2*u[3])-1 )*(1-1.9*u[1]),
3.78*u[3]*(1-u[3])+0.2*u[2]]

@inline function manual_jacobian(x)
    return @SMatrix [3.8*(1 - 2x[1]) -0.05*(1-2x[3]) 0.1*(x[2] + 0.35);
    -0.19((x[2] + 0.35)*(1-2x[3]) - 1)  0.1*(1-2x[3])*(1-1.9x[1])  -0.2*(x[2] + 0.35)*(1-1.9x[1]);
    0.0  0.2  3.78(1-2x[3]) ]
end

@inline function manual_jacobian!(J, x)
  J[1,1] = 3.8*(1 - 2x[1]); J[1,2] = -0.05*(1-2x[3]); J[1,3] = 0.1*(x[2] + 0.35);
  J[2,1] = 0.19((x[2] + 0.35)*(1-2x[3]) - 1); J[2,2] = 0.1*(1-2x[3])*(1-1.9x[1])
  J[2,3] = -0.2*(x[2] + 0.35)*(1-1.9x[1])
  J[3,1] = 0.0; J[3,2] = 0.2; J[3,3] = 3.78(1-2x[3])
  return J
end

@inline forward_jacobian(u) = ForwardDiff.jacobian(eom, u)

u = SVector{3}(rand(3)); J = rand(3,3)

# sanity check
@assert manual_jacobian(u) == forward_jacobian(u)

println("---- eom(::SVector) benchmark results ----")
display(@benchmark eom($u))

println("---- manual_jacobian(::SVector) benchmark results ----")
display(@benchmark manual_jacobian($u))

println("---- forward_jacobian(::SVector) benchmark results ----")
display(@benchmark forward_jacobian($u))

# sanity check
@assert manual_jacobian(u) == forward_jacobian(u)

@inline forward_jacobian!(J, u) = ForwardDiff.jacobian!(J, eom, u)

println("---- manual_jacobian_inplace(::SVector) benchmark results ----")
display(@benchmark manual_jacobian!($J, $u))

println("---- forward_jacobian_inplace(::SVector) benchmark results ----")
display(@benchmark forward_jacobian!($J, $u))


# Using Base Arrays #
#-------------------#

@inline function eom!(du, x)
  x1, x2, x3 = x[1], x[2], x[3]
  du[1] = 3.8*x1*(1-x1)-0.05*(x2+0.35)*(1-2*x3)
  du[2] = 0.1*((x2+0.35)*(1-2*x3)-1 )*(1-1.9*x1)
  du[3] = 3.78*x3*(1-x3)+0.2*x2
end

@inline function manual_jacobian!(J, x)
  J[1,1] = 3.8*(1 - 2x[1]); J[1,2] = -0.05*(1-2x[3]); J[1,3] = 0.1*(x[2] + 0.35);
  J[2,1] = 0.19((x[2] + 0.35)*(1-2x[3]) - 1); J[2,2] = 0.1*(1-2x[3])*(1-1.9x[1])
  J[2,3] = -0.2*(x[2] + 0.35)*(1-1.9x[1])
  J[3,1] = 0.0; J[3,2] = 0.2; J[3,3] = 3.78(1-2x[3])
  return J
end

@inline forward_jacobian!(J, du, u, cfg) = ForwardDiff.jacobian!(J, eom!, du, u, cfg)

u = rand(3)
du = zeros(3)
J = zeros(3, 3)
cfg = ForwardDiff.JacobianConfig(eom!, du, u)

println("---- eom!(::Vector, ::Vector) benchmark results ----")
display(@benchmark eom!($du, $u))

println("---- manual_jacobian!(::Matrix, ::Vector) benchmark results ----")
display(@benchmark manual_jacobian!($J, $u))

println("---- forward_jacobian!(::Matrix, ::Vector, ::Vector, ::JacobianConfig) benchmark results ----")
display(@benchmark forward_jacobian!($J, $du, $u, $cfg))
