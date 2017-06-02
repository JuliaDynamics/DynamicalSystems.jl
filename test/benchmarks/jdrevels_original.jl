using ForwardDiff, BenchmarkTools, StaticArrays

# Using StaticArrays #
#--------------------#

@inline eom(u) = @SVector [10.0(u[2]-u[1]), u[1]*(28.0-u[3]) - u[2], u[1]*u[2] - (8/3)*u[3]]

@inline function manual_jacobian(u)
    i = one(eltype(u))
    return @SMatrix [       -10i    10i     0i;
                     (28i - u[3])    -i  -u[1];
                             u[2]  u[1]  -8i/3]
end

@inline forward_jacobian(u) = ForwardDiff.jacobian(eom, u)

u = SVector{3}(rand(3))

# sanity check
@assert manual_jacobian(u) == forward_jacobian(u)

println("---- eom(::SVector) benchmark results ----")
display(@benchmark eom($u))

println("---- manual_jacobian(::SVector) benchmark results ----")
display(@benchmark manual_jacobian($u))

println("---- forward_jacobian(::SVector) benchmark results ----")
display(@benchmark forward_jacobian($u))

# Using Base Arrays #
#-------------------#

@inline function eom!(du, u)
    du[1] = 10.0(u[2]-u[1])
    du[2] = u[1]*(28.0-u[3]) - u[2]
    du[3] = u[1]*u[2] - (8/3)*u[3]
end

@inline function manual_jacobian!(J, u)
    i = one(eltype(J))
    J[1,1] = -10i
    J[1,2] = 10i
    J[1,3] = 0i
    J[2,1] = 28i - u[3]
    J[2,2] = -i
    J[2,3] = -u[1]
    J[3,1] = u[2]
    J[3,2] = u[1]
    J[3,3] = -8i/3
    return J
end

@inline forward_jacobian!(J, du, u, cfg) = ForwardDiff.jacobian!(J, eom!, du, u, cfg)

u = rand(3)
du = zeros(3)
J = zeros(3, 3)
cfg = ForwardDiff.JacobianConfig(eom!, du, u)

# sanity check
@assert manual_jacobian!(J, u) == forward_jacobian!(J, du, u, cfg)

println("---- eom!(::Vector, ::Vector) benchmark results ----")
display(@benchmark eom!($du, $u))

println("---- manual_jacobian!(::Matrix, ::Vector) benchmark results ----")
display(@benchmark manual_jacobian!($J, $u))

println("---- forward_jacobian!(::Matrix, ::Vector, ::Vector, ::JacobianConfig) benchmark results ----")
display(@benchmark forward_jacobian!($J, $du, $u, $cfg))
