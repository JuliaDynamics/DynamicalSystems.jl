#= Conclusions:
Use SVectors and assignment form everywhere! :D
=#
println("Testing speed of Jacobian calculation for different types of e.o.m.")

using ForwardDiff, BenchmarkTools, StaticArrays

println("\n\nUsing StaticArray assignment form")
#--------------------#

# @inline eom(u) = @SVector [10.0(u[2]-u[1]), u[1]*(28.0-u[3]) - u[2], u[1]*u[2] - (8/3)*u[3]]
@inline eom(u) = SVector{3}(10.0(u[2]-u[1]), u[1]*(28.0-u[3]) - u[2], u[1]*u[2] - (8/3)*u[3])

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

println("\n\nUsing MVector assignment form")
#--------------------#

# @inline eom(u) = @MVector [10.0(u[2]-u[1]), u[1]*(28.0-u[3]) - u[2], u[1]*u[2] - (8/3)*u[3]]
@inline eom_mm(u) = MVector{3}(10.0(u[2]-u[1]), u[1]*(28.0-u[3]) - u[2], u[1]*u[2] - (8/3)*u[3])

@inline function manual_jacobian_mm(u)
    i = one(eltype(u))
    return @MMatrix [       -10i    10i     0i;
                     (28i - u[3])    -i  -u[1];
                             u[2]  u[1]  -8i/3]
end

@inline forward_jacobian_mm(u) = ForwardDiff.jacobian(eom, u)

u = MVector{3}(rand(3))

# sanity check
@assert manual_jacobian_mm(u) == forward_jacobian_mm(u)

println("---- eom(::MVector) benchmark results ----")
display(@benchmark eom($u))

println("---- manual_jacobian(::MVector) benchmark results ----")
display(@benchmark manual_jacobian($u))

println("---- forward_jacobian(::MVector) benchmark results ----")
display(@benchmark forward_jacobian($u))



println("\n\nUsing Base Arrays (in-place form)")
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


println("\n\nUsing MArrays & in-place")
#------------------------------------------------#

@inline function eom_m!(du::MVector, u::MVector)
    du[1] = 10.0(u[2]-u[1])
    du[2] = u[1]*(28.0-u[3]) - u[2]
    du[3] = u[1]*u[2] - (8/3)*u[3]
end

@inline function manual_jacobian!(J::MMatrix, u::MVector)
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

@inline forward_jacobian!(J, du, u, cfg) = ForwardDiff.jacobian!(J, eom_m!, du, u, cfg)

u = MVector{3}(rand(3))
du = MVector{3}(rand(3))
J = MMatrix{3,3}(rand(3,3))
cfg = ForwardDiff.JacobianConfig(eom_m!, du, u)

# sanity check
@assert manual_jacobian!(J, u) == forward_jacobian!(J, du, u, cfg)

println("---- eom_m!(::MVector, ::MVector) benchmark results ----")
display(@benchmark eom!($du, $u))

println("---- manual_jacobian!(::MMatrix, ::MVector) benchmark results ----")
display(@benchmark manual_jacobian!($J, $u))

println("---- forward_jacobian!(::MMatrix, ::MVector, ::MVector, ::JacobianConfig) benchmark results ----")
display(@benchmark forward_jacobian!($J, $du, $u, $cfg))
