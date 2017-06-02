using StaticArrays, OrdinaryDiffEq, DiffEqBase, BenchmarkTools
tspan = (0.0,1.0)
println("Testing speed of solving ODEs\n")
println("---Using StaticArrays: f(t, u) = du::SVector")
#--------------------#
@inline eom(t, u) = SVector{3}(10.0(u[2]-u[1]), u[1]*(28.0-u[3]) - u[2], u[1]*u[2] - (8/3)*u[3])

# println("ODEProblem creation:")
# display(@benchmark ODEProblem(eom,SVector{3}(rand(3)),$tspan))
prob = ODEProblem(eom,SVector{3}(10rand(3)),tspan)
println("Solution (for total time 1)")
sol = solve(prob,Tsit5(),reltol=1e-8,abstol=1e-8)
display(@benchmark solve($prob,Tsit5(),reltol=1e-8,abstol=1e-8))
sleep(0.1)

println("\n---Using Base Arrays: f(t, u) = du::Vector")
#--------------------#
@inline eom(t, u) = [10.0(u[2]-u[1]), u[1]*(28.0-u[3]) - u[2], u[1]*u[2] - (8/3)*u[3]]

# println("ODEProblem creation:")
# display(@benchmark ODEProblem(eom,rand(3),$tspan))
prob = ODEProblem(eom,10rand(3),tspan)
println("Solution (for total time 1)")
sol = solve(prob,Tsit5(),reltol=1e-8,abstol=1e-8)
display(@benchmark solve($prob,Tsit5(),reltol=1e-8,abstol=1e-8))
sleep(0.1)



println("\n---Using Base & in-place: f!(t, u, du)::Vector")

@inline function eom!(t, u, du)
    du[1] = 10.0(u[2]-u[1])
    du[2] = u[1]*(28.0-u[3]) - u[2]
    du[3] = u[1]*u[2] - (8/3)*u[3]
end

# println("ODEProblem creation:")
# display(@benchmark ODEProblem(eom!,rand(3),$tspan))
prob = ODEProblem(eom!,10rand(3),tspan)
println("Solution (for total time 1)")
sol = solve(prob,Tsit5(),reltol=1e-8,abstol=1e-8)
display(@benchmark solve($prob,Tsit5(),reltol=1e-8,abstol=1e-8))
sleep(0.1)


#=
println("\n---Using MVector & in-place: f!(t, u, du)::MVector")
@inline function eomM!(t, u::MVector, du::MVector)
    du[1] = 10.0(u[2]-u[1])
    du[2] = u[1]*(28.0-u[3]) - u[2]
    du[3] = u[1]*u[2] - (8/3)*u[3]
end

println("ODEProblem creation:")
probM = ODEProblem(eomM!,MVector{3}(10rand(3)),tspan)
display(@benchmark ODEProblem(eomM!,MVector{3}(10rand(3)),$tspan))
println("Solution (for total time 1)")
sol = solve(probM,Tsit5(),reltol=1e-8,abstol=1e-8)
display(@benchmark solve($probM,Tsit5(),reltol=1e-8,abstol=1e-8))
=# #this doesnt work currently
