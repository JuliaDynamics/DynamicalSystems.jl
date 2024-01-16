using DynamicalSystems, GLMakie, ModelingToolkit

@parameters σ=28.0 ρ=10.0 β=8/3
@variables t x(t)=5.0 y(t)=0.0 z(t)=0.0 w(t)
D = Differential(t)

eqs = [D(x) ~ σ * (y - x),
    D(y) ~ x * (ρ - z) - y,
    D(z) ~ w - β * z,
    w ~ x*y]

@named lorenz = ODESystem(eqs)
sys = structural_simplify(lorenz)

tspan = (0.0, 100.0)
prob = ODEProblem(sys, nothing, (0.0, 100.0))
ds = CoupledODEs(prob)

parameter_sliders = Dict(
    σ => 5:0.1:50,
    ρ => 0:0.1:200,
    β => 0:0.01:5,
)

f(u) = u[1]*u[3]

# observables = [1, z, w, f]
observables = [1, 2]

fig, dsobs = interactive_trajectory_timeseries(ds, observables;
    parameter_sliders
)

display(fig)