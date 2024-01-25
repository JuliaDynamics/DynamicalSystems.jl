using DynamicalSystems, GLMakie, ModelingToolkit

@variables t
D = Differential(t)
@parameters begin
    a = 0.29
    b = 0.14
    c = 4.52
end
@variables begin
    x(t) = 10.0
    y(t) = 0.0
    z(t) = 1.0
    nlt(t) # nonlinear term
end

eqs = [
    D(x) ~ -y - z,
    D(y) ~ x + a*y,
    D(z) ~ b + nlt - z*c,
    nlt ~ z*x, # observed variable
]

@named roessler = ODESystem(eqs)
model = structural_simplify(roessler)

tspan = (0.0, 100.0)
prob = ODEProblem(sys)
ds = CoupledODEs(prob)
model = referrenced_sciml_model(ds)

parameter_sliders = Dict(
    # can use integer indexing
    1 => 0:0.01:1,
    # the global scope symbol
    b => 0:0.01:1,
    # or the symbol obtained from the MTK model
    model.c => 0:0.01:10,
)

norm(u) = sqrt(u[1]*u[1] + u[2]*u[2])

observables = [
    1,         # can use integer indexing,
    z,         # MTK state variable
    model.nlt, # MTK observed variable
    norm,      # or arbitrary function of the state
]

# same as above, any indexing works:
idxs = Any[1, y, 3]

u0s = [
    # no fancy indexing here yet; numbers must correspond to state variables
    [-4.0, -4, 0.1],
    [4.0, 4, 0.2],
    [5.72, 0.28, 0.21],
    [-5.72, 0.0, 0.0],
]

update_theme!(fontsize = 14)

fig, dsobs = interactive_trajectory_timeseries(ds, observables, u0s;
    parameter_sliders, statespace_axis = true, Δt = 0.01,
    tail = 1000, idxs,
    figure = (size = (1100, 650),)
)

display(fig)