using DynamicalSystems, GLMakie, ModelingToolkit
# Import canonical time and time-derivative from MTK,
# however use the unitless versions as we don't need units here
using ModelingToolkit: t_nounits as t, D_nounits as D

# Define the variables and parameters in symbolic format
@parameters begin
    a = 0.29
    b = 0.14
    c = 4.52
    d = 1.0
end
@variables begin
    x(t) = 10.0
    y(t) = 0.0
    z(t) = 1.0
    nlt(t) # nonlinear term
end

# Create the equations of the model
eqs = [
    D(x) ~ -y - z,
    D(y) ~ x + a*y,
    D(z) ~ b + nlt - z*c,
    nlt ~ d*z*x, # observed variable
]

# Create the model via ModelingToolkit
@named roessler = ODESystem(eqs, t)
# Do not split parameters so that integer indexing can be used as well
model = structural_simplify(roessler; split = false)
# Cast it into an `ODEProblem` and then into a `DynamicalSystem`.
# Due to low-dimensionality it is preferred to cast into out of place
prob = ODEProblem{false}(model, nothing, (0.0, Inf); u0_constructor = x->SVector(x...))
ds = CoupledODEs(prob)
# If you have "lost" the model, use:
model = referrenced_sciml_model(ds)

# Define which parameters will be interactive during the simulation
parameter_sliders = Dict(
    # can use integer indexing
    1 => 0:0.01:1,
    # the global scope symbol
    b => 0:0.01:1,
    # the symbol obtained from the MTK model
    model.c => 0:0.01:10,
    # or a `Symbol` with same name as the parameter
    # (which is the easiest and recommended way)
    :d => 0.8:0.01:1.2,
)

# Define what variables will be visualized as timeseries
norm(u) = sqrt(u[1]*u[1] + u[2]*u[2])
observables = [
    1,         # can use integer indexing,
    z,         # MTK state variable (unknown)
    model.nlt, # MTK observed variable
    :y,        # `Symbol` instance with same name
    norm,      # or arbitrary function of the state
]

# Define what variables will be visualized as state space trajectory
# same as above, any indexing works, but ensure to make the vector `Any`
# so that integers are not converted to symbolic variables
idxs = Any[1, y, 3]

u0s = [
    # we can specify dictionaries, each mapping the variable to its value
    # un-specified variables get the value they currently have in `ds`
    Dict(:x => -4, :y => -4, :z => 0.1),
    Dict(:x => 4, :y => 3, :z => 0.1),
    Dict(:x => -5.72),
    Dict(:x => 5.72, :y => 0.28, :z => 0.21),
]

update_theme!(fontsize = 14)
tail = 1000

fig, dsobs = interactive_trajectory_timeseries(ds, observables, u0s;
    parameter_sliders, Δt = 0.01, tail, idxs,
    figure = (size = (1100, 650),)
)

step!(dsobs, 2tail)

display(fig)