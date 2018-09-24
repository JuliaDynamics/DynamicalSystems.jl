# Local Modeling

!!! tip "Reconstruction parameters"
    Don't forget that **DynamicalSystems.jl** also has functions for estimating good parameters for delay embedding: [`estimate_delay`](@ref) and [`estimate_dimension`](@ref).

## Local Model Prediction
```@docs
localmodel_tsp
AbstractLocalModel
```
## Single Timeseries Example
We will predict the future of a (relatively simple) timeseries:
```@example tspred
using TimeseriesPrediction # Re-exports DynamicalSystemsBase

ds = Systems.roessler(0.1ones(3))
dt = 0.1
data = trajectory(ds, 1000; dt=dt)
N_train = 6001
s_train = data[1:N_train, 1]
s_test  = data[N_train:end,1]

ntype = FixedMassNeighborhood(3)

p = 500
s_pred = localmodel_tsp(s_train, 4, 15, p; ntype=ntype)

using PyPlot
figure()
plot(550:dt:600, s_train[5501:end], label = "training (trunc.)", color = "C1")
plot(600:dt:(600+p*dt), s_test[1:p+1], color = "C3", label = "actual signal")
plot(600:dt:(600+p*dt), s_pred, color = "C0", ls="--", label="predicted")

title("Pool of points: $(N_train), predicted points: $(p)")
xlabel("\$t\$"); ylabel("\$x\$")
legend(loc="upper left")
tight_layout()
savefig("single_tspred.png"); nothing # hide
```
![](single_tspred.png)


## Multiple Timeseries Example
Predicting multivariate timeseries works the same as with scalar timeseries.
```@example tspred
using TimeseriesPrediction

ds = Systems.roessler(ones(3))
dt = 0.1
data = trajectory(ds, 1000; dt=dt)
N_train = 1501
s_train = data[1:N_train, SVector(1,2)]
#Identical to data[1:N_train, 1:2] but much faster
s_test  = data[N_train:end, SVector(1,2)]

p = 100; stepsize = 5
s_pred_10 = localmodel_tsp(s_train, 3, 15, p;  stepsize = stepsize)

using PyPlot; figure(figsize=(12,6))

idx_prev = 200 # how many previous points to show
tf = Int((N_train - 1)*dt) # final time of test set

# Plot real x-coordinate
plot((tf - idx_prev*dt):dt:tf, s_train[N_train-idx_prev:end,1],
    label = "real x", color = "C1")
plot(tf:dt:(tf+p*dt*stepsize), s_test[1:p*stepsize+1,1], color = "C1")
# Plot predicted x-coordinate
plot(tf:dt*stepsize:(tf+p*dt*stepsize), s_pred_10[:,1], color = "C3",
lw=0.5, marker="s", ms = 4.0, label="pred. x", alpha = 0.75)

# Plot real y-coordinate
plot((tf - idx_prev*dt):dt:tf, s_train[N_train-idx_prev:end,2],
    label = "real y", color = "C0")
plot(tf:dt:(tf+p*dt*stepsize), s_test[1:p*stepsize+1,2], color = "C0")
# Plot predicted y-coordinate
plot(tf:dt*stepsize:(tf+p*dt*stepsize), s_pred_10[:,2], color = "C9",
lw=0.5, marker="s", ms = 4.0, label="pred. y", alpha = 0.75)

# Plot start of prediction
plot([tf,tf],[-12,12], "--", color="black", alpha = 0.5)

title("Pool of points: $(N_train), predicted points: $(p)")
xlabel("\$t\$"); ylabel("\$x, y\$")
legend(loc="upper left")
tight_layout()
savefig("multiple_tspred.png"); nothing # hide
```
![](multiple_tspred.png)

## Error Measures
Being able to evaluate model performance without looking at plots can be very helpful when trying to quantify its error as well as finding good parameters in the first place.

```@docs
MSEp
```

Here is an example function that employs `MSEp` to find good parameters. It takes in
a timeseries `s` and ranges for the dimensions, delays and number of nearest neighbors to
 try. Keyword arguments are `valid_len`, which is the number of
prediction steps, and `num_tries` the number of different starting points to choose.

It then calculates `MSEp` for all parameter combinations and returns the best parameter
set.
```@example tspred
function estimate_param(s::AbstractVector,
    dims, delay, K; valid_len=100, num_tries=50)
    result = Dict{NamedTuple,Float64}()
    step = 1
    for D ∈ dims, τ ∈ delay
        s_train = @view s[1:end-D*τ-valid_len-num_tries-50]
        s_test = @view s[end-(D-1)*τ-valid_len-num_tries:end]
        R = reconstruct(s_train,D,τ)
        R_test = reconstruct(s_test,D,τ)
        tree = KDTree(R[1:end-1])
        for k ∈ K
            ntype = FixedMassNeighborhood(k)
            result[(D=D,τ=τ,k=k)] =
            MSEp(R, tree, R_test, valid_len; ntype=ntype)
        end
    end
    best_param = collect(keys(result))[findmin(collect(values(result)))[2]]
    return best_param
end

ds = Systems.roessler(0.1ones(3))
dt = 0.1
data = trajectory(ds, 1000; dt=dt)
N_train = 6001
s_train = data[1:N_train, 1]
s_test  = data[N_train:end,1]

D, τ, k = estimate_param(s_train, 1:4, [10, 15, 30], 2:4)
```

## Cool animation!
This is an animation of timeseries prediction of the `z` variable of the Roessler system.
On the left you can see the time evolution of the whole system with the
chaotic attractor indicated in gray. The right side is a plot of the `z` component of the
system. The actual values are displayed in green. In red you can see the iteratively
predicted version. As training set it used part of the attractor shown in gray on the left.

<video controls="controls">
<source src="https://github.com/JuliaDynamics/JuliaDynamicsDocumentation.jl/blob/master/animations/tspred/tspred_animation_zRossler.mp4?raw=true" type="video/mp4">
</video>

You can find the script that produced this animation in
`DynamicalSystems/coolanimations/roessler_Z_tspred.jl`.
