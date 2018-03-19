# Local Modeling

Local Modeling predicts timeseries using a delay embedded state space reconstruction.
It finds the nearest neighbors of a query point within this reconstructed space and applies
a local model to make a prediction. "Local" model refers to the fact that the images (future points) of the [`neighborhood`](@ref) of a point are the only component used
to make a prediction.

This Local Modeling is proven to be a very effective tool for timeseries of low-dimensional chaotic attractors.


!!! tip "Reconstruction parameters"
    Don't forget that **DynamicalSystems.jl** also has functions for estimating
    good parameters for
    a [`Reconstruction`](@ref): [`estimate_delay`](@ref) and [`estimate_dimension`](@ref).

## Local Model Prediction
```@docs
localmodel_tsp
AbstractLocalModel
```
## 1D Example
We will predict the future of a (relatively simple) timeseries:
```julia
using TimeseriesPrediction

ds = Systems.roessler(ones(3))
dt = 0.1
data = trajectory(ds, 1000; dt=dt)
N_train = 6001
s_train = data[1:N_train, 1]
s_test  = data[N_train:end,1]

method = AverageLocalModel(2)
ntype = FixedMassNeighborhood(2)

p = 500
s_pred = localmodel_tsp(s_train, 3, 15, p; method=method, ntype=ntype)

s_pred_10 = localmodel_tsp(s_train, 3, 15, p÷10;
    method=method, ntype=ntype, stepsize = 10)

using PyPlot
figure()
plot(550:dt:600, s_train[5501:end], label = "training (trunc.)", color = "C1")
plot(600:dt:(600+p*dt), s_test[1:p+1], color = "C3", label = "actual signal")
plot(600:dt:(600+p*dt), s_pred, color = "C0", ls="--", label="predicted")
plot(600:dt*10:(600+p*dt), s_pred_10, color = "C2",
    lw=0, marker="s", label="pred. step=10")
title("AverageLocalModel, Training points: $(N_train), attempted prediction: $(p)",
size = 18)
xlabel("\$t\$"); ylabel("\$x\$")
legend(loc="upper left")
tight_layout()
```
![Average local model prediction](https://i.imgur.com/VJSjHMI.png)

## 2D Example
Predicting multivariate timeseries works the same as with scalar timeseries.
```julia
using DynamicalSystemsBase
using TimeseriesPrediction
using StaticArrays: SVector

ds = Systems.roessler(ones(3))
dt = 0.1
data = trajectory(ds, 1000; dt=dt)
N_train = 1501
s_train = data[1:N_train, SVector(1,2)]
#Identical to data[1:N_train, 1:2] but much faster
s_test  = data[N_train:end, SVector(1,2)]

p = 100
stepsize = 5
s_pred_10 = localmodel_tsp(s_train, 3, 15, p;  stepsize = stepsize)

using PyPlot; figure(figsize=(12,6))

idx_prev = 200 # how many previous points to show
tf = Int((N_train - 1)*dt) # final time of test set

# Plot real x-coordinate
plot((tf - idx_prev*dt):dt:tf, s_train[N_train-idx_prev:end,1],
    label = "real x", color = "C1")
plot(tf:dt:(tf+p*dt*stepsize), s_test[1:p*stepsize+1,1], color = "C1")
# Plot predicted x-coordinate
plot(tf:dt*stepsize:(tf+p*dt*stepsize), s_pred_10[:,1], color = "C0",
lw=0.5, marker="s", ms = 4.0, label="pred. x")

# Plot real y-coordinate
plot((tf - idx_prev*dt):dt:tf, s_train[N_train-idx_prev:end,2],
    label = "real y", color = "C2")
plot(tf:dt:(tf+p*dt*stepsize), s_test[1:p*stepsize+1,2], color = "C1")
# Plot predicted y-coordinate
plot(tf:dt*stepsize:(tf+p*dt*stepsize), s_pred_10[:,2], color = "C4",
lw=0.5, marker="s", ms = 4.0, label="pred. y")

# Plot separatrix
plot([tf,tf],[-12,12], "--", color="black", alpha = 0.5)

title("AverageLocalModel, Training points: $(N_train), attempted prediction: $(p), step=$(stepsize)", size = 14)
xlabel("\$t\$"); ylabel("\$x\$")
legend(loc="upper left")
tight_layout()
```
![2D Average local model prediction](https://i.imgur.com/yjPSvx9.png)

## Error Measures
Being able to evaluate model performance without looking at plots can be very helpful
when trying to quantify its error as well as finding good parameters in the first place.

```@docs
MSEp
```

Here is an example function that employs `MSEp` to find good parameters. It takes in
a timeseries `s` and ranges for the dimensions, delays, number of nearest neighbors,
and weighting degrees to try. Keyword arguments are `valid_len`, which is the number of
prediction steps, and `num_tries` the number of different starting points to choose.

It then calculates `MSEp` for all parameter combinations and returns the best parameter
set.
```julia
function estimate_param(s::AbstractVector,
    dims, delay, K, N; valid_len=100, num_tries=50)
    Result = Dict{SVector{4,Int},Float64}()
    step = 1
    for n ∈ N
        method = LocalAverageModel(n)
        for D ∈ dims, τ ∈ delay
            s_train = @view s[1:end-D*τ-valid_len-num_tries-50]
            s_test = @view s[end-(D-1)*τ-valid_len-num_tries:end]
            R = Reconstruction(s_train,D,τ)
            R_test = Reconstruction(s_test,D,τ)
            tree = KDTree(R[1:end-1])
            for k ∈ K
                ntype = FixedMassNeighborhood(k)
                Result[@SVector([D,τ,k,n])] =
                MSEp(R, tree, R_test, valid_len; method=method, ntype=ntype,
                    stepsize=stepsize)
            end
        end
    end
    best_param = collect(keys(Result))[findmin(values(Result))[2]]
    return best_param
end
```
