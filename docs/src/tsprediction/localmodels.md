# Local Modeling

Local Modeling predicts timeseries using a delay embedded state space reconstruction.
It finds the nearest neighbors of a query within this reconstructed space and applies
a local model to make a prediction.

Local Modeling is an effective tool for timeseries of low-dimensional chaotic attractors.


!!! tip "Reconstruction parameters"
    We also offer functions that estimate good parameters for
    a [`Reconstruction`](@ref). These can be found [here](definition/reconstruction/#estimating-reconstruction-parameters).

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
Predicting multivariate timeseries works the same as with skalar timeseries.
```julia
using DynamicalSystemsBase
using TimeseriesPrediction
using StaticArrays: SVector

ds = Systems.roessler(ones(3))
dt = 0.1
data = trajectory(ds, 1000; dt=dt)
N_train = 1001
s_train = data[1:N_train, SVector(1,2)]  
#Identical to data[1:N_train, 1:2] but much faster
s_test  = data[N_train:end,SVector(1,2)]

p = 200
s_pred_10 = localmodel_tsp(s_train, 3, 15, p÷10;  stepsize = 10)

using PyPlot
begin
    figure(figsize=(8,4))
    plot(90:dt:100, s_train[901:end,1], label = "x (trunc.)", color = "C1")
    plot(100:dt:(100+p*dt), s_test[1:p+1,1], color = "C1")
    plot(100:dt*10:(100+p*dt), s_pred_10[:,1], color = "blue",
    lw=0, marker="s", label="pred. x step=10")
    plot(100:dt:(100+p*dt), s_test[1:p+1,2], color = "green")
    plot(90:dt:100, s_train[901:end,2], label = "y (trunc.)", color = "green")
    plot(100:dt*10:(100+p*dt), s_pred_10[:,2], color = "red",
        lw=0, marker="s", label="pred. y step=10")
    plot([100,100],[-12,12], "--", color="black")
    title("AverageLocalModel, Training points: $(N_train),
     attempted prediction: $(p)", size = 14)
    xlabel("\$t\$"); ylabel("\$x\$")
    legend(loc="upper left")
    tight_layout()
end

```
![2D Average local model prediction](https://i.imgur.com/Dm10TGO.png)

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
