# Local Modeling

Local Modeling predicts timeseries using a delay embedded state space reconstruction.
It finds the nearest neighbors of a query within this reconstructed space and applies
a local model to make a prediction.

Local Modeling is an effective tool for timeseries of low-dimensional chaotic attractors.

As a part of the `DynamicalSystems` package, the delay embedding is natively supported and
does not need to be done beforehand. Instead embedding parameters can be passed directly
to the prediction function `localmodel_tsp`.

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

ds = Systems.roessler(ones(3))
dt = 0.1
data = trajectory(ds, 1000; dt=dt)
N_train = 6001
s_train = data[1:N_train, 1:2]
s_test  = data[N_train:end, 1:2]

p = 200

s_pred = localmodel_tsp(s_train, 3, 15, p)
s_pred_10 = localmodel_tsp(s_train, 3, 15, p÷10;  stepsize = 10)

using PyPlot
figure(figsize=(8,4))
plot(580:dt:600, s_train[5801:end,1], label = "training x (trunc.)", color = "C1")
plot(580:dt:600, s_train[5801:end,2], label = "training y (trunc.)", color = "green")
plot(600:dt:(600+p*dt), s_test[1:p+1,1], color = "C3", label = "actual x")
plot(600:dt:(600+p*dt), s_test[1:p+1,2], color = "C3", label = "actual y")
plot(600:dt:(600+p*dt), s_pred[:,1], color = "C0", ls="--", label="predicted x")
plot(600:dt:(600+p*dt), s_pred[:,1], color = "C2", ls="--", label="predicted y")
plot(600:dt*10:(600+p*dt), s_pred_10[:,1], color = "C0",
    lw=0, marker="s", label="pred. x step=10")
plot(600:dt*10:(600+p*dt), s_pred_10[:,2], color = "C2",
    lw=0, marker="s", label="pred. y step=10")
title("AverageLocalModel, Training points: $(N_train), attempted prediction: $(p)",
size = 14)
xlabel("\$t\$"); ylabel("\$x\$")
legend(loc="top left")
tight_layout()
```
![2D Average local model prediction](https://i.imgur.com/4gtzGFG.png)

## Error Measures
Here you write why this exists and what it offers.
```@docs
MSE1
MSEp
```
