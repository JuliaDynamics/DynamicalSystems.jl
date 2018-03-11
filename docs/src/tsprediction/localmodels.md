# Local Modelling

Here is a paragraph that says describes local modelling in 2 sentences. (it uses nearest neighbors).

Here is another one that says when is local modelling good.

You also say that it reconstructs.

!!! tip "Reconstruction parameters"
    Don't forget that we also offer functions that estimate good parameters for
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
To predict many timeseries you have to reconstruct them manually, and get the result
out manually as well:

Here you put an example where you reconstruct 2D timeseries of roessler,
predict both of them and then plot both predictions versus both real signals.


## Error Measures
Here you write why this exists and what it offers.
```@docs
MSE1
MSEp
```
