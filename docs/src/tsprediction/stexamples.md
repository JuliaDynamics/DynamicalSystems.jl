```@meta
EditURL = "https://github.com/TRAVIS_REPO_SLUG/blob/master/"
```

# Spatio-Temporal Prediction Examples
In this page we are simply running files from the
`examples` folder of the `TimeseriesPrediction` package.

This is how you can (programmatically) find this folder:

```@example stexamples
using TimeseriesPrediction
exdir = dirname(dirname(pathof(TimeseriesPrediction)))*"/examples"
```

## Temporal Prediction: Kuramoto-Sivashinsky
*(this requires `FFTW` to be installed)*

This example predicts the temporal evolution of a one-dimensional field
U, along with a time vector T,
which has to be represented as vectors of vectors.
Where the field comes from does not matter, but to make the example
runnable we load one of the test systems of `TimeseriesPrediction`.

In this example we use the solution of Kuramoto Sivashinsky equation.

Importantly, the results are compared with the "real" evolution of the
system.

In the plots, the x axis is space and y axis is time.

### Produce field U (Kuramoto Sivashinsky)

```@example stexamples
using PyPlot
using TimeseriesPrediction

testdir = dirname(dirname(pathof(TimeseriesPrediction)))*"/test"
@assert isdir(testdir)
include(testdir*"/ks_solver.jl")

Ntrain = 10000
p = 100
N = Ntrain + p

U, T = KuramotoSivashinsky(64, 22, N÷4, 0.25)
```

### Temporal prediction of field U

```@example stexamples
Q = length(U[1]) # spatial length
pool = U[1:Ntrain]
test = U[Ntrain:N]

D = 10
τ = 1
B = 10
k = 1
ntype = FixedMassNeighborhood(4)
method = AverageLocalModel()

em = cubic_shell_embedding(pool, D,τ,B,k,PeriodicBoundary())
pcaem= PCAEmbedding(pool,em)

@time pred = temporalprediction(pool,pcaem, p;ntype=ntype, method=method, progress = false)

err = [abs.(test[i]-pred[i]) for i=1:p+1]
println("Maximum error: ", maximum(maximum(e) for e in err))
```

### Plot the result

Deduce field extremal values

```@example stexamples; continued = true
vmax = max(maximum(maximum(s) for s in test),
           maximum(maximum(s) for s in pred))
vmin = min(minimum(minimum(s) for s in test),
           minimum(minimum(s) for s in pred))
```

Transform data for imshow

```@example stexamples
ptest = cat(test..., dims = 2)
ppred = cat(pred..., dims = 2)
perr = cat(err..., dims = 2)
```

plot plot plot

```@example stexamples
fig = figure(figsize=(8,8))
ax1 = subplot2grid((3,1), (0,0))
ax2 = subplot2grid((3,1), (1,0))
ax3 = subplot2grid((3,1), (2,0));

im1 = ax1[:imshow](ppred, cmap="viridis", vmin = vmin, vmax = vmax,
aspect = "auto", extent = (T[Ntrain], T[N], 1, Q))
im2 = ax2[:imshow](ptest, cmap="viridis", vmin = vmin, vmax = vmax,
aspect = "auto", extent = (T[Ntrain], T[N], 1, Q))
im3 = ax3[:imshow](perr, cmap="inferno", vmin = 0, vmax = vmax-vmin,
aspect = "auto", extent = (T[Ntrain], T[N], 1, Q))

for (j, (im, ax)) in enumerate(zip([im1,im2,im3], [ax1,ax2,ax3]))
    colorbar(im, ax = ax, fraction=0.04, pad=0.01)# format="%.1f")
    if j < 3
        ax[:set_xticklabels]([])
    end
end
ax1[:set_title]("Prediction")
ax2[:set_title]("Real evolution")
ax3[:set_title]("Absolute error")

ax2[:set_ylabel]("space")
ax3[:set_xlabel]("time")
tight_layout(w_pad=0.1, h_pad=0.00001)
```

## Cross Prediction: Barkley Model

This example cross-predicts a field U from a field V.
Both fields have to be represented as vectors of matrices.
Where the fields come from does not matter, but to make the example
runnable we load one of the test systems of TimeseriesPrediction.

This example uses cubic shell embedding and a linear Barkley model.

Importantly, the results are compared with the "real" evolution of the
system.

### Simulate a test system

```@example stexamples
using PyPlot
using TimeseriesPrediction

testdir = dirname(dirname(pathof(TimeseriesPrediction)))*"/test"
@assert isdir(testdir)
include(testdir*"/system_defs.jl")

Nx = 50
Ny = 50
Tskip = 100
Ttrain = 500
Ttest = 10
T = Tskip +Ttrain + Ttest

U, V = barkley_const_boundary(T, Nx, Ny)
```

### Cross predict field U from field V

```@example stexamples
D = 2
τ = 1
B = 2
k = 1
bc = ConstantBoundary(20.0)

source_train = V[Tskip + 1:Tskip + Ttrain]
target_train = U[Tskip + 1:Tskip + Ttrain]
source_pred  = V[Tskip + Ttrain  - D*τ + 1:  T]
target_test  = U[Tskip + Ttrain  - D*τ + 1:  T]

em = cubic_shell_embedding(source_train, D,τ,B,k,bc)
pcaem = PCAEmbedding(source_train, em) # PCA speeds things up!

@time target_pred = crossprediction(source_train, target_train, source_pred, em;
progress = false)

err = [abs.(target_test[i]-target_pred[i]) for i=1:Ttest]

println("Maximum error: ", maximum(maximum(e) for e in err))
```

### Plot prediction

Deduce field extremal values

```@example stexamples; continued = true
source_max = maximum(maximum(s) for s in source_pred)
target_max = max(maximum(maximum(s) for s in target_test),
                 maximum(maximum(s) for s in target_pred))
source_min = minimum(minimum(s) for s in source_pred)
target_min = min(minimum(minimum(s) for s in target_test),
                 minimum(minimum(s) for s in target_pred))
```

Plot various predicted frames (only the last one shown here)

```@example stexamples; continued = true
for i in [1, length(err)÷2, length(err)]

    fig = figure(figsize=(10,10))
    ax1 = subplot2grid((2,2), (0,0))
    ax2 = subplot2grid((2,2), (0,1))
    ax3 = subplot2grid((2,2), (1,0))
    ax4 = subplot2grid((2,2), (1,1))
    im1 = ax1[:imshow](source_pred[i], cmap="viridis", vmin = source_min, vmax = source_max)
    im2 = ax2[:imshow](target_test[i], cmap="cividis", vmin = target_min, vmax = target_max)
    im3 = ax3[:imshow](target_pred[i], cmap="cividis", vmin = target_min, vmax = target_max)
    im4 = ax4[:imshow](err[i], cmap="inferno", vmin = 0, vmax = target_max - target_min)
    for (im, ax) in zip([im1,im2,im3,im4], [ax1,ax2,ax3,ax4])
        ax[:get_xaxis]()[:set_ticks]([])
        ax[:get_yaxis]()[:set_ticks]([])
        colorbar(im, ax = ax, fraction=0.046, pad=0.04)#, format="%.1f")
```

ax[:minorticks_off]()

```@example stexamples
    end
    ax1[:set_title]("Source")
    ax2[:set_title]("Target Test")
    ax3[:set_title]("Target Cross-Pred.")
    ax4[:set_title]("absolute error")
    tight_layout(w_pad=0.6, h_pad=0.00001)
    suptitle("frame $i")
end
```

## Temporal Prediction: Periodic Nonlinear Barkley Model

This example predicts the temporal evolution of a field U,
which has to be represented as vectors of matrices.
Where the field comes from does not matter, but to make the example
runnable we load one of the test systems of TimeseriesPrediction.

This example uses light cone embedding and a nonlinear Barkley model.

Importantly, the results are compared with the "real" evolution of the
system.

### Simulate a test system

```@example stexamples
using PyPlot
using TimeseriesPrediction

testdir = dirname(dirname(pathof(TimeseriesPrediction)))*"/test"
@assert isdir(testdir)
include(testdir*"/system_defs.jl")

Nx = 50
Ny = 50
Tskip = 100
Ttrain = 600
Ttest = 10
T = Tskip + Ttrain + Ttest

U, V = barkley_periodic_boundary_nonlin(T, Nx, Ny)
```

### Temporal prediction of field U

```@example stexamples
D = 2
τ = 1
r₀ = 1
c = 1
bc = PeriodicBoundary()

pool = U[Tskip + 1 : Tskip + Ttrain]
test  = U[Tskip + Ttrain : T]

em = light_cone_embedding(pool, D,τ,r₀,c,bc)
pcaem = PCAEmbedding(pool, em) # PCA speeds things up!

@time pred = temporalprediction(pool, em, Ttest; progress = false)

err = [abs.(test[i]-pred[i]) for i=1:Ttest+1]

println("Maximum error: ", maximum(maximum(e) for e in err))
```

### Plot prediction

Deduce field maximum values

```@example stexamples; continued = true
vmax = max(maximum(maximum(s) for s in test),
           maximum(maximum(s) for s in pred))
vmin = min(minimum(minimum(s) for s in test),
           minimum(minimum(s) for s in pred))
```

plot plot plot

```@example stexamples; continued = true
for i in [1, length(err)÷2, length(err)]

    fig = figure(figsize=(10,3))
    ax1 = subplot2grid((1,3), (0,0))
    ax2 = subplot2grid((1,3), (0,1))
    ax3 = subplot2grid((1,3), (0,2))

    im1 = ax1[:imshow](pred[i], cmap="viridis", vmin = vmin, vmax = vmax)
    im2 = ax2[:imshow](test[i], cmap="viridis", vmin = vmin, vmax = vmax)
    im3 = ax3[:imshow](err[i], cmap="inferno", vmin = 0, vmax = vmax-vmin)
    for (im, ax) in zip([im1,im2,im3], [ax1,ax2,ax3])
        ax[:get_xaxis]()[:set_ticks]([])
        ax[:get_yaxis]()[:set_ticks]([])
        colorbar(im, ax = ax, fraction=0.046, pad=0.04)#, format="%.1f")
```

ax[:minorticks_off]()

```@example stexamples
    end
    ax1[:set_title]("Prediction")
    ax2[:set_title]("Real evolution")
    ax3[:set_title]("Absolute error")
    suptitle("frame $i")
    tight_layout(w_pad=0.6, h_pad=0.00001)
    subplots_adjust(top=0.75)
end
```

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

