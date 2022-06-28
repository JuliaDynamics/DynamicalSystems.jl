# Windowed RQA

In some cases, specially with very long time series, it may be suitable to perform the analysis at different points, considering only a limited *window* of data around each observation. The macro [`@windowed`](@ref) modifies the behaviour of the basic functions to calculate RQA parameters in that fashion. For instance, if `rmat` is a 10<sup>4</sup>&times;10<sup>4</sup> recurrence matrix, then
```julia
@windowed determinism(rmat, theiler=2, lmin=3) width=1000 step=100
```
will return a 91-element vector, such that each value is the determinism associated to a 1000-point fragment, starting at every 100 points (i.e. at `1`, `101`, &hellip; `9001`).

The general syntax of that macro is:
```julia
@windowed expr w                 #1
@windowed expr width=w step=s    #2
```
where:

* `expr` is an expression used to calculate RQA parameters
* `w` is the width of the window for relevant data around each point.
* `s` is the step or distance between points where the calculations are done (starting in the first point).

To prevent syntax failures in the expansion of the macro, identify the RQA function (`rqa`, `recurrencerate`, `determinism`,...) directly by its name (avoid aliases), and use simple variable names (not complex expressions) for the arguments. On the other hand, the windowing options `width` and `step` can be given in any order. If `step` is omitted, the calculations are done at every point, and the keyword `width` may be omitted. (However, using `step=1` may be computationally very expensive, and that will provide just overly redundant results around each point, so it is advisable to set `step` a relatively big fraction of the window `width`.)

The value returned by the macro will normally be a vector with the same type of numbers as expected by `expr`. In the case of `@windowed rqa(...) ...`, it will return a NamedTuple with a similar structure as in the default `rqa` function, but replacing scalar values by vectors.

The macro `@windowed` can also be applied to the functions that calculate recurrence matrices (`RecurrenceMatrix`, `CrossRecurrenceMatrix`, `JointRecurrenceMatrix`). That creates a sparse matrix with the same size as if the macro was not used, but only containing valid values for pairs of points that belong to the `w` first main diagonals (i.e. the separation in time from one point to the other is `w` or smaller). The &lsquo;step&rsquo; parameter `s` has no effect on those functions. Such &lsquo;windowed&rsquo; matrices can be used as the input arguments to calculate windowed RQA parameters, obtaining the same results as if the complete matrix was used (under certain conditions, see below). For instance, the following calculations are equivalent:

```julia
# Using complete matrix
rmat = RecurrenceMatrix(x, 1.5)
d = @windowed determinism(rmat) width=1000 step=250

# Using windowed matrix
rmatw = @windowed RecurrenceMatrix(x, 1.5) 1000
d = @windowed determinism(rmatw) width=1000 step=250
```

The main difference between the two alternatives is that the second one will be faster and consume less memory. To ensure the equivalence between both approaches, the window width used to create the matrix must be greater than the one used to calculate the RQA parameters. Otherwise, the computation of RQA parameters might involve data points whose value is not well defined. Besides, the threshold to identify recurrences should be referred to a fixed scale. For instance:

```julia
rmat  =           RecurrenceMatrix(x, 0.1, scale=maximum)
rmatw = @windowed RecurrenceMatrix(x, 0.1, scale=maximum) 1000
rmat[1:1000,1:1000] == rmatw[1:1000,1:1000] # FALSE!!!
```
In this example, the `1000×1000` blocks of both matrices differ, because the threshold `0.1` is scaled with respect to the maximum distance between all points of `x` in `rmat`, but in the case of `rmatw` the scale changes between subsets of points. Something similar may happen if the recurrence matrix is calculated for a fixed recurrence rate (with the option `fixedrate=true`).

## Docstring
```@docs
@windowed
```

## Alternative syntax for `@windowed`

The following ways of using the macro `@windowed` are equivalent:

```julia
y = @windowed f(x,...) w
@windowed y=f(x,...) w
y = @windowed(f(x,...), w)
@windowed(y=f(x,...), w)
```

In all four cases, the width parameter `w` might have been qualified with a keyword as `width=w`. If the step parameter is added, the keyword qualification is mandatory.
