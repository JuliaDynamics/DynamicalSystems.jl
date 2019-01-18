# Quantification & Analysis functions
A [`RecurrenceMatrix`](@ref) can be analyzed in several ways to yield information about the dynamics of the trajectory. All these various *measures* and functions are collectively called "Recurrence Quantification Analysis" (RQA).

To understand how each measure can be useful, we suggest to see the review articles listed in our documentation strings, namely:

1. N. Marwan *et al.*, "Recurrence plots for the analysis of complex systems", *Phys. Reports 438*(5-6), 237-329 (2007).
2. N. Marwan & C.L. Webber, "Mathematical and computational foundations of recurrence quantifications", in: Webber, C.L. & N. Marwan (eds.), *Recurrence Quantification Analysis. Theory and Best Practices*, Springer, pp. 3-43 (2015).

You can also check the wikipedia page for [Recurrence quantification analysis](https://en.wikipedia.org/wiki/Recurrence_quantification_analysis).

The functions described in this page all accept a recurrence matrix (`x`), see the [previous page](rplots).

## RQA Measures
### Classical RQA Measures
```@docs
recurrencerate
determinism
dl_average
dl_max
dl_entropy
divergence
trend
```
---

### Extended RQA Measures
```@docs
laminarity
trappingtime
vl_average
vl_max
vl_entropy
```
---

### Recurrence Time Measures
```@docs
meanrecurrencetime
nmprt
```
---

## Keyword table
Since most of the above functions can be fined tuned with keyword arguments, here is a table summarizing them that could be of use:


| Argument  | Default   | Functions | Description |
| --------  | --------  | --------- | -----------
| `scale`   | 1         | `RecurrenceMatrix`<br/>`CrossRecurrenceMatrix`<br/>`JointRecurrenceMatrix` | Function or fixed number to scale the threshold or radius that is used to identify recurrences. Use `maximum` if the threshold is to be taken as a fraction of the maximum distance, `mean` if it is a fraction of the mean distance, etc., and `1` (identity scale, applied by default) to keep the threshold without scaling. |
| `theiler` | 0         | `recurrencerate`<br/>`determinism`<br/>`_average`<br/>`_max`<br/>`divergence`<br/>`_entropy`<br/>`trend`<br/>`laminarity`<br/>`trappingtime`<br/> `meanrecurrencetime`<br/>`nmprt` | Theiler window: number of diagonals around the LOI **excluded** from the analysis. The value `0` means that the LOI is _included_ in the analysis. Use `1` to exclude the LOI. |
| `lmin`    | 2         | `determinism`<br/>`dl_average`<br/>`_max`<br/>`divergence`<br/>`_entropy`<br/>`laminarity`<br/>`trappingtime`<br/> `meanrecurrencetime`<br/>`nmprt` | Minimum length of the recurrent structures (diagonal or vertical) considered in the analysis. |
| `border`  | 10        | `trend`  | Number of diagonals excluded from the analysis near the border of the matrix. |

## All-in-one Bundle
In case you need all of the aforementioned functions and you don't want to write 10 lines of code to compute them all (since they are so many) we also got an all-in-one function that computes all of them and returns a dictionary of the result!
```@docs
rqa
```
---

!!! note "Return values for empty histograms"
    It may be the case that for a given recurrence matrix some structures do not exist at all. For example there are recurrence matrices that have no vertical lengths (or no vertical lengths with length less than `lmin`). In such cases the behavior of our RQA pipeline is the following:

    1. Quantities that represent maximum or average values are `0.0`.
    2. Quantities that represent entropies are `NaN`.

---

See also the [`@windowed`](@ref) macro for a windowed version of [`rqa`](@def)!

## Recurrence Structures Histograms
The functions that we list in this page internally compute histograms of some recurrence structures, like e.g. the vertical lengths.
You can access these values directly with the following function:
```@docs
recurrencestructures
```
