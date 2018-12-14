# Quantification & Analysis functions
A [`recurrencematrix`](@ref) can be analyzed in several ways to yield information about the dynamics of the timeseries. All these various quantities and functions listed in this section are also listed in the wikipedia page for [Recurrence quantification analysis](https://en.wikipedia.org/wiki/Recurrence_quantification_analysis). We also suggest the review articles listed in our documentation strings, namely:

1. N. Marwan *et al.*, "Recurrence plots for the analysis of complex systems", *Phys. Reports 438*(5-6), 237-329 (2007).
2. N. Marwan & C.L. Webber, "Mathematical and computational foundations of recurrence quantifications", in: Webber, C.L. & N. Marwan (eds.), *Recurrence Quantification Analysis. Theory and Best Practices*, Springer, pp. 3-43 (2015).

```@docs
recurrencerate
determinism
avgdiag
maxdiag
divergence
trend
laminarity
trappingtime
maxvert
RecurrenceAnalysis.entropy
```
Since most of the above functions are fined tuned with keyword arguments, here is a table summarizing them that could be of use


| Argument  | Default   | Functions | Description |
| --------  | --------  | --------- | -----------
| `scale`   | 1         | `recurrencematrix`<br/>`crossrecurrencematrix`<br/>`jointrecurrencematrix` | Function or fixed number to scale the threshold or radius that is used to identify recurrences. Use `maximum` if the threshold is to be taken as a fraction of the maximum distance, `mean` if it is a fraction of the mean distance, etc., and `1` (identity scale, applied by default) to keep the threshold without scaling. |
| `theiler` | 0         | `recurrencerate`<br/>`determinism`<br/>`avgdiag`<br/>`maxdiag`<br/>`divergence`<br/>`entropy`<br/>`trend`<br/>`laminarity`<br/>`trappingtime`<br/>`maxvert` | 'Theiler' window: number of diagonals around the LOI excluded from the analysis. |
| `lmin`    | 2         | `determinism`<br/>`avgdiag`<br/>`maxdiag`<br/>`divergence`<br/>`entropy`<br/>`laminarity`<br/>`trappingtime`<br/>`maxvert` | Minimum length of the recurrent structures (diagonal or vertical) considered in the analysis. |
| `border`  | 10        | `trend`  | Number of diagonals excluded from the analysis near the border of the matrix. |

## All-in-one Bundle
In case you need all of the aforementioned functions and you don't want to write 10 lines of code to compute them all (since they are so many) we also got an all-in-one function that computes all of them and returns a dictionary of the result!
```@docs
rqa
```
