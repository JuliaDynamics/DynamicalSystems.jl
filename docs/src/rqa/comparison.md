# Comparison of RQA software
There are several software packages that can be used to perform Recurrence Quantification Analysis (RQA), in different programming languages. Many of them can be found at <http://www.recurrence-plot.tk/programmes.php>. This document compares the features of the package [**RecurrenceAnalysis**](https://www.github.com/JuliaDynamics/RecurrenceAnalysis.jl) (version 1.0-alpha) for Julia with other packages that have been developed for popular cross-platform programming languages, namely:

* [**crqa**](https://cran.r-project.org/package=crqa) version 1.0.7 for R (Coco & Dale, 2014).
* [**CRP toolbox**](http://tocsy.pik-potsdam.de/CRPtoolbox/) version 5.22 for Matlab (Marwan et al., 2007)
* [**pyunicorn**](http://www.pik-potsdam.de/~donges/pyunicorn/) version 0.5.1 for Python (Donges et al., 2015).

The state of the art of RQA for Python is the package [PyRQA](https://pypi.org/project/PyRQA/), based on parallel computing and capable of using efficient hardware resources (Rawald, Sips & Marwan, 2017), but it has not used for this comparison.

## Definition of Recurrence Plots

The first step of RQA is the calculation of a *recurrence plot* (RP) &mdash; or a *cross-recurrence plot* (CRP) or a *joint recurrence plot* (JRP) when two time series are analyzed. Those plots are numerically defined as Boolean (and usually sparse) matrices that can be called (*cross-*, *joint-*) *recurrence matrices*. There is not a unique algorithm for the calculation of those matrices, and each software package provides different options that are described next.

### Normalization of the time series

The original time series can be normalized prior to their analysis. Some packages provide this operation as part of the RQA:

* **RecurrenceAnalysis**: No normalization method included.
* **crqa**: Normalization in the in the [0&ndash;1] interval or as z-score (zero mean and unit standard deviation of the sample).
* **CRP Toolbox**: Z-score normalization.
* **pyunicorn**: Z-score normalization.

### Metric of the distances

The distances between the points of the embedded time series can be calculated by different types of norms. The metrics for the norms supported by the different packages are:

* **RecurrenceAnalysis**: Euclidean, Infinity or Manhattan distance.
* **crqa**: Euclidean distance.
* **CRP Toolbox**: Euclidean, Infinity, Manhattan, or normalized Euclidean (angle).
* **pyunicorn**: Euclidean, Inifinity or Manhattan distance.

### Scale of the threshold

The distances between the points of the embedded time series can be scaled before comparing them with the threshold to evaluate recurrences. The scales supported by the different packages are:

* **RecurrenceAnalysis**: Any fixed, absolute scale, or a function (e.g. maximum, mean, etc.) of the observed distances.
* **crqa**: No scaling, or scaling by the maximum or mean of the observed distances.
* **CRP Toolbox**: No scaling.
* **pyunicorn**: No scaling.

### Fixed recurrence rate:

Instead of using a threshold for the value of the distances between points of the embedded time series, it is possible to define the recurrence plots by setting a fixed rate or number of recurrent points. This can be done globally or locally (same number or rate of recurrences for every single point of the time series).

* **RecurrenceAnalysis**: Fixed recurrence rate supported.
* **crqa**: No direct method to fix the recurrence rate provided.
* **CRP Toolbox**: Fixed recurrence rate (global) or fixed amount of nearest neighbours (local) supported.
* **pyunicorn**: Fixed global or local recurrence rate supported.

### Visualization

**crqa** and the **CRP Toolbox** have their own methods to visualize the RPs. The recurrence matrices created by **pyunicorn** can be visualized with the function `matshow` from the **matplotlib** package in Python. The matrices created by **RecurrenceAnalysis** are shown in the console as text-based representations of sparse matrices, and there is a function to convert them into numeric matrices plottable by other packages.

## Classical RQA

The classical RQA consists in the calculation of the recurrence rate (RR), and some parameters related to the distribution of diagonal recurrent structures that are longer than a minimum length, which is usually 2 points. Table 1 shows the parameters that are directly calculated by each software package.

Table 1. Parameters of classical RQA

| Parameter                     | RecurrenceAnalysis | crqa | CRP Toolbox | pyunicorn |
| ---------                     |:------------------:|:----:|:-----------:|:---------:|
| Recurrence rate (RR)          |   X                |  X   |   X         |  X        |
| Determinism (DET)             |   X                |  X   |   X         |  X        |
| Average diagonal (L)          |   X                |  X   |   X         |  X        |
| Maximum diagonal (Lmax)       |   X                |  X   |   X         |  X        |
| Divergence (1/Lmax)           |   X                |      |             |           |
| Entropy (ENTR)                |   X                |  X   |   X         |  X        |
| Number of diagonal lines (NL) |   X[*]             |  X   |             |           |
| Relative ENTR (ENTR/NL)       |                    |  X   |             |           |
| Trend (TND)                   |   X                |      |             |           |

[*]: #table1_note
<a name="table2_note">(*)</a> The number of lines is not directly given as output of any function by RecurrenceAnalysis, but can be easily obtained from the histograms of recurrence structures, which can be calculated.

Other differences between the packages are commented on next.

### Scale of ratios (RR, DET)

**crqa** defines RR and DET as a percentage, whereas the other packages present them as ratios in the [0&ndash;1] interval.

### Exclusion of the LOI and Theiler window

All the points of the line of identity (LOI), i.e. the main diagonal of RPs and JRPs, are always recurrent, and when a finite, non-zero threshhold is used to evaluate if two points are recurrent, so are most of the points near to the LOI. For that reason the LOI and the diagonals within a certain distance from the LOI &mdash; the so-called Theiler window &mdash; are often excluded from the analysis. The different packages give different options for this:

* **RecurrenceAnalysis**: LOI excluded by default in RPs and JRPs. Theiler window defined by the user.
* **crqa**: No diagonal excluded by default. Theiler window defined by the user.
* **CRP Toolbox**: LOI excluded by default in RPs and JRPs. Theiler window defined by the user.
* **pyunicorn**: No diagonal excluded. No support for Theiler window.

### Adjustment of RR when diagonals are excluded

The number of points in (cross, joint-) recurrence matrices is greater than the number of possible recurrent points when the LOI is excluded or a Theiler window is applied. This affects the calculation of the recurrence rate, which is handled in a different way by each package:

* **RecurrenceAnalysis**: The RR is adjusted taking into account the actual number of potential recurrent points for RPs and JRPs. In CRPs the number of recurrent points is divided by the full size of the matrix.
* **crqa**: The number of recurrent points is divided by the full size of the matrix.
* **CRP Toolbox**: The RR is adjusted taking into account the actual number of potential recurrent points.
* **pyunicorn**: There is no possible exclusion of recurrent points, so this issue does not apply.

## Extended RQA

The advanced RQA consists in the calculation of parameters that describe the distribution of vertical recurrent structures. This can be further extended by the analysis of the periods between vertical recurrent structures (recurrence time based measures). The parameters calculated by the different packages are shown in table 2.


Table 2. Parameters of extended RQA

| Parameter                   | RecurrenceAnalysis | crqa | CRP Toolbox | pyunicorn |
| ---------                   |:------------------:|:----:|:-----------:|:---------:|
| Laminarity (LAM)            |   X                |  X   |   X         | X         |
| Trapping time (TT)          |   X                |  X   |   X         | X         |
| Maximum vertical (Vmax)     |   X                |      |   X         | X         |
| Average recurrence time     |   X                |      |   X         | X         |
| Entropy of recurrence times |   X                |      |   X         | X         |
| Maximum recurrence time     |   X                |      |             | X         |
| Entropy of recurrence times |   X                |      |             | X         |

### Estimation of recurrence times

There are different possibilites of characterizing and estimating recurrence times (cf. Marwan and Webber, 2015), and each package supports different definitions (excluding crqa, which does not have parameters based on recurrence times):

* **RecurrenceAnalyis**: recurrence times are estimated as the distance between the central points of consecutive vertical structures of the recurrence plot.
* **CRP Toolbox**: the two types of recurrence times defined by Gao and Cai (2000) are calculated: first type recurrence times are calculated as the distance between each consecutive pair of recurrent points in a column of the matrix; second type recurrence times are estimated as the distance between the start of consecutive vertical structures.
* **pyunicorn**: recurrence times are estimated as the length of "white lines" between vertical structures of the plot.

### Scale of ratios LAM

**crqa** defines LAM as a percentage, whereas the other packages present it as a ratio in the [0&ndash;1] interval.

### Exclusion of the LOI and Theiler window

The points in the LOI and a Theiler window around the LOI can be excluded in extended RQA for the same reasons as in the classical RQA, although the effect is slightly different: while diagonal structures are fully ruled out, long vertical structures can be "split" into pairs of shorter structures. This is handled in a different way by each package.

* **RecurrenceAnalysis**: Different Theiler windows can be set for the calculation of classical and extended RQA.
* **crqa**: The same Theiler window is applied for the calculation of all parameters.
* **CRP Toolbox**: The same Theiler window is applied for the calculation of all parameters.
* **pyunicorn**: No support for Theiler windows.

## Computation times

Figure 1 shows a benchmark of the times taken by the four packages that are compared to compute the RQA of a time series (the first coordinate of a Rossler system with parameters *a*=0.25, *b*=0.25, and *c*=4, starting at (0,0,0) excluding the first 1000 points), taking samples of increasing length from *N*=250 to *N*=3000. The *x*-coordinate of the trajectories was embedded in three dimensions with a fixed delay of 6 samples; the recurrence matrix of the embedded time series was calculated using a fixed threshold ε=1.2, and all RQA parameters were calculated. 

The computations have been done in Linux and Windows machines, although the CRP Toolbox (Matlab) has only been tested in the Windows machine, and pyunicorn (Python) in the Linux machine. The Linux machine was an Intel(R) Core(TM) i3-4100M CPU @ 2.50GHz, 8GiB SODIMM DDR3 1600 MHz, running Ubuntu 16.04.5-LTS. The Windows machine was an Intel(R) Core(TM) i3-3240 CPU @ 3.40GHz CPU, 4GB DDR3 RAM 1333 MHz, running Windows 10 Enterprise. The software versions were:

* **RecurrenceAnalysis** 0.3.0 running on Julia 1.0.0
* **crqa** 1.0.7 running on R 3.5.1
* **CRP toolbox** 5.22 running on Matlab R2013b
* **pyunicorn**0.5.1 running on Python 2.7.12.

The first panel shows the times (milliseconds) in a natural scale, and the second in a logarithmic scale to see better the differences between Julia and Python. Coco and Dale (2014) already demonstrated the better performance of crqa in R with respect to the CRP Toolbox in Matlab, with a similar benchmark study using binary time series, although the differences shown here with the Rossler series is more striking.

Both packages are outperformed by RecurrenceAnalysis and pyunicorn in Julia and Python, respectively, which show a similar speed in a natural scale, although pyunicorn is significantly faster, as clearly seen in the logarithmic scale: it is is about 5 times faster than RecurrenceAnalysis, between 20 and 50 times faster than crqa, and hundreds of times faster than the CRP Toolbox.

[[/imgs/comparison_rqa_noscale.png|Speed comparison (no scale)]]
[[/imgs/comparison_rqa_logscale.png|Speed comparison (logarithmic scale)]]

Figure 1. Median computation speeds of RQA by the different packages.

An advantage of RecurrenceAnalysis with respect to all the other packages is that is is entirely written in Julia (without C or C++ code), so it is easy to inspect, extend and improve. Future versions will include features that are currently missing, and possibly improve the performance of the code to gain further computation speed.

## References

Coco, M., Dale, R. (2014). Cross-recurrence quantification analysis of categorical and continuous time series: an R package. *Frontiers in Psychology* 5, 510. DOI: 10.3389/fpsyg.2014.00510

Donges, J.F., Heitzig, J., Beronov, B., Wiedermann, M., Runge, J., Feng, Q.-Y., Tupikina, L., Stolbova, V., Donner, R.V., Marwan, N., Dijkstra, H.A., Kurths, J. (2015). Unified functional network and nonlinear time series analysis for complex systems science: The pyunicorn package, *Chaos* 25, 113101, DOI: 10.1063/1.4934554

Gao, J. Cai, H. (2000). On the structures and quantification of recurrence plots, *Physics Letters A*, 270, 75–87. DOI:10.1016/S0375-9601(00)00304-2

Marwan, N., Romano, M.C., Thiel, M., Kurths, J. (2007). Recurrence Plots for the Analysis of Complex Systems, *Physics Reports*, 438(5-6), 237-329. DOI: 10.1016/j.physrep.2006.11.001

Marwan, N., Webber, C.L. (2015) "Mathematical and computational foundations of recurrence quantifications", in: Webber, C.L. & Marwan, N. (eds.), *Recurrence Quantification Analysis. Theory and Best Practices*, Springer, pp. 3-43.

Rawald, T., Sips, M., Marwan, N. (2017). PyRQA—Conducting recurrence quantification analysis on very long time series efficiently. *Computers & Geosciences*, 104(C), 101-108. DOI: 10.1016/j.cageo.2016.11.016
