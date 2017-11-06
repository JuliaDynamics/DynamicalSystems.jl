* Added the Generalized Alignment Index method: `gali`
* Added GALI for continuous systems.
* Added the Henon Helies system in famous_systems. TODO is still to add
  callbacks that conserve energy.

# v0.6.0:
* Many minor bugfixes and many corrections to the documentation examples.
* `Dataset` now always converts the input into vector of `SVectors` for better
  interaction with NearestNeighbors.jl
* **[BREAKING]** : `timeseries` function was renamed to `trajectory`.
    * Now the documentation is clear: `trajectory` means a set of one-dimensional
      timeseries (or a *multi-dimensional* timeseries), while `timeseries` always means a one-dimensional timeseries
* Removed function `perform_non0hist`, as it had no reason to exist.


# v0.5.0:
* Bugfix where `periodicorbits` was not considering `disttol` keyword.
* Removed `rescale` argument from `lyapunov`. This is technically breaking.
* Bugfixes in `lyapunov`.
* Improved default parameters in `lyapunov` as well as type-stability.

# v0.4.5:
* Major bugfix in the `lyapunov` function where the time averaging was actually
  not happening correctly for discrete systems, resulting in slight inaccuracy.
* Added option to choose whether you want the convergence timeseries for `lyapunov`,
  for both discrete and continuous systems.
* Added option to supply your own functions to rescale and initialize the test
  state in `lyapunov` for both continuous and discrete systems.
* Bugfix related to `estimate_delay`.

# v0.4.0:
* Increased speed of Standard Map system (by using `while` instead of `mod`).
* Added a huuuge method due to Schmelcher & Diakonos which finds stable and
  unstable fixed points of any order for discrete maps!
* Added `append!` method for `Dataset`.

# v0.3.0:
* Added method to estimate optimal delay time `estimate_delay`.
* Added abstract `AbstractDataset` supertype.
* Completely changed the way `Reconstruction` is done.
* Added more methods, interaction and smart indexing related to the (new)
  `Reconstruction` object. See the new docstring of `reconstruct` for usage.
* Added a new method that enables the calculation of the maximum Lyapunov
  exponent from a numeric timeseries (like an experimental recording).
  This method (called `numericallyapunov`) operates on a reconstruction, **using
  4 different algorithms**. See the docstring of `numericallyapunov` for
  more info!
* Added `double_pendulum` to famous systems.

# v0.2.0: Major Improvements and new features
+ Improved algorithm for deducing best boxsizes (using `min_pairwise_distance`)
+ All numerical data is now structured and expected to be in the `Dataset` Type,
  which by itself defines a new interface (identical to matrix).
+ Major speed improvements in lyapunov function calculation.
+ Major speed improvements in histogram computation.
+ A new qr method for small matrices, faster than `Base.qr` (see `qr_sq`).
+ Reconstruction interface.
+ New and better docs.
+ Better integration of `SVectors` in the package all around.
+ Better default keywords for attractor dimension and lyapunov exponents.
+ Equations of motion on all famous systems.

# v0.1.0: First release
Normal changelog will be kept from this point onwards.
