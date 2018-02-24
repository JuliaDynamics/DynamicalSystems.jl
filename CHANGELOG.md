# v0.8.0

DynamicalSystems.jl was separated into 3 packages, to prepare for the arrival
of more packages into the ecosystem.

The packages now are:
1. This one.
2. DynamicalSystemsBase.jl
3. ChaosTools.jl

Now DynamicalSystems.jl is simply a package coordinator that loads and exports
all packages of the ecosystem, as well as documentation host.

From now on, the changes will be logged at each individual package's CHANGELOG.md
file for clarity.

## Bugfixes and Enchancements
* `geneneralized_dim` allows one to pass the box sizes to be used.
  *  Also improved the default values for `estimate_boxsizes`.

# v0.7.0
## BREAKING
* Changed the definition of all pre-defined systems. Now the documentation
  also suggests users to use Functors for system definition.
* The only mutable system types are `DiscreteDS` and `DiscereteDS1D`. The others
  are immutable.
* Completely removed the possibility to return the convergence timeseries for
  functions like `lyapunov`. This tripled the necessary amount of code and
  makes maintaining and expanding much harder. Since our source code is very
  clear and concise, if users really want those convergence timeseries they can
  create their own modified version of the functions.
* Removed the "name" field (and keyword) from all types. It was an unneccessary
  complication since the eom function/type name is displayed anyway.
* `Reconstruction` and `reconstruct` immediately create the data upon call.
  Now `Reconstruction` is almost identical to a dataset but it has one more
  type-parameter `τ`.
  * This made the source code much cleaner and more concise!
  * The constructor `Reconstruction` is the function to be used. The function
  `reconstruct` is only used internally.
  * Using generated functions creating a `Reconstruction` is crazy fast!
* `AbstractDataset` now has 2 parameters, `{D, T}` (for dimension and type).
  * All-round improvement of method definitions. Everything is now done through
  the `Abstract` type, allowing more flexibility.

## New Additions
* New visualization routines that are compiled "on-demand" when users type `using PyPlot`:
    * `phasespace` that plots phasespaces of discrete maps.
    * New documentation page about visualizations.
* Now users can initialize an empty `Dataset` to append datasets later: `Dataset{D,T}()`.

## Bugfixes and Enchancements
* Created a dedicated method that integrates the variational equations
  of motion for a system. This allows all functions to simply call
  this one instead of each defining a new integrator. This will make DynamicalSystems
  more extendable as now it is much easier to add e.g. new methods for
  chaos detection.
* Tripled the speed of `non0hist`!!! Note to self: DO NOT use broadcasting on basic operations like `+. *, /` etc. when using `SVector`s!
* Much better documentation for system definition.

# v0.6.1
## New Additions
* Generalized Alignment Index method: `gali`
* GALI for continuous systems.
* Documentation & examples for `gali`.

# v0.6.0:
## BREAKING
* `timeseries` function was renamed to `trajectory`.
    * Now the documentation is clear: `trajectory` means a set of one-dimensional
      timeseries (or a *multi-dimensional* timeseries), while `timeseries`
      always means a one-dimensional timeseries
* The `jacob` field of `ContinuousDS` has been changed to `jacob!`. Now the function
  operates in-place `jacob!(J, state)` instead of returning an `SMatrix`. This
  is similar to how `BigDiscreteDS` operate.

## New Additions
* Hypermassive Docs restructuring + upgrade. Now all documentation strings have dedicated
  sections for calling the function, describing the algorithims, etc. etc.
* Automatic Jacobian for continuous systems!
* `BigDiscreteDS` type: a new system type for discrete systems with very high dimensionality.
* Added automated Jacobian computation for `BigDiscreteDS`.
* Henon Heiles system in famous_systems.
* Nonlinearly coupled standard maps in famous systems.
* Lorenz96 in famous_systems.


## Bugfixes and Enhancements
* **Major** improvements in the speed of all `lyapunovs` functions. Especially the `DiscreteDS` one is now 3 times faster due to the `qr` provided by StaticArrays.jl.
* Change the pretty-printing of all systems, by adding an additional
  `name` field everywhere that can be used for convenience.
* Bugfix in `lyapunovs` where `u_modified!(integ, true)` was not used.
* Bugfix of the printing of `Dataset` on Juno where the set was always printed
  on the console.
* Major bugfix in `lyapunov` for `DiscreteDS` where the test state was initialized
  before the transient iterations!!!
* Many minor bugfixes and many corrections to the documentation examples.
* `Dataset` now always converts the input into vector of `SVectors` for better
  interaction with NearestNeighbors.jl
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
