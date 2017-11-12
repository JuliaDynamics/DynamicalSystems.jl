using DynamicalSystems

ti = time()

# Mathematics:
include("math_tests.jl")
# System Evolution:
include(joinpath("systems", "discrete_systems.jl"))
include(joinpath("systems", "continuous_systems.jl"))
# lyapunov Exponents:
include(joinpath("lyapunovs", "discrete_lyapunov.jl"))
include(joinpath("lyapunovs", "continuous_lyapunov.jl"))
# Entropies (and attractor dimensions)
include("entropy_dimension.jl")
# Nonlinear Timeseries Analysis:
include("nlts_tests.jl")
# Periodicity:
include("periodicity_tests.jl")
# Chaos Detection:
include("chaos_detection_tests.jl")

ti = time() - ti
println("\nTest took total time of:")
println(round(ti, 3), " seconds or ", round(ti/60, 3), " minutes")
