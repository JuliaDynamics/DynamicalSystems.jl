ti = time()
include("discrete_systems.jl")
include("discrete_lyapunov.jl")
include("continuous_systems.jl")
ti = time() - ti
println("Test took total time of:")
println(round(ti, 3), " seconds or ", round(ti/60, 3), " minutes")
