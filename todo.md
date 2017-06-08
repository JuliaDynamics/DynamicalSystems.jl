1. If I am solving ODEs with stiff problems, and using e.g. the Rosenbrock23
  then you can use the integrator interface of the ODE package to get the QR decompoisition
  directly from the solver (because it uses the QR to do it). Implement a method
  for lyapunov for stiff problems.

1. Add much more systems in the `famous_systems.jl`

1. add progress meter: make it possible for user to choose between Juno progress
  and REPL progress.
