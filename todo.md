1. add tests for solving with many different kwargs using ODE

1. If I am solving ODEs with stiff problems, and using e.g. the Rosenbrock23
  then you can use the integrator interface of the ODE package to get the QR
  decompoisition
  directly from the solver (because it uses the QR to do it). Implement a method
  for lyapunov for stiff problems.

1. add progress meter: make it possible for user to choose between Juno progress
  and REPL progress.

1. Add much more systems in the `famous_systems.jl`
  1. standard map
  1. Van der Pol oscillator
  1. Arlonds Cat map
  1. Bakers map
  1. Henon map
  1. Rossler (continuous)
  1. Heiles
  1. Lotka Voltera
  1. Kaplan Yorke system

1. create external documentation repo
