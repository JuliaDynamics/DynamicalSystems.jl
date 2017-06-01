using DynamicalSystems, SymPy
#######################################################################################
#                                     Constructors                                    #
#######################################################################################
"""
    DiscreteDynamicalSystem <: DynamicalSystem
# Fields:
* `u::AbstractVector` : Current state-vector of the system (initialized as the
  initial conditions).
* `J::AbstractMatrix` : Jacobian matrix of the e.o.m. at the current `u`.
* `eom!::Function` : The function representing the equations of motion in an
  **in-place form with only one argument**: `eom!(x)` with `x` the current state.
  The function updates `x` in-place with the next state.
* `jacobian!::Function` : An in-place function `jacobian!(J, u)` that given a
  state-vector calculates the corresponding Jacobian of the e.o.m.
  and writes it in-place for `J`.
# Constructors:
1. `DiscreteDynamicalSystem(u0, eom!::Function, jacobian!::Function)`
  creates a system with user-provided functions for the equations of motion and the
  Jacobian of them (most efficient). *Always* use this if you know the Jacobian of
  your system.
2. `DiscreteDynamicalSystem(u0, eom!::Function)`
  uses the package `ForwardDiff` for automatic (numeric) forward
  differentiation to calculate the `jacobian!` function.
"""
struct DiscreteDynamicalSystem <: DynamicalSystem
  u::AbstractVector
  J::AbstractArray
  eom!::Function
  jacobian!::Function
end

function DiscreteDynamicalSystem(u0, eom!::Function, jacob!::Function)
  J = similar(u0, (length(u0), length(u0)))
  jacob!(J, u0)
  DiscreteDynamicalSystem(u0, J, eom!, jacob!)
end

function DiscreteDynamicalSystem(u0, eom!::Function)
  J = similar(u0, (length(u0), length(u0)))
  L = length(u0)
  fakef = (u) -> (un = copy(u); eom!(un); return un)
  jacob! = (J, x) -> ForwardDiff.jacobian!(J, fakef, x)
  DiscreteDynamicalSystem(u0, J, eom!, jacob!)
end

# Use SymPy for the method:
function DiscreteDynamicalSystem(u0, f::Array{SymPy.Sym,1}, variables::Array{SymPy.Sym,1})
