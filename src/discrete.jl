####################################
#          Discrete                #
####################################
### REWOOOOOOOOOOOOOOOOOOORK###
# to be the similar as continuous
struct DiscreteDynamicalSystem <: DynamicalSystem
  eom!::Function
  jacobian!::Function
  u::AbstractVector
  J::AbstractArray
end
function DiscreteDynamicalSystem(u0, eom!::Function, jacob!::Function)
  x = rand(dim)
  y = eom!(x)
  if y !== x
    erst = "Equations of motion are not of the expected type for `DiscreteDynamicalSystem`!"
    erst*= "\nIt is expected that if `y = eom!(x)` then y==x, so that the"
    erst*= "\ne.o.m. are in-place and also return the result!"
    error(erst)
  end
  xx = rand(dim,dim)
  try
    jacob!(xx, x)
  catch er
    erst = "Jacobian of the equations of motion is not of the expectedt form for"
    erst*= "\n`DiscreteDynamicalSystem`. It is expected that the jacobian is in-place,"
    erst*= "\nso that jacob!(J, x) writes the result in J if x is the state-vector!"
    erst*= "\nWhile trying to access jacob!(J, x), got error: $er"
    error(erst)
  end
  DiscreteDynamicalSystem(dim, eom!, jacob!)
end
function DiscreteDynamicalSystem(dim::Int, eom!::Function)
  x = rand(dim)
  y = eom!(x)
  if y !== x
    erst = "Equations of motion are not of the expected type for `DiscreteDynamicalSystem`!"
    erst*= "\nIt is expected that if `y = eom!(x)` then y==x, so that the"
    erst*= "\ne.o.m. are in-place and also return the result!"
    error(erst)
  end
  jacob! = (J, x) -> ForwardDiff.jacobian!(J, eom!, x)
  DiscreteDynamicalSystem(dim, eom!, jacob!)
end
