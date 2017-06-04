# v0.0.1: Initial Commit
The fundamental types (and their structure) of this package are decided:

* `DiscreteDS` : **immutable** type representing a discrete dynamical system. The fields are the system state (which is a **SVector**) `state`, the equations of motion `eom` and the Jacobian function `jacob`. The e.o.m. is a function that returns
  an SVector given some vector `eom(u) -> un::SVector`. Therefore it is *not* an in-place function!. The Jacobian function does the same as well: `jacob(u) -> J::SMatrix`
* `ContinuousDS` :

### Comments
Making a `DiscreteDynamicalSystem` with subfield `DiscreteDS` was a really bad idea!
