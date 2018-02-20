infomsg = """
~~~ DynamicalSystems.jl update message ~~~

Our library has been reworked from the ground-up, and is now better than ever!

Most top level functions have slightly altered signature, read the docstrings!

All `DynamicalSystem` instances are immutable, and eight (8) system
combinations are now possible:

* Continuous or Discrete.
* In-place or out-of-place (large versus small systems).
* Auto-differentiated or not (for the Jacobian function).

Please see the documentation strings of the system types you are using!

Also, be sure to visit the updated documentation here:
https://juliadynamics.github.io/DynamicalSystems.jl/latest/

"""

info(infomsg)
