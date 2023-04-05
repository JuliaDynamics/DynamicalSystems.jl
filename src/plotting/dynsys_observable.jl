# Unlike Agents.jl, here we don't need a wrapper struct.
# The `DynamicalSystem` itself contains all information, and hence
# it can be used directly as an observable
"""
    step!(ds::Observable{<:DynamicalSystem}, args...)

Call `step!` on the system of the observable and then update the observable.
"""
function SciMLBase.step!(ds::Observable{<:DynamicalSystem}, args...)
    step!(ds[], args...)
    update!(ds)
    return
end
