# This is a very light struct that contains the trajectory end state
# and the trajectory tail. Note that for convenience it always contains the state
# as a vector of observables, each observable containing each of the
# parallel states of the dynamical system
struct DynamicalSystemObservable
    pds::ParallelDynamicalSystem # reference to the dynamical system
    state_observables::Vector{Observable}
    tail_observables::Vector{Observable}
end


"""
    step!(ds::Observable{<:DynamicalSystem}, args...)

Call `step!` on the system of the observable and then update the observable.
"""
function SciMLBase.step!(ds::Observable{<:DynamicalSystem}, args...)
    step!(ds[], args...)
    update!(ds)
    return
end
