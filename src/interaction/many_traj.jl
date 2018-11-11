using DynamicalSystems, PyPlot

# Have an "Events" function that takes a state and returns true
# or false and if true a horizontal line is plotted.

ds = Systems.henonheiles()
D = dimension(ds)

potential(x, y) = 0.5(x^2 + y^2) + (x^2*y - (y^3)/3)
energy(x,y,px,py) = 0.5(px^2 + py^2) + potential(x,y)
const E = energy(get_state(ds)...)
function complete(y, py, x)
    V = potential(x, y)
    Ky = 0.5*(py^2)
    Ky + V ≥ E && error("Point has more energy!")
    px = sqrt(2(E - V - Ky))
    return ic = [x, y, px, -py]
end

# %%
tr = trajectory(ds, 500.0, dt = 0.1, u0 = get_state(ds))
tvec = 0:0.1:500.0

function event(t, state)
    abs(state[1]) < 0.01 && state[3] > 0
end

cols = columns(tr)
mini, maxi = minmaxima(tr)

# find event lines
evns = Int[]
for j in eachindex(tvec)
    event(tvec[j], tr[j]) && push!(evns, j)
end

# plot everything
fig, axs = subplots(D, 1, sharex = true)
for i in 1:D
    sca(axs[i])
    PyPlot.plot(tvec,cols[i])
    for v in evns
        PyPlot.plot((tvec[v], tvec[v]), (mini[i], maxi[i]), color = (0,0,0,0.25))
    end
    if i == 1
        PyPlot.plot(tvec, zeros(length(tvec)))
    end
end
