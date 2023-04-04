using DynamicalSystems
using OrdinaryDiffEq
using GLMakie
using DataStructures: CircularBuffer

# double pendulum system
@inbounds function doublependulum_rule(u, p, t)
    G, L1, L2, M1, M2 = p
    du1 = u[2]
    φ = u[3] - u[1]
    Δ = (M1 + M2) - M2*cos(φ)*cos(φ)
    du2 = (M2*L1*u[2]*u[2]*sin(φ)*cos(φ) +
               M2*G*sin(u[3])*cos(φ) +
               M2*L2*u[4]*u[4]*sin(φ) -
               (M1 + M2)*G*sin(u[1]))/(L1*Δ)

    du3 = u[4]
    du4 = (-M2*L2*u[4]*u[4]*sin(φ)*cos(φ) +
               (M1 + M2)*G*sin(u[1])*cos(φ) -
               (M1 + M2)*L1*u[2]*u[2]*sin(φ) -
               (M1 + M2)*G*sin(u[3]))/(L2*Δ)

    return SVector{4}(du1, du2, du3, du4)
end

const L1 = 1.0
const L2 = 1.0
p0 = (G=10.0, L1, L2, M1 = 1.0, M2 = 1.0)
u0 = [π/3, 0, 3π/4, -2]
# Solve diffeq with constant step for smoother curves
diffeq = (alg = Vern9(), adaptive = false, dt = 0.01, abstol = 1e-9, reltol = 1e-9)
dp = CoupledODEs(doublependulum_rule, u0, p0; diffeq)









# Extract xy coordinates from state
function xycoords(state)
    θ1 = state[1]
    θ2 = state[3]
    x1 = L1 * sin(θ1)
    y1 = -L1 * cos(θ1)
    x2 = x1 + L2 * sin(θ2)
    y2 = y1 - L2 * cos(θ2)
    return x1,x2,y1,y2
end

function progress_for_one_step!(ds)
    step!(ds)
    u = current_state(ds)
    return xycoords(u)
end


# Initialize observables
x1,x2,y1,y2 = xycoords(u0)
rod   = Observable([Point2f(0, 0), Point2f(x1, y1), Point2f(x2, y2)])
balls = Observable([Point2f(x1, y1), Point2f(x2, y2)])
tail = 300 # length of plotted trajectory, in units of `dt`
traj = CircularBuffer{Point2f}(tail)
fill!(traj, Point2f(x2, y2))
traj = Observable(traj)

# %% Initialize figure
fig = Figure(resolution = (1200, 800)); display(fig)
ax = Axis(fig[1,1])

# Plot observables
c = to_color(:purple)
tailcol = [RGBAf(c.r, c.g, c.b, (i/tail)^2) for i in 1:tail]
lines!(ax, traj; linewidth = 3, color = tailcol)

lines!(ax, rod; linewidth = 4, color = :purple)

scatter!(ax, balls; marker = :circle, strokewidth = 2,
    strokecolor = :purple,
    color = :black, markersize = [8, 12]
)