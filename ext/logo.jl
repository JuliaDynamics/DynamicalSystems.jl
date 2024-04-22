using DynamicalSystems
using OrdinaryDiffEq
using CairoMakie
using DataStructures: CircularBuffer

# double pendulum dynamical system
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

    return SVector(du1, du2, du3, du4)
end

const L1 = 1.0
const L2 = 1.0
p0 = (G=10.0, L1, L2, M1 = 1.0, M2 = 1.0)
u0 = [π/3, 0, 3π/4, -2]

# Solve diffeq with constant step for smoother curves
diffeq = (alg = Vern9(), adaptive = false, dt = 0.005)
dp = CoupledODEs(doublependulum_rule, u0, p0; diffeq)

# map state to x-y coordinates in 2D space
function xycoords(state)
    θ1 = state[1]
    θ2 = state[3]
    x1 = L1 * sin(θ1)
    y1 = -L1 * cos(θ1)
    x2 = x1 + L2 * sin(θ2)
    y2 = y1 - L2 * cos(θ2)
    return SVector(x1,x2,y1,y2)
end

# Initialize observables
x1,x2,y1,y2 = xycoords(current_state(dp))
rod   = Observable([Point2f(0, 0), Point2f(x1, y1), Point2f(x2, y2)])
balls = Observable([Point2f(0, 0), Point2f(x1, y1), Point2f(x2, y2)])
tail = 3000 # length of plotted trajectory, in units of `dt`
traj = CircularBuffer{Point2f}(tail)
fill!(traj, Point2f(x2, y2))
traj = Observable(traj)

# %% Initialize figure
fig = Figure(size = (800, 800), backgroundcolor = :transparent)
ax = Axis(fig[1,1]; backgroundcolor = :transparent)
ax.limits = ((-L1-L2-0.1, L1 + L2+0.1), (-L1-L2-0.1, L1 + L2 + 0.1))
hidedecorations!(ax)
hidespines!(ax)

# Plot observables
# from https://github.com/JuliaGraphics/Luxor.jl/blob/e4fe3eb50e14fbe113dcfcf2e6a0d8ed8ad613da/src/juliagraphics.jl
julia_blue    = RGBf(0.251, 0.388, 0.847)
julia_purple  = RGBf(0.584, 0.345, 0.698)
julia_green   = RGBf(0.22, 0.596, 0.149)
julia_red     = RGBf(0.796, 0.235, 0.2)

lighter_green  = RGBf(0.376, 0.678, 0.318)
lighter_red    = RGBf(0.835, 0.388, 0.361)
lighter_blue   = RGBf(0.4, 0.51, 0.878)
lighter_purple = RGBf(0.667, 0.475, 0.757)

c = julia_blue
tailcol = [RGBAf(c.r, c.g, c.b, (i/tail)^(1.2)) for i in 1:tail]
# trajwidth = [4*(i/tail) for i in 1:tail]

trajline = lines!(ax, traj; color = tailcol, linewidth = 4)

lines!(ax, balls; linewidth = 8, color = :black)

scatter!(ax, balls; marker = :circle, strokewidth = 10,
    strokecolor = [julia_green, julia_red, julia_purple],
    color = [lighter_green, lighter_red, lighter_purple],
    markersize = 160,
)

function animstep!(dp)
    step!(dp)
    x1,x2,y1,y2 = xycoords(current_state(dp))
    rod[] = [Point2f(0, 0), Point2f(x1, y1), Point2f(x2, y2)]
    balls[] = [Point2f(0, 0), Point2f(x1, y1), Point2f(x2, y2)]
    push!(traj[], Point2f(x2, y2))
    notify(traj)
end
function animstep!(dp, t)
    t0 = current_time(dp)
    while current_time(dp) < t0 + t
        animstep!(dp)
    end
end


fig

# %% initial condition that leads to nice Julia logo, fantastic!
u0 = [π/2, +0.1, 3π/4, -3]
reinit!(dp, u0)

animstep!(dp, 6.57)

fig

# %% find initial condition that leads to nice Julia logo:
# u0 = [π/2 + 0.1, +0.03, 0.3, -3.78]
u0 = [π/2 + 0.1, +2.03, 0.3, +5.412]
# u0 = [π/2 + 0.1, +2.03, 0.3, +5.42]
reinit!(dp, u0)


dmin = 0.5
function distance_from_optimal(u)
    logocoords = SVector{2}[[-cosd(60), -sind(60)], [cosd(60), -sind(60)]]
    x1,x2,y1,y2 = xycoords(u)
    pos = SVector{2}[(x1, y1), (x2, y2)]
    d = euclidean(logocoords[1], pos[1]) + euclidean(logocoords[2], pos[2])
    return d/2
end

d = distance_from_optimal(current_state(dp))


# scatter!(ax, logocoords, markersize = 30)

animstep!(dp, 6.56)

while d > 0.01
    animstep!(dp)
    d = distance_from_optimal(current_state(dp))
    if current_time(dp) > 1000.0
        break
    end
end

d
current_time(dp)

fig

save(desktop("logo_transparent.png"), fig; px_per_unit = 4)
trajline.visible = false
save(desktop("logo_no_tail.png"), fig; px_per_unit = 4)
trajline.visible = true


# %% step, animate

animstep!(dp)

fig


# %% Search for a good initial condition
using DynamicalSystems.StateSpaceSets.Distances
function good_ic(u)
    X, t = trajectory(dp, 1000.0, u; Δt = 0.25)
    X = map(u -> xycoords(u)[3:4], X) # position of last ball
    N = length(X)
    D = fill(Inf, N, N)
    for i in 1:N
        for j in (i+1):N
            D[i, j] = Euclidean()(X[i], X[j])
        end
    end

    dmin, i = findmin(D)
    t1, t2 = t[i[1]], t[i[2]]

    # D = pairwise(Euclidean(), vec(X))
end

reinit!(dp, u)
animstep!(dp, 304.5)
fig
animstep!(dp, 690.75 - 304.5)
fig
