using DynamicalSystems
using OrdinaryDiffEq
using CairoMakie
using DataStructures: CircularBuffer
desktop() = joinpath(homedir(), "Desktop")
desktop(args...) = joinpath(desktop(), args...)

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
dt = 0.005
diffeq = (alg = Vern9(), adaptive = false, dt)
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
tail = 2000 # length of plotted trajectory, in units of `dt`
traj = CircularBuffer{Point2f}(tail)
fill!(traj, Point2f(x2, y2))
traj = Observable(traj)

# %% Initialize figure
fig = Figure(size = (800, 800), backgroundcolor = :transparent)
ax = Axis(fig[1,1]; backgroundcolor = :transparent, aspect = DataAspect() )
# this is maximum possible limits:
# ax.limits = ((-L1-L2-0.1, L1 + L2+0.1), (-L1-L2-0.1, L1 + L2 + 0.1))
# this is reduced size in height:
ax.limits = ((-L1-L2-0.1, L1 + L2+0.1), (-L1-L2-0.1, (L1 + L2 + 0.1)/2))
ylims!(-0.8660254037844386/2 - 2.1, -0.8660254037844386/2 + 2.1) # center y of axis for equal aspect ratio
# ax.aspect = 1
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

lighter_green  = julia_green
lighter_red    = julia_red
lighter_blue   = julia_blue
lighter_purple = julia_purple

# color of the trajectory (fading out)
c = julia_blue
tailcoltransp = [RGBAf(c.r, c.g, c.b, (i/tail)^(1.2)) for i in 1:tail]
tailcol = [RGBf(c.r, c.g, c.b) for i in 1:tail]
trajline = lines!(ax, traj; color = tailcoltransp, linewidth = 4)


# rods that connect the pendulum
lines!(ax, balls; linewidth = 12, color = :black)

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


# %% find initial condition that leads to nice Julia logo:
using DynamicalSystems.StateSpaceSets.Distances
u0 = [π/2 + 0.1, +2.03, 0.3, +5.412]
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

animstep!(dp, 6.56)

while d > 0.01
    animstep!(dp)
    d = distance_from_optimal(current_state(dp))
    if current_time(dp) > 1000.0
        break
    end
end

d
tf = current_time(dp)

# reported final time of evolution from given initial condition:
# 168.38499999991936

# okay, save high quality version:
ax.backgroundcolor = :transparent
CairoMakie.save(desktop("juliadynamics_logo.png"), fig; px_per_unit = 4)
# and one without tail
trajline.visible = false
save(desktop("juliadynamics_logo_no_tail.png"), fig; px_per_unit = 4)
trajline.visible = true
# and one more with white background
ax.backgroundcolor = :white
CairoMakie.save(desktop("juliadynamics_logo_white.png"), fig; px_per_unit = 4)
# and a dark background
ax.backgroundcolor = "#1e1e20"
CairoMakie.save(desktop("juliadynamics_logo_dark.png"), fig; px_per_unit = 4)
ax.backgroundcolor = :transparent


# save(desktop("logo_transparent.png"), fig; px_per_unit = 4)
# trajline.visible = false
# save(desktop("logo_no_tail.png"), fig; px_per_unit = 4)
# trajline.visible = true

fig

# %% perform video animation animate from some t start to tf
ax.backgroundcolor = :transparent
reinit!(dp, u0)
resize!(fig, 800, 800)

span = tail*dt
ts = tf - span
# ts += 20dt # for whatever reason we have to do this correction
animstep!(dp, ts) # initial state

# fig.backgroundcolor = :
dtrecord = dt*5
frames = range(0, span; length = Int(span ÷ dtrecord))
@show length(frames)

record(fig, desktop("juliadynamics_logo_anim.mp4"), frames; framerate = 30) do i # i = frame number
    animstep!(dp, dtrecord)
    fig
end

fig
