using DynamicalSystems
using OrdinaryDiffEqVerner
using CairoMakie
using DataStructures: CircularBuffer
desktop() = joinpath(homedir(), "Desktop")
desktop(args...) = joinpath(desktop(), args...)
desktop() = @__DIR__

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
lima = 0.1
ax.limits = ((-L1-L2-lima, L1 + L2+lima), (-L1-L2-lima, (L1 + L2 + lima)/2))
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
lighter_blue = RGBf(0.4, 0.51, 0.878)
lighter_purple = RGBf(0.667, 0.475, 0.757)

function lighten(c, f = 1.2)
    c = to_color(c)
    hsl = Makie.HSLA(c)
    neg = Makie.RGBAf(Makie.HSLA(hsl.h, hsl.s, clamp(hsl.l*f, 0.0, 1.0), hsl.alpha))
    neg = Makie.RGBf(Makie.HSL(hsl.h, hsl.s, clamp(hsl.l*f, 0.0, 1.0)))
    return neg
end

lighter_blue_x = lighten(lighter_blue)

lighter_green  = julia_green
lighter_red    = julia_red
lighter_blue   = julia_blue
lighter_purple = julia_purple

# color of the trajectory (fading out)
c = lighter_blue_x
tailcoltransplight = [RGBAf(c.r, c.g, c.b, (i/tail)^(1.2)) for i in 1:tail]
c = julia_blue
tailcoltransp = [RGBAf(c.r, c.g, c.b, (i/tail)^(1.2)) for i in 1:tail]
tailcol = [RGBf(c.r, c.g, c.b) for i in 1:tail]
trajline = lines!(ax, traj; color = tailcoltransp, linewidth = 4)


# rods that connect the pendulum
rodlines = lines!(ax, balls; linewidth = 12, color = :black)

juliaballs = scatter!(ax, balls; marker = :circle, strokewidth = 10,
    strokecolor = [julia_green, julia_red, julia_purple],
    color = [lighter_green, lighter_red, lighter_purple],
    markersize = 160,
)

function animstep!(dp)
    step!(dp)
    update!(current_state(dp))
end

function update!(u)
    x1,x2,y1,y2 = xycoords(u)
    rod[] = [Point2f(0, 0), Point2f(x1, y1), Point2f(x2, y2)]
    balls[] = [Point2f(0, 0), Point2f(x1, y1), Point2f(x2, y2)]
    push!(traj[], Point2f(x2, y2))
    notify(traj)
end

function animstep!(dp::DynamicalSystem, t)
    t0 = current_time(dp)
    while current_time(dp) < t0 + t
        animstep!(dp)
    end
end

fig


# %% find initial condition that leads to nice Julia logo:
using DynamicalSystems.StateSpaceSets.Distances
u0 = SVector{4}([π/2 + 0.1, +2.03, 0.3, +5.412])
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
uf = current_state(dp)

# reported final time of evolution from given initial condition:
# 168.38499999991936

# Final state
# 5.757958811321035
# -5.924686393536801
# 58.11993236370379
# -2.839077350331716

fig


# %% okay, save high quality version:
ax.backgroundcolor = :transparent
CairoMakie.save(desktop("juliadynamics_logo.png"), fig; px_per_unit = 4)
# and one without tail
trajline.visible = false
save(desktop("juliadynamics_logo_no_tail.png"), fig; px_per_unit = 4)
trajline.visible = true
# and one more with white background
ax.backgroundcolor = :white
fig.scene.backgroundcolor = to_color(:white)
CairoMakie.save(desktop("juliadynamics_logo_white.png"), fig; px_per_unit = 4)
# and a dark background
rodlines.color = :white
trajline.color = tailcoltransp
ax.backgroundcolor = "#1e1e20"
fig.scene.backgroundcolor = to_color("#1e1e20")
CairoMakie.save(desktop("juliadynamics_logo_dark.png"), fig; px_per_unit = 4)
ax.backgroundcolor = :transparent
fig.scene.backgroundcolor = to_color(:transparent)
CairoMakie.save(desktop("juliadynamics_logo_dark_transp.png"), fig; px_per_unit = 4)

# %% add text axis
resize!(fig, 1600, 800)
texax = Axis(fig[:, 0]; backgroundcolor = "#1e1e20")

resize!(fig, 1600, 600)

lima = 0.1
ax.limits = ((-L1-L2-lima, L1 + L2+lima), (-L1-L2-lima, (L1 + L2 + lima)/2))

fig

# %% Add DynamicalSystems.jl font

empty!(texax)
hidedecorations!(texax)
hidespines!(texax)

dstext = text!(texax, 1.0, 0.5;
    text = "Dynamical\nSystems.jl", font = joinpath(@__DIR__, "montserrat", "Montserrat-Medium.ttf"),
    align = (:right, :center), justification = :right, space = :relative,
    color = :white, fontsize = 130, offset = (-35.0, 0),
)

scatter!(
    [0.69, 0.885], [0.725, 0.445]; markersize = 40,
    color = julia_blue
)

# add vertical line
vertline = lines!(texax, [0.99, 0.99], [0.2, 0.8]; color = :white, linewidth = 10,)
texax.limits = ((0, 1), (0, 1))
fig


# %% save full logo!
# dark version - solid
trajline.color = tailcoltransplight
ax.backgroundcolor = "#1e1e20"
texax.backgroundcolor = "#1e1e20"
fig.scene.backgroundcolor = to_color("#1e1e20")
rodlines.color = :white
CairoMakie.save(desktop("juliadynamics_full_logo_dark.png"), fig; px_per_unit = 4)

# dark version - transparent
texax.backgroundcolor = :transparent
ax.backgroundcolor = :transparent
fig.scene.backgroundcolor = to_color(:transparent)

CairoMakie.save(desktop("juliadynamics_full_logo_dark_transparent.png"), fig; px_per_unit = 4)

# light version - transparent
trajline.color = tailcoltransp
rodlines.color = :black
dstext.color = :black
vertline.color = :black

CairoMakie.save(desktop("juliadynamics_full_logo_light_transparent.png"), fig; px_per_unit = 4)

# light version - solid
texax.backgroundcolor = :white
ax.backgroundcolor = :white
fig.scene.backgroundcolor = to_color(:white)

CairoMakie.save(desktop("juliadynamics_full_logo_light.png"), fig; px_per_unit = 4)


# %% reset back to solid dark for animation
rodlines.color = :white
trajline.color = tailcoltransplight
dstext.color = :white
vertline.color = :white

ax.backgroundcolor = "#1e1e20"
texax.backgroundcolor = "#1e1e20"
fig.scene.backgroundcolor = to_color("#1e1e20")

# %% perform video animation animate from some t start to tf
# What I want to do here is the following:
# animate the pendulum motion so that it slows down as it reaches the final
# state, and it stops once it reaches the final state.

# For this, it is better to have a full ODE solution

# we use solve here, it is so much simpler for adjusting time
prob = ODEProblem((u,p,t) -> -doublependulum_rule(u,p,t), uf, (0.0, 100.0), p0)

sol = solve(prob; alg = Vern9(), dt, dense = true, adaptive = false)

fig

# %%
# We need to adjust limits of the pendulum axis and ball size because
# the pendulum goes outside the axis range
lima = 0.25
ax.limits = ((-L1-L2-lima, L1 + L2+lima), (-L1-L2-lima, (L1 + L2 + lima)/2))
juliaballs.markersize = 140

# time
span = tail*dt # total time to record for
freq = 10 # record every 10th step
times = 2span:-dt:0 # having 2 span allows us to bring the tail to the starting point

record(fig, desktop("juliadynamics_logo_anim.mp4"); framerate = 30) do io
    for (i, t) in enumerate(times)
        update!(sol(t))
        if t ≤ span # we only record the final `span`
            i % freq == 0 && recordframe!(io)
        end
    end
end

fig
