# Trajectory Highlighter
## Docstrings
```@docs
trajectory_highlighter
```
## Function Video
```julia
using InteractiveChaos, Makie

ds = Systems.henonheiles()

# Grid of initial conditions at given energy:
energy(x,y,px,py) = 0.5(px^2 + py^2) + potential(x,y)
potential(x, y) = 0.5(x^2 + y^2) + (x^2*y - (y^3)/3)
function generate_ics(E, n)
    ys = range(-0.4, stop = 1.0, length = n)
    pys = range(-0.5, stop = 0.5, length = n)
    ics = Vector{Vector{Float64}}()
    for y in ys
        V = potential(0.0, y)
        V ≥ E && continue
        for py in pys
            Ky = 0.5*(py^2)
            Ky + V ≥ E && continue
            px = sqrt(2(E - V - Ky))
            ic = [0.0, y, px, py]
            push!(ics, [0.0, y, px, py])
        end
    end
    return ics
end

density = 15
tfinal = 2000.0
tgali = 1000.0
E = 0.13
ics = generate_ics(E, density)

tinteg = tangent_integrator(ds, 4)

regularity = Float64[]; psos = Dataset{2, Float64}[]
trs = Dataset{3, Float64}[]
@time for u in ics
    # compute gali (using advanced usage)
    reinit!(tinteg, u, orthonormal(4,4))
    push!(regularity, gali(tinteg, tgali, 1, 1e-12)[2][end]/tgali)
    push!(psos, poincaresos(ds, (1, 0.0), 2000.0; u0 = u, idxs = [2, 4]))
    tr = trajectory(ds, 200.0, u)[:, [1, 2, 4]]
    push!(trs, tr)
end

# %%
# 2D version:
trajectory_highlighter(psos, regularity; α = 0.05, hname = "regularity")
# 3D version:
trajectory_highlighter(trs[1:10:end], regularity[1:10:end];
nbins = 10, α = 0.05, linewidth = 4.0, hname = "regularity")

```
2D Version:

<video width="100%" height="auto" controls autoplay loop>
<source src="https://raw.githubusercontent.com/JuliaDynamics/JuliaDynamics/master/videos/interact/highlighter2D.mp4?raw=true" type="video/mp4">
</video>

3D Version:

<video width="100%" height="auto" controls autoplay loop>
<source src="https://raw.githubusercontent.com/JuliaDynamics/JuliaDynamics/master/videos/interact/highlighter3D.mp4?raw=true" type="video/mp4">
</video>
