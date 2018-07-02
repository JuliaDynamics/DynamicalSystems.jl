using DynamicalSystems, PyPlot

cd()
mkpath("hhanim")

ds = Systems.henonheiles()

# Grid of initial conditions at given energy:

energy(x,y,px,py) = 0.5(px^2 + py^2) + potential(x,y)
potential(x, y) = 0.5(x^2 + y^2) + (x^2*y - (y^3)/3)
function generate_ics(E, n)
    ys = linspace(-0.4, 1.0, n)
    pys = linspace(-0.5, 0.5, n)
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

Emin = 0.1
Emax = 1/6
Es = linspace(Emin, Emax, 100)

density = 20
tfinal = 2000.0
tgali = 1000.0
tinteg = tangent_integrator(ds, 4)

# Initialize required stuff for plotting
mpl = PyPlot.matplotlib
norm = mpl[:colors][:Normalize](0.0, 1.0)
c_m = mpl[:cm][:Spectral]
s_m = mpl[:cm][:ScalarMappable](cmap=c_m, norm=norm)
s_m[:set_array]([])


# %%
ioff()
for j in 1:length(Es)
    j==1 && (tim = time())
    E = Es[j]
    println("Energy $(j)/$(length(Es))")
    figure(figsize = (12, 8))
    ax = gca()

    ics = generate_ics(E, density)

    for u in ics
        # compute section:
        psos = poincaresos(ds, (1, 0.0), tfinal; u0 = u)

        # compute gali (using advanced usage)
        set_state!(tinteg, u)
        set_deviations!(tinteg, orthonormal(4,4))
        reinit!(tinteg, tinteg.u)
        regularity = ChaosTools._gali(tinteg, tgali, 1, 1e-12)[2][end]/tgali

        # Plot PSOS with color corresponding to regularity
        y = psos[:, 2]; py = psos[:, 4]
        plot(y, py, ls = "None", ms = 0.8, marker = "o",
        color = s_m[:to_rgba](regularity), alpha = 0.75)
    end
    cb = colorbar(s_m, ticks=[])
    # cb = colorbar(s_m, ticks=[lowerend, 1.0])
    # cb[:ax][:set_yticklabels](["chaotic", "periodic"], rotation="vertical")
    cb[:set_label]("regularity")
    xlim(-0.55, 1.05)
    ylim(-0.65, 0.65)
    xlabel("\$y\$")
    ylabel("\$p_y\$")
    text(0.7, 0.9, "\$E = $(round(E, 5))\$", transform = ax[:transAxes])
    tight_layout()
    savefig("hhanim/galipsos_$j.png", dpi = 200)
    if j == 1
        ttim = time()
        println("Average time per energy: $(round(ttim - tim, 3)) secs")
    end
end
close("all")
ion()

# %%
framerate = 10
savename = "hhanim/galipsos"
anim = `ffmpeg -y -framerate $(framerate) -start_number 1 -i $(savename)_%d.png
-c:v libx264 -pix_fmt yuv420p -preset veryslow -profile:v baseline
-filter:v scale=1200:-1 -b:v 2048k -level 3.0 $(savename).mp4`
run(anim)

if false
    for i in 1:length(Es)
        rm(savename*"_$(i).png")
    end
end
