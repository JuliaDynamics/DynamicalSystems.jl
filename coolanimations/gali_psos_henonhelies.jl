using DynamicalSystems, PyPlot

cd(@__DIR__)
mkpath("hhanim")

ds = Systems.henonheiles()

# Random initial conditions at given energy:
potential(x, y) = 0.5(x^2 + y^2) + (x^2*y - (y^3)/3)
function randomic(E)
    φ = 2π*rand()
    y = 1.4rand() - 0.4
    V = potential(0, y)
    while V ≥ E
        y = rand() - 0.4
        V = potential(0, y)
    end
    p = sqrt(2*(E - V))
    return [0, y, p*cos(φ), p*sin(φ)]
end

Emin = 0.1
Emax = 1/6
Es = linspace(Emin, Emax, 100)

total = 100
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
    E = Es[j]
    println("Energy $(j)/$(length(Es))")
    figure(figsize = (12, 8))
    ax = gca()
    for i in 1:total
        u = randomic(E)

        # compute section:
        psos = poincaresos(ds, (1, 0.0), tfinal; u0 = u)

        # compute gali (using advanced usage)
        set_state!(tinteg, u)
        set_deviations!(tinteg, orthonormal(4,4))
        reinit!(tinteg, tinteg.u)
        regularity = ChaosTools._gali(tinteg, tgali, 1, 1e-12)[2][end]/tgali

        # Plot PSOS with color corresponding to regularity
        y = psos[:, 2]; py = psos[:, 4]
        plot(y, py, ls = "None", ms = 1.0, marker = "o",
        color = s_m[:to_rgba](regularity))
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
    savefig("hhanim/galipsos_$j.png")
end
close("all")
ion()
