function DynamicalSystems.interactive_poincaresos_scan(A::AbstractStateSpaceSet, j::Int; kwargs...)
    interactive_poincaresos_scan([A], j; kwargs...)
end

function DynamicalSystems.interactive_poincaresos_scan(
        As::Vector{<:AbstractStateSpaceSet}, j::Int;
        linekw = (), scatterkw = (), direction = -1,
        colors = [COLORS[i] for i in 1:length(As)]
    )

    for A in As; @assert dimension(A) == 3; end
    @assert j ∈ 1:3
    mi, ma = total_minmaxima(As)

    otheridxs = SVector(setdiff(1:3, j)...)

    figure = Figure(size = (1600, 800))
    display(figure)
    ax = figure[1, 1] = Axis3(figure)
    axp = figure[1, 2] = Axis(figure)
    sg = SliderGrid(figure[2,:],
        (label = "$(('x':'z')[j]) =", range =  range(mi[j], ma[j]; length = 1001),
        startvalue = (ma[j]+mi[j])/2)
    )
    y = sg.sliders[1].value

    for i in 1:length(As)
        A = As[i]
        # plot 3D trajectories
        lines!(ax, vec(A); color = colors[i], transparency = false, linekw...)

        # Poincare sos
        psos = lift(y) do y
            poincaresos(A, (j, y); direction, warning = false)
        end
        psos2d = lift(p -> vec(p[:, otheridxs]), psos)
        psos3d = lift(p -> vec(p), psos)
        Makie.scatter!(axp, psos2d; color = colors[i], scatterkw...)
        Makie.scatter!(ax, psos3d; color = colors[i], markersize = 5, scatterkw...)
    end

    xlims!(axp, mi[otheridxs[1]], ma[otheridxs[1]])
    ylims!(axp, mi[otheridxs[2]], ma[otheridxs[2]])

    # plot transparent plane
    ss = [mi...]
    ws = [(ma - mi)...]
    ws[j] = 0

    p = lift(y) do y
        ss[j] = y
        o = Makie.Point3f(ss...)
        w = Makie.Point3f(ws...)
        Makie.FRect3D(o, w)
    end

    a = RGBAf(0,0,0,0)
    c = RGBAf(0.2, 0.2, 0.25, 1.0)
    img = Makie.ImagePattern([c a; a c])
    Makie.mesh!(ax, p; color = img)

    return figure
end

function total_minmaxima(As::Vector{<:AbstractStateSpaceSet})
    mi, ma = minmaxima(As[1])
    mi = Vector(mi); ma = Vector(ma)
    for j in 2:length(As)
        mi2, ma2 = minmaxima(As[j])
        for i in 1:size(As[j], 2)
            mi[i] = min(mi2[i], mi[i])
            ma[i] = max(ma2[i], ma[i])
        end
    end
    return mi, ma
end

