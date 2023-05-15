export brainscan_poincaresos

"""
    brainscan_poincaresos(A::Dataset, j::Int; kwargs...)
    brainscan_poincaresos(As::Vector{Dataset}, j::Int; kwargs...)
Launch an interactive application for scanning a Poincare surface of section of `A`
like a "brain scan", where the plane that defines the section can be arbitrarily
moved around via a slider. Return `figure, ax3D, ax2D`.

The input dataset must be 3 dimensional, and here the crossing plane is always
chosen to be when the `j`-th variable of the dataset crosses a predefined value.
The slider automatically gets all possible values the `j`-th variable can obtain.

If given multiple datasets, the keyword `colors` attributes a color to each one, e.g.
`colors = [JULIADYNAMICS_COLORS[mod1(i, 6)] for i in 1:length(As)]`.

The keywords `linekw, scatterkw` are named tuples that are propagated as keyword arguments
to the line and scatter plot respectively, while the keyword `direction = -1` is propagated
to the function `DyamicalSystems.poincaresos`.
"""
function brainscan_poincaresos(A::DynamicalSystems.AbstractDataset, j::Int; kwargs...)
    brainscan_poincaresos([A], j; kwargs...)
end

function brainscan_poincaresos(
        As::Vector{<:DynamicalSystems.AbstractDataset}, j::Int;
        linekw = (), scatterkw = (), direction = -1,
        colors = [CYCLIC_COLORS[i] for i in 1:length(As)]
    )

    for A in As; @assert size(A, 2) == 3; end
    @assert j ∈ 1:3
    mi, ma = total_minmaxima(As)

    otheridxs = DynamicalSystems.SVector(setdiff(1:3, j)...)

    figure = Figure(resolution = (1600, 800))
    display(figure)
    ax = figure[1, 1] = Axis3(figure)
    axp = figure[1, 2] = Axis(figure)
    # old
    # sll = labelslider!(
    #     figure, "$(('x':'z')[j]) =", range(mi[j], ma[j]; length = 100);
    #     sliderkw = Dict(:startvalue => (ma[j]+mi[j])/2)
    # )
    # figure[2, :] = sll.layout
    # y = sll.slider.value
    # new
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
            DynamicalSystems.poincaresos(A, (j, y); direction, warning = false)
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
    img = Makie.ImagePattern([c a; a c]);
    Makie.mesh!(ax, p; color = img);

    return figure, ax, axp
end

function total_minmaxima(As::Vector{<:DynamicalSystems.AbstractDataset})
    mi, ma = DynamicalSystems.minmaxima(As[1])
    mi = Vector(mi); ma = Vector(ma)
    for j in 2:length(As)
        mi2, ma2 = DynamicalSystems.minmaxima(As[j])
        for i in 1:size(As[j], 2)
            mi[i] = min(mi2[i], mi[i])
            ma[i] = max(ma2[i], ma[i])
        end
    end
    return mi, ma
end

