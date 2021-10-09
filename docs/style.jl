struct CyclicContainer <: AbstractVector{String}
    c::Vector{String}
    n::Int
end
CyclicContainer(c) = CyclicContainer(c, 0)

Base.length(c::CyclicContainer) = typemax(Int)
Base.size(c::CyclicContainer) = size(c.c)
Base.iterate(c::CyclicContainer, state=1) = Base.iterate(c.c, state)
Base.getindex(c::CyclicContainer, i) = c.c[(i-1)%length(c.c) + 1]
function Base.getindex(c::CyclicContainer)
    c.n += 1
    c[c.n]
end
Base.iterate(c::CyclicContainer, i = 1) = iterate(c.c, i)

COLORSCHEME = [
    "#6D44D0",
    "#2CB3BF",
    "#1B1B1B",
    "#DA5210",
    "#03502A",
    "#866373",
]

COLORS = CyclicContainer(COLORSCHEME)
LINESTYLES = CyclicContainer(["-", ":", "--", "-."])

using PyPlot
using3D()

PyPlot.rc("lines", lw = 2)
PyPlot.rc("axes", grid = true)
PyPlot.rc("grid", color = "0.75", alpha = 0.75)

PyPlot.rc("font", size = 22) # set default fontsize
PyPlot.rc("xtick", labelsize = 20)
PyPlot.rc("ytick", labelsize = 20)
PyPlot.rc("axes", labelsize = 24)
PyPlot.rc("legend", fontsize = 20)
# PyPlot.rc("font", family = "Times New Roman") # Serif main font
PyPlot.rc("font", family = "DejaVu Sans") # sans main font
# PyPlot.rc("mathtext", rm = "sanserif", fontset="dejavusans") # sans math font
PyPlot.rc("mathtext", rm = "serif", fontset="cm") # serif math font

for z in ("x", "y")
    PyPlot.rc("$(z)tick.major", size = 6, width = 1.2)
    PyPlot.rc("$(z)tick.minor", size = 3, visible = false)
end

figx = 9
figy = 7
PyPlot.rc("figure", figsize = (figx, figy))
PyPlot.rc("savefig", dpi = 300, transparent = false, format = "png")

# set default color cycle
PyPlot.rc("axes", prop_cycle = matplotlib.cycler(color=COLORS.c))
