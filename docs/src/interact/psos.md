# Interactive Poincaré Surface of Section
## Docstrings
```@docs
interactive_poincaresos
```

## Function Video
```julia
using InteractiveChaos, Makie

ds = Systems.henonheiles()

potential(x, y) = 0.5(x^2 + y^2) + (x^2*y - (y^3)/3)
energy(x,y,px,py) = 0.5(px^2 + py^2) + potential(x,y)
const E = energy(get_state(ds)...)

function complete(y, py, x)
    V = potential(x, y)
    Ky = 0.5*(py^2)
    Ky + V ≥ E && error("Point has more energy!")
    px = sqrt(2(E - V - Ky))
    ic = [x, y, px, py]
    return ic
end

chaotic = get_state(ds)
stable = [0., 0.1, 0.5, 0.]

plane = (1, 0.0)

psos = interactive_poincaresos(ds, plane, (2, 4), complete; markersizes = (-5, -1))
```

<video width="100%" height="auto" controls autoplay loop>
<source src="https://raw.githubusercontent.com/JuliaDynamics/JuliaDynamics/master/videos/interact/interactive_psos.mp4?raw=true" type="video/mp4">
</video>

## Video Tutorial
<iframe width="560" height="315" src="https://www.youtube.com/watch?v=SozXxa7blJs" frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
