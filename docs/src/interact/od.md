# Interactive Orbit Diagram
## Docstrings
```@docs
interactive_orbitdiagram
scaleod
```
## Function Video
```julia
using InteractiveChaos, Makie

i = 1
p_index = 1

systems = [(Systems.logistic(), 3.0, 4.0, "r"),
           (Systems.henon(), 0.8, 1.4, "a"),
           (Systems.standardmap(), 0.6, 1.2, "k")]

ds, p_min, p_max, parname = systems[1]

oddata = interactive_orbitdiagram(
           ds, i, p_index, p_min, p_max;
           parname = parname
         )
```

<video width="100%" height="auto" controls autoplay loop>
<source src="https://raw.githubusercontent.com/JuliaDynamics/JuliaDynamics/master/videos/interact/orbitdiagram.mp4?raw=true" type="video/mp4">
</video>

## Video Tutorial
<iframe width="560" height="315" src="https://www.youtube.com/watch?v=XECr8iyqfDc" frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
