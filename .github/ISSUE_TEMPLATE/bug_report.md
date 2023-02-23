---
name: Bug report
about: Create a report to help us improve
title: ''
labels: ''
assignees: ''

---

**Describe the bug**
A clear and concise description of what the bug is.

**Minimal Working Example**
Please provide a piece of code that leads to the bug you encounter.

If the code is **runnable**, it will help us identify the problem faster.

**Package versions**

Please provide the versions you use. To do this, run the code:
```julia
using Pkg
Pkg.status([
    "DynamicalSystems",
    "StateSpaceSets", "DynamicalSystemsBase", "RecurrenceAnalysis", "FractalDimensions", "DelayEmbeddings", "ComplexityMeasures", "TimeseriesSurrogates", "PredefinedDynamicalSystems", "Attractors", "ChaosTools"
    ];
    mode = PKGMODE_MANIFEST
)
```
