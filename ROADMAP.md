# Project roadmap

## Scope

The `UnderwaterAcoustics.jl` project aspires to be the go-to tool for underwater acoustic propagation modeling. In particular, the following topics are within the scope of the project:

- Basic tools for computing various quantities related to underwater acoustics.
- Various underwater acoustic propagation models based on ray theory, normal modes, parabolic equations, fast field programs, finite element methods, etc.
- Sediment and seabed acoustics, to the extent it affects acoustics in the water.
- Physical oceanography, to the extent it affects acoustics in the water.
- Bubble acoustics, since bubbles in water contribute significantly acoustic propagation and ambient noise in water.
- Replay and stochastic acoustic channel models that allow signals to be transmitted through a simulated underwater channel.
- Data-driven acoustic modeling combining machine learning techniques with the physics of underwater acoustics.
- Modeling of acoustic channel variability.
- Modeling of ocean acoustic noise.

The base package aims to be lightweight and so avoids too many heavy dependencies. More complex models may be implemented as related packages (such as `AcousticsToolbox.jl`, `AcousticRayTracers.jl`, etc) that add functionality to the main package.

## Out of scope

Specifically, the following are out of scope for the project:

- Seismic and in-air acoustics, beyond its impact on underwater acoustics.
- Physical oceanography, beyond its impact on underwater acoustics.
- Acoustic signal processing. See [SignalAnalysis.jl](https://github.com/org-arl/SignalAnalysis.jl), which is a closely related package that provides acoustic signal processing tools. `UnderwaterAcoustics.jl` leverages on tools provided by this package.

## Features & enhancements

For a list of planned features and enhancements, see [issues](https://github.com/org-arl/UnderwaterAcoustics.jl/issues).
