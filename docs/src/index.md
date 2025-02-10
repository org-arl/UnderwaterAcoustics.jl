# UnderwaterAcoustics.jl
### Julia toolbox for underwater acoustic modeling

```@meta
CurrentModule = UnderwaterAcoustics
```

## Highlights

- Underwater acoustic propagation modeling with pluggable models
- 2D/3D underwater acoustic simulation tools
- Differentiable and probabilistic underwater acoustic modeling
- Underwater acoustics utility functions

## Installation

```julia-repl
julia>]
pkg> add UnderwaterAcoustics
```

!!! warning "API change"

    The API has changed significantly in `UnderwaterAcoustics v0.4`. If you have code that depends
    on the old API, you may wish to refer to the [Porting guide](@ref) before you upgrade.

## Getting started

- [Propagation & channel modeling](@ref)
- [Channel replay](@ref)

## Documentation

- [Underwater acoustics](@ref)
- [Propagation & channel modeling API](@ref)
- [Developing your own propagation or channel models](@ref)
- [General utilities](@ref)
- [Porting guide](@ref) from `v0.3` to `v0.4` API

## Talks & publications

- Mandar Chitre, "[Underwater Acoustics in the age of differentiable and probabilistic programming](https://www.facebook.com/watch/live/?v=2473971036238315)", UComms 2020 webinar, 3 December 2020.
- Mandar Chitre, "[Differentiable Ocean Acoustic Propagation Modeling](https://arl.nus.edu.sg/wp-content/uploads/2023/04/Chitre_Differentiable-Ocean-Acoustic-Propagation-Modeling.pdf)," in OCEANS 2023 IEEE/MTS – Limerick, 5-8 June 2023.
