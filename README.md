[![CI](https://github.com/org-arl/UnderwaterAcoustics.jl/workflows/CI/badge.svg)](https://github.com/org-arl/UnderwaterAcoustics.jl/actions)
[![doc-stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://org-arl.github.io/UnderwaterAcoustics.jl/stable)
[![doc-dev](https://img.shields.io/badge/docs-latest-blue.svg)](https://org-arl.github.io/UnderwaterAcoustics.jl/dev)
[![Codecov](https://codecov.io/gh/org-arl/UnderwaterAcoustics.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/org-arl/UnderwaterAcoustics.jl)

# UnderwaterAcoustics.jl

### Julia toolbox for underwater acoustic modeling

![](https://org-arl.github.io/UnderwaterAcoustics.jl/dev/images/txloss1.png)

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

## Related packages

**NOTE: In version 0.2, `RaySolver` and `Bellhop` models have been moved out to separate packages.**

- Install [`AcousticRayTracers.jl`](https://github.com/org-arl/AcousticRayTracers.jl) for `RaySolver` model
- Install [`AcousticsToolbox.jl`](https://github.com/org-arl/AcousticsToolbox.jl) for `Bellhop` and `Kraken` models

## Getting started

- Propagation modeling toolkit -- [quickstart guide](https://org-arl.github.io/UnderwaterAcoustics.jl/stable/pm_basic.html)
- Probabilistic propagation modeling -- [tutorial](https://org-arl.github.io/UnderwaterAcoustics.jl/stable/tut_turing.html)
- Differentiable propagation modeling -- [tutorial](https://org-arl.github.io/UnderwaterAcoustics.jl/stable/tut_autodiff.html)

## Contributing

Contributions in the form of bug reports, feature requests, ideas/suggestions, bug fixes, code enhancements, and documentation updates are most welcome. Please read [contribution guidelines](CONTRIBUTING.md) if you wish to start contributing.

## Talks & publications

- Mandar Chitre, "[Underwater Acoustics in the age of differentiable and probabilistic programming](https://www.facebook.com/watch/live/?v=2473971036238315)", UComms 2020 webinar, 3 December 2020.
