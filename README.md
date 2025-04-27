[![doc-stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://org-arl.github.io/UnderwaterAcoustics.jl)
[![CI](https://github.com/org-arl/UnderwaterAcoustics.jl/workflows/CI/badge.svg)](https://github.com/org-arl/UnderwaterAcoustics.jl/actions)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![Codecov](https://codecov.io/gh/org-arl/UnderwaterAcoustics.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/org-arl/UnderwaterAcoustics.jl)
[![ColPrac](https://img.shields.io/badge/ColPrac-contributing-blueviolet)](CONTRIBUTING.md)

# UnderwaterAcoustics.jl

**Julia toolbox for underwater acoustic modeling**

## Overview

The `UnderwaterAcoustics.jl` ecosystem provides a set of tools for modeling and simulating underwater acoustic propagation. It defines a core application programming interface (API), and provides a set of core differentiable propagation models and utilities. It also provides support for replay channels, where measurements from the ocean are used to empirically construct channel models. The package is designed to be extensible, allowing other packages to add models. Several packages (`AcousticsToolbox.jl`, `AcousticRayTracers.jl`, `VirtualAcousticOcean.jl`, etc) add more propagation models and related tools, thus forming a rich ecosystem of models and tools for underwater acoustic modeling and simulation.

For more information, see [documentation](https://org-arl.github.io/UnderwaterAcoustics.jl).

## Highlights

- Underwater acoustic propagation modeling API with pluggable models
- Differentiable and probabilistic underwater acoustic models
- 2D/3D underwater acoustic simulation tools
- Replay channel and noise models
- Underwater acoustics utility functions

> [!IMPORTANT]
> The API has changed significantly in `UnderwaterAcoustics v0.4`. If you have code that depends
> on the old API, you may wish to refer to the [Porting guide](https://org-arl.github.io/UnderwaterAcoustics.jl/porting.html) before you upgrade.

## Contributing

Contributions in the form of bug reports, feature requests, ideas/suggestions, bug fixes, code enhancements, and documentation updates are most welcome. Please read [contribution guidelines](https://github.com/org-arl/UnderwaterAcoustics.jl/blob/master/CONTRIBUTING.md) if you wish to start contributing.

## Talks & publications

- Mandar Chitre, "[Differentiable Ocean Acoustic Propagation Modeling](https://arl.nus.edu.sg/wp-content/uploads/2023/04/Chitre_Differentiable-Ocean-Acoustic-Propagation-Modeling.pdf)," in OCEANS 2023 IEEE/MTS – Limerick, 5-8 June 2023. [[doi]](https://doi.org/10.1109/OCEANSLimerick52467.2023.10244307)
- Mandar Chitre, "[Underwater Acoustics in the age of differentiable and probabilistic programming](https://www.facebook.com/watch/live/?v=2473971036238315)", UComms 2020 webinar, 3 December 2020.

## Citing

If you use `UnderwaterAcoustics.jl` in your work or are influenced by its ideas, please cite:
```
@inproceedings{chitre2023ua,
  author={Chitre, Mandar},
  booktitle={OCEANS 2023 - Limerick},
  title={Differentiable Ocean Acoustic Propagation Modeling},
  year={2023},
  pages={1-8},
  doi={10.1109/OCEANSLimerick52467.2023.10244307}
}
```
