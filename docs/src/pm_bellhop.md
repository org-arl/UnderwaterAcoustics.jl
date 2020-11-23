# Bellhop

```@meta
CurrentModule = UnderwaterAcoustics
```

!!! note

    To use Bellhop, you first need to install [OALIB Bellhop](http://oalib.hlsresearch.com/AcousticsToolbox/) and ensure you have `bellhop.exe` available on your `PATH`.

Bellhop is an interface to the [OALIB Bellhop](http://oalib.hlsresearch.com/AcousticsToolbox/) 2D Gaussian beam tracer. If `UnderwaterAcoustics` can find `bellhop.exe`, the Bellhop model will be available:

```julia-repl
julia> models()
3-element Array{Any,1}:
 PekerisRayModel
 RaySolver
 Bellhop
```

Additional options available with [`Bellhop`](@ref):

- `nbeams` -- number of beams used for ray tracing (default: auto)
- `minangle` -- minimum beam angle in degrees (default: auto)
- `maxangle` -- maximum beam angle in degrees (default: auto)
- `debug` -- if `true`, intermediate Bellhop files are made available for user inspection (default: `false`)

**Example:**

```julia
using UnderwaterAcoustics
using Plots

env = UnderwaterEnvironment(
  seasurface = Vacuum,
  seabed = SandyClay,
  ssp = SampledSSP(0.0:20.0:40.0, [1540.0, 1510.0, 1520.0], :smooth),
  bathymetry = SampledDepth(0.0:50.0:100.0, [40.0, 35.0, 38.0], :linear)
)
pm = Bellhop(env)
tx = AcousticSource(0.0, -5.0, 1000.0)
r = rays(pm, tx, -60°:2°:60°, 100.0)
plot(env; sources=[tx], rays=r)
```

![](images/rays2.png)
