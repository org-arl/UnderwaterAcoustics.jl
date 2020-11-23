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

env = UnderwaterEnvironment(seasurface=Vacuum)
pm = Bellhop(env; nbeams=180)
tx = AcousticSource(0.0, -5.0, 1000.0)
rx = AcousticReceiver(100.0, -10.0)
loss = transmissionloss(pm, tx, rx)
```

