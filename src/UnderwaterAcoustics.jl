module UnderwaterAcoustics

using DSP
using LinearAlgebra
using SignalAnalysis

# basic underwater acoustics
include("basic.jl")

# propagation modeling
include("pm_core.jl")
include("pm_basic.jl")
include("pm_isoray.jl")

end # module
