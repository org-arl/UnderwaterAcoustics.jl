module UnderwaterAcoustics

using Requires

using DSP
using LinearAlgebra
using SignalAnalysis
using Interpolations

# basic underwater acoustics
include("basic.jl")

# propagation modeling
include("pm_core.jl")
include("pm_basic.jl")
include("pm_pekeris.jl")

# plot recipes
function __init__()
  @require Plots="91a5bcdd-55d7-5caf-9e0b-520d859cae80" begin
    include("plot.jl")
  end
end

end # module
