module UnderwaterAcoustics

using Requires

using DSP: amp2db, db2amp
# using LinearAlgebra
# using SignalAnalysis

# basic underwater acoustics
include("basic.jl")

# propagation modeling
include("pm_core.jl")
include("pm_basic.jl")
include("pm_pekeris.jl")
include("pm_bellhop.jl")
include("pm_all.jl")

# plot recipes
function __init__()
  @require Plots="91a5bcdd-55d7-5caf-9e0b-520d859cae80" begin
    include("plot.jl")
  end
end

end # module
