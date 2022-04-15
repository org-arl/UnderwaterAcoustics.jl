module UnderwaterAcoustics

using Requires
using DocStringExtensions

using DSP: amp2db, db2amp, pow2db, db2pow
using Printf

export °

const ° = π/180.0

# basic underwater acoustics
include("uw_basic.jl")

# propagation modeling
include("pm_core.jl")
include("pm_basic.jl")
include("pm_pekeris.jl")
include("pm_all.jl")

# plot recipes
function __init__()
  @require Plots="91a5bcdd-55d7-5caf-9e0b-520d859cae80" begin
    include("plot.jl")
  end
end

end # module
