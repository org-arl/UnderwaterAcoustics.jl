module UnderwaterAcoustics

# common utilities
include("utils.jl")

# underwater acoustics
include("uw_basic.jl")

# propagation modeling
include("pm_api.jl")
include("pm_stdlib.jl")
include("pm_pekeris.jl")

end # module
