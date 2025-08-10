module UnderwaterAcoustics

# common utilities
include("utils.jl")

# underwater acoustics
include("basic.jl")

# propagation modeling
include("api.jl")
include("stdlib.jl")
include("pekeris.jl")
include("adiabatic.jl")
include("replay.jl")

end # module
