module WormModel

export std_pars
export calc_fst
export model

include("fst.jl")
include("multiloci.jl")
include("std_pars.jl")
include("utils.jl")
include("model.jl")

end # module
