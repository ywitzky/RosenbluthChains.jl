module RosenbluthChains

using StaticArrays, Mmap

include("./Helper.jl")
include("./VectorTypes.jl")
include("./SimulationData.jl")
include("./AbstractTypes.jl")
include("./HelperFuncs.jl")
include("./Measurements.jl")
include("./MMapOutput.jl")
include("./RandomNumberFuncs.jl")
include("./RosenbluthMethod.jl")
include("./SpecificModels.jl")
include("./PERM.jl")

include("./Models/Cosine_BondAngels.jl")
include("./Models/LJ_Repulsion.jl")
include("./Models/Calvados2/Calvados2.jl")

end # module RosenbluthChains
