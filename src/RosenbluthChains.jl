module RosenbluthChains

using StaticArrays, Mmap

include("./VectorTypes.jl")
include("./SimulationData.jl")
include("./AbstractTypes.jl")
include("./HelperFuncs.jl")
include("./Measurements.jl")
include("./MMapOutput.jl")
include("./RandomNumberFuncs.jl")
include("./RosenbluthMethod.jl")
include("./SpecificModels.jl")

include("./Models/Cosine_BondAngels.jl")


end # module RosenbluthChains
