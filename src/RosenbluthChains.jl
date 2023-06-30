module RosenbluthChains

using StaticArrays
using Mmap

include("./VectorTypes.jl")
include("./SimulationData.jl")
include("./AbstractTypes.jl")
include("./HelperFuncs.jl")
include("./Measurements.jl")
include("./MMapOutput.jl")
include("./RandomNumberFuncs.jl")
include("./RosenbluthMethod.jl")
include("./SpecificModels.jl")


end # module RosenbluthChains
