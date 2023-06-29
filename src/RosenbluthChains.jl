module RosenbluthChains

using StaticArrays
using Mmap

include("./AbstractTypes.jl")
include("./HelperFuncs.jl")
include("./Measurements.jl")
include("./MMapOutput.jl")
include("./Models.jl")
include("./RandomNumberFuns.jl")
include("./RosenbluthMethod.jl")
include("./SimulatinData.jl")
include("./SpecificModels.jl")
include("./VectorTypes.jl")

end # module RosenbluthChains
