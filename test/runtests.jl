using Test, RosenbluthChains
mkpath("./tmp/")

include("./NoMeasurement.jl")
#include("./SpecificModels_Test.jl")
#include("./BoltzmannFaktor_Test.jl")
#include("./RandomNumberTest.jl")
#include("./Cosine_BondAngles_Test.jl")
#include("./LJ_Repulsion_Test.jl")

include("./Calvados2/Calvados2_Test.jl")


#rm("./tmp/")
