using Test, RosenbluthChains
mkpath("./tmp/")

VERBOSE=true
DOPLOTS=true

include("./Helper.jl")
include("./NoMeasurement.jl")

@testset "All tests" begin
    #include("./SpecificModels_Test.jl")
    #include("./Angle_Tests.jl")
    #include("./BoltzmannFaktor_Test.jl")
    #include("./Cosine_BondAngles_Test.jl")
    include("./LJ_Repulsion_Test.jl")
    #include("./PERM_Test.jl")
    #include("./Calvados2/Calvados2_Test.jl")

    include("./RandomNumberTest.jl")### has to be last, currently has side effects

end
#rm("./tmp/")
