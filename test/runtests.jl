using Test, RosenbluthChains, Aqua, Scratch

TestPath = Scratch.get_scratch!(RosenbluthChains, "test")
PkgPath="/"*joinpath(split(pathof(RosenbluthChains),"/")[1:end-2])


mkpath("$TestPath/tmp/")

println("Visualisation of the test can be found at $TestPath/tmp/ if DOPLOTS is set to true.")

VERBOSE=true
DOPLOTS=true
EXTENSIVE = false ### not suitable for daily testing since test arent robust

include("./Helper.jl")
include("./NoMeasurement.jl")

@testset "Aqua" begin
    Aqua.test_all(RosenbluthChains; deps_compat=(ignore=[:Test,:Mmap, :Printf, :Random, :Scratch], ), project_extras=true,piracies = false )
    ### ignore standard libraries, not sure how to deal with them im PackageCompatUI/add compats manually
    ### Test fails in project extras since they dont get excluded their normaly
end


@testset "All tests" begin
    include("./SpecificModels_Test.jl")
    include("./Angle_Tests.jl")
    include("./BoltzmannFaktor_Test.jl")
    include("./LJ_Repulsion_Test.jl")
    include("./Calvados2/Calvados2_Test.jl")
    include("./RandomNumberTest.jl")
    include("./Cosine_BondAngles_Test.jl")
    include("./PERM_Test.jl")
end

