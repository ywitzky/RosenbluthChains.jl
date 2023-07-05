using Test, RosenbluthChains
mkpath("./tmp/")

include("./NoMeasurement.jl")
include("./SpecificModels_Test.jl")
#println("Code so great, no test implemented yet.") 


rm("./tmp/")