

N= 100
θ=collect(range(0.01,π/4, N-2))
#θ=reverse(collect(range(0.01,π/4, N-2)))

θ=ones(N-2)*π/8
φ=collect(range(0.01,π, N-2))
φ=ones(N-2)*π/6
φ[1:10] .= ones(10)*π/8
r=3.8
Data = SimData("./tmp/", 1.0, N, 1      , 1, 1)
KP = SimulationParameters( FixedBondParameters(r), FixedBondAngles(θ), FixedTorsionAngles(φ))
Meas = RunSim(Data,KP, NoMeasurement());

@testset "fixed bond length" begin
    @test all(r .≈ [ norm(Data.xyz[i+1]-Data.xyz[i]) for i in 1:N-1   ])
end;

@testset "fixed bond angles" begin
    @test all((θ .≈ [ angle(Data.xyz[i+1]-Data.xyz[i],Data.xyz[i+2]-Data.xyz[i+1] ) for i in 1:N-2   ]) )
end;


@testset "fixed torsion angles" begin
    #https://www.math.fsu.edu/~quine/MB_10/6_torsion.pdf
    rel  = [Data.xyz[i+1]-Data.xyz[i] for i in 1:N-1]
    tors3 = [torsion(rel[i:i+2]...) for i in 1:N-3]
    @test all(φ[2:end] .≈ tors3) 
end;
