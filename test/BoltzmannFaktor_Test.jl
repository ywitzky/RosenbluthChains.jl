
N_Test=1_000_000

N_Trial=8
N= 100
θ=collect(range(0.01,π/4, N-2))
φ=collect(range(0.01,π, N-2))
r=3.8
Data = SimData("./tmp/", 1.0, N, N_Trial      , 1, 1)
KP = SimulationParameters( FixedBondParameters(r), FixedBondAngles(θ), FixedTorsionAngles(φ))
Meas = RunSim(Data,KP, NoMeasurement());


LogDist=[0.1, 0.2,0.3, 0.452, 0.5923, 0.534, 0.134, 0.234]

### copy of FixedBondParameters
struct TestBondParameters{T<: Real} <: AbstractBondParam
    ### Bond Length
    BondLength::T
end 

@inline function InitSimParam(data::SimData,param::TestBondParameters ) 
    fill!(data.trial_radius, param.BondLength)
end


@inline function SetTrialRadius(data::SimData, param::TestBondParameters) 
    ### is set once in set up and doesnt change afterwards
end

@inline function GetTrialBoltzmannWeight(data::SimData,param::TestBondParameters) nothing end


@testset "Rosenbluth Weights & Boltzmann factors" begin
    Dists= [[9.0, 1.0, 0.0, 0.0, 0.0,0.0,0.0, 0.0], [0.05,0.05,0.2,0.4,0.05,0.05, 0.1,0.1]]
    for Dist in Dists
        function RosenbluthChains.GetTrialBoltzmannWeight(data::SimData,param::TestBondParameters) 
            data.LogBoltzmannFaktor .+= log.(Dist)
            nothing
        end
        Data.BoltzmannFaktor .= Dist
        result = [begin Data.id=4; ChooseTrialPosition(Data, KP);   Data.tid;  end for _ in 1:N_Test]
        hist = [ count(x->x==(a), result) for a in 1:N_Trial]./N_Test

        @test all( abs.(Dist/sum(Dist) .-hist).<5*10^-3)
    end
    function RosenbluthChains.GetTrialBoltzmannWeight(data::SimData,param::TestBondParameters) 
        data.LogBoltzmannFaktor .+= LogDist
        nothing
    end
    
    RosenbluthChains.GetTrialBoltzmannWeight(Data, KP)
    @test all(Data.BoltzmannFaktor.*Data.NTrials .≈ exp.(LogDist))
end;
