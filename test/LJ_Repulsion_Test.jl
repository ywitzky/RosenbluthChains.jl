using Plots, LsqFit

mutable struct AMeasurement{T<:Number} <: AbstractMeasurement 
    SumRosenbluthWeights::T
    REE::Vector{T}
    Weight::Vector{T}
    AMeasurement(NBeads::I, NSteps::I) where {I<:Integer} = new{Float64}(zero(Float64 ), zeros(Float64,NSteps), zeros(Float64, NSteps)) #where {T::BigFloat}
end


function RosenbluthChains.InitMeasurement(data::SimData, param::SimulationParameters, Tmp::AMeasurement) 
    return AMeasurement(data.NBeads, data.BatchSize)
end 

function RosenbluthChains.MeasureAfterBatch(data::SimData, param::SimulationParameters,Measurement::AMeasurement) 
end

function RosenbluthChains.MeasureAfterChainGrowth(data::SimData, param::SimulationParameters, Measurement::AMeasurement) 
    #println(data.xyz)
    Measurement.REE[data.id_in_batch]    = norm(data.xyz[1]-data.xyz[end])
    #println(norm(data.xyz[1]-data.xyz[end]))
    #println(Measurement.REE[data.id_in_batch] )
    Measurement.Weight[data.id_in_batch] = getRosenbluthWeigth(data, param)
    #println(data.RosenbluthWeight)
    Measurement.SumRosenbluthWeights += data.RosenbluthWeight
end
    

N= 100
N_Trial=100
θ=collect(range(0.01,π/4, N-2))
#θ=reverse(collect(range(0.01,π/4, N-2)))

θ=ones(N-2)*π/8
φ=collect(range(0.01,π, N-2))
φ=ones(N-2)*π/6
φ[1:10] .= ones(10)*π/8
r=3.8

#Meas = RunSim(Data,KP, NoMeasurement());

LJ(r,σ, ϵ) =  r < 2^(1.0/6.0)*σ ? 4*ϵ*((σ/r)^12-(σ/r)^6)  + ϵ : 0.0

REE_Fit(N, x) = @. x[1]*(N-1)^x[2]

@testset "LJ Potential" begin
    N_Trial=100
    N=100

    Data = SimData("./tmp/", 1.0, N, N_Trial     , 1, 1)
    KP = SimulationParameters( FixedBondParameters(r), RandBondAngle(), RandTorsion(), LJ_Repulsion(ones(Int64,N),Dict(1 => r), 3.0 ))
    Data.xyz .= [Vector3([0.0,0.0,0.0]) for _ in 1:N]
    x = collect(LinRange(0.8,1.6, N_Trial))*r
    Data.id=3
    Data.trial_positions = [Vector3([1.0,1.0,1.0])*1/sqrt(3)*x[i] for i in 1:N_Trial]
    RosenbluthChains.GetTrialBoltzmannWeight(Data,KP.SAWParam)
    fig=plot(x, -Data.LogBoltzmannFaktor, label="Sim result")
    plot!(x, LJ.(x, r, 3.0), label="theory", linestyle=:dot, ylim=(-0.1, 40))
    savefig(fig, "./tmp/LJ_Pot.pdf")
    @test Data.LogBoltzmannFaktor ≈  -LJ.(x, r, 3.0)

   
    N_Trial=8
    N_Batch=10_000
    #N=100
    N_Vals = collect(10:25:210)
    μ_val= zeros(length(N_Vals))
    for (runID,N) in enumerate(N_Vals)
        Data = SimData("./tmp/", 1.0, N, N_Trial, N_Batch, 1)
        KP = SimulationParameters( FixedBondParameters(r), RandBondAngle(), RandTorsion(), LJ_Repulsion(ones(Int64,N),Dict(1 => r), 1.0 ))
        Meas = RunSim(Data,KP, AMeasurement(N, N_Batch));
        #println(Meas.REE)
        #println(Meas.Weight)
        #println(Meas.SumRosenbluthWeights)
        μ = sum(Meas.REE.*Meas.Weight)/Meas.SumRosenbluthWeights
        σ   = sqrt(sum( Meas.Weight.*(μ.-Meas.REE).^2)/Meas.SumRosenbluthWeights)
        Δμ = σ/sqrt((N_Batch-1))

        μ_val[runID] = μ
    end

    fit = curve_fit(REE_Fit, N_Vals, μ_val, [r, 0.58])

    fig = plot(N_Vals, μ_val, yaxis=:log, xaxis=:log, label= "data");
    
    plot!(N_Vals, REE_Fit(N_Vals, fit.param), label="fit");
    savefig(fig, "./tmp/Scaling.pdf");

    (ΔN, Δν) = standard_errors(fit)
    (N_fit, ν_fit) = fit.param


    @test 0.58 < ν_fit+3.0*Δν && 0.58 > ν_fit-3.0*Δν && Δν<0.10

    println("ν= $(ν_fit)±$(Δν) vs. 0.58...")
    
end;




