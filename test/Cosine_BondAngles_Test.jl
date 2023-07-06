using Plots

mutable struct AMeasurement{T<:Number} <: AbstractMeasurement 
    SumRosenbluthWeights::T
    AvgCosBondAngle::Vector{T}

    AMeasurement(NBeads::I) where {I<:Integer} = new{Float64}(zero(Float64), zeros(Float64, NBeads-2))
end


function RosenbluthChains.InitMeasurement(data::SimData, param::SimulationParameters, Tmp::AMeasurement) 
    return AMeasurement(data.NBeads)
end 

function RosenbluthChains.MeasureAfterBatch(data::SimData, param::SimulationParameters,Measurement::AMeasurement) 

end

function RosenbluthChains.MeasureAfterChainGrowth(data::SimData, param::SimulationParameters, Measurement::AMeasurement) 
    rel = [data.xyz[i+1]-data.xyz[i] for i in 1:data.NBeads-1]
    angles=[angle( rel[i],rel[i+1]) for i in 1:data.NBeads-2]
    Measurement.AvgCosBondAngle .+= getRosenbluthWeigth(data, param).*cos.(angles).^2
    #print(data.id_in_batch,"  ")
    println(data.id_in_batch,"  ", getRosenbluthWeigth(data, param))#,"   ", data.LogRosenbluthWeight, "   ",  param.BondAngleParam.K[data.tid]*(1+cos(data.trial_angle[data.tid])))
    #println( getRosenbluthWeigth(data, param).* [angle( rel[i],rel[i+1]) for i in 1:N-2])
    Measurement.SumRosenbluthWeights += getRosenbluthWeigth(data, param)#data.RosenbluthWeight
end
    


N_Trial=8
N= 10
N_cut=15
N_Test=100#0
φ=collect(range(0.01,π, N-2))
r=3.8
Data = SimData("./tmp/", 1.0, N, N_Trial      , N_Test, 1)

μ=100.0
σ=1.0
KP = SimulationParameters( FixedBondParameters(r), GaussianLp_Cosine_BondAngles(μ, σ), FixedTorsionAngles(φ), IdealChain())
@time Meas = RunSim(Data,KP, AMeasurement(N));

println(Meas.SumRosenbluthWeights)
println(extrema(Meas.AvgCosBondAngle))
AvgCos= Meas.AvgCosBondAngle[3:end]./Meas.SumRosenbluthWeights
sort!(AvgCos)
println(extrema(AvgCos))
lp = @. -1.0*r/log(AvgCos)#[N_cut:end-N_cut])
println(extrema(lp))
μ_test = sum(lp)/(N-2)
σ_test = sqrt(sum((lp.-μ_test).^2)/((N-2)*(N-3)))

println("σ:$μ vs: $μ_test")
println("σ:$σ vs: $σ_test")


lp = @. -1.0*r/log(AvgCos)

fig = scatter(AvgCos, lp, xlabel="<cos(x)>", ylabel="lp", yaxis=:log, label="")
savefig(fig,"./tmp.png")
