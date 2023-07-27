using Printf

mutable struct AMeasurement{T<:Number} <: AbstractMeasurement 
    SumRosenbluthWeights::T
    AvgCosBondAngle::Vector{T}
    AvgLp::Vector{T}
    OldWeight::T
    AMeasurement(NBeads::I) where {I<:Integer} = new{Float32}(zero(Float32 ), zeros(Float32, NBeads-2), zeros(Float32, NBeads-2), zero(Float32)) #where {T::BigFloat}
end


function RosenbluthChains.InitMeasurement(data::SimData, param::SimulationParameters, Tmp::AMeasurement) 
    return AMeasurement(data.NBeads)
end 

function RosenbluthChains.MeasureAfterBatch(data::SimData, param::SimulationParameters,Measurement::AMeasurement) 
BatchWeigth= Measurement.SumRosenbluthWeights- Measurement.OldWeight
     Measurement.AvgCosBondAngle ./=  BatchWeigth
     Measurement.AvgLp .+= -3.8.\log.(Measurement.AvgCosBondAngle) .* (BatchWeigth)
     Measurement.OldWeight = Measurement.SumRosenbluthWeights
     fill!(Measurement.AvgCosBondAngle, 0.0)
end

function RosenbluthChains.MeasureAfterChainGrowth(data::SimData, param::SimulationParameters, Measurement::AMeasurement) 
    rel = [data.xyz[i+1]-data.xyz[i] for i in 1:data.NBeads-1]
    angles=[angle( rel[i],rel[i+1]) for i in 1:data.NBeads-2]
    Measurement.AvgCosBondAngle .+= getRosenbluthWeigth(data, param).*cos.(angles)#.*sin.(angles)
    #print(data.id_in_batch,"  ")
    #=
    println(getRosenbluthWeigth(data, param))
    println(data.xyz[1:10])
    println(rel[1:10])
    println(angles[1:10])
    =#
    #println(data.id_in_batch,"  ", getRosenbluthWeigth(data, param))#,"   ", data.LogRosenbluthWeight, "   ",  param.BondAngleParam.K[data.tid]*(1+cos(data.trial_angle[data.tid])))
    #println( getRosenbluthWeigth(data, param).* [angle( rel[i],rel[i+1]) for i in 1:N-2])
    Measurement.SumRosenbluthWeights += data.RosenbluthWeight #getRosenbluthWeigth(data, param)#
end
    
angles=collect(range(0,π/2.0,10_000))

N_Trial=8
N= 60
N_cut=15
N_Test=500_000
φ=ones(N)*2#collect(range(0.01,π, N-2))
r=3.8
μ=100.0
σ=1.0
Data = SimData("./tmp/", 1.0, N, N_Trial      , N_Test, 1)
KP = SimulationParameters( FixedBondParameters(r), Cosine_BondAngles(ones(Float64, N)*μ), FixedTorsionAngles(φ), IdealChain())
@time Meas = RunSim(Data,KP, AMeasurement(N));

f(K,x) = @. cos(x)*sin(x)*exp(-1.0*K*(1.0-cos(x)))
h(K, y)=@. exp(-1.0*K*(1.0-cos(y))) ### P(Angle)
norm(x) = x./sum(x)


N_Manual=10_000_000


@testset "Cosine_BondAngles & GaussianLp_Cosine_BondAngles " begin
    #println(RosenbluthChains.AvgCos(2.0))
    
    @test (RosenbluthChains.AvgCos(2.0) - 0.537315) <10^-5

    Data.trial_angle .= collect(range(0,π,N_Trial))
    RosenbluthChains.CompTrigonometricTrialBondAngles(Data)
    RosenbluthChains.CompTrigonometricTrialTorsionAngles(Data)

    RosenbluthChains.ResetSim(Data, KP)
    RosenbluthChains.GetTrialBoltzmannWeight(Data,KP)
    @time choose_ids = [begin ChooseTrialPosition(Data, KP); Data.tid; end for _ in 1:N_Manual]
    @time hist = [ count(x->x==(a), choose_ids) for a in 1:N_Trial]./N_Manual
    println(hist)
    println(norm(h(KP.BondAngleParam.K[1], Data.trial_angle)))
    @test all( abs.(hist .-norm(h(KP.BondAngleParam.K[1], Data.trial_angle))) .< 10^-3)

    cos_x = hist.*cos.(Data.trial_angle).*sin.(Data.trial_angle)
    @test all(abs.(norm(cos_x) .-norm(f.(KP.BondAngleParam.K[1], Data.trial_angle))) .< 10^-3)

    
    AvgCos= Meas.AvgCosBondAngle[1:end]./Meas.SumRosenbluthWeights
    #sort!(AvgCos)
    println(Meas.SumRosenbluthWeights)
    println(RosenbluthChains.AvgCos(KP.BondAngleParam.K[1]))
    println(extrema(AvgCos))
    #println(AvgCos)
    println(sum(AvgCos)/(N-2),"  ",  -3.8/log(sum(AvgCos)/(N-2)))
    μ_test = -r/log(sum(AvgCos)/(N-2))
    lp = @. -r/log(AvgCos)#[N_cut:end-N_cut])

    σ_test = sqrt(sum((lp.-μ_test).^2)/((N-2)*(N-3)))
    println("σ:$μ vs: $μ_test")
    println("σ:$σ vs: $σ_test")
    println(" $(extrema(lp)) vs $(μ)")
    println(-3.8/log(RosenbluthChains.AvgCos(KP.BondAngleParam.K[1])))
    
    #fig = scatter(AvgCos, lp, xlabel="<cos(x)>", ylabel="lp", yaxis=:log, label="")
    #savefig(fig,"./tmp.png")

end;


#### Inverse Persistence length
#=
N_Trial=8
N= 15
N_cut=15
N_Test=500#_000
N_Batches=10
φ=collect(range(0.01,π, N-2))
r=3.8
μ=0.01
σ=0.001
Data = SimData("./tmp/", 1.0, N, N_Trial      , N_Test, N_Batches)
KP = SimulationParameters( FixedBondParameters(r), GaussianInvLp_Cosine_BondAngles(μ,σ), FixedTorsionAngles(φ), IdealChain())
@time Meas = RunSim(Data,KP, AMeasurement(N));

f(K,x) = @. cos(x)*sin(x)*exp(-1.0*K*(1.0-cos(x)))
h(K, y)=@. exp(-1.0*K*(1.0-cos(y))) ### P(Angle)
norm(x) = x./sum(x)


N_Manual=10_000_000
    
rand_num = rand(KP.BondAngleParam.InvPersistenceDist, 100_000) #1-0
bla =sum(rand_num)/100_000
println("Check Rand Dist: μ:$(bla), σ:$(sqrt(sum((rand_num.-bla).^2)/(99_999)))")

### maybe error i first 3
AvgCos= Meas.AvgCosBondAngle[3:end]./Meas.SumRosenbluthWeights
AvgLp= Meas.AvgLp[3:end]./Meas.SumRosenbluthWeights
#sort!(AvgCos)
println("Rosenbluth weight: ", Meas.SumRosenbluthWeights)
println("K=", KP.BondAngleParam.K[1], " <cos(K)>", RosenbluthChains.AvgCos(KP.BondAngleParam.K[1]),"  lp(K)=", -3.8./log.(RosenbluthChains.AvgCos.(KP.BondAngleParam.K[1])))
println("lp=", AvgLp)
println("Min/Max (<cos(x)>)", extrema(AvgCos))
#println(AvgCos)
println("AvgAvgCos $(sum(AvgCos)/(length(AvgCos)))   lp=$(-3.8/log( sum(AvgCos)/(length(AvgCos))))")
Invlp = @. -log(AvgCos)/r#[N_cut:end-N_cut])
μ_test = (sum(Invlp)/(length(Invlp)))
σ_test = sqrt(sum((Invlp.-μ_test).^2)/(length(Invlp)*(length(Invlp)-1)))
@printf "μ: %3.8f  vs: %3.8f \n" μ μ_test

@printf "σ: %3.8f vs: %3.8f \n" σ σ_test

println(" $(extrema(Invlp)) vs $(μ)")
=#