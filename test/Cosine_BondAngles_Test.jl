using Printf, Plots, StatsBase

function ComputeKullbackLeiblerDivergence(P,Q)
    return sum([p*log(p/q) for (p,q) in zip(P,Q) if p>0 ])
end

mutable struct AMeasurement{T<:Number} <: AbstractMeasurement 
    SumRosenbluthWeights::T
    AvgCosBondAngle::Array{T}
    Weight::Vector{T}
    AMeasurement(NBeads::I, NSteps::I) where {I<:Integer} = new{BigFloat}(zero(BigFloat ), zeros(BigFloat,(NSteps,  NBeads-2)), zeros(BigFloat, NSteps)) #where {T::BigFloat}
end


function RosenbluthChains.InitMeasurement(data::SimData, param::SimulationParameters, Tmp::AMeasurement) 
    return AMeasurement(data.NBeads, data.BatchSize)
end 

function RosenbluthChains.MeasureAfterBatch(data::SimData, param::SimulationParameters,Measurement::AMeasurement) 
    nothing
end

function RosenbluthChains.MeasureAfterChainGrowth(data::SimData, param::SimulationParameters, Measurement::AMeasurement) 
    rel = [data.xyz[i+1]-data.xyz[i] for i in 1:data.NBeads-1]
    angles=[angle( rel[i],rel[i+1]) for i in 1:data.NBeads-2]

    Measurement.AvgCosBondAngle[data.id_in_batch,:] .= BigFloat.(cos.(angles))
    Measurement.Weight[data.id_in_batch] = getRosenbluthWeigth(data, param)
    Measurement.SumRosenbluthWeights += data.RosenbluthWeight #getRosenbluthWeigth(data, param)#
end
    

N_Trial=8
N_Test=10_000
N= 100
N_cut=15
φ=ones(N)*1#collect(range(0.01,π, N-2))
r=3.8
μ=3.0
σ=0.1
Data = SimData("./tmp/", 1.0, N, N_Trial      , N_Test, 1)


N_Manual=10_000_000


@testset "Cosine_BondAngles & GaussianLp_Cosine_BondAngles " begin

N_Rand = 10_000_000
fac=1.05    

ψ = collect(0:0.01:π)
bin_center = (ψ[2:end].+ψ[1:end-1])/2.0
Normal(μ, σ,x) = 1/sqrt(2pi*σ^2) * exp(-(μ-x)^2/(2*σ^2))

for (μ, σ) in [(50.0, 5.0),(50.0, 10.0), (3.0, 0.1), (3.0, 0.2), (3.0, 1.0),(3.0, 1.5), (50.0, 1.0), (500.0, 5.0)]
    ### K function first
    K_target= RosenbluthChains.solveRecursive.(RosenbluthChains.lpToCosAngle.(μ),RosenbluthChains.AvgCos , 0.01,200.0 )
    KP = SimulationParameters( FixedBondParameters(r), GaussianK_Cosine_BondAngles(K_target, σ), FixedTorsionAngles(φ), IdealChain())

    tuple = [RosenbluthChains.give_rand(KP.BondAngleParam.Sampler) for i in 1:N_Rand]
    K = getindex.(tuple, 1)

    μ_K = sum(K)/N_Rand
    σ_K = sqrt(sum( (μ_K.-K).^2)/(N_Rand-1))
    Δμ_K = σ_K/sqrt(N_Rand)

    angles = [ x for (k,x) in tuple if (k>K_target && k<K_target*fac)]

    hist =  fit(Histogram, angles, ψ) ### computes histogram
    StatsBase.normalize(hist, mode=:density)
    new_hist =hist.weights ./ sum(hist.weights*(ψ[2]-ψ[1]))
    val= RosenbluthChains.Cosine_BondAngle_func.(K_target*(1.0+fac)/2.0, bin_center)
    theory_hist=val./sum(val*0.01)

    @test ComputeKullbackLeiblerDivergence(new_hist,theory_hist) <0.1 ### close to 0 means very similar distributions

    fig = histogram(angles, normalized=true, label="histogram")
    plot!(bin_center, new_hist, label="Statsbase Hist")
    plot!(bin_center, theory_hist, label="theory")
    savefig(fig, "./tmp/$(μ)_$(σ)_angle_test.pdf")







    max_width =min(5.0,K_target/(σ*1.0001))
    x = LinRange(K_target-max_width*σ, K_target+max_width*σ, 1000)
    vals = Normal.(K_target, σ, x)
    cutOffGaussian= vals./(sum(vals)*(2.0*max_width*σ/1000.0))
    σ_th = sqrt(sum( cutOffGaussian .*(K_target.- x).^2)/(sum(cutOffGaussian)))
    #println("wanted σ: $(σ) estimated: $(σ_K), theory prediction for cutoff  $(σ_th) testmeasure: $((σ_th-σ_K)/σ_th)")




    fig=histogram(K, normalized=true, label="hist")
    plot!(x, cutOffGaussian, label="cutoff Gaussian")
    savefig(fig, "./tmp/$(μ)_$(σ)_K_hist.pdf")


    @test K_target<μ_K + 3.0*Δμ_K && K_target>μ_K - 3.0*Δμ_K
    @test (σ_th-σ_K)/σ_th<0.0025

    
    
    Meas = RunSim(Data,KP, AMeasurement(N, N_Test));
    AvgCos=sum([Meas.AvgCosBondAngle[i,:].*Meas.Weight[i] for i in 1:N_Test])./Meas.SumRosenbluthWeights#[1:end]
    lp = AvgCos


    lp = -3.8./log.(AvgCos)
    μ_test = sum(lp)/(N-2)
    σ_test = sqrt(sum((lp.-μ_test).^2)/((N-3)))

    println("wanted μ: $(μ) estimated: $(μ_test), σ  $(σ_test)")
    if max_width>=5 ### test doesnt work if K is cutoff
        @test μ<μ_test + 3.0*σ_test && μ>μ_test - 3.0*σ_test
    end
                

    ### lp func
    KP = SimulationParameters( FixedBondParameters(r), GaussianLp_Cosine_BondAngles(μ, σ), FixedTorsionAngles(φ), IdealChain())

    tuple = [RosenbluthChains.give_rand(KP.BondAngleParam.Sampler) for i in 1:N_Rand]
    lp = getindex.(tuple, 1)

    K_target = RosenbluthChains.solveRecursive.(RosenbluthChains.lpToCosAngle.(μ),RosenbluthChains.AvgCos , 0.0001,200.0 )


    angles = [ x for (k,x) in tuple if (k>μ && k<μ*fac)]
    angles2 = getindex.(tuple, 2)


    hist =  fit(Histogram, angles, ψ) ### computes histogram
    StatsBase.normalize(hist, mode=:density)
    new_hist =hist.weights ./ sum(hist.weights*(ψ[2]-ψ[1]))
    val= RosenbluthChains.Cosine_BondAngle_func.(K_target*(1.0+fac)/2.0, bin_center)
    theory_hist=val./sum(val*0.01)

    μ_lp = sum(lp)/N_Rand
    σ_lp = sqrt(sum( (μ_lp.-lp).^2)/(N_Rand-1))
    Δμ_lp = σ_lp/(N_Rand)

    fig = histogram(lp, normalized=true, label="histogram")
    max_width =min(5.0,μ/(σ*1.0001))
    x = LinRange(μ-max_width*σ, μ+max_width*σ, 1000)
    vals = Normal.(μ, σ, x)
    cutOffGaussian= vals./(sum(vals)*(2.0*max_width*σ/1000.0))
    σ_th = sqrt(sum( cutOffGaussian .*(μ.- x).^2)/(sum(cutOffGaussian)))

    plot!(x,cutOffGaussian, label="Theory. renormed\n cutoff Gaussian")
    #plot!(bin_center, theory_hist, label="theory")
    savefig(fig, "./tmp/$(μ)_$(σ)_lp_test.pdf")


    @test abs(μ-μ_lp)/μ<10^-3
    @test (σ_th-σ_lp)/σ_th<0.001 
    #println("wanted σ: $(σ) estimated: $(σ_lp), theory prediction for cutoff  $(σ_th) testmeasure: $((σ_th-σ_lp)/σ)")
    fig = histogram(angles, normalized=true, label="histogram")
    histogram!(angles2, normalized=true, label="histogram all angles" )
    plot!(bin_center, new_hist, label="Statsbase Hist")
    plot!(bin_center, theory_hist, label="theory")
    savefig(fig, "./tmp/$(μ)_$(σ)_lp_angle_test.pdf")

    @test ComputeKullbackLeiblerDivergence(new_hist,theory_hist) <0.1



                
    Meas = RunSim(Data,KP, AMeasurement(N, N_Test));
    AvgCos=sum([Meas.AvgCosBondAngle[i,:].*Meas.Weight[i] for i in 1:N_Test])./Meas.SumRosenbluthWeights
    lp = AvgCos

    lp = -3.8./log.(AvgCos)
    μ_test = sum(lp)/(N-2)
    σ_test = sqrt(sum((lp.-μ_test).^2)/((N-3)))
    #println("μ:$μ vs: $μ_test")
    #println("σ: vs: $σ_test")
    #println(" $(extrema(lp)) vs $(μ)")
    println("wanted μ: $(μ) estimated: $(μ_test), σ  $(σ_test)")

    @test μ<μ_test + 3.0*σ_test && μ>μ_test - 3.0*σ_test

end
end