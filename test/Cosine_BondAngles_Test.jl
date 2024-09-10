using Printf, Plots, StatsBase,Distributions



N_Trial=8
N_Test=40_000
N= 300
N_cut=15
φ=ones(N)*1#collect(range(0.01,π, N-2))
r=3.8
μ=3.0
σ=0.1
Data = SimData("./tmp/", 1.0, N, N_Trial      , N_Test, 1)



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
    



N_Manual=10_000_000


@testset "Cosine_BondAngles & GaussianLp_Cosine_BondAngles " begin

N_Rand = 10_000_000
fac=1.05    

ψ = collect(0:0.01:π)
bin_center = (ψ[2:end].+ψ[1:end-1])/2.0
Normal(μ, σ,x) = 1/sqrt(2pi*σ^2) * exp(-(μ-x)^2/(2*σ^2))

pers(κ, b=3.8) =  -b/log(coth(κ)-1.0/κ) ### lp(κ)


for (μ, σ) in [ (3.0, 0.1), (3.0, 0.2),(3.0, 1.0),(3.0, 1.5), (50.0, 1.0), (500.0, 50.0)]

    ### test whether we can solve for prefactor in the energy potential to get the accurate persistence length on average
    ### K function first
    K_target= RosenbluthChains.solveRecursive.(RosenbluthChains.lpToCosAngle.(μ),RosenbluthChains.AvgCos , 0.01,200.0 ) ### solve for prefactor
    KP = SimulationParameters( FixedBondParameters(r), GaussianK_Cosine_BondAngles(K_target, σ), FixedTorsionAngles(φ), IdealChain()) ### run simulation and measure

    tuple = [RosenbluthChains.give_rand(KP.BondAngleParam.Sampler) for i in 1:N_Rand] ### draw random number according to prefactor K
    K = getindex.(tuple, 1) ### 
    
    μ_K_est, Δμ_K_est, σ_K_est, Δσ_K_est = ComputeMeanAndVariance(K)

    angles = [ x for (k,x) in tuple if (k>K_target && k<K_target*fac)] ### compute angles in simulation

    hist =  fit(Histogram, angles, ψ) ### computes histogram of angles
    StatsBase.normalize(hist, mode=:density)
    new_hist =hist.weights ./ sum(hist.weights*(ψ[2]-ψ[1])) ### reweight
    val= RosenbluthChains.Cosine_BondAngle_func.(K_target*(1.0+fac)/2.0, bin_center) ### analytical solution to histogram
    theory_hist=val./sum(val*0.01) ### reweight

    @test ComputeKullbackLeiblerDivergence(new_hist,theory_hist) <0.1 ### close to 0 means very similar distributions

    ### visualise histogram
    if DOPLOTS
        fig = histogram(angles, normalized=true, label="histogram") 
        plot!(bin_center, new_hist, label="Statsbase Hist")
        plot!(bin_center, theory_hist, label="theory")
        savefig(fig, "./tmp/$(μ)_$(σ)_angle_test.pdf")
    end

    ### set limits to the sampling range of K: here 5 σ are choosen
    max_width =min(5.0,K_target/(σ*1.0001))
    x = LinRange(K_target-max_width*σ, K_target+max_width*σ, 1000) ### get values in [μ-5σ, μ+5σ]
    vals = Normal.(K_target, σ, x) ### compute normal distribution
    cutOffGaussian= vals./(sum(vals)*(2.0*max_width*σ/1000.0)) ### renormalise since we cut off values outside [μ-5σ, μ+5σ]
    σ_th = sqrt(sum( cutOffGaussian .*(K_target.- x).^2)/(sum(cutOffGaussian))) ### compute variance according to this distribution


    ### visualise histogram of K values vs theory histogram
    if DOPLOTS
        fig=histogram(K, normalized=true, label="hist") 
        plot!(x, cutOffGaussian, label="cutoff Gaussian")
        savefig(fig, "./tmp/$(μ)_$(σ)_K_hist.pdf")
    end

    ### compare theory and simulation for K and and its variance
    @test K_target<μ_K_est + 3.0*Δμ_K_est && K_target>μ_K_est - 3.0*Δμ_K_est 
    @test (σ_th<σ_K_est + 3.0*Δσ_K_est && σ_th>σ_K_est - 3.0*Δσ_K_est) || (σ_th-σ_K_est)/σ_th < 10^-2

    
    if VERBOSE
        println("input: μ $μ, σ $σ")
        println("K-theory: K: $μ_K_est±$(Δμ_K_est), K_th $(K_target) σ_th $(σ_th) σ_K $(σ_K_est)±$(Δσ_K_est)")
        println("cos(φ(Κ)): $( RosenbluthChains.lpToCosAngle.(μ)), cos(φ_test): $(RosenbluthChains.AvgCos(K_target))")
        println("lp(cos(φ)): $(-3.8/log(RosenbluthChains.AvgCos(K_target)))vs. $(μ)")
    end
    
    
    ### NOW Gaussian LP Distributions
    KP = SimulationParameters( FixedBondParameters(r), GaussianLp_Cosine_BondAngles(μ, σ, N), FixedTorsionAngles(φ), IdealChain()) ### redefine force field

    ### test the gaussian initialisation of the κ values
    NT = 10_000
    κ_left, κ_right = RosenbluthChains.getKappaBorders(μ,σ)
    left = pers(κ_left)
    right = pers(κ_right)

    kappas = GaussianLp_Cosine_BondAngles(μ,σ, NT).κ_vec[1:end-2]
    lps = pers.(kappas, 3.8)


    xrange = collect(LinRange( left,right, 50))

    hist =  fit(Histogram, lps, xrange) ### computes histogram
    StatsBase.normalize(hist, mode=:density)
    hist =hist.weights ./ sum(hist.weights*(xrange[2]-xrange[1]))

    xrange = (xrange[1:end-1]+xrange[2:end])./2.0

    theory = Distributions.pdf.(Distributions.Normal(μ, σ),xrange )
    theory ./= sum(theory)*(maximum(xrange)-minimum(xrange))/length(xrange)
    σ_th = sqrt(sum( theory .*(μ.- xrange).^2)/(sum(theory))) ### compute variance according to the cutoff normal distribution
    

    μ_est, Δμ_est, σ_est, Δσ_est = ComputeMeanAndVariance(lps)



    if DOPLOTS
        fig = histogram(pers.(kappas, 3.8), normalized=true, label="histogram", xlabel="lp", ylabel="P(lp)")
        plot!(xrange, hist, ylabel="hist2")
        plot!(xrange,theory, label="Theory")
        savefig(fig, "./tmp/fixed_kappa_$(μ)_$(σ)_gaussian_lp_test.pdf")
    end

    @test μ<μ_est + 3.0*Δμ_est && μ>μ_est - 3.0*Δμ_est
    @test σ_th<σ_est + 3.0*Δσ_est && σ_th>σ_est - 3.0*Δσ_est

    @test ComputeKullbackLeiblerDivergence(hist,theory) < 0.5 ### could be decreased if NT is larger ⇒ long test times



    ### Run simulation to test results
    Meas = RunSim(Data,KP, AMeasurement(N, N_Test));
    AvgCos=sum([Meas.AvgCosBondAngle[i,:].*Meas.Weight[i] for i in 1:N_Test])./Meas.SumRosenbluthWeights

    ### compute mean and standard error of average persistence lengths
    lp = -3.8./log.(AvgCos)

    μ_est, Δμ_est, σ_est, Δσ_est = ComputeMeanAndVariance(lp)
    
    if VERBOSE
        println("LP: wanted μ: $(Float32(μ)),  estimated: μ=$(Float32(μ_est))±$(Float32(Δμ_est)), ")
        println("LP: wanted σ: $(Float32(σ)), possible $(σ_th), estimated σ=$(Float32(σ_est))±$(Float32(Δσ_est))")
        println("\n\n\n")
    end
   
    @test μ<μ_est + 3.0*Δμ_est && μ>μ_est - 3.0*Δμ_est
    @test σ_th<σ_est + 3.0*Δσ_est && σ_th>σ_est - 3.0*Δσ_est
end
  
end
