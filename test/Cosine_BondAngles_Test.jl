using Printf, Plots, StatsBase



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

    ### test whether we can solve for prefactor in the energy potential to get the accurate persistence length on average
    ### K function first
    K_target= RosenbluthChains.solveRecursive.(RosenbluthChains.lpToCosAngle.(μ),RosenbluthChains.AvgCos , 0.01,200.0 ) ### solve for prefactor
    KP = SimulationParameters( FixedBondParameters(r), GaussianK_Cosine_BondAngles(K_target, σ), FixedTorsionAngles(φ), IdealChain()) ### run simulation and measure

    tuple = [RosenbluthChains.give_rand(KP.BondAngleParam.Sampler) for i in 1:N_Rand] ### draw random number according to prefactor K
    K = getindex.(tuple, 1) ### 

    μ_K = sum(K)/N_Rand ### compute average
    σ_K = sqrt(sum( (μ_K.-K).^2)/(N_Rand-1)) ### compute the variance
    Δμ_K = σ_K/sqrt(N_Rand) ### compute error

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
    @test K_target<μ_K + 3.0*Δμ_K && K_target>μ_K - 3.0*Δμ_K 
    @test (σ_th-σ_K)/σ_th<0.0025

    if VERBOSE
        println("input: μ $μ, σ $σ")
        println("K-theory: K: $μ_K, K_th $(K_target) σ_th $(σ_th) σ_K $(σ_K)")
        println("cos(φ(Κ)): $( RosenbluthChains.lpToCosAngle.(μ)), cos(φ_test): $(RosenbluthChains.AvgCos(K_target))")
        println("lp(cos(φ)): $(-3.8/log(RosenbluthChains.AvgCos(K_target)))vs. $(μ)")
    end

    RosenbluthChains.lpToCosAngle.(μ),RosenbluthChains.AvgCos 
    ### simulations to mesure persistence length with gaussian K
    Meas = RunSim(Data,KP, AMeasurement(N, N_Test));
    AvgCos=sum([Meas.AvgCosBondAngle[i,:].*Meas.Weight[i] for i in 1:N_Test])./Meas.SumRosenbluthWeights#[1:end]
    lp = AvgCos



    lp = -3.8./log.(AvgCos)
    μ_test = sum(lp)/(N-2) ### estimate mean 
    σ_test = sqrt(sum((lp.-μ_test).^2)/((N-3))) ### estimate variance
    Δμ_test = σ_test/sqrt(N-3) ### estimate error of mean

    ### estimate error of variance
    μ4 =  sum((lp.-μ_test).^4)/(N-2) #https://stats.stackexchange.com/questions/156518/what-is-the-standard-error-of-the-sample-standard-deviation
    n =Float64(N-3)
    Δσ_test = sqrt(1/n*(μ4-(n-3.0)/(n-1.0)*σ_test^4))/(2.0*σ_test)

    if VERBOSE
        println("smaller than 0: $(sum([1 for x in AvgCos if x<0]))/$(length(AvgCos))")
        println("K: wanted μ(lp=$(Float32(μ))) estimated: $(Float32(μ_test))±$(Float32(Δμ_test)), σ(μ(lp))  $(Float32(σ_test))±$(Float32(Δσ_test))")
    end
    if max_width>=5 ### test doesnt work if K is cutoff
        @test μ<μ_test + 3.0*Δμ_test && μ>μ_test - 3.0*Δμ_test
    end



    ### NOW Gaussian LP Distributions
    KP = SimulationParameters( FixedBondParameters(r), GaussianLp_Cosine_BondAngles(μ, σ), FixedTorsionAngles(φ), IdealChain()) ### redefine force field

    ### draw random numbers as in simulation
    tuple = [RosenbluthChains.give_rand(KP.BondAngleParam.Sampler) for i in 1:N_Rand] 
    lp = getindex.(tuple, 1)

    K_target = RosenbluthChains.solveRecursive.(RosenbluthChains.lpToCosAngle.(μ),RosenbluthChains.AvgCos , 0.0001,200.0 ) ### compute target value of K for wanted persistence length μ


    ### compute angles
    angles = [ x for (k,x) in tuple if (k>μ && k<μ*fac)]  ### only angles for which k is close to the mean.
    angles2 = getindex.(tuple, 2) ### take all angles


    hist =  fit(Histogram, angles, ψ) ### computes histogram
    StatsBase.normalize(hist, mode=:density)
    new_hist =hist.weights ./ sum(hist.weights*(ψ[2]-ψ[1]))
    val= RosenbluthChains.Cosine_BondAngle_func.(K_target*(1.0+fac)/2.0, bin_center)
    theory_hist=val./sum(val*0.01)

    μ_lp = sum(lp)/N_Rand ### mean persistence length
    σ_lp = sqrt(sum( (μ_lp.-lp).^2)/(N_Rand-1)) ### variance persistence length
    Δμ_lp = σ_lp/(N_Rand) ### standard error of persistence length

    max_width =min(5.0,μ/(σ*1.0001)) ### compute cutoff width
    x = LinRange(μ-max_width*σ, μ+max_width*σ, 1000) ### valus in [μ-w*σ, μ+w*σ]
    vals = Normal.(μ, σ, x) ### compute gaussian
    cutOffGaussian= vals./(sum(vals)*(2.0*max_width*σ/1000.0))### reweigh with cutoff effects
    σ_th = sqrt(sum( cutOffGaussian .*(μ.- x).^2)/(sum(cutOffGaussian))) ### compute theoretial variance

    if DOPLOTS
        fig = histogram(lp, normalized=true, label="histogram")
        plot!(x,cutOffGaussian, label="Theory. renormed\n cutoff Gaussian")
        savefig(fig, "./tmp/$(μ)_$(σ)_lp_test.pdf")
    end

    @test abs(μ-μ_lp)/μ<10^-3
    @test (σ_th-σ_lp)/σ_th<0.001 

    if DOPLOTS
        fig = histogram(angles, normalized=true, label="histogram")
        histogram!(angles2, normalized=true, label="histogram all angles" )
        plot!(bin_center, new_hist, label="Statsbase Hist")
        plot!(bin_center, theory_hist, label="theory")
        savefig(fig, "./tmp/$(μ)_$(σ)_lp_angle_test.pdf")
    end

    @test ComputeKullbackLeiblerDivergence(new_hist,theory_hist) <0.1



    ### Run simulation to test results
    Meas = RunSim(Data,KP, AMeasurement(N, N_Test));
    AvgCos=sum([Meas.AvgCosBondAngle[i,:].*Meas.Weight[i] for i in 1:N_Test])./Meas.SumRosenbluthWeights

    ### compute mean and standard error of average persistence lengths
    lp = -3.8./log.(AvgCos)
    μ_test = sum(lp)/(N-2)
    σ_test = sqrt(sum((lp.-μ_test).^2)/((N-3)))
    Δμ_test = σ_test/sqrt(N-3)

    ### compute mean and standard error of average variances of persistence lengths
    ### cannot be done since some angles of the potentials cause cos(angle)<0 ⇒ -3.8/log(cos(angle)); is not defined comparison to distribution has to be sufficient
    #=
    μs = ones(N-2)*μ_test
    println("μs: $(Float32.(μs))")
    println("data: $(Float32.(Meas.AvgCosBondAngle[1:30]))")
    σ_arr =sum([(-3.8./log.(Meas.AvgCosBondAngle[i,:]).-μs).^2 .*Meas.Weight[i] for i in 1:N_Test])./Meas.SumRosenbluthWeights ### estimate σ^2
    println(Float32.(σ_arr))

    σ_arr = sqrt.(σ_arr) ### estimate of σ

    μ_σ = sum(σ_arr)/(N-2)  ### mean of σ's
    σ_σ = sqrt(sum((σ_arr.-μ_σ).^2)/((N-3))) ###estimate variance of σ's
    Δμ_σ = σ_test/sqrt(N-3) ### standard error of mean of σ's
 =#
    if VERBOSE
        println("LP: wanted μ: $(Float32(μ)), σ: $(Float32(σ_th)) estimated: μ=$(Float32(μ_test))±$(Float32(Δμ_test))")
        #println(" σ =$(Float32(μ_σ))±$(Float32(Δμ_σ))")
        println("\n\n\n")
    end
   

    @test μ<μ_test + 3.0*Δμ_test && μ>μ_test - 3.0*Δμ_test
    

end
if VERBOSE
    println("The GaussianLP_Cosine_BondAngles tend to fail for broad distributions.") ### @TODO: Fix GaussianLP_Cosine_BondAngles
end

end

