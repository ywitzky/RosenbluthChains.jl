using StatsBase,Plots



@testset "random bond angles" begin
    N= 200_000
    θ=collect(0.01:0.01:π)
    η=collect(0.01:0.01:π)

    ϕ=collect(0.01:0.01:2*π)

    x = collect(0:0.01:1)

    Data = SimData("$TestPath/tmp/", 1.0, N, 1      , 1, 1)
    KP = SimulationParameters( FixedBondParameters(1.0), RandBondAngle(), RandTorsion(), IdealChain())
    Result = RunSim(Data,KP, NoMeasurement());

    vecs = [Data.xyz[i+1]-Data.xyz[i] for i in 1:N-1]
    angles = [ angle(vecs[i],-1*vecs[i+1]) for i in 1:N-2   ]
    cross = [ ×(vecs[i],vecs[i+1]) for i in 1:N-2]

    orthoangle = [angle(vecs[i], vecs[i+1]-vecs[i]*(vecs[i+1]*vecs[i]))  for i in 1:N-2 ]
    orthoangle2 = [angle(cross[i], vecs[i+1]-vecs[i]*(vecs[i+1]*vecs[i]))  for i in 1:N-2 ]
    torsion_vec = [ angle(cross[i],cross[i+1]) for i in 1:N-3]
    torsion2 = [torsion(vecs[i], vecs[i+1], vecs[i+2]) for i in 1:N-3 ]


    angles_hist = fit(Histogram, angles, θ)
    StatsBase.normalize(angles_hist, mode=:density)
    angles_hist =angles_hist.weights ./ sum(angles_hist.weights*(θ[2]-θ[1]))
    θ = (θ[2:end].+θ[1:end-1])/2.0

    theory = 0.5.*sin.(θ)
    theory ./= (sum(theory)*(θ[2]-θ[1]))
    if DOPLOTS
        fig = Plots.bar(θ, angles_hist, label="rand bond angle")
        plot!(θ, theory, label="theory")
        savefig(fig, "$TestPath/tmp/RandomBondAngle.pdf");
    end

    
    @test ComputeKullbackLeiblerDivergence(angles_hist,theory) <0.5 ### we can lower that value if we take more samples


    torsion_hist = fit(Histogram, torsion_vec, η)
    StatsBase.normalize(torsion_hist, mode=:density)
    torsion_hist =torsion_hist.weights ./ sum(torsion_hist.weights*(η[2]-η[1]))
    torsion2_hist = fit(Histogram, torsion2, η)
    StatsBase.normalize(torsion2_hist, mode=:density)
    torsion2_hist =torsion2_hist.weights ./ sum(torsion2_hist.weights*(η[2]-η[1]))

    η = (η[2:end].+η[1:end-1])/2.0

    theory = η.*0.0.+1.0./(π)

    if DOPLOTS
        fig = Plots.bar(η, torsion_hist, label="rand cross bond angle")
        plot!(η, torsion2_hist, label="torsion2")
        plot!(η, theory, label="theory")
        savefig(fig, "$TestPath/tmp/RandomTorsionAngle.pdf");

        fig = Plots.scatter(angles[1:10_000], torsion2[1:10_000], xlabel="polar angle", ylabel="azimuth angle\ntorsion angle", label="")
        savefig(fig, "$TestPath/tmp/RandomBondAngleTorsionAngleCorrelation.pdf");
    end

    @test ComputeKullbackLeiblerDivergence(torsion_hist,theory) <0.5

    ### test angles build orthogonal vectors as intended
    @test all(orthoangle .≈ π/2)
    @test all(orthoangle2 .≈ π/2) 
end;

using Plots

### compute root mean square RG and REE values for ideal chains and compare to theory
@testset "ideal chain Scalings" begin
for (b, N) in zip([1.0, 1.0,0.5, 3.0], [20, 100, 200, 300])
    N_Batch=50_000
    Data = SimData("../tmp/", 1.0, N, 8, N_Batch, 1)

    Model = SimulationParameters( FixedBondParameters(b), RandBondAngle(), RandTorsion(), IdealChain())

    Result = RunSim(Data, Model, RosenbluthChains.RG_Measurement(Data.FolderPath, N_Batch))

    ### numerical stability
    Result.Weights ./= maximum(Result.Weights)
    Total_Weights = sum(Result.Weights)

    RG_avg, RG_err = ComputeSqrtMeanError(Result.RGs, Result.Weights ; NIntervals =100)
    REE_avg, REE_err = ComputeSqrtMeanError(Result.REEs, Result.Weights; NIntervals =100)

    RG_exp = sqrt(b^2*(N)/6.0)
    REE_exp = sqrt(b^2*(N-1))

    Result.RGs .= sqrt.( Result.RGs)
    Result.REEs .= sqrt.( Result.REEs)
    bins = collect(0:1:ceil(Int32,maximum(Result.REEs)))

    REE_hist = fit(Histogram, Result.REEs, weights(Result.Weights), bins)
    StatsBase.normalize(REE_hist, mode=:density)
    REE_hist =REE_hist.weights ./ sum(REE_hist.weights*(bins[2]-bins[1]))
    bins = (bins[2:end].+bins[1:end-1])/2.0
    
    theory_hist = @. sqrt(3/(2*π*N*b^2))^3*exp(-3.0*bins^2/(2.0*N*b^2))*bins^2
    theory_hist ./= sum(theory_hist)*(bins[2]-bins[1])

    if DOPLOTS
        fig=plot(bins, REE_hist, label="Sim result", xlabel="REE", ylabel="P(REE)")
        plot!(bins,theory_hist, label="theory", linestyle=:dot)
        savefig(fig, "$TestPath/tmp/IdealChainGyrationradius_$(b)_$(N).pdf")
    end

    @test ComputeKullbackLeiblerDivergence(REE_hist,theory_hist) <0.1

    if VERBOSE
        println("RG $(RG_avg) ± $(RG_err): theory $(RG_exp)  rel :$(RG_avg/RG_exp)")
        println("REE $(REE_avg) ± $(REE_err): theory $(REE_exp) rel :$(REE_avg/REE_exp)")
    end

    @test abs(RG_exp - RG_avg)<3.0*RG_err && RG_err<0.1
    @test abs(REE_exp - REE_avg)<3.0*REE_err && REE_err<0.5

    RosenbluthChains.close(Result)

end
end