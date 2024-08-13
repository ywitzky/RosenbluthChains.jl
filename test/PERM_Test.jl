
using StatsBase,Plots,LsqFit

using Base.Iterators: partition


#=
@testset "PERM: ideal chain Scalings" begin
    for (b, N) in zip([1.0,0.5, 3.0], [100, 200, 300])
        N_Batch=50_000
        Data = SimData("../tmp/", 1.0, N, 8, N_Batch, 1)

        Model = SimulationParameters( FixedBondParameters(b), RandBondAngle(), RandTorsion(), IdealChain())

        PERM = PermData(N; cMax=2.0, cMin=0.2,PreAverage=N_Batch)
        Result = RunSim(Data, Model, RG_Measurement(N_Batch, ones(N)), PERM)

        ### numerical stability
        Result.Weights ./= maximum(Result.Weights)
        Total_Weights = sum(Result.Weights)

        RG_avg, RG_err = ComputeSqrtMeanError(Result.RGs, Result.Weights)
        REE_avg, REE_err = ComputeSqrtMeanError(Result.REEs, Result.Weights)

        RG_exp = sqrt(b^2*(N-1)/6.0)
        REE_exp = sqrt(b^2*(N-1))

        Result.REEs .= sqrt.( Result.REEs)
        bins = collect(0:1:ceil(Int32,maximum(Result.REEs)))

        REE_hist = fit(Histogram, Result.REEs, weights(Result.Weights), bins)
        StatsBase.normalize(REE_hist, mode=:density)
        REE_hist =REE_hist.weights ./ sum(REE_hist.weights*(bins[2]-bins[1]))
        bins = (bins[2:end].+bins[1:end-1])/2.0
        
        theory_hist = @. sqrt(3/(2*π*N*b^2))^3*exp(-3.0*bins^2/(2.0*N*b^2))*bins^2
        theory_hist ./= sum(theory_hist)*(bins[2]-bins[1])

        if DOPLOTs
            fig=Plots.bar(bins, REE_hist, label="Sim result", xlabel="REE", ylabel="P(REE)")
            plot!(bins,theory_hist, label="theory", linestyle=:dot)
            savefig(fig, "./tmp/PERM_IdealChainGyrationradius_$(b)_$(N).pdf")
            @test ComputeKullbackLeiblerDivergence(REE_hist,theory_hist) <0.1
        end

        @test (RG_exp - RG_avg)<3.0*RG_err
        @test (REE_exp - REE_avg)<3.0*REE_err
    end
end
=#


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

@testset "PERM: LJ Potential SAW Scaling" begin
    N_Trial=32
    N_Batch=25_000 #100_000
    #N=100
    N_Vals = collect(30:15:150)
    μ_val= zeros(length(N_Vals))
    Δμ_val= zeros(length(N_Vals))
    for (runID,N) in enumerate(N_Vals[1:end])
        Data = SimData("./tmp/", 1.0, N, N_Trial, N_Batch, 1)
        KP = SimulationParameters( FixedBondParameters(r), RandBondAngle(), RandTorsion(), LJ_Repulsion(ones(Int64,N),Dict(1 => r), 1.0 ))

        val = 1.5
        ratio = 200.0
        PERM = PermData(N; cMax=val, cMin=val/ratio, PreAverage=5_000)

        Result = RunSim(Data,KP, RG_Measurement(N_Batch, ones(N)), PERM);

        PlotsWeightHistorgam(Result.Weights,"./tmp/PERM_LJ_Weights_$(N)",collect(-2:0.1:3.0))

        if VERBOSE
            println("N $N")
            println("unique $(PERM.Total_Grown) vs. $(N_Batch)")
            println("prune $(PERM.TotalPrune) vs. $(PERM.TotalEnrich)")
            println("\n\n\n\n")
        end

        RG_avg, RG_err = ComputeSqrtMeanError(Result.RGs, Result.Weights)
        REE_avg, REE_err = ComputeSqrtMeanError(Result.REEs, Result.Weights)


        μ_val[runID] = REE_avg
        Δμ_val[runID] = REE_avg

    end

    fit = curve_fit(REE_Fit, N_Vals, μ_val,1.0./Δμ_val, [r, 0.58])

    if DOPLOTS
        fig = plot(N_Vals, μ_val,yerr=Δμ_val, yaxis=:log, xaxis=:log, label= "data");
        plot!(N_Vals, REE_Fit(N_Vals, fit.param), label="fit");
        savefig(fig, "./tmp/PERM_LJ_Scaling.pdf");
    end

    (ΔN, Δν) = stderror(fit)
    (N_fit, ν_fit) = fit.param


    @test 0.58 < ν_fit+3.0*Δν && 0.58 > ν_fit-3.0*Δν && Δν<0.10

    println("ν= $(ν_fit)±$(Δν) vs. 0.58...")
end;




