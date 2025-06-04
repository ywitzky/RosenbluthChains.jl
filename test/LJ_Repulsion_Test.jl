using Plots, LsqFit


N= 100
N_Trial=100
θ=ones(N-2)*π/8
φ=ones(N-2)*π/6
φ[1:10] .= ones(10)*π/8
r=3.8


LJ(r,σ, ϵ) =  r < 2^(1.0/6.0)*σ ? 4*ϵ*((σ/r)^12-(σ/r)^6)  + ϵ : 0.0

REE_Fit(N, x) = @. x[1]*(N-1)^x[2]

@testset "LJ Potential" begin
    N_Trial=100
    N=100

    ### Replicate cut off Lennard Jones potential with package intrinsics
    Data = SimData("$TestPath/tmp/", 1.0, N, N_Trial     , 1, 1)
    KP = SimulationParameters( FixedBondParameters(r), RandBondAngle(), RandTorsion(), LJ_Repulsion(ones(Int64,N),Dict(Int64(1) => r), 3.0 ))
    Data.xyz .= [Vector3([0.0,0.0,0.0]) for _ in 1:N]
    x = collect(LinRange(0.8,1.6, N_Trial))*r
    Data.id=3
    Data.trial_positions = [Vector3([1.0,1.0,1.0])*1/sqrt(3)*x[i] for i in 1:N_Trial]
    RosenbluthChains.GetTrialBoltzmannWeight(Data,KP.SAWParam)

    if DOPLOTS
        fig=plot(x, -Data.LogBoltzmannFaktor, label="Sim result")
        plot!(x, LJ.(x, r, 3.0), label="theory", linestyle=:dot, ylim=(-0.1, 40))
        savefig(fig, "$TestPath/tmp/LJ_Pot.pdf")
    end
    @test Data.LogBoltzmannFaktor ≈  -LJ.(x, r, 3.0)

   
    ### run simulations for multiple N to determine scaling exponent ν
    N_Trial=16
    N_Batch=10_000
    N_Vals = collect(25:25:250)
    μ_val= zeros(length(N_Vals))
    Δμ_val= zeros(length(N_Vals))
    for (runID,N) in enumerate(N_Vals)
        Data = SimData("$TestPath/tmp/", 1.0, N, N_Trial, N_Batch, 1)
        KP = SimulationParameters( FixedBondParameters(r), RandBondAngle(), RandTorsion(), LJ_Repulsion(ones(Int64,N),Dict(Int64(1) => r), 1.0 ))
        Meas = RunSim(Data,KP, RosenbluthChains.RG_Measurement(Data.FolderPath, N_Batch));

        REE_avg, REE_err = ComputeSqrtMeanError(Meas.REEs, Meas.Weights)

        μ_val[runID] =  REE_avg #μ
        Δμ_val[runID] = REE_err
    end

        fit = curve_fit(REE_Fit, N_Vals, μ_val,1.0./Δμ_val, [r, 0.58])

    if DOPLOTS
        fig = plot(N_Vals.-1, μ_val, yerr=Δμ_val, yaxis=:log, xaxis=:log, label= "data");
        plot!(N_Vals, REE_Fit(N_Vals, fit.param), label="fit");
        savefig(fig, "$TestPath/tmp/Scaling.pdf");
    end

    (ΔN, Δν) = stderror(fit)
    (N_fit, ν_fit) = fit.param


    @test 0.58 < ν_fit+3.0*Δν && 0.58 > ν_fit-3.0*Δν && Δν<0.33 

    if VERBOSE
        println("ν= $(ν_fit)±$(Δν) vs. 0.58...")
    end
end;




