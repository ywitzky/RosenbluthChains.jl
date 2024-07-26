using Plots

function Debye(r, D, q1,q2, ϵr,relT, rcut)
    ### T needs to be relative to 298K: T = 0.8 means 298*0.8 Degrees
    kt = 4.11 # 10^-21 J at 298K need to rescale if temperature changes
    e = 1.602 # 10-19 C
    ϵ0 = 8.854 # 10^-12  As/(Vm)
    if r < rcut
        return q1*q2*e^2/(4*π*ϵr*ϵ0*(kt*relT))*exp(-r/D)/r*10^5 ### float 
    else
        return 0
    end
end

function er(T)
    return 5321.0/T+233.76-0.9297*T+1.417*10.0^-3*T^2-8.292*10^-7*T^3
end

function LJ(r, σ, ϵ)
    return 4*ϵ*((σ/r)^12-(σ/r)^6)
end

function AshbaughHatch(r, σ1, σ2,λ1, λ2, ϵ, rcut, rel_T)
    σ=(σ1+σ2)/2.0
    λ=(λ1+λ2)/2.0
    ϵ /= rel_T
    if r > rcut
        return 0
    else
        if r< 2^(1/6)*σ
            return LJ(r, σ, ϵ)-λ*LJ(rcut, σ, ϵ)+ϵ*(1-λ)
        else
            return λ*(LJ(r, σ, ϵ)-LJ(rcut, σ, ϵ))
        end
    end
end

N = 3
N_Trial = 1000


### Take calvados2 values
OneToLambda = Dict('D'=>0.0416, 'E'=>0.000693, 'H'=>0.4663, 'K'=>0.1790, 'F'=>0.8672, 'W'=>0.9893, 'Y'=>0.9774, 'A'=>0.27432, 'R'=>0.7307, 'N'=>0.42558, 'C'=>0.5615, 'Q'=>0.3934, 'G'=>0.7058, 'I'=>0.5423, 'L'=>0.6440, 'M'=>0.5308, 'P'=>0.3593, 'S'=>0.4625, 'T'=>0.3713, 'V'=>0.2083, '#'=>0.4625, '&'=>0.3713, '@'=>0.189)

### https://github.com/KULL-Centre/papers/blob/main/2022/CG-cutoffs-Tesei-et-al/SC/05_24_20/residues.csv only for test purpose
OneToLambda_05_24_20 = Dict('D'=>0.05517, 'E'=>0.001712, 'H'=>0.4488, 'K'=>0.1785, 'F'=>0.8823, 'W'=>0.9926, 'Y'=>0.9318, 'A'=>0.2985, 'R'=>0.7307, 'N'=>0.4255, 'C'=>0.5571, 'Q'=>0.3618, 'G'=>0.6645, 'I'=>0.5362, 'L'=>0.6051, 'M'=>0.5480, 'P'=>0.3452, 'S'=>0.5616, 'T'=>0.3135, 'V'=>0.2805, '#'=>0.4625, '&'=>0.3713)

OneToSigma= Dict('A' => 5.04, 'R'=> 6.56, 'N'=> 5.68 , 'D'=> 5.58, 'C'=> 5.48, 'Q'=> 6.02, 'E'=> 5.92, 'G'=> 4.5, 'H'=> 6.08 , 'I'=> 6.18, 'L'=> 6.18, 'K'=> 6.36, 'M'=> 6.18, 'F'=> 6.36, 'P'=> 5.56, 'S'=>5.18 , 'T'=> 5.62, 'W'=> 6.78, 'Y'=> 6.46, 'V'=> 5.86, '#'=>6.36, '&'=> 6.62, '@'=> 7.38)

OneToCharge=Dict('D'=>-1.0, 'E'=>-1.0, 'H'=>0.0, 'K'=>1.0, 'F'=>0.0, 'W'=>0.0, 'Y'=>0.0, 'A'=>0.0, 'R'=>1.0, 'N'=>0.0, 'C'=>0.0, 'Q'=>0.0, 'G'=>0.0, 'I'=>0.0, 'L'=>0.0, 'M'=>0.0, 'P'=>0.0, 'S'=>0.0, 'T'=>0.0, 'V'=>0.0, 'X'=>0.0, '#'=>-2., '&'=>-2., '@'=>-2. )
#
rc_ah = 20.0 ## AA
rc_deb = 40.0 # AA
ϵ_ah = 0.8368 # kJ/mol
ϵ_ah /= 2.479 # now in kT at 298 K

NA = 6.02214086 * 10.0^23 
kt = 4.11 # 10^-21 J at 298K need to rescale if temperature changes
e = 1.602 # 10-19 C
ϵ0 = 8.854 # 10^-12  As/(Vm)
### 1-2 interactions are turned of, only amino acids 1 and 3 matter
@testset "Debye_Ashbaught_Hatch Potential" begin
    for (cnt,(Sequence, pH, rel_T, cs_deb)) in enumerate(zip(["FFF","HHH", "MMA", "KEE"], [6.0, 6.5, 7., 8.0], [0.8, 0.9, 1.0, 1.1], [0.15, 0.2, 0.25, 0.3]))
    #println("\n")
    #println(Sequence)
    OneToCharge['H'] = 1. / ( 1.0 + 10.0^(pH-6) ) ### HIS is pH-Dependent for calvados2
    ϵr_deb = er(298.0*rel_T)
    #println("er out: $(ϵr_deb)")
    B = e^2/(4*π*ϵ0*ϵr_deb*kt*rel_T)#*10^-5m ### in Angstroem
    D = sqrt(1.0/(8.0*π*B*cs_deb*NA)) *10^1 *10^10### convert to Angstroem

    #println("λ_out $((OneToLambda[Sequence[1]]+OneToLambda[Sequence[3]])/2.0)")
    #println("D out: $(D)")
    ϵ_ah_l = deepcopy(ϵ_ah)

    Model = SimulationParameters( FixedBondParameters(3.8), RandBondAngle(), RandTorsion(), Debye_Ashbaugh_Hatch(deepcopy(Sequence), OneToSigma, OneToLambda, ϵ_ah_l, rc_ah, OneToCharge, rc_deb, cs_deb, rel_T))

    σ = [OneToSigma[c] for c in Sequence]

    ### recreate each time to ensure that everything is filled with 0
    Data = SimData("../tmp/", 1.0, N, N_Trial, 1, 1)
    Data.xyz .= [Vector3([0.0,0.0,0.0]) for _ in 1:N]
    x = collect(LinRange( 0.85*minimum(σ),rc_deb, N_Trial))
    Data.id=3
    ### prepare independent trial positions that sample force field
    Data.trial_positions = [Vector3([1.0,1.0,1.0])*1/sqrt(3)*x[i] for i in 1:N_Trial]


    RosenbluthChains.GetTrialBoltzmannWeight(Data,Model.SAWParam)
    fig=plot(x, -Data.LogBoltzmannFaktor, label="Sim result")
    theory = (AshbaughHatch.(x, OneToSigma[Sequence[1]], OneToSigma[Sequence[3]],OneToLambda[Sequence[1]], OneToLambda[Sequence[3]], ϵ_ah_l, rc_ah, rel_T) .+ Debye.(x, D,OneToCharge[Sequence[1]],OneToCharge[Sequence[3]], ϵr_deb,rel_T, rc_deb))

    plot!(x,theory, label="theory", linestyle=:dot)
    savefig(fig, "./tmp/Debye_Ashbaugh_Hatch_Pot_$(Sequence)_pH$(pH)_$(298.0*rel_T)K_$(cs_deb).pdf")
    @test all( -Data.LogBoltzmannFaktor .≈ theory)
    end
end

