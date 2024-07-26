export Calvados2

include("./HarmonicBonds.jl")
include("./Debye_Ashbaugh_Hatch.jl")

OneToLambda = Dict('D'=>0.0416, 'E'=>0.000693, 'H'=>0.4663, 'K'=>0.1790, 'F'=>0.8672, 'W'=>0.9893, 'Y'=>0.9774, 'A'=>0.27432, 'R'=>0.7307, 'N'=>0.42558, 'C'=>0.5615, 'Q'=>0.3934, 'G'=>0.7058, 'I'=>0.5423, 'L'=>0.6440, 'M'=>0.5308, 'P'=>0.3593, 'S'=>0.4625, 'T'=>0.3713, 'V'=>0.2083, '#'=>0.4625, '&'=>0.3713, '@'=>0.189)

OneToSigma= Dict('A' => 5.04, 'R'=> 6.56, 'N'=> 5.68 , 'D'=> 5.58, 'C'=> 5.48, 'Q'=> 6.02, 'E'=> 5.92, 'G'=> 4.5, 'H'=> 6.08 , 'I'=> 6.18, 'L'=> 6.18, 'K'=> 6.36, 'M'=> 6.18, 'F'=> 6.36, 'P'=> 5.56, 'S'=>5.18 , 'T'=> 5.62, 'W'=> 6.78, 'Y'=> 6.46, 'V'=> 5.86, '#'=>6.36, '&'=> 6.62, '@'=> 7.38)

OneToCharge=Dict('D'=>-1.0, 'E'=>-1.0, 'H'=>0.0, 'K'=>1.0, 'F'=>0.0, 'W'=>0.0, 'Y'=>0.0, 'A'=>0.0, 'R'=>1.0, 'N'=>0.0, 'C'=>0.0, 'Q'=>0.0, 'G'=>0.0, 'I'=>0.0, 'L'=>0.0, 'M'=>0.0, 'P'=>0.0, 'S'=>0.0, 'T'=>0.0, 'V'=>0.0, 'X'=>0.0, '#'=>-2., '&'=>-2., '@'=>-2. )


function Calvados2(Sequence::String, Temperature::R; OneToLambda=OneToLambda, OneToSigma=OneToSigma,OneToCharge=OneToCharge,  SaltConcentration=-1.0, pH=7.0, rc_ah=20.0, rc_deb=40) where {R<:Real}
    ### drop powers of 10 for numerical stability, compute by hand
    kt = 4.11 # 10^-21 J at 298K need to rescale if temperature changes
    NA = 6.02214086 #* 10.0^23 

    rel_T = Temperature/298.0
    kbt = kt * rel_T ### fit kbT to the temperature
    if SaltConcentration==-1.0
        SaltConcentration=0.15
    end

    OneToCharge['H'] = 1. / ( 1.0 + 10.0^(pH-6) ) 

    return SimulationParameters( HarmonicBondLength(803.3/(NA*kbt), 3.8, collect(0.1:0.1:7.0)), RandBondAngle(), RandTorsion(), Debye_Ashbaugh_Hatch(Sequence, OneToSigma, OneToLambda, 0.8368*10/(NA*kbt), 20.0, OneToCharge, 40.0, SaltConcentration, rel_T))
end