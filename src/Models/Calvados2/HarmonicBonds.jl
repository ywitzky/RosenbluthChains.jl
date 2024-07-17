
mutable struct HarmonicBondLength{T<:Real} <: RosenbluthChains.AbstractBondParam 
    k  ::T ### bond strength for all particles
    r0 ::T ### minima of the bond length potential
    Δr ::T ### sampling width=> take values from interval [r0-Δr, r0+Δr]
end

@inline function SetTrialRadius(data::SimData, param::HarmonicBondLength)
    rand!(data.trial_radius, eltype(data.TType)) ### created NTrial many random number ∈ [0,1]
    data.trial_radius .*= 2.0*param.Δr ### ∈ [0, 2Δr]
    data.trial_radius .+= param.r0-param.Δr ### ∈ [r0-Δr, r0+Δr]
    return nothing
 end

 @inline function GetTrialBoltzmannWeight(data::SimData,param::HarmonicBondLength) 
    fill!(data.tmp4, 0.0) ### temporary array of length NTrials that i will use to store the energy

    data.tmp4  .= @. 1.0/2.0*(data.trial_radius-r0)^2 ### compute energies

    ### here comes a numerical trick to stabilise, we add the energies of al potentials first and take the exponential in "GetTrialBoltzmannWeight" before we have to choose the positions in "ChooseTrialPosition"
    # -= since we have to add all energies and take exp(-kt*E); kt=1
    data.LogBoltzmannFaktor .-= data.tmp4
    return nothing
end