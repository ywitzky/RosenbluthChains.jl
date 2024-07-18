using SpecialFunctions

export HarmonicBondLength_Slow, HarmonicBondLength

mutable struct HarmonicBondLength_Slow{T<:Real} <: RosenbluthChains.AbstractBondParam 
    k  ::T ### bond strength for all particles
    r0 ::T ### minima of the bond length potential
    Δr ::T ### sampling width=> take values from interval [r0-Δr, r0+Δr]
    HarmonicBondLength_Slow(k ::T, r0::T, Δr::T) where {T<:Real} = new{T}(k, r0, Δr) ### constructs a Measurement 
end




@inline function SetTrialRadius(data::SimData, param::HarmonicBondLength_Slow)
    rand!(data.trial_radius, eltype(data.TType)) ### created NTrial many random number ∈ [0,1]
    data.trial_radius .*= 2.0*param.Δr ### ∈ [0, 2Δr]
    data.trial_radius .+= param.r0-param.Δr ### ∈ [r0-Δr, r0+Δr]
    return nothing
 end

 @inline function GetTrialBoltzmannWeight(data::SimData,param::HarmonicBondLength_Slow) 
    fill!(data.tmp4, 0.0) ### temporary array of length NTrials that i will use to store the energy

    data.tmp4  .= @. 1.0/2.0*param.k*(data.trial_radius-param.r0)^2 ### compute energies

    ### here comes a numerical trick to stabilise, we add the energies of al potentials first and take the exponential in "GetTrialBoltzmannWeight" before we have to choose the positions in "ChooseTrialPosition"
    # -= since we have to add all energies and take exp(-kt*E); kt=1
    data.LogBoltzmannFaktor .-= data.tmp4
    return nothing
end

function RosenbluthChains.InitSimParam(data::SimData{R, I}, param::HarmonicBondLength_Slow{R}) where {R<:Real, I<: Integer}
    nothing 
end


### cumulative probability density function of the energy term 0.5*k*(x-r0)^2
cdf(x,k,r0) = sqrt(π/(2*k))*(1+erf(sqrt(k/2)*(x-r0)))

mutable struct HarmonicBondLength{T<:Real} <: RosenbluthChains.AbstractBondParam 
    k  ::T ### bond strength for all particles
    r0 ::T ### minima of the bond length potential
    inv_cdf::Interpolations.Extrapolation{T, 1, Interpolations.GriddedInterpolation{T, 1, Vector{T}, Gridded{Linear{Throw{OnGrid}}}, Tuple{Vector{T}}}, Gridded{Linear{Throw{OnGrid}}}, Line{Nothing}} ### interpolated approximation of the inverse cumulative distribution function
    
    ### xgrids are the interpolation points used to construct the inverse cumulative distribution function
    HarmonicBondLength(k ::T, r0::T,  xgrids::Vector{T}) where {T<:Real} = begin
        y = cdf.(xgrids,k,r0)
        y ./= maximum(y)
        Interpolations.deduplicate_knots!(y, move_knots=false) ### gets rid of multiples of the same y values
        inv_cdf = linear_interpolation(y, xgrids,extrapolation_bc=Line())

        new{T}(k, r0, inv_cdf) 
    end
end


@inline function SetTrialRadius(data::SimData, param::HarmonicBondLength)
    data.trial_radius .= param.inv_cdf.(rand(data.NTrials)) 
    return nothing
 end

 @inline function GetTrialBoltzmannWeight(data::SimData,param::HarmonicBondLength) 
    return nothing
end

function RosenbluthChains.InitSimParam(data::SimData{R, I}, param::HarmonicBondLength{R}) where {R<:Real, I<: Integer}
    nothing  ### setup done in constructor of HarmonicBondLength
end