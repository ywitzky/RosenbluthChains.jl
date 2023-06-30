
export RandTorsion, FixedBondParameters, FixedBondAngles, GaussLpBondAngles, IdealChain, SAWParameters, SimulationParameters
using Distributions 


struct RandTorsion <: AbstractTosionAngleParam end

struct FixedBondParameters{T<: Real} <: AbstractBondParam
    ### Bond Length
    BondLength::T
end 

struct FixedBondAngles{T<:Real} <: AbstractBondAngleParam
    BondAngles::Vector{T} 
 end

struct GaussLpBondAngles{T<:Real} <: AbstractBondAngleParam
    PersistenceDist::Normal{T}
    GaussLpBondAngles(Dist::Normal{T}) where {T<:Real} = new{T}(Dist)
end

struct IdealChain <: AbstractSelfAvoidanceParameters end

struct SAWParameters{T<:Real, I<:Integer} <: AbstractSelfAvoidanceParameters
    ### Self avoiding Parameters
    σ_Ex::T
    ϵ_Ex::T
    offset_Ex::T
end

@inline function SetTrialTorsionAngle(data::SimData, param::RandTorsion)
    data.trial_torsion_angle .= _rand(eltype(data.TType)(2*π), data.NTrials)
end

@inline function SetTrialRadius(data::SimData, param::FixedBondParameters) 
    ### is set once in set up and doesnt change afterwards
end

@inline function SetTrialBondAngle(data::SimData,param::FixedBondAngles)
    fill!(data.trial_angle, param.BondAngles[data.id-2])
end

@inline function SetTrialBondAngle(data::SimData,param::GaussLpBondAngles)
    data.trial_angle.= lpToAngle.(rand(param.PersistenceDist, data.NTrials))
end





@inline function InitSimParam(data::SimData,param::FixedBondParameters ) 
    fill!(data.trial_radius, param.BondLength)
end

@inline function InitSimParam(data::SimData,param::GaussLpBondAngles ) end

@inline function InitSimParam(data::SimData,param::FixedBondAngles ) end

@inline function InitSimParam(data::SimData,param::SAWParameters ) end

@inline function InitSimParam(data::SimData,param::IdealChain ) end



SimulationParameters( Bond::FixedBondParameters{T},Angle::Union{GaussLpBondAngles{T},FixedBondAngles{T}})  where {T<:Real} = SimulationParameters{T,Int64}(Bond,Angle,RandTorsion(),IdealChain())