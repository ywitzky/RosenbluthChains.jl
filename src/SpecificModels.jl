
export RandTorsion, FixedBondParameters, FixedBondAngles, FixedTorsionAngles, GaussLpBondAngles, IdealChain, SAWParameters, SimulationParameters
using Distributions 


struct RandTorsion <: AbstractTorsionAngleParam end

struct FixedBondParameters{T<: Real} <: AbstractBondParam
    ### Bond Length
    BondLength::T
end 

struct FixedBondAngles{T<:Real} <: AbstractBondAngleParam
    BondAngles::Vector{T} 
    FixedBondAngles(BondAngles::Vector{T}) where {T<:Real} = new{T}(BondAngles)
 end

struct FixedTorsionAngles{T<:Real} <: AbstractTorsionAngleParam
    TorsionAngles::Vector{T} 
    FixedTorsionAngles(TorsionAngles::Vector{T}) where {T<:Real} = new{T}(TorsionAngles)
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
    # data.trial_torsion_angle .= 
   rand!(data.trial_torsion_angle, eltype(data.TType))
   data.trial_torsion_angle.*= (eltype(data.TType))(2.0*π)
   CompTrigonometricTrialTorsionAngles(data)
end

@inline function SetTrialTorsionAngle(data::SimData, param::FixedTorsionAngles)
    fill!(data.trial_torsion_angle, param.TorsionAngles[data.id-2])
    CompTrigonometricTrialTorsionAngles(data)
end

@inline function SetTrialRadius(data::SimData, param::FixedBondParameters) 
    ### is set once in set up and doesnt change afterwards
end

@inline function SetTrialBondAngle(data::SimData,param::FixedBondAngles)
    fill!(data.trial_angle, param.BondAngles[data.id-2])
    CompTrigonometricTrialBondAngles(data)
end

@inline function SetTrialBondAngle(data::SimData,param::GaussLpBondAngles)
    data.trial_angle.= lpToRigidAngle.(rand(param.PersistenceDist, data.NTrials))
end

function ChooseTrialPosition(data::SimData,param::IdealChain) ### Assume no external potential
    data.tid  = rand(1:data.NTrials)
end

@inline function GetTrialBoltzmannWeight(data::SimData,param::FixedBondParameters) nothing end

@inline function GetTrialBoltzmannWeight(data::SimData,param::FixedBondAngles) nothing end

@inline function GetTrialBoltzmannWeight(data::SimData,param::FixedTorsionAngles) nothing end



@inline function InitSimParam(data::SimData,param::FixedBondParameters ) 
    fill!(data.trial_radius, param.BondLength)
end

@inline function InitSimParam(data::SimData,param::GaussLpBondAngles ) end

@inline function InitSimParam(data::SimData,param::FixedBondAngles ) end

@inline function InitSimParam(data::SimData,param::FixedTorsionAngles ) end

@inline function InitSimParam(data::SimData,param::SAWParameters ) end

@inline function InitSimParam(data::SimData,param::IdealChain ) end

@inline function InitSimParam(data::SimData,param::RandTorsion ) end




SimulationParameters( Bond::FixedBondParameters{T},Angle::Union{GaussLpBondAngles{T},FixedBondAngles{T}})  where {T<:Real} = SimulationParameters{T,Int64}(Bond,Angle,RandTorsion(),IdealChain())

SimulationParameters( Bond::FixedBondParameters{T},Angle::Union{GaussLpBondAngles{T},FixedBondAngles{T}, },  Torsion ::AbstractTorsionAngleParam)  where {T<:Real} = SimulationParameters{T,Int64}(Bond,Angle,Torsion,IdealChain())