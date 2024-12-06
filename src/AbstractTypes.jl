export SimulationParameters , InitMeasurement, MeasureAfterChainGrowthMeasureAfterBatch, SaveMeasurements, AbstractMeasurement, ChooseTrialPosition, clear

abstract type AbstractPERMData end

abstract type AbstractBondParam end

abstract type AbstractBondAngleParam end

abstract type AbstractTorsionAngleParam end

abstract type AbstractSelfAvoidanceParameters end

abstract type AbstractMeasurement end

struct NoPERM <:AbstractPERMData ### turns perm off 
    Measure::Bool
    NoPERM() = new(false)
end 


struct SimulationParameters{T<:Real, I<:Integer}
    BondParam::AbstractBondParam
    BondAngleParam::AbstractBondAngleParam
    TorsionAngleParam::AbstractTorsionAngleParam
    SAWParam::AbstractSelfAvoidanceParameters
end


### functions only compute the Log of the boltzmann faktor and finale evaluation of the exponent is done together
function GetTrialBoltzmannWeight(data::SimData,param::SimulationParameters)
        fill!(data.LogBoltzmannFaktor, 0.0)

        ### @TODO: Automatise so that nothing can be forgotten.
        GetTrialBoltzmannWeight(data, param.BondParam)
        GetTrialBoltzmannWeight(data, param.BondAngleParam)
        GetTrialBoltzmannWeight(data, param.TorsionAngleParam)
        GetTrialBoltzmannWeight(data, param.SAWParam)

        #map(field-> GetTrialBoltzmannWeight(data,  getfield(param, field)), fieldnames(typeof(param)))
        data.BoltzmannFaktor .= @. exp(data.LogBoltzmannFaktor)/data.NTrials
        nothing
end


@inline function SetTrialRadius(data::SimData,param::SimulationParameters)
    SetTrialRadius(data, param.BondParam)
end

@inline function SetTrialBondAngle(data::SimData,param::SimulationParameters)
    SetTrialBondAngle(data, param.BondAngleParam)
end

@inline function SetTrialTorsionAngle(data::SimData,param::SimulationParameters)
    SetTrialTorsionAngle(data, param.TorsionAngleParam)
end


@inline function GetUnboundInteractions(data::SimData,param::SimulationParameters)
    GetUnboundInteractions(data, param.SAWParam)
end

@inline function InitSimParam(data::SimData,param::SimulationParameters)
    for field in fieldnames(typeof(param))
        InitSimParam(data, getfield(param, field))
    end
    data.max_interaction_length = get_max_bondlength(data, param.BondParam) + get_non_bonded_cutoff(data,param.SAWParam)
end

function GetTrialBoltzmannWeight(data::SimData,param::AbstractBondParam) nothing end

function GetTrialBoltzmannWeight(data::SimData,param::AbstractBondAngleParam) nothing end

function GetTrialBoltzmannWeight(data::SimData,param::AbstractTorsionAngleParam) nothing end

function GetTrialBoltzmannWeight(data::SimData,param::AbstractSelfAvoidanceParameters) nothing end


#const small_diff = 0.999999

function ChooseTrialPosition(data::SimData{T,I},param::SimulationParameters) where {T<:Real, I<:Integer}
    GetTrialBoltzmannWeight(data,param)

    ### avoid energies that are so small that they crash the floating point precision and turn to "nothing"
    data.mask = isnothing.(data.BoltzmannFaktor)
    data.BoltzmannFaktor[data.mask] .= zero(T)

    cumsum!(data.tmp4, data.BoltzmannFaktor)

    rnd_num=rand(T)*data.tmp4[end]
    data.tid = findfirst(x->x>=rnd_num, data.tmp4) 
    #data.tid   = typeof(data.tid_tmp)==Nothing ? data.NTrials : data.tid_tmp
    data.btmp .= exp.(data.LogBoltzmannFaktor)
    data.RosenbluthWeight *= sum(data.btmp)/data.NTrials

    data.xyz[data.id] .= data.trial_positions[data.tid]
    data.x[data.id] = data.trial_positions[data.tid][1]
    data.y[data.id] = data.trial_positions[data.tid][2]
    data.z[data.id] = data.trial_positions[data.tid][3]

    nothing
end

function InitMeasurement(data::SimData, param::SimulationParameters, _::AbstractMeasurement) end 

function MeasureAfterChainGrowth(data::SimData, param::SimulationParameters, Measurement::AbstractMeasurement) end

function MeasureAfterBatch(data::SimData, param::SimulationParameters,Measurement::AbstractMeasurement) end

function SaveMeasurements(data::SimData, param::SimulationParameters,Measurement::AbstractMeasurement) end


SimulationParameters( Bond::AbstractBondParam,Angle::AbstractBondAngleParam, Torsion::AbstractTorsionAngleParam, SAW::AbstractSelfAvoidanceParameters)  = SimulationParameters{Float64,Int64}(Bond,Angle,Torsion,SAW)

function clear(data::SimData, param::AbstractSelfAvoidanceParameters) nothing end


function get_max_bondlength(data::SimData,Bond::AbstractBondParam) 
    return eltype(data.TType)(Inf) 
end

function get_non_bonded_cutoff(data::SimData,NonBonded::AbstractSelfAvoidanceParameters)
    return eltype(data.TType)(Inf) 
end