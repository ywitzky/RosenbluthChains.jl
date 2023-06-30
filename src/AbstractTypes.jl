export SimulationParameters , InitMeasurement, MeasureAfterChainGrowthMeasureAfterBatch, SaveMeasurements, AbstractMeasurement

abstract type AbstractBondParam end

abstract type AbstractBondAngleParam end

abstract type AbstractTosionAngleParam end

abstract type AbstractSelfAvoidanceParameters end

abstract type AbstractMeasurement end


struct SimulationParameters{T<:Real, I<:Integer}
    BondParam::AbstractBondParam
    BondAngleParam::AbstractBondAngleParam
    TorsionAngleParam::AbstractTosionAngleParam
    SAWParam::AbstractSelfAvoidanceParameters
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
end

function GetTrialBoltzmannWeight(data::SimData,param::SimulationParameters)
    #  printstyled("No BoltzmannWeights implemented yet.\n" ; color=:red)
  end

function ChooseTrialPosition(data::SimData,param::SimulationParameters)
    ChooseTrialPosition(data, param.SAWParam)
end

function InitMeasurement(data::SimData, param::SimulationParameters, _::AbstractMeasurement) end 

function MeasureAfterChainGrowth(data::SimData, param::SimulationParameters, Measurement::AbstractMeasurement) end

function MeasureAfterBatch(data::SimData, param::SimulationParameters,Measurement::AbstractMeasurement) end

function SaveMeasurements(data::SimData, param::SimulationParameters,Measurement::AbstractMeasurement) end