struct NoMeasurement{T<:Number} <: AbstractMeasurement 
    NothingImportant::T
    NoMeasurement() = new{Float64}()
end


function RosenbluthChains.InitMeasurement(data::SimData, param::SimulationParameters, Tmp::NoMeasurement) 
    return NoMeasurement()
end 

function RosenbluthChains.MeasureAfterBatch(data::SimData, param::SimulationParameters,Measurement::NoMeasurement) 
end

function RosenbluthChains.MeasureAfterChainGrowth(data::SimData, param::SimulationParameters, Measurement::NoMeasurement) 
end
    
