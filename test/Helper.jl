using Base.Iterators: partition

function ComputeKullbackLeiblerDivergence(P,Q)
    return sum([p*log(p/q) for (p,q) in zip(P,Q) if p>0 ])
end

function ComputeMeanError(data::Vector{R}, Weights::Vector{R}, NIntervals = 10) where {R<:Real}
    subinvertals  = collect(partition(1:length(data), ceil(Int32,length(data)/NIntervals)))
    avg = sqrt(sum(data.*Weights)/sum(Weights))

    avg_array = [ sqrt(sum(data[interval].*Weights[interval])/sum(Weights[interval]))  for interval in subinvertals]

    error = sqrt(sum((avg.-avg_array).^2))/NIntervals

    return avg, error
end



mutable struct RG_Measurement{T<:Real} <: AbstractMeasurement 
    RGs::Vector{T}
    REEs::Vector{T}
    Weights::Vector{T}
    Masses::Vector{T}
    M_Total::T
    RG_Measurement(NProteins::Int, Masses::Vector{Float64}) = new{Float64}(zeros(NProteins), zeros(NProteins), zeros(NProteins), Masses, sum(Masses)) 
end


### initialise measurments with empty vectors
RosenbluthChains.InitMeasurement(data::SimData, param::SimulationParameters, Tmp::RG_Measurement) =  RG_Measurement(data.BatchSize, Tmp.Masses)

function RosenbluthChains.MeasureAfterChainGrowth(data::SimData, param::SimulationParameters,Measurement::RG_Measurement) 
    com = sum(Measurement.Masses.*data.xyz)./Measurement.M_Total

    Measurement.RGs[data.id_in_batch] = sum(Measurement.Masses.*sqr_norm.([pos-com for pos in data.xyz]))/Measurement.M_Total
    Measurement.REEs[data.id_in_batch] = sqr_norm(data.xyz[end]-data.xyz[1])
    Measurement.Weights[data.id_in_batch] =  getRosenbluthWeigth(data, param)
end

RosenbluthChains.MeasureAfterBatch(data::SimData, param::SimulationParameters, Measurement::RG_Measurement) = nothing