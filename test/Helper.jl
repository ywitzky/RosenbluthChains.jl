using Base.Iterators: partition

function ComputeKullbackLeiblerDivergence(P,Q)
    return sum([p*log(p/q) for (p,q) in zip(P,Q) if p>0 ])
end

function ComputeMeanError(data::Vector{R}, Weights::Vector{R}; NIntervals = 10) where {R<:Real}
    subinvertals  = collect(partition(1:length(data), floor(Int32,length(data)/Float32(NIntervals))))[1:end-1]
    avg = sum(data.*Weights)/sum(Weights)

    sub_Weights = [sum(Weights[interval]) for interval in subinvertals]

    avg_array = [ sum(data[interval].*Weights[interval])/sub_Weights[id]  for (id,interval) in enumerate(subinvertals)]

    neff = sum(sub_Weights)^2/(sum(sub_Weights.^2))

    error = sqrt(sum((avg.-avg_array).^2 .*sub_Weights)/(sum(sub_Weights)*(neff-1)))

    return avg, error
end

function ComputeSqrtMeanError(data::Vector{R}, Weights::Vector{R}; NIntervals = 10) where {R<:Real}
    subinvertals  = collect(partition(1:length(data), floor(Int32,length(data)/Float32(NIntervals))))[1:end-1]
    avg = sqrt(sum(data.*Weights)/sum(Weights))

    sub_Weights = [sum(Weights[interval]) for interval in subinvertals]

    avg_array = [ sqrt(sum(data[interval].*Weights[interval])/sub_Weights[id])  for (id,interval) in enumerate(subinvertals)]

    neff = sum(sub_Weights)^2/(sum(sub_Weights.^2))

    error = sqrt(sum((avg.-avg_array).^2 .*sub_Weights)/(sum(sub_Weights)*(neff-1)))

    return avg, error
end

function PlotsWeightHistorgam(Weights, path::String, limits=-8:0.1:0)
    limits=collect(limits)
    b = 10.0 .^ (limits)
    fig = Plots.histogram(collect(Weights); bin=b, xscale=:log10, yscale=:identity,  xlim=extrema(b), ylabel="Weights")
    Plots.savefig(fig, path)
end


#### measurement for radius of gyration and end-to-end distance length
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






#### measurement for persistence length
mutable struct AMeasurement{T<:Number} <: AbstractMeasurement 
    SumRosenbluthWeights::T
    AvgCosBondAngle::Array{T}
    Weight::Vector{T}
    AMeasurement(NBeads::I, NSteps::I) where {I<:Integer} = new{BigFloat}(zero(BigFloat ), zeros(BigFloat,(NSteps,  NBeads-2)), zeros(BigFloat, NSteps)) #where {T::BigFloat}
end


function RosenbluthChains.InitMeasurement(data::SimData, param::SimulationParameters, Tmp::AMeasurement) 
    return AMeasurement(data.NBeads, data.BatchSize)
end 

function RosenbluthChains.MeasureAfterBatch(data::SimData, param::SimulationParameters,Measurement::AMeasurement) 
    nothing
end

function RosenbluthChains.MeasureAfterChainGrowth(data::SimData, param::SimulationParameters, Measurement::AMeasurement) 
    rel = [data.xyz[i+1]-data.xyz[i] for i in 1:data.NBeads-1]
    angles=[angle( rel[i],rel[i+1]) for i in 1:data.NBeads-2]

    Measurement.AvgCosBondAngle[data.id_in_batch,:] .= BigFloat.(cos.(angles))
    Measurement.Weight[data.id_in_batch] = getRosenbluthWeigth(data, param)
    Measurement.SumRosenbluthWeights += data.RosenbluthWeight #getRosenbluthWeigth(data, param)#
end
    