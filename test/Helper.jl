using Base.Iterators: partition

function ComputeMeanAndVariance(data::Vector{T}) where {T<:Real}
    N = length(data)
    μ_est = sum(data)/(N)
    σ_est = sqrt(sum((data.-μ_est).^2)/((N-1)))
    Δμ_est = σ_est/sqrt(N-1)


    ### estimate error of variance
    μ4 =  sum((data.-μ_est).^4)/(N) #https://stats.stackexchange.com/questions/156518/what-is-the-standard-error-of-the-sample-standard-deviation
    n =Float64(N-1)
    Δσ_est = sqrt(1/n*(μ4-(n-3.0)/(n-1.0)*σ_est^4))/(2.0*σ_est)

    return μ_est, Δμ_est, σ_est, Δσ_est
end


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









