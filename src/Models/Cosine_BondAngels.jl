export Cosine_BondAngles, SetTrialBondAngle, GaussianLp_Cosine_BondAngles

using Bessels, Distributions, NumericalIntegration

struct Cosine_BondAngles{T<:Real} <: AbstractBondAngleParam
    PersistenzLengths::Vector{T} 
    K::Vector{T} 
    Cosine_BondAngles(PersistenzLengths::Vector{T}) where {T<:Real} = new{T}(PersistenzLengths, zeros(T, length(PersistenzLengths)))
end

### is not optimised for very stiff polymers, where some values will never occur
@inline function SetTrialBondAngle(data::SimData,Cosine_BondAngles)
    data.trial_angle.= rand(eltype(data.TType), data.NTrials) * eltype(data.TType)(π)
end

function GetTrialBoltzmannWeight(data::SimData,param::Cosine_BondAngles)
    data.LogBoltzmannFaktor .+= @. param.K[data.id]*(1+cos(data.trial_angle))
    nothing
end

### Init function doesnt allow <cos> values larger 1000 and smaller than ≈1/3
@inline function InitSimParam(data::SimData,param::Cosine_BondAngles ) 
    param.K .= -1.0.*solveRecursive.(lpToCosAngle.(param.PersistenzLengths), AvgCos, 0.1,500.0 )
end

### Analytical solution to <cos(x)>=integral cos(x)exp(+K(1+cos(x)))/Z; assuming K is negativ to have potential
#function AvgCos(K)
#    return besseli1(K)/besseli0(K)
#end

function AvgCos(K, x)
    f(y)=@. cos(y)^2*exp(-1.0*K*(1.0-cos(y)))
    return integrate(x, f(x))/besseli1(K)
end


### Same potential but, K's are drawn such that lp is Gaussian distributed
mutable struct GaussianLp_Cosine_BondAngles{T<:Real, I<:Integer} <: AbstractBondAngleParam
    PersistenceDist::Normal{T}
    N_Hashes::Int64
    K_Hashfunc::Function
    K_HashTable::Vector{T}
    K::Vector{T} 
    μ::T
    σ::T
    GaussianLp_Cosine_BondAngles(μ::T, σ::T; N_Hashes::I=200, K_Hashfunc=x-> trunc(Int32,(x-(μ-5.0*σ))*(N_Hashes)/(10.000000001*σ)+1)) where {T<:Real, I<:Integer} = new{T,I}(Normal(μ, σ), N_Hashes,K_Hashfunc, zeros(T,N_Hashes), zeros(T, 1), μ, σ)
end

### Init function doesnt allow <cos> values larger 1000 and smaller than ≈1/3
function InitSimParam(data::SimData,param::GaussianLp_Cosine_BondAngles ) 
    ### Init K_HashTable
    angles=collect(range(0,π,10_000))
    for lp in range(param.μ-5.0*param.σ,param.μ+5.0*param.σ, param.N_Hashes)
        param.K_HashTable[param.K_Hashfunc(lp)] = -1.0.*solveRecursive.(lpToCosAngle.(lp), K->AvgCos(K, angles), 10^-6,10.0 )
    end
    #println(param.K_HashTable)
    param.K=zeros(eltype(data.TType), data.NTrials)
end

### is not optimised for very stiff polymers, where some values will never occur
@inline function SetTrialBondAngle(data::SimData,param::GaussianLp_Cosine_BondAngles)
    data.trial_angle.= rand(eltype(data.TType), data.NTrials) * eltype(data.TType)(π)
end

function GetTrialBoltzmannWeight(data::SimData,param::GaussianLp_Cosine_BondAngles)
    ### get hashed K values according to Gaus distribution of lp   
    data.tmp4 .= rand(param.PersistenceDist, 1)#data.NTrials)
    for (id,val) in  enumerate(data.tmp4)
         param.K[id] = get(param.K_HashTable, param.K_Hashfunc(val) , -1.0.*solveRecursive.(lpToCosAngle(val), y->AvgCos(y, collect(range(0,π,1000))), 10^-6,10.0 ))
    end
    println(data.tmp4, " ", param.K[1])
    data.LogBoltzmannFaktor .+= @. param.K[1]*(1-cos(data.trial_angle))
    nothing
end
