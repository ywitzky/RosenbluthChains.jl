export Cosine_BondAngles, SetTrialBondAngle, GaussianLp_Cosine_BondAngles, GaussianInvLp_Cosine_BondAngles, GaussianK_Cosine_BondAngles

using  Distributions, Interpolations



function integrate_Cosine_BondAngle_toA(K::T, A::T) where {T<:Real}
    return (1.0-exp(K*(cos(A)-1.0)))/K
end

function Cosine_BondAngle_func(K::T, x::T) where {T<:Real}
    return exp(-K*(1.0-cos(x)))*sin(x)
end

function maxima_pos(K::T) where {T<:Real}
    return 2.0*atan(sqrt(sqrt(4.0*K^2+1.0)-2.0*K))
end



struct Cosine_BondAngle_Sampler{I<:Int, T<: Real} <: Sampleable{Univariate,Continuous}
    K::T
    N::I
    borders::Vector{T}
    width::T
    likelyhood::Vector{T}
    maxima::Vector{T}
    work_array::Array{Bool}
    Cosine_BondAngle_Sampler(K::T, N=30::Int64) where {T<:Real} = begin
        borders = LinRange(0,π, N+1)

        max_vals = [max( Cosine_BondAngle_func(K, borders[i]), Cosine_BondAngle_func(K, borders[i+1])) for i in 1:N]
        ind =  findfirst(maxima_pos(K).>borders )
        max_vals[ind] = max(max_vals[ind] , Cosine_BondAngle_func(K,maxima_pos(K)))
        new{Int64, T}( K, N, borders, π/Float64(N), integrate_Cosine_BondAngle_toA.(K,borders[2:end]), max_vals, zeros(Bool, N))
    end
end

function give_rand( s::Cosine_BondAngle_Sampler) 
    #s.work_array.=(rand(eltype(s.likelyhood))*s.likelyhood[end]).<= s.likelyhood
    bla = rand(eltype(s.likelyhood))*s.likelyhood[end]
    ind = findfirst(x->x>=bla, s.likelyhood)
    return rand_help(ind, s)
end

function rand_help(ind::Int64,s::Cosine_BondAngle_Sampler)
    rand_pos = _rand_off(s.width, s.borders[ind])
    while rand() > Cosine_BondAngle_func(s.K, rand_pos)/s.maxima[ind]
        rand_pos = _rand_off(s.width, s.borders[ind])
    end
    return rand_pos
end



struct Cosine_BondAngles{T<:Real} <: AbstractBondAngleParam
    PersistenzLengths::Vector{T} 
    AngleGenerator::Vector{Cosine_BondAngle_Sampler{Int64, T}}
    Cosine_BondAngles(PersistenzLengths::Vector{T}) where {T<:Real} = new{T}(PersistenzLengths, Cosine_BondAngle_Sampler.(solveRecursive.(lpToCosAngle.(PersistenzLengths),AvgCos , 0.1,200.0 ), 30))
end

### is not optimised for very stiff polymers, where some values will never occur
@inline function SetTrialBondAngle(data::SimData,param::Cosine_BondAngles)
    for i in 1:data.NTrials
        data.trial_angle[i] =  give_rand(param.AngleGenerator[data.id])
    end
    CompTrigonometricTrialBondAngles(data)
end

function GetTrialBoltzmannWeight(data::SimData,param::Cosine_BondAngles)
    nothing
end

### Init function doesnt lp values larger ≈700A and smaller than 2A
@inline function InitSimParam(data::SimData,param::Cosine_BondAngles ) 
    nothing
end



function AvgCos(K) #### coth(K)-1/K
    return  (coth(K)-1.0/K) #1.0/K*(exp(-2.0*K)*(exp(2.0*K)*(K-1.0)+K+1.0))/(1.0-exp(-2.0*K))
end




struct Cosine_BondAngle_Prefactor_Sampler{I<:Int, T<: Real} <: Sampleable{Univariate,Continuous}
    μ::T
    σ::T
    K_borders::Vector{T}
    K_cdf::Vector{T}
    K_max::Vector{}
    K_width::T
    max_val::T
    P::Any
    Cosine_BondAngle_Prefactor_Sampler(μ::T, σ::T, P, K; NK=10::Int64, NA=10::Int64) where {T<:Real} = begin
        K_borders = collect(LinRange(μ-5.0*σ,μ+5.0*σ, NK+1))
        K_width= K_borders[2]-K_borders[1]

        K_mean = (K_borders[2:end].+K_borders[1:end-1])/2.0
        K_cdf = cdf.(Normal(μ, σ), K_borders[2:end])
        K_max  = [max(pdf(Normal(μ, σ),K_borders[i]), pdf(Normal(μ, σ),K_borders[i+1]) ) for i in 1:NK]
        K_ind =  findfirst(μ.<=K_borders) -1
        K_max[K_ind] = max(pdf(Normal(μ, σ),μ), K_max[K_ind])
        max_val = maximum([P(k,maxima_pos(K(k))) for k in  collect(LinRange(μ-1.0*σ,μ+1.0*σ, 100))]) ### safety margin
        new{Int64, T}(μ, σ,K_borders,K_cdf,K_max, K_width,max_val, P)
    end
end


function give_rand( s::Cosine_BondAngle_Prefactor_Sampler)
    ### draw gaussian for K or lp
    rand_K  = _rand_off(10.0*s.σ, s.μ-5.0*s.σ)
    rnd_num = rand(eltype(s.K_cdf))*s.K_cdf[end]
    ind = findfirst(x->x>=rnd_num, s.K_cdf)
    rand_K = _rand_off(s.K_width, s.K_borders[ind])
    while rand()*s.K_max[ind] > pdf(Normal(s.μ, s.σ), rand_K)
        rand_K = _rand_off(s.K_width,  s.K_borders[ind])
    end

    ### draw angle accordingly
    rand_ϕ = rand(eltype(rand_K))*π
    while rand()*s.max_val >  s.P(rand_K,rand_ϕ)
        rand_ϕ = rand(eltype(rand_K))*π
    end
    return (rand_K, rand_ϕ)
end



### Same potential but, K's are drawn such that lp is Gaussian distributed
mutable struct GaussianK_Cosine_BondAngles{T<:Real} <: AbstractBondAngleParam
    μ::T
    σ::T
    Sampler::Cosine_BondAngle_Prefactor_Sampler
    GaussianK_Cosine_BondAngles(μ::T, σ::T) where {T<:Real} = begin 
        new{T}(μ, σ, Cosine_BondAngle_Prefactor_Sampler(μ, σ, Cosine_BondAngle_func, x->x))
    end
end

### Init function doesnt allow <cos> values larger 1000 and smaller than ≈1/3
function InitSimParam(data::SimData,param::GaussianK_Cosine_BondAngles ) 
    nothing
end

### is not optimised for very stiff polymers, where some values will never occur
@inline function SetTrialBondAngle(data::SimData,param::GaussianK_Cosine_BondAngles)
    for i in 1:data.NTrials
        data.trial_angle[i] = give_rand(param.Sampler)[2]
    end
    CompTrigonometricTrialBondAngles(data)
    nothing
end

@inline function GetTrialBoltzmannWeight(data::SimData,param::GaussianK_Cosine_BondAngles)
    nothing
end


##############################

### Same potential but, K's are drawn such that lp is Gaussian distributed
mutable struct GaussianLp_Cosine_BondAngles{T<:Real} <: AbstractBondAngleParam
    μ::T
    σ::T
    Sampler::Cosine_BondAngle_Prefactor_Sampler
    GaussianLp_Cosine_BondAngles(μ::T, σ::T) where {T<:Real} = begin
        ### make linear interpolation
        lp = collect(LinRange(μ-6.0*σ, μ+6.0*σ, 500))
        K_val= solveRecursive.(lpToCosAngle.(lp),AvgCos , 0.0001,10000.0 )
        K = linear_interpolation(lp, K_val)
        H(x1,x2) = Cosine_BondAngle_func(K(x1), x2)
        new{T}(μ, σ, Cosine_BondAngle_Prefactor_Sampler(μ, σ, H, K))
    end
end

### Init function doesnt allow <cos> values larger 1000 and smaller than ≈1/3
function InitSimParam(data::SimData,param::GaussianLp_Cosine_BondAngles ) 
    nothing
end

### is not optimised for very stiff polymers, where some values will never occur
@inline function SetTrialBondAngle(data::SimData,param::GaussianLp_Cosine_BondAngles)
    for i in 1:data.NTrials
        data.trial_angle[i] = give_rand(param.Sampler)[2]
    end
    CompTrigonometricTrialBondAngles(data)
    nothing
end

function GetTrialBoltzmannWeight(data::SimData,param::GaussianLp_Cosine_BondAngles)
    nothing
end