using Random, Distributions

export _rand, _rand2π, _rand_off
### Random numbers
_rand(len::T)     where {T<:Real} = rand(T) * len
_rand_off(len::T) where {T<:Real} = rand(T) * len - len / 2.0
_rand_off(width::T,A::T) where {T<:Real} = rand(T) * width + A
_rand2π(len::T)   where {T<:Real} = _rand(2.0*T(π))

#=
### @TODO Implement more general, do proper implementation into Distributions

function integrate_Cosine_BondAngle_toA(K::T, A::T) where {T<:Real}
    return (1.0-exp(K*(cos(A)-1.0)))/K
end

function Cosine_BondAngle_func(K::T, x::T) where {T<:Real}
    return exp(-K*(1.0-cos(x)))*sin(x)
end

function maxima_pos(K::T) where {T<:Real}
    return 2*atan(sqrt(sqrt(4.0*K^2+1)-2.0*K))
end

struct Cosine_BondAngle_Sampler{I<:Int, T<: Real} <: Sampleable{Univariate,Continuous}
    K::T
    N::I
    borders::Vector{T}
    width::T
    likelyhood::Vector{T}
    maxima::Vector{T}# = [max(maxima(K, Cosine_BondAngle_func(K, borders[i]), Cosine_BondAngle_func(K, borders[i+1]))) for i in 1:N]

    Cosine_BondAngle_Sampler(K::T, N=30::Int64) where {T<:Real} = begin
        borders = LinRange(0,π, N+1)
        #println([max(maxima(K), Cosine_BondAngle_func(K, borders[i]), Cosine_BondAngle_func(K, borders[i+1])) for i in 1:N])
        max_vals = [max( Cosine_BondAngle_func(K, borders[i]), Cosine_BondAngle_func(K, borders[i+1])) for i in 1:N]
        ind =  findfirst(maxima_pos(K).>borders )
        max_vals[ind] = max(max_vals[ind] , Cosine_BondAngle_func(K,maxima_pos(K)))
        new{Int64, T}( K, N, borders, π/Float64(N), integrate_Cosine_BondAngle_toA.(K,borders[2:end]), max_vals)
    end
end

function give_rand( s::Cosine_BondAngle_Sampler) 
    ind = findfirst((rand(eltype(s.likelyhood))*s.likelyhood[end]).<= s.likelyhood)
    return rand_help(ind, s)
end

function rand_help(ind::Int64,s::Cosine_BondAngle_Sampler)
    rand_pos = _rand_off(s.width, s.borders[ind])
    while rand() > Cosine_BondAngle_func(s.K, rand_pos)/s.maxima[ind]
        rand_pos = _rand_off(s.width, s.borders[ind])
    end
    return rand_pos
end
=#