export Cosine_BondAngles, SetTrialBondAngle, GaussianLp_Cosine_BondAngles,  GaussianK_Cosine_BondAngles, GaussianFixedK_Cosine_BondAngles, κ_from_lp

#GaussianInvLp_Cosine_BondAngles,
using  Distributions, Interpolations, Distributions, HCubature



function integrate_Cosine_BondAngle_toA(K::T, A::T) where {T<:Real}
    return (1.0-exp(K*(cos(A)-1.0)))/K
end

### P(θ|κ)
function Cosine_BondAngle_func(K::T, x::T) where {T<:Real}
    return exp(-K*(1.0-cos(x)))*sin(x)
end

function maxima_pos(K::T) where {T<:Real}
    return 2.0*atan(sqrt(sqrt(4.0*K^2+1.0)-2.0*K))
end

κ_from_lp(κ) = solveRecursive(lpToCosAngle(κ),AvgCos , eltype(κ)( 10^-20), eltype(κ)(1.0*10^7)) 




struct Cosine_BondAngle_Sampler{I<:Integer, T<: Real} <: Sampleable{Univariate,Continuous}
    K::T
    N::I
    borders::Vector{T}
    width::T
    likelyhood::Vector{T}
    maxima::Vector{T}
    work_array::Array{Bool}
    Cosine_BondAngle_Sampler(K::T, N=30::Int64) where {T<:Real} = begin
        borders = LinRange{T}(0,π, N+1)

        max_vals = [max( Cosine_BondAngle_func(K, borders[i]), Cosine_BondAngle_func(K, borders[i+1])) for i in 1:N]
        ind =  findfirst(maxima_pos(K).>borders )
        max_vals[ind] = max(max_vals[ind] , Cosine_BondAngle_func(K,maxima_pos(K)))
        new{Int64, T}( K, N, borders, π/Float64(N), integrate_Cosine_BondAngle_toA.(K,borders[2:end]), max_vals, zeros(Bool, N))
    end
end

function give_rand( s::Cosine_BondAngle_Sampler{I,T}, data::SimData{T,I}) where {I<:Integer, T<:Real}
    #s.work_array.=(rand(eltype(s.likelyhood))*s.likelyhood[end]).<= s.likelyhood
    bla = rand(T)*s.likelyhood[end]
    for i in axes(s.likelyhood,1)
        if @inbounds s.likelyhood[i] > bla
            return rand_help(i, s, data)
        end
    end
    #ind = findfirst(x->x>=bla, s.likelyhood)
end

function rand_help(ind::Int64,s::Cosine_BondAngle_Sampler, data::SimData)
    @inbounds data.rand_val = _rand_off(s.width, s.borders[ind])
    @inbounds while rand()*s.maxima[ind] > Cosine_BondAngle_func(s.K, data.rand_val)
        @inbounds data.rand_val = _rand_off(s.width, s.borders[ind])
    end
    #return data.rand_pos
    nothing
end



struct Cosine_BondAngles{T<:Real} <: AbstractBondAngleParam
    PersistenzLengths::Vector{T} 
    AngleGenerator::Vector{Cosine_BondAngle_Sampler{Int64, T}}
    Cosine_BondAngles(PersistenzLengths::Vector{T}) where {T<:Real} = new{T}(PersistenzLengths, Cosine_BondAngle_Sampler.(solveRecursive.(lpToCosAngle.(PersistenzLengths),AvgCos , 0.1,200.0 ), 30))
end

### is not optimised for very stiff polymers, where some values will never occur
@inline function SetTrialBondAngle(data::SimData,param::Cosine_BondAngles)
    for i in 1:data.NTrials
        give_rand(param.AngleGenerator[data.id], data)
        data.trial_angle[i] =  data.rand_val
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


### Z[ k_] := Integrate[ Exp[-k*(1 - Cos[\[Theta]])]*Sin[\[Theta]], {\[Theta], 0, Pi}]
### P[\[Theta]_, k_] := Exp[-k*(1 - Cos[\[Theta]])]*Sin[\[Theta]]/Z[ k]
### AvgCos[k_] :=  Integrate[P[\[Theta], k]*Cos[\[Theta]], {\[Theta], 0, Pi}]
function AvgCos(K) #### coth(K)-1/K
    return  (coth(K)-1.0/K)
end



### sampler for random number of θ of gaussian K
#=struct Cosine_BondAngle_Prefactor_Sampler{I<:Integer, T<: Real} <: Sampleable{Univariate,Continuous}
    μ::T    ### normal distribution center 
    σ::T    ### normal distribution variance
    K_borders::Vector{T}### bins of the subdivisions
    K_cdf::Vector{T}    ### cumulated probability of bins
    K_max::Vector{} ### maximum of probability distribution in bin
    K_width::T  ### bin width
    max_val::T  ### maximum of the probability function of θ 
    P::Any      ### probability function of θ
    width::T    ### sampling width in κ
    offset::T   ### samling offset in κ
    Cosine_BondAngle_Prefactor_Sampler(μ::T, σ::T, P, K, width; NK=10::Int64, NA=10::Int64) where {T<:Real} = begin
        
        K_borders = collect(LinRange(μ-width*σ,μ+width*σ, NK+1)) ### get bins for K in [μ-width, μ+width]
        K_width= K_borders[2]-K_borders[1] ### bin width

        K_cdf = cdf.(Normal(μ, σ), K_borders[2:end]).-cdf(Normal(μ,σ), μ-width*σ) ### get likely hood for each bin assuming normal distributed κ
        K_max  = [max(pdf(Normal(μ, σ),K_borders[i]), pdf(Normal(μ, σ),K_borders[i+1]) ) for i in 1:NK] ### get maximal value of likelyhood for each bin 
        K_ind =  findfirst(μ.<=K_borders) -1 ### get bin with center of normal distribution
        K_max[K_ind] = max(pdf(Normal(μ, σ),μ), K_max[K_ind]) ### update maximal likely hood for the extrema of the function
        max_val = maximum([P(k,maxima_pos(K(k))) for k in  collect(LinRange(μ-1.0*σ,μ+1.0*σ, 100))]) ### maximum of the probability function of θ 
        new{Int64, T}(μ, σ,K_borders,K_cdf,K_max, K_width,max_val, P,width*2.0*σ,  μ-width*σ)
    end
end
=#

### sampler for random number of θ of gaussian K
struct Cosine_BondAngle_Prefactor_Sampler{T<: Real} <: Sampleable{Univariate,Continuous}
    #μ::T    ### normal distribution center 
    #σ::T    ### normal distribution variance
    K_borders::Vector{T}### bins of the subdivisions
    K_cdf::Vector{T}    ### cumulated probability of bins
    K_max::Vector{} ### maximum of probability distribution in bin
    K_width::T  ### bin width
    max_val_Pθ::T  ### maximum of the probability function of θ 
    P_Θ::Any      ### probability function of θ
    P_κ::Any      ### probability function of κ
    width::T    ### sampling width in κ
    offset::T   ### samling offset in κ
end


function give_rand( s::Cosine_BondAngle_Prefactor_Sampler)
    ### draw gaussian for K or lp
    rand_K  = _rand_off(s.width, s.offset) ### x*σ, μ-x*σ
    rnd_num = rand(eltype(s.K_cdf))*s.K_cdf[end]
    ind = 1
    @inbounds for i in axes(s.K_cdf,1) ### replicate findfirst(x->x>=rnd_num, s.K_cdf)
        @inbounds if s.K_cdf[i] > rnd_num
            ind = i
            break
        end
    end

    #ind = findfirst(x->x>=rnd_num, s.K_cdf)
    @inbounds rand_K = _rand_off(s.K_width, s.K_borders[ind])
    @inbounds while rand()*s.K_max[ind] > s.P_κ(rand_K)
        rand_K = _rand_off(s.K_width,  s.K_borders[ind])
    end

    ### draw angle accordingly
    rand_ϕ = rand(eltype(rand_K))*π
    @inbounds while rand()*s.max_val_Pθ >  s.P_Θ(rand_K,rand_ϕ)
        rand_ϕ = rand(eltype(rand_K))*π
    end
    return (rand_K, rand_ϕ)
end



### Same potential but, K's are drawn such that lp is Gaussian distributed
mutable struct GaussianK_Cosine_BondAngles{T<:Real} <: AbstractBondAngleParam
    μ::T
    σ::T
    Sampler::Cosine_BondAngle_Prefactor_Sampler
    GaussianK_Cosine_BondAngles(μ::T, σ::T; NK=10) where {T<:Real} = begin 
        max_width =min(5.0,μ/(σ*1.0001))

        K_borders = collect(LinRange(μ-max_width*σ,μ+max_width*σ, NK+1)) ### get bins for K in [μ-width, μ+width]
        K_width= K_borders[2]-K_borders[1] ### bin width

        K_cdf = Distributions.cdf.(Distributions.Normal(μ, σ), K_borders[2:end]).-Distributions.cdf(Distributions.Normal(μ,σ), μ-max_width*σ) ### get likely hood for each bin assuming normal distributed κ
        K_max  = [max(Distributions.pdf(Distributions.Normal(μ, σ),K_borders[i]), Distributions.pdf(Distributions.Normal(μ, σ),K_borders[i+1]) ) for i in 1:NK] ### get maximal value of likelyhood for each bin 
        K_ind =  findfirst(μ.<=K_borders) -1 ### get bin with center of normal distribution
        K_max[K_ind] = max(Distributions.pdf(Distributions.Normal(μ, σ),μ), K_max[K_ind]) ### update maximal likely hood for the extrema of the function
      

        P_k(κ) = Distributions.pdf(Distributions.Normal(μ, σ), κ) ### probability of κ
        ### P_θ(θ) = Cosine_BondAngle_func, P(θ|κ)

        max_val_Pθ = maximum([Cosine_BondAngle_func(k,maxima_pos(k)) for k in  collect(LinRange(μ-1.0*σ,μ+1.0*σ, 100))]) ### maximum of the probability function of θ 

        new{T}(μ, σ, Cosine_BondAngle_Prefactor_Sampler{T}(K_borders, K_cdf, K_max, K_width, max_val_Pθ,Cosine_BondAngle_func,  P_k,  max_width*2.0*σ,μ-max_width*σ ))
    end
end

### Init function doesnt allow <cos> values larger 1000 and smaller than ≈1/3
function InitSimParam(data::SimData,param::GaussianK_Cosine_BondAngles ) 
    nothing
end

### is not optimised for very stiff polymers, where some values will never occur
@inline function SetTrialBondAngle(data::SimData,param::GaussianK_Cosine_BondAngles)
    for i in 1:data.NTrials
        @inbounds data.trial_angle[i] = give_rand(param.Samplers)[2]
    end
    CompTrigonometricTrialBondAngles(data)
    nothing
end

@inline function GetTrialBoltzmannWeight(data::SimData,param::GaussianK_Cosine_BondAngles)
    nothing
end



##############################
### !!! doesnt really work since we integrate cos(θ) κ and θ simulataneously in the analysis
### Same potential but, K's are drawn such that lp is Gaussian distributed
#=
mutable struct GaussianLp_Cosine_BondAngles_Old{T<:Real} <: AbstractBondAngleParam
    μ::T
    σ::T
    Sampler::Cosine_BondAngle_Prefactor_Sampler
    GaussianLp_Cosine_BondAngles(μ::T, σ::T; NK=10) where {T<:Real} = begin
        κ_center=κ_from_lp(μ)# solveRecursive(lpToCosAngle(μ),AvgCos , 0.0001,10000.0 ) ### solve for κ for which lp=μ
       
        max_width =min(5.0,μ/(σ*1.0001)) ### max_width in lp space
        κ_max_width = κ_from_lp(max_width) #solveRecursive(lpToCosAngle(max_width),AvgCos , 0.0001,10000.0 )
        κ_max_width = min(κ_center, κ_max_width) ### ensure κ>=0 and restrict sampling to values in which lp>0 and lp in [μ-5*σ, μ+5*σ]

        println("κ_cen : $(κ_center), w: $(κ_max_width)")

        K_borders = collect(LinRange(κ_center-κ_max_width,κ_center+κ_max_width, NK+1)) ### distribution has longer tail for κ→0 ; get bins for K in [μ-width, μ+width]
        K_width= K_borders[2]-K_borders[1] ### bin width

        ### cdf not known analytically ⇒ determine numerically
        pers(κ, b) =  -b/log(coth(κ)-1.0/κ) ### lp(κ)
        P(κ,μ, σ, b=3.8) = Distributions.pdf.(Distributions.Normal(μ, σ),pers.(κ,b)) ### PDF of κ, not normalized
        P_κ(κ)=  P(κ,μ, σ, 3.8) #(P_κ)
        #normal.(pers.(κ,b),μ,σ)
        ### P_θ(θ) = Cosine_BondAngle_func, P(θ|κ)

        (norm_val, error) = HCubature.hcubature(P_κ,  [0.00000001], [κ_from_lp(μ+σ*max_width)*1.2]; rtol=10^-6, atol=10^-7, maxevals=1000, initdiv=100) ### compute normalisation


        ### compute cumulated probability of bin numerically
        K_cdf = [ HCubature.hcubature(P_κ, [left], [right]; rtol=10^-6, atol=10^-7, maxevals=1000, initdiv=100)[1][1] for (left, right) in zip(K_borders[1:end-1], K_borders[2:end])] 
        K_cdf ./=norm_val[1] ### norm cdf to 1

        K_max = [maximum(P_κ(LinRange(left, right,100))) for (left, right) in zip(K_borders[1:end-1], K_borders[2:end])] .*1.05 ### approximate maximum in bin and multiply by safety margin

        max_val_Pθ  = maximum([Cosine_BondAngle_func(k, θ) for k in LinRange(κ_center-κ_max_width,κ_center+κ_max_width, NK*11) for θ in LinRange(0,π, 21)])*1.10 ### approximation of  maximum of P(θ|κ) ∀ κ, θ with 10% safety margin
        #=
        ### make linear interpolation
        lp = collect(LinRange(μ-(max_width)*σ, μ+(max_width)*σ, 500)) ### get lp in [μ-w, μ+w]
        K_val= solveRecursive.(lpToCosAngle.(lp),AvgCos , 0.0001,10000.0 ) ### compute κ for the lp [μ-w, μ+w]
        K = linear_interpolation(lp, K_val) ### construct spline interpolator for lp ⇒ κ
        H(x1,x2) = Cosine_BondAngle_func(K(x1), x2) ### H(lp, θ) = P(θ|κ(lp))
       
        new{T}(μ, σ,  Cosine_BondAngle_Prefactor_Sampler(μ, σ, H, K, max_width)) 
         =#
        new{T}(μ, σ, Cosine_BondAngle_Prefactor_Sampler{T}(K_borders, K_cdf, K_max, K_width, max_val_Pθ , Cosine_BondAngle_func, P_κ, 2*κ_max_width, κ_center-κ_max_width))
    end
end
=#

function getKappaBorders(μ::T, σ::T; b=3.8 ) where {T<:Real}
    κ_center=κ_from_lp(μ)### solve for κ for which lp=μ
       
    max_width =min(5.0,μ/(σ*1.000000001)) ### max_width in lp space
    κ_max_width = κ_center-κ_from_lp(μ-max_width*σ)-10^-2 ### sampling length in κ, numerical stability gone if κ≈0
    κ_max_width = min(κ_center, κ_max_width) ### ensure κ>=0 and restrict sampling to values in which lp>0 and lp in [μ-5*σ, μ+5*σ]

    ### need asymetric borders in κ to have symmetric borders in lp
    ### distribution has longer tail for κ→0 ; get bins for K ∈ [left, right] such that lp ∈ [μ-w*σ, μ+w*σ]
    left = κ_center-κ_max_width
    right = μ+ (μ- -b/log(AvgCos(κ_center-κ_max_width))) ### in units of lp
    right = κ_from_lp(right)

    return left, right, max_width
end


mutable struct GaussianLp_Cosine_BondAngles{T<:Real} <: AbstractBondAngleParam
    μ::T
    σ::T
    κ_vec::Vector{T}
    Sampler::Vector{Cosine_BondAngle_Sampler}
    GaussianLp_Cosine_BondAngles(μ::T, σ::T,NBeads::I; NK=10, b=3.8) where {I<:Integer,T<:Real} = begin

        left, right, max_width  = getKappaBorders(μ, σ; b=b)
        lp_left = -b/log(AvgCos(left))
        lp_right= -b/log(AvgCos(right))

        ### make linear interpolation
        lp = collect(LinRange(lp_left, lp_right, 500)) ### get lp in [μ-w, μ+w]
        K_val= κ_from_lp.(lp) ### compute κ for the lp [μ-w, μ+w]
        K_func = linear_interpolation(lp, K_val) ### construct spline interpolator for lp ⇒ κ

        #println(" $(lp_left), $(lp_right)")

        ### linear interpolation of derivative of κ(lp)
        lp_mean = (lp[2:end].+lp[1:end-1])./2.0
        κ_der = (K_func.(lp[2:end])-K_func.(lp[1:end-1]))./2.0
        κ_der_func = linear_interpolation(lp_mean, κ_der,extrapolation_bc=Line()) ### construct spline interpolator for derivative of lp(κ) ⇒ κ

        pers(κ, b) =  -b/log(coth(κ)-1.0/κ) ### lp(κ)
        P(κ,μ, σ, b=3.8) = Distributions.pdf.(Distributions.Normal(μ, σ),pers.(κ,b)) * (1.0/κ_der_func(pers(κ,b))) ### PDF of κ, not normalized

        #P(κ,μ, σ, b=T(3.8)) = 1/(sqrt(2*π)*σ^2)*exp(-(μ-pers(κ,b))^2/(2*σ^2))  ### PDF of κ, not normalized

        P_κ(κ)=  P(κ,μ, σ, 3.8) #(P_κ)
        max_val_Pκ  = maximum([P_κ(k) for k in LinRange{T}(left*1.01, right*0.99, NK*101)])*1.10 ### approximation of  maximum of P(κ) ∀ κ with 10% safety margin


        κ_total_width=right-left
        ### draw κ for each triplet according to distribution of κ
        κ_vec = zeros(T, NBeads)
        for i in 1:NBeads-2
            rand_κ= rand(T)*κ_total_width+left  ### draw random κ in [left,right]
            while P_κ(rand_κ)< rand(T)*max_val_Pκ ### accept according to prob. of κ
                rand_κ= rand(T)*κ_total_width+left ### draw random κ in [center-width, center+width]
            end
            κ_vec[i] = rand_κ 
        end

        new{T}(μ, σ,κ_vec, [Cosine_BondAngle_Sampler(κ) for κ in κ_vec])
    end
end

mutable struct GaussianFixedK_Cosine_BondAngles{T<:Real} <: AbstractBondAngleParam
    μ::T
    σ::T
    κ_vec::Vector{T}
    Sampler::Vector{Cosine_BondAngle_Sampler}
    GaussianFixedK_Cosine_BondAngles(μ::T, σ::T,NBeads::I; NK=10, b=3.8) where {I<:Integer,T<:Real} = begin
        P_κ= Distributions.Normal(μ, σ)
        κ_vec = rand(P_κ, NBeads) 
        new{T}(μ, σ,κ_vec, [Cosine_BondAngle_Sampler(κ) for κ in κ_vec])
    end
end

### Init function doesnt allow <cos> values larger 1000 and smaller than ≈1/3
function InitSimParam(data::SimData,param::GaussianLp_Cosine_BondAngles ) 
    nothing
end

function InitSimParam(data::SimData,param::GaussianFixedK_Cosine_BondAngles ) 
    nothing
end

### is not optimised for very stiff polymers, where some values will never occur
@inline function SetTrialBondAngle(data::SimData,param::Union{GaussianLp_Cosine_BondAngles,GaussianFixedK_Cosine_BondAngles})
    for i in 1:data.NTrials
        give_rand(param.Sampler[data.id-2], data) ### writes into data.rand_Val
        data.trial_angle[i] = data.rand_val
    end
    CompTrigonometricTrialBondAngles(data)
    nothing
end

function GetTrialBoltzmannWeight(data::SimData,param::GaussianLp_Cosine_BondAngles)
    nothing
end
