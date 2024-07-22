export Debye_Ashbaugh_Hatch

### combines Debye potential and ashbaugh hatch potential
struct Debye_Ashbaugh_Hatch{T<:Real, I<:Integer} <: AbstractSelfAvoidanceParameters
    
    ah_ϵ::T #### saved as 4*epsilon
    ah_λ::Array{T}
    ah_sqr_cutoffs::Array{T} ### 2^(1/6)*σ
    ah_σ_sqr_mean::Array{T}
    ah_types::Vector{I}
    ah_rcut_sqr::T
    ah_λu_cut::Array{T}
    ah_shift::Array{T}
    deb_pref::T
    deb_rcut_sqrs::T
    deb_qq::Array{T}
    deb_D::T
    LJ_Repulsion(Sequence::String , TypeToSigma::Dict{Char, T}, TypeToLambda::Dict{Char, T}, ϵ_ah::T, TypeToCharge::Dict{Char, T}, ϵ_rdeb::T, rc_ah::T) where {T<:Real, I<:Integer} =begin
        Types = Vector{Char}(Sequence)

        Unique_Types = Set(Types)
        NTypes = length(Unique_Types)

        TypeToID = Dict(Dict([(y,x) for (x,y) in enumerate(Unique_Types)]))

        λ = [TypeToLambda[c] for c in Sequence]

        ah_BeadType = [TypeToID[c] for c in Sequence]

        σ_mean = zeros(T, (length(Sequence),length(Sequence)))
        λ_mean = zeros(T, (length(Sequence),length(Sequence)))
        qq = zeros(T, (length(Sequence),length(Sequence)))
        for (id,type) in enumerate(Types)
            for (id,type2) in enumerate(Types)
                σ_mean[id, id2] = (TypeToSigma[type]+TypeToSigma[type2])/2.0
                λ_mean[id, id2] = (TypeToLambda[type]+TypeToLambda[type2])/2.0
                qq[id, id2] = TypeToCharge[type]*TypeToCharge[type2]
            end
        end

        ulj(ϵ,σ,r) = 4*ϵ*((σ/r)^12-(σ/r)^6)

        λu_cut  = λ_mean.*ulj.(ϵ_ah, σ_mean, rc_ah)
        ah_shift = @.ϵ_ah*(1-λ_mean)

        sqr_cutoffs = 2.0^(1.0/6.0).*σ_mean
        return new{T,I}(4*ϵ,sqr_cutoffs.^2, σ_mean.^2, CellLinkedList(maximum(sqr_cutoffs)) )
    end
end


function InitSimParam(data::SimData,param::Debye_Ashbaugh_Hatch ) 
    nothing
end

### "naiv implementation" (n^2 loop) seemed to be the fastest for reasonable lenght polymers/proteins when developing LJ_Repulsion
@inline function GetTrialBoltzmannWeight(data::SimData,param::Debye_Ashbaugh_Hatch)
    fill!(data.tmp4, 0.0) ### saves energy for every test position  

    for id in 1:data.id-2
        for tid in 1:data.NTrials
            data.btmp[tid] = sqr_norm(data.trial_positions[tid]-data.xyz[id]) ### square distance
        end
        data.btmp3 .= param.deb_rcut_sqr.>data.btmp ### mask for debye cutoff

        ### compute the debye potential
        if param.deb_qq[id, data.id] != 0.0
            data.tmp5 .= sqrt.(data.btmp) # r
            data.tmp4 .+=  data.btmp3 * (param.deb_pref .* param.deb_qq[id, data.id] .*exp.(-data.tmp5./param.deb_D)./data.tmp5 ) 
        end

        ### compute the ashbaugh hatch potential
        if any(data.btmp3) ### skip all computations if no distance is smaller than largest cutof
            data.tmp5 .= param.sqr_cutoffs[id, data.id].>data.btmp ### mask for cutoff 2^(1/6)*σ, 0 if larger than cutoff
            data.btmp2 .= param.ah_rcut_sqr.>data.btmp.-data.tmp5 ### mask for cutoff r in [2^(1/^)*σ, ah_rc], 0 else
     
            data.btmp .= param.σ_sqr_mean[id, data.id]./data.btmp ### ^2
            data.btmp .= data.btmp.*data.btmp.*data.btmp ### ^6
            data.btmp .= param.ϵ .*(data.btmp.*data.btmp .-data.btmp) ### U_LJ(r), ϵ is actually 4*ϵ

            data.tmp4 .+= data.tmp5.* (data.btmp.- ah_λu_cut[id, data.id] .+ ah_shift) ### add everything smaller r<2^(1/6)σ
            
            data.tmp4 .+=  data.btmp2 * (data.btmp.*param.ah_λ[id, data.id]  .- .- ah_λu_cut[id, data.id] ) ### add everything smaller in [2^(1/^)*σ, ah_rc],
        end
    end
    data.LogBoltzmannFaktor .-=  data.tmp4
    nothing
end

