export Debye_Ashbaugh_Hatch

### combines Debye potential and ashbaugh hatch potential
struct Debye_Ashbaugh_Hatch{T<:Real} <: AbstractSelfAvoidanceParameters
    ah_ϵ::T #### saved as 4*epsilon
    ah_λ::Array{T}
    ah_sqr_cutoffs::Array{T} ### 2^(1/6)*σ
    ah_σ_sqr_mean::Array{T}
    ah_rcut_sqr::T ### 20AA usually
    ah_λu_cut::Array{T}
    ah_shift::Array{T}
    deb_rcut_sqr::T
    deb_qq::Array{T}
    deb_D::T
   
    Debye_Ashbaugh_Hatch(Sequence::String , TypeToSigma::Dict{Char, T}, TypeToLambda::Dict{Char, T}, ϵ_ah::T, rc_ah::T, TypeToCharge::Dict{Char, T}, rc_deb::T, cs_deb::T, rel_Temp::T) where {T<:Real} =begin
        Types = Vector{Char}(Sequence)

        σ_mean = zeros(T, (length(Sequence),length(Sequence)))
        λ_mean = zeros(T, (length(Sequence),length(Sequence)))
        qq = zeros(T, (length(Sequence),length(Sequence)))
        for (id,type) in enumerate(Types)
            for (id2,type2) in enumerate(Types)
                σ_mean[id, id2] = (TypeToSigma[type]+TypeToSigma[type2])/2.0
                λ_mean[id, id2] = (TypeToLambda[type]+TypeToLambda[type2])/2.0
                qq[id, id2] = TypeToCharge[type]*TypeToCharge[type2]
            end
        end

        ϵ_ah_l =ϵ_ah/rel_Temp 
        ulj(ϵ,σ,r) = 4*ϵ*((σ/r)^12-(σ/r)^6)

        λu_cut  = λ_mean.*ulj.(ϵ_ah_l, σ_mean, rc_ah)
        ah_shift = @.ϵ_ah_l*(1-λ_mean)
        sqr_cutoffs = 2.0^(1.0/6.0).*σ_mean

        ### compute exponents by hand to avoid float effects
        Temp = 298.0*rel_Temp

        kt = 4.11 # 10^-21 J at 298K need to rescale if temperature changes
        e = 1.602 # 10-19 C
        ϵ0 = 8.854 # 10^-12  As/(Vm)
        NA = 6.02214086 #* 10.0^23 ## * 10.0^23 

        ϵ_rdeb = 5321.0/Temp+233.76-0.9297*Temp+1.417*10.0^-3*Temp^2-8.292*10.0^-7*Temp^3
        B = e^2/(4*π*ϵ0*ϵ_rdeb*kt*rel_Temp) ### 10^-5 m 
        D = sqrt(1.0/(8.0*π*B*cs_deb*NA*1000.0))*10 ### convert to Angstroem, paper uses ionic strengt per M, so per Liter = 1000 mol/m^3

        qq .*= e^2/(4*π*ϵ0*ϵ_rdeb*(kt*rel_Temp))*10^5 ### correct to Angstroem

        return new{T}(4*ϵ_ah_l,λ_mean,sqr_cutoffs.^2, σ_mean.^2, rc_ah^2, λu_cut, ah_shift, rc_deb^2, qq, D)
    end
end


function InitSimParam(data::SimData,param::Debye_Ashbaugh_Hatch ) 
    nothing
end

### "naiv implementation" (n^2 loop) seemed to be the fastest for reasonable lenght polymers/proteins when developing LJ_Repulsion
### we aim at SIMD operations for all trial positions
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
            data.tmp4 .+= data.btmp3 .* (param.deb_qq[id, data.id] .*exp.(-data.tmp5./param.deb_D)./data.tmp5 ) 
        end

        ### compute the ashbaugh hatch potential
        if any(data.btmp3.!=0.0) ### skip all computations if no distance is smaller than largest cutof
            data.tmp5 .= param.ah_sqr_cutoffs[id, data.id].>data.btmp ### mask for cutoff 2^(1/6)*σ, 0 if larger than cutoff
            data.btmp2 .= (param.ah_rcut_sqr.> data.btmp).-data.tmp5 ### mask for cutoff r in [2^(1/^)*σ, ah_rc], 0 else
     
            data.btmp .= param.ah_σ_sqr_mean[id, data.id]./data.btmp ### ^2
            data.btmp .= data.btmp.*data.btmp.*data.btmp ### ^6
            data.btmp .= param.ah_ϵ .*(data.btmp.*data.btmp .-data.btmp) ### U_LJ(r), ϵ is actually 4*ϵ

            data.tmp4 .+= data.tmp5.* (data.btmp.- param.ah_λu_cut[id, data.id] .+ param.ah_shift[id, data.id]) ### add everything smaller r<2^(1/6)σ
            
            data.tmp4 .+=  data.btmp2 .* (data.btmp.*param.ah_λ[id, data.id] .- param.ah_λu_cut[id, data.id] ) ### add everything in [2^(1/^)*σ, ah_rc],
        end
        
    end
    data.LogBoltzmannFaktor .-=  data.tmp4
    nothing
end

