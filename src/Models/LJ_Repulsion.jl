export LJ_Repulsion

struct LJ_Repulsion{T<:Real} <: AbstractSelfAvoidanceParameters
    ϵ::T #### saved as 4*epsilon
    sqr_cutoffs::Array{T}
    σ_sqr_mean::Array{T}
    LJ_Repulsion(BeadType::Vector{I}, TypeToSigma::Dict{I, T}, ϵ::T) where {T<:Real, I<:Integer} =begin
        σ = zeros(T, length(BeadType))
        σ_mean = zeros(T, (length(BeadType),length(BeadType)))
        for (id, type) in enumerate(BeadType)
            σ[id] = TypeToSigma[type]
        end


        for (id, type) in enumerate(BeadType)
            for (id2, type2) in enumerate(BeadType)
                σ_mean[id, id2] = (σ[id]+σ[id2])/2.0
            end
        end
        sqr_cutoffs = 2^(1.0/6.0).*σ_mean
        return new{T}(4*ϵ,sqr_cutoffs.^2, σ_mean.^2 )
    end
end

function InitSimParam(data::SimData,param::LJ_Repulsion ) 
    nothing
end

###naiv implementation
@inline function GetTrialBoltzmannWeight(data::SimData,param::LJ_Repulsion)
    fill!(data.tmp4, 0.0) ### saves energy for every test position  
    for id in 1:data.id-1
        for tid in 1:data.NTrials
        data.btmp[tid] = sqr_norm(data.trial_positions[tid]-data.xyz[id])
        end ### square distance
        #println("id: $id , $(data.btmp)")
        data.tmp5 .= param.sqr_cutoffs[id, data.id].>data.btmp ### mask for cutoff, 0 if larger than cutoff
        data.btmp .= param.σ_sqr_mean[id, data.id]./data.btmp ### ^2
        data.btmp .= data.btmp.*data.btmp.*data.btmp ### ^6
        data.tmp4 .+= data.tmp5.* param.ϵ .*(data.btmp.*data.btmp .-data.btmp.+0.25)### +0.25 is offset at cutoff, ϵ is actually 4*ϵ
    end
    #println("energy: $(data.tmp4) ")
    data.LogBoltzmannFaktor .-=  data.tmp4
    nothing
end

@inline function Get3DHash(Position::Vector3{T}, Spacing::T) where {T<: Real}
    return Int32.(floor.(Position./Spacing))
end


### naiv implementation, probably possibilities to speed up
mutable struct CellList{I <: Integer}
    dict::Dict{Vector3{I}, Vector{I}}
    HashFunc::Function
    offsets::Vector{Vector3{I}}
    neighbourcells::Vector{Vector3{I}}
    CellList(GridSpacing::T) where{T<:Real} = new{Int32}(Dict{Vector3{Int32}, Vector{Int32}}(),     x-> Get3DHash(x, GridSpacing), [ [x,y,z] for x in [-1,0,1] for y in [-1,0,1] for z in [-1,0,1]], [[0,0,0] for i in 1:27])

    add(CL, pos, id) = begin
        ind = HashFunc(pos)
        if ind in keys(CL.dict)
            append!(CL.dict[ind], id)
        else
            CL.dict[ind]=[id]
        end
    end
    neighbours(CL, pos) = begin
        ind = CL.HashFunc(pos)
        CL.neighbourcells .= ind .+ CL.offsets
        return [getkey(CL.dict, i, empty) for i in CL.neighbourcells if i in keys(CL.dict)]
    end
end

function clear(CL::CellList{I})where{I<:Integer} 
    CL.dict = Dict{Vector3{I}, Vector{I}}()
end