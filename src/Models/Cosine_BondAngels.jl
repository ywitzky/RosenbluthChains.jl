struct Cosine_BondAngles{T<:Real} <: AbstractBondAngleParam
    PersistenzLengths::Vector{T} 
    K::Vector{T} 
    FixedBondAngles(BondAngles::Vector{T}) where {T<:Real} = new{T}(BondAngles)
 end
