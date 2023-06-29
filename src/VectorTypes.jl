import Base: +,-,*,/,show, angle
const Vector3{T} = MVector{3, T}
const Matrix3{T} = MArray{Tuple{3,3},T}

@inline sqr_norm(vec::Vector3{T})  where {T<:Real} = vec[1]*vec[1]+vec[2]*vec[2]+vec[3]*vec[3]
@inline norm(vec::Vector3{T}) where {T<:Real} = sqrt(vec[1]*vec[1]+vec[2]*vec[2]+vec[3]*vec[3])

Vector{Vector3{T}}(N::I) where {I<:Integer, T<:Number} = [Vector3{T}(0,0,0) for _ in 1:N]

Matrix3{T}() where {T<:Real} =  Matrix3{T}(zeros(T,3,3))
    
### overload for speed up
#@inline function *(a::Matrix3{T}, b::Vector3{T}) where{T<:Real} a.*b end

@inline function -( a::Vector3{T},b::Vector3{T}, c::Vector3{T}) where{T<:Real}
    c[1] =a[1]-b[1]
    c[2] =a[2]-b[2]
    c[3] =a[3]-b[3]
    nothing
end

### overload for speed up
@inline function *(a::Matrix3{T}, b::Vector3{T}, c::Vector3{T}) where{T<:Real}
    @inbounds c[1] = a[1]*b[1]+a[2]*b[2]+a[3]*b[3]
    @inbounds c[2] = a[4]*b[1]+a[5]*b[2]+a[6]*b[3]
    @inbounds c[3] = a[7]*b[1]+a[8]*b[2]+a[9]*b[3]
    nothing
end


@inline function ×( a::Vector3{T},b::Vector3{T}, c::Vector3{T}) where{T<:Real}
    c[1] =a[2]*b[3]-a[3]*b[2]
    c[2] =a[3]*b[1]-a[1]*b[3]
    c[3] =a[1]*b[2]-a[2]*b[1]
    nothing
end

function angle(a::Vector3{T},b::Vector3{T}) where {T<:Number} 
   return acos(sum(a.*b)/(norm(a)*norm(b)))
end

function show(a::Matrix3{Any}) #where {T<:Real} 
    "$(a[1]) $(a[2]) $(a[3])\n$(a[4]) $(a[5]) $(a[6])\n$(a[7]) $(a[8]) $(a[9])"
end
