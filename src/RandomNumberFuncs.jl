using Random

export _rand, _rand2π, _rand_off
### Random numbers
_rand(len::T)     where {T<:Real} = rand(T) * len
_rand_off(len::T) where {T<:Real} = rand(T) * len - len / 2.0
_rand_off(len::T,A::T) where {T<:Real} = rand(T) * len + A
_rand2π(len::T)   where {T<:Real} = _rand(2.0*T(π))