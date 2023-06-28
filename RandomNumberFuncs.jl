### Random numbers
_rand(len::T)     where {T<:Real} = rand(T) * len
_rand_off(len::T) where {T<:Real} = rand(T) * len - len / 2.0
_rand_off(A::T, len::T) where {T<:Real} = () -> rand(T) * len - len / 2.0
_rand2π(len::T)   where {T<:Real} = _rand(2.0*T(π))