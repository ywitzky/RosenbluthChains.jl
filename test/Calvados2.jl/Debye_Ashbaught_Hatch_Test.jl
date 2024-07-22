
e_charge = 0.0
function Debye(r, D, q1,q2,ϵ0, ϵr)
    return q1*q2*e_charge^2/(4*π*ϵr*ϵ0)*exp(-r/D)/r
end