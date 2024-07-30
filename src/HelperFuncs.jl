export lpToRigidAngle,solveRecursive

function lpToRigidAngle(lps, r0=3.8)
    return @. acos(exp(-r0/lps))
end

function lpToCosAngle(lps, r0=3.8)
    return @. exp(-r0/lps)
end

function InvlpToCosAngle(InvLps, r0=3.8)
    return @. exp(-r0*InvLps)
end

### Solve Recursive for Value that solves Func(Value)=Target
function solveRecursive(Target::T, Func::Function, Left::T, Right::T)where {T<:Real}
    Middle = (Left+Right)/2.0
    #println("Left: $Left, Right: $Right, F(middle)=$(Func(Middle)), Target: $Target")
    if abs(Func(Middle)-Target)<10^-14 return Middle end
    if abs(Middle-Left)<10^-14 || abs(Middle-Right)<10^-14 return Middle end
    if Func(Middle)<Target
        return solveRecursive(Target, Func, Middle, Right)
    else
        return solveRecursive(Target, Func, Left, Middle)
    end
end

@inline function CompTrigonometricTrialBondAngles(data::SimData)
    data.cos_trial_angle .= cos.(data.trial_angle)
    data.sin_trial_angle .=  @. sign(data.trial_angle-π).*  √(1-data.cos_trial_angle^2) ## #only works on interval [0, 2π]@.  √(1-data.cos_trial_angle^2)
    nothing
end

@inline function CompTrigonometricTrialTorsionAngles(data::SimData)
    data.cos_trial_torsion_angle .= cos.(data.trial_torsion_angle)
    data.sin_trial_torsion_angle .= @. sign(data.trial_torsion_angle-π).* √(1-data.cos_trial_torsion_angle^2) #only works on interval [0, 2π]  #@. √(1-data.cos_trial_torsion_angle^2)
    nothing
end