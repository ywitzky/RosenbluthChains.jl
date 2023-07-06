export lpToRigidAngle,solveRecursive

function lpToRigidAngle(lps, r0=3.8)
    return @. acos(exp(-r0/lps))
end

function lpToCosAngle(lps, r0=3.8)
    return @. exp(-r0/lps)
end

### Solve Recursive for Value that solves Func(Value)=Target
function solveRecursive(Target::T, Func::Function, Left::T, Right::T)where {T<:Real}
    
    Middle = (Left+Right)/2.0
    #println("Left: $Left, Right: $Right, F(middle)=$(Func(Middle)), Target: $Target")
    if abs(Func(Middle)-Target)<10^-10 return Middle end
    if abs(Middle-Left)<10^-10 || abs(Middle-Right)<10^-10 return Middle end
    if Func(Middle)<Target
        return solveRecursive(Target, Func, Left, Middle)
    else
        return solveRecursive(Target, Func, Middle, Right)
    end
end
