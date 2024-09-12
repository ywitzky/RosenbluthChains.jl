export lpToRigidAngle,solveRecursive

using LoopVectorization

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
    data.sin_trial_angle .= @. sign(π.-data.trial_angle).*  √(1-data.cos_trial_angle^2) ## #only works on interval [0, 2π]@.  √(1-data.cos_trial_angle^2)
    nothing
end

@inline function CompTrigonometricTrialTorsionAngles(data::SimData)
    data.cos_trial_torsion_angle .= cos.(data.trial_torsion_angle)
    data.sin_trial_torsion_angle .= @. sign(π.-data.trial_torsion_angle).* √(1-data.cos_trial_torsion_angle^2) #only works on interval [0, 2π]  #@. √(1-data.cos_trial_torsion_angle^2)
    nothing
end

function getPotentialNeighbors(data::SimData)
    ### 2d or 3d hashing are equally fast for low number of beads => no need to not do 3d hashing 


    max_bondlength=data.max_interaction_length

    @inbounds x_max = floor.(Int32, (data.xyz[data.id-1][1]+max_bondlength)) 
    @inbounds x_min = floor.(Int32, (data.xyz[data.id-1][1]-max_bondlength)) 

    #LoopVectorization.@avx data.id_arr[1:data.id-2] .=  floor.(Int32, (data.x[1:data.id-2]))

    ### not using the above version removes the allocations
    LoopVectorization.@avx data.id_arr .=  floor.(Int32, (data.x))

    #@inbounds floor_f(data.id_arr[1:data.id-2] ,data.x[1:data.id-2], Int32)

    @inbounds y_max = floor.(Int32, (data.xyz[data.id-1][2]+max_bondlength)) 
    @inbounds y_min = floor.(Int32, (data.xyz[data.id-1][2]-max_bondlength)) 

    #LoopVectorization.@avx  data.id_arr2[1:data.id-2] .=  floor.(Int32, (data.y[1:data.id-2])) 

    LoopVectorization.@avx data.id_arr2 .=  floor.(Int32, (data.y))


    @inbounds z_max = floor.(Int32, (data.xyz[data.id-1][3]+max_bondlength)) 
    @inbounds z_min = floor.(Int32, (data.xyz[data.id-1][3]-max_bondlength)) 


    #LoopVectorization.@avx data.id_arr3[1:data.id-2] .=  floor.(Int32, (data.z[1:data.id-2]))

    LoopVectorization.@avx data.id_arr3 .=  floor.(Int32, (data.z))


    cnt= 0
    @inbounds for i in 1:data.id-2
        @inbounds if data.id_arr[i]>=x_min && data.id_arr[i]<=x_max
            @inbounds if data.id_arr2[i]>=y_min && data.id_arr2[i]<=y_max 
                @inbounds if data.id_arr3[i]>=z_min && data.id_arr3[i]<=z_max
                    cnt += 1
                    @inbounds data.id_arr4[cnt] =  i
                end
            end
        end
    end
    return @view data.id_arr4[1:cnt]

    #return  [id for (id, (xid, yid, zid)) in enumerate(zip(data.id_arr[1:data.id-2], data.id_arr2[1:data.id-2], data.id_arr3[1:data.id-2])) if xid>=x_min && xid<=x_max &&  yid>=y_min && yid<=y_max  && zid>=z_min && zid<=z_max]
end
