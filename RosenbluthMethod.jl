function RunSim(data::SimData, param::SimulationParameters)
### Init
LA_Obj=InitLammpsAnalysisScript(data, param)
InitSimParam(data, param)
Measurements = InitMeasurementStorage(data, param)

### Computational Part
mainLoop(data, param, Measurements, LA_Obj)

### Write Out
MeasurementsToPickle(data, Measurements)
return Measurements
end

function CalcRotationMatrix(out_vec::Vector3{T},angle::T, data::SimData, param::SimulationParameters) where {T<:Real}
    #https://en.wikipedia.org/wiki/Rotation_matrix
    #section Rotation matrix from axis and angle
    c = cos(angle)
    f = 1-c
    s = sin(angle)
    length = norm(out_vec)
    data.tmp1 .=  out_vec./length

    data.RotationMatrix[1]=c+data.tmp1[1]*data.tmp1[1]*f
    data.RotationMatrix[2]=  data.tmp1[1]*data.tmp1[2]*f-data.tmp1[3]*s
    data.RotationMatrix[3]=  data.tmp1[1]*data.tmp1[3]*f+data.tmp1[2]*s
    
    data.RotationMatrix[4]=  data.tmp1[2]*data.tmp1[1]*f+data.tmp1[3]*s
    data.RotationMatrix[5]=c+data.tmp1[2]*data.tmp1[2]*f
    data.RotationMatrix[6]=  data.tmp1[2]*data.tmp1[3]*f-data.tmp1[1]*s

    data.RotationMatrix[7]=  data.tmp1[3]*data.tmp1[1]*f-data.tmp1[2]*s
    data.RotationMatrix[8]=  data.tmp1[3]*data.tmp1[2]*f+data.tmp1[1]*s
    data.RotationMatrix[9]=c+data.tmp1[3]*data.tmp1[3]*f
end

function ComputeTrialPositions(data::SimData, param::SimulationParameters)
    for n in 1:data.NTrials
        CalcRotationMatrix(data.crossproduct,data.trial_angle[n] , data, param)
        *(data.RotationMatrix,data.current , data.tmp2 )
        
        CalcRotationMatrix(data.current ,_rand(eltype(data.TType)(2*π) ) , data, param )
        *(data.RotationMatrix,data.tmp2, data.new_vec)
        
        length = norm(data.new_vec)

        data.tmp3 .= data.new_vec.* (data.trial_radius[n]/length)

        data.trial_positions[n][:] .= (data.xyz[data.id-1].+ data.tmp3)
    end
end

function mainLoop( data::SimData, param::SimulationParameters, Measurements::Dict{String, Array{T}}, LA_Obj::LA.SimData{T2, N}) where {T<:Real, N<:Integer, T2<:Real}
    for I in 0:data.NumberOfBatches-1
        for i in 1:data.BatchSize
            ResetSim(data, param)
            SetFirstThreeBeads(data, param)
            ComputeBeadsIteratively(data,param)
            CopyXYZToLA(data,param, LA_Obj,i)
        end
        DoBatchMeasurements(data, param, Measurements, LA_Obj, I)
    end
end

function ResetSim(data::SimData, param::SimulationParameters)
    data.id=1
    data.tid=1
end

function ComputeBeadsIteratively(data::SimData, param::SimulationParameters)
    while(data.id<=data.NBeads)
        SetTrialRadius(data, param)
        SetTrialBondAngle(data, param)
        ComputeTrialPositions(data,param)

        #CalculateExternalPotential(data, param);

        ChooseTrialPosition(data, param)
        
        data.xyz[data.id] .= data.trial_positions[data.tid]
        
        data.tmp1 .= data.xyz[data.id].-data.xyz[data.id-1]
        ×( data.xyz[data.id-2]-data.xyz[data.id-1], data.tmp1,data.crossproduct)
        data.current .= data.tmp1

        data.id+=1
    end
end