using ProgressMeter

export RunSim, getRosenbluthWeigth

function RunSim(data::SimData, param::SimulationParameters, TmpMeas::AbstractMeasurement, PERM=NoPERM()) 
    InitSimParam(data, param)
    Measurements = InitMeasurement(data, param, TmpMeas)

    ### Computational Part
    mainLoop(data, param, Measurements, PERM)

    ### Write Out
    SaveMeasurements(data, param, Measurements)
    return Measurements
end

function mainLoop( data::SimData, param::SimulationParameters, Measurement::AbstractMeasurement, perm::NoPERM)

    for data.batch_id in 0:data.NumberOfBatches-1
        #println("Batch $(data.batch_id+1) /$(data.NumberOfBatches)")
        @showprogress 1 "Computing batch $(data.batch_id+1)" for data.id_in_batch in 1:data.BatchSize#, printing_delay=0.1
            ResetSim(data, param)
            SetFirstThreeBeads(data, param)
            ComputeBeadsIteratively(data,param, perm)
            MeasureAfterChainGrowth(data, param,Measurement)
        end
        MeasureAfterBatch(data, param, Measurement)
    end
end

function ComputeTrialPositions(data::SimData, param::SimulationParameters)
    for n in 1:data.NTrials       
        u1=data.crossproduct/norm(data.crossproduct)
        u2=Vector3{Float64}(0,0,0)
        ×( data.current/norm(data.current),-1.0 .*u1, u2)

        data.new_vec.=(
            (u2.*data.cos_trial_torsion_angle[n]
            .+ u1.*data.sin_trial_torsion_angle[n]).* data.sin_trial_angle[n]
            .+data.current./norm(data.current).*data.cos_trial_angle[n])          .*data.trial_radius[n]

        data.trial_positions[n][:] .= (data.xyz[data.id-1].+ data.new_vec)
    end
end

function getRosenbluthWeigth(data::SimData, param::SimulationParameters)
    return deepcopy(data.RosenbluthWeight) # = exp(data.LogRosenbluthWeight) # exp(BigFloat(data.LogRosenbluthWeight))
end

function ResetSim(data::SimData, param::SimulationParameters)
    data.id=1
    data.tid=1
    data.RosenbluthWeight=1.0
    data.LogRosenbluthWeight=0.0
    clear(data, param.SAWParam)
end

function SetFirstThreeBeads(data::SimData, param::SimulationParameters)
    ### First bead always at (0,0,0)
    SetTrialRadius(data, param)
    GetTrialBoltzmannWeight(data, param)
    ChooseTrialPosition(data, param)
    
    ### place second bead somewhere in spherical coordinates
    φ = _rand(eltype(data.TType)(2*π) )
    θ = acos(_rand_off(2.0, -1.0))   ### polar angle to dont oversample the poles.
    
    data.xyz[2][1]=  data.trial_radius[data.tid]*sin(θ)*cos(φ)
    data.xyz[2][2]=  data.trial_radius[data.tid]*sin(θ)*sin(φ)
    data.xyz[2][3]=  data.trial_radius[data.tid]*cos(θ)

    ### go for third bead.
    data.id=3

    SetTrialRadius(data, param)
    SetTrialBondAngle(data,param)
    SetTrialTorsionAngle(data, param)

    data.crossproduct[1] = -sin(φ)
    data.crossproduct[2] = cos(φ)
    data.crossproduct[3] = 0.0 
    
    data.current .= data.xyz[2]

    ComputeTrialPositions(data,param)
    ChooseTrialPosition(data, param)
    data.xyz[3] .= data.trial_positions[data.tid]

    data.tmp1 .= data.xyz[3].-data.xyz[2]
    ×(data.xyz[2]-data.xyz[1],data.tmp1,data.crossproduct)
    data.current .= data.tmp1
    data.id+=1
end

function ComputeBeadsIteratively(data::SimData, param::SimulationParameters, perm::AbstractPERMData)
    while(data.id<=data.NBeads)
        SetTrialRadius(data, param)
        SetTrialBondAngle(data, param)
        SetTrialTorsionAngle(data, param)
        ComputeTrialPositions(data,param)

        #CalculateExternalPotential(data, param);

        ChooseTrialPosition(data, param, perm)
        

        if any(isnan.(data.xyz[data.id]))
            println(data.tid, "  ", data.id)
            println("trial angle $(data.trial_angle), \ncos $(data.cos_trial_angle) \nsin $(data.sin_trial_angle) \n test:$(data.sin_trial_angle.^2 .+ data.cos_trial_angle.^2)")
            println("trial  torsion angle $(data.trial_torsion_angle), \ncos $(data.cos_trial_torsion_angle) \nsin $(data.sin_trial_torsion_angle)  \n test:$(data.sin_trial_torsion_angle.^2 .+ data.cos_trial_torsion_angle.^2)")

            println("pos ", data.trial_positions)
            break
        end
        data.tmp1 .= data.xyz[data.id].-data.xyz[data.id-1]
        ×( data.xyz[data.id-1]-data.xyz[data.id-2], data.tmp1,data.crossproduct)
        data.current .= data.tmp1
        data.id+=1

    end
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

