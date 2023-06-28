@inline function SetTrialRadius(data::SimData, param::FixedBondParameters) 
    ### is set once in set up and doesnt change afterwards
end

@inline function SetTrialBondAngle(data::SimData,param::FixedBondAngles)
    fill!(data.trial_angle, param.BondAngles[data.id-2])
end

@inline function SetTrialBondAngle(data::SimData,param::GaussLpBondAngles)
    data.trial_angle.= lpToAngle.(rand(param.PersistenceDist, data.NTrials))
end
