@inline function SetTrialRadius(data::SimData, param::FixedBondParameters) 
    ### is set once in set up and doesnt change afterwards
end

@inline function SetTrialBondAngle(data::SimData,param::FixedBondAngles)
    fill!(data.trial_angle, param.BondAngles[data.id-2])
end

@inline function SetTrialBondAngle(data::SimData,param::GaussLpBondAngles)
    data.trial_angle.= lpToAngle.(rand(param.PersistenceDist, data.NTrials))
end




@inline function InitSimParam(data::SimData,param::FixedBondParameters ) 
    fill!(data.trial_radius, param.BondLength)
end

@inline function InitSimParam(data::SimData,param::GaussLpBondAngles ) end

@inline function InitSimParam(data::SimData,param::FixedBondAngles ) end

@inline function InitSimParam(data::SimData,param::SAWParameters ) end

@inline function InitSimParam(data::SimData,param::IdealChain ) end