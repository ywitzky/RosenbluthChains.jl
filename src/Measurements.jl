function InitLammpsAnalysisScript(data::SimData, param::SimulationParameters)
    LA_Obj = LA.initReducedData_NoConfStorage(data.FolderPath, data.BatchSize, 1, [1], [data.NBeads])
    LA_Obj
end

function CopyXYZToLA(data::SimData, param::SimulationParameters, LA_Obj,i)
    ### change from Array of Sturct to Struct of Array
    LA_Obj.x[:, i] .= getindex.(data.xyz, 1)
    LA_Obj.y[:, i] .= getindex.(data.xyz, 2)
    LA_Obj.z[:, i] .= getindex.(data.xyz, 3)
end

function DoBatchMeasurements(data::SimData, param::SimulationParameters, Measurements::Dict{String, Array{T}}, LA_Obj::LA.SimData{T2, N2}, I::N) where {T<:Real, T2<:Real,N<:Integer,N2<:Integer}
    MeasureRange= I*data.BatchSize+1:(I+1)*data.BatchSize
    LA.computeCOM!(LA_Obj)
    LA.computeInertiaTensor(LA_Obj)
    Measurements["λ"][:,:,MeasureRange] .= LA_Obj.InertiaTensorEigVals[:, :,:]

    LA.computeBondAngles(LA_Obj)
    Measurements["l_p"][:,I+1] .= LA_Obj.LocalPersistenceLength[:]
    #Measurements["BondAngle"][:,I+1] .= LA_Obj.BondAngles[:]
    
end

function InitMeasurementStorage(data::SimData, param::SimulationParameters)
    TotalStep = data.NumberOfBatches*data.BatchSize
    Measurements = Dict{String, Array{eltype(data.trial_radius)}}()
    
    Measurements["λ"]   = zeros(eltype(data.trial_radius), (3,1,TotalStep) )
    Measurements["l_p"] = zeros(eltype(data.trial_radius), (data.NBeads-2, data.NumberOfBatches))
    #Measurements["BondAngle"] = zeros(eltype(data.trial_radius), (data.NBeads-2, data.NumberOfBatches))
    return Measurements
end

function MeasurementsToPickle(data::SimData,Measurements::Dict{String, Array{T}}) where {T<:Real}
    lp = Measurements["l_p"] 
   # Measurements["l_p"] = reshape(sum(lp, dims=2)/data.NumberOfBatches, data.NBeads-2)
   # Measurements["Δl_p"] = reshape(sum((lp.-Measurements["l_p"]).^2, dims=2)/(data.NumberOfBatches-1), data.NBeads-2)

    JLD2.save(data.FolderPath*"Data/Measurements.jld2", Measurements)
end