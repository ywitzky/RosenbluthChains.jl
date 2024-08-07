using ProgressMeter


export Noperm, PermData


### oriented on THE paper of Peter Grassberger
mutable struct PermData{T<:Real,I<:Integer} <: AbstractPERMData
    alternatives_ids::Vector{I}
    alternatives_weights::Vector{T}
    N_alloc::I
    N_alternatives::I ### tracks number of replicas left to generate from same ancestor
    Total_Grown::I ### number of different "tours" that have been started
    TotalPrune::I
    TotalEnrich::I
    MaxWeight::Vector{T}
    MinWeight::Vector{T}
    cMaxWeight::T
    cMinWeight::T
    Z::Vector{T} ### current estimate of partition sum
    N_meas::Vector{T} ### current estimate of partition sum

    Measure::Bool
    PreAverage::I
    PermData(NBeads;TType=Float64, IType=Int32, N=1000, cMax=3.0,cMin=0.3, PreAverage=1000) = new{TType,IType}(zeros(IType, N), zeros(TType, N), N, zero(IType) , zero(IType), zero(IType), zero(IType), TType(Inf).*ones(TType, NBeads),zeros(TType, NBeads),cMax, cMin,zeros(TType, NBeads),zeros(TType, NBeads), true, PreAverage)
end



function mainLoop( data::SimData, param::SimulationParameters, Measurement::AbstractMeasurement, perm::PermData)

    for data.batch_id in 0:data.NumberOfBatches-1
        #@showprogress 1 "Computing batch $(data.batch_id+1)" 
        data.id_in_batch = 1
        while data.id_in_batch <= data.BatchSize
            if perm.N_alternatives == 0 ### no alternative replicas left with same ancestor
                ResetSim(data, param, perm)
                SetFirstThreeBeads(data, param)
            else
                pickNextReplica(data,param, perm)
            end
            ComputeBeadsIteratively(data,param, perm)
            if perm.Measure
                MeasureAfterChainGrowth(data, param,Measurement)
                UpdateAcceptanceWeights(data, perm)

                data.id_in_batch += 1

            end
        end
        MeasureAfterBatch(data, param, Measurement)
    end
end

function ResetSim(data::SimData, param::SimulationParameters, perm::PermData)
    ResetSim(data, param) ### do normal clean up 
    if perm.Measure ### track all the tries where at least one replica was succefully grown
        perm.Total_Grown += 1 
    end
    perm.Measure=true
end

@inline function pickNextReplica(data::SimData, param::SimulationParameters, perm::PermData)
    ### take over the weights and ids from the last enrichment step
    data.id = perm.alternatives_ids[perm.N_alternatives]
    data.RosenbluthWeight = perm.alternatives_weights[perm.N_alternatives]
    #println("pick replica weight $( data.RosenbluthWeight ) vs avg. $(perm.Z/(data.batch_id * data.BatchSize + data.id_in_batch)), max $(perm.MaxWeight) min $(perm.MinWeight)")

    data.tmp1 .= data.xyz[data.id].-data.xyz[data.id-1]
    ×( data.xyz[data.id-1].-data.xyz[data.id-2],data.tmp1,data.crossproduct) ### set crossproduct to old values
    data.current .= data.tmp1 ### set bond direction to old values

    perm.N_alternatives -= 1
    perm.Measure = true
    nothing
end

@inline function prune(data::SimData, param::SimulationParameters, perm::PermData)
    perm.TotalPrune += 1
    if perm.N_alternatives > 0
        pickNextReplica(data, param, perm)
    else
        data.id = data.NBeads ### breaks loop in ComputeBeadsIteratively
        perm.Measure = false ### avoid structures that are half created through pruning
    end
    nothing
end

@inline function enrich(data::SimData{T,I}, param::SimulationParameters, perm::PermData) where {T<:Real,I<:Integer}
    ### allocate more memory if its too little
    if perm.N_alternatives>perm.N_alloc
        tmp = deepcopy(data.perm.alternatives_weights)
        perm.alternatives_weights = zeros(T, 2*perm.N_alloc)
        perm.alternatives_weights[1:perm.N_alloc] .= tmp

        tmp = deepcopy(data.perm.alternatives_ids)
        perm.alternatives_ids = zeros(I, 2*perm.N_alloc)
        perm.alternatives_ids[1:perm.N_alloc] .= tmp
    end

    ### half the weight, since we create a second copy
    data.RosenbluthWeight /= 2.0

    perm.N_alternatives +=1
    perm.alternatives_ids[perm.N_alternatives] = data.id
    perm.alternatives_weights[perm.N_alternatives] = data.RosenbluthWeight

    perm.TotalEnrich += 1

    nothing
end

@inline function ChooseTrialPosition(data::SimData,param::SimulationParameters, perm::NoPERM)
    ChooseTrialPosition(data::SimData,param::SimulationParameters)
    nothing
end

@inline function ChooseTrialPosition(data::SimData,param::SimulationParameters, perm::PermData)
    ChooseTrialPosition(data, param) ### normaly choose position out of trial positions

    #println("weight $(data.RosenbluthWeight ) $(data.id)")

    rel_weight  = data.RosenbluthWeight/data.id### weight per Residue
    if rel_weight< perm.MinWeight[data.id]
        data.RosenbluthWeight *= 2.0
        if rand() < 0.5
            #println("prune")

            prune(data, param, perm)
        end
    else
        if rel_weight > perm.MaxWeight[data.id]
            enrich(data, param, perm)
        end
    end
    perm.Z[data.id] += data.RosenbluthWeight
    perm.N_meas[data.id] += 1

    nothing
end


@inline function UpdateAcceptanceWeights(data::SimData{T,I}, perm::PermData) where {T<:Real,I<:Integer}
    if perm.Measure 
        NMeasured = data.batch_id * data.BatchSize + data.id_in_batch #-1 ### the main loop updated too early for this use case
        #perm.Z += data.RosenbluthWeight
        if NMeasured >= perm.PreAverage ### start updating after one chain is grown
            #perm.Z += data.RosenbluthWeight# (data.RosenbluthWeight+ T(NMeasured-1)*perm.Z)/T(NMeasured) ### stream computation of mean of rosenbluth weights

            perm.MaxWeight = @. perm.Z*perm.cMaxWeight/perm.N_meas
            perm.MinWeight = @. perm.Z*perm.cMinWeight/perm.N_meas
            #=
            println("pre: $(data.RosenbluthWeight)")
            println("Z: $(perm.Z), cmin: $(perm.cMinWeight), cmax: $(perm.cMaxWeight)")
            println("min $(perm.MinWeight), max: $(perm.MaxWeight)")
            println()
            =#
        end
    end
    nothing
end