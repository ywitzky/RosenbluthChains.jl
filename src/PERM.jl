using ProgressMeter


export Noperm, PermData


### oriented on THE paper of Peter Grassberger
mutable struct PermData{T<:Real,I<:Integer} <: AbstractPERMData
    K::I
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
    μ::T
    σ_sqr::T

    Measure::Bool
    PreAverage::I
    PermData(NBeads;K=1, TType=Float64, IType=Int32, N=1000, cMax=3.0,cMin=0.3, PreAverage=1000) = new{TType,IType}(K, zeros(IType, N), zeros(TType, N), N, zero(IType) , zero(IType), zero(IType), zero(IType), TType(Inf).*ones(TType, NBeads),zeros(TType, NBeads),cMax, cMin,zeros(TType, NBeads),zeros(TType, NBeads), zero(TType), zero(TType), true, PreAverage)
end



function mainLoop( data::SimData, param::SimulationParameters, Measurement::AbstractMeasurement, perm::PermData)

    for data.batch_id in 0:data.NumberOfBatches-1
        #@showprogress 1 "Computing batch $(data.batch_id+1)" 
        data.id_in_batch = 1
        while data.id_in_batch <= data.BatchSize
            if perm.N_alternatives == 0 ### no alternative replicas left with same ancestor
                UpdateAcceptanceWeights(data, perm)### do not update within tour

                ResetSim(data, param, perm)
                SetFirstThreeBeads(data, param)
                data.RosenbluthWeight
            else
                pickNextReplica(data,param, perm)
            end
            ComputeBeadsIteratively(data,param, perm)
            if perm.Measure
                MeasureAfterChainGrowth(data, param,Measurement)

                #=
                if data.id_in_batch > perm.PreAverage-10
                    println(data.RosenbluthWeight)
                end=#
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
    data.id = deepcopy(perm.alternatives_ids[perm.N_alternatives])
    data.RosenbluthWeight = deepcopy(perm.alternatives_weights[perm.N_alternatives])

    data.tmp1 .= data.xyz[data.id].-data.xyz[data.id-1]
    ×( data.xyz[data.id-1].-data.xyz[data.id-2],data.tmp1,data.crossproduct) ### set crossproduct to old values
    data.current .= data.tmp1 ### set bond direction to old values


    if isinf(data.RosenbluthWeight)
        println("data.id $(data.id), batch_id $(data.id_in_batch), perm.N_alternatives $(perm.N_alternatives) \n\n\n")
        data.id_in_batch = data.BatchSize ### just break
    end

    perm.N_alternatives -= 1
    perm.Measure = true
    nothing
end

@inline function prune(data::SimData, param::SimulationParameters, perm::PermData)
    #println("prune")

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
    #println("enrich")
    ### allocate more memory if its too little
    if perm.N_alternatives+perm.K>perm.N_alloc
        tmp = deepcopy(perm.alternatives_weights)
        perm.alternatives_weights = zeros(T, 2*perm.N_alloc)
        perm.alternatives_weights[1:perm.N_alloc] .= tmp

        tmp = deepcopy(perm.alternatives_ids)
        perm.alternatives_ids = zeros(I, 2*perm.N_alloc)
        perm.alternatives_ids[1:perm.N_alloc] .= tmp
        perm.N_alloc *= 2
    end

    ### half the weight, since we create a second copy
    #println("pre $( data.RosenbluthWeight)")
    data.RosenbluthWeight = data.RosenbluthWeight / T(perm.K+1) 
    #println("post $( data.RosenbluthWeight)")

    for i in 1:perm.K
        perm.N_alternatives += 1

        perm.alternatives_ids[perm.N_alternatives] = deepcopy(data.id)
        perm.alternatives_weights[perm.N_alternatives] = deepcopy(data.RosenbluthWeight)
    end

    perm.TotalEnrich += perm.K

    nothing
end

@inline function ChooseTrialPosition(data::SimData,param::SimulationParameters, perm::NoPERM)
    ChooseTrialPosition(data::SimData,param::SimulationParameters)
    nothing
end

@inline function ChooseTrialPosition(data::SimData,param::SimulationParameters, perm::PermData)

    ChooseTrialPosition(data, param) ### normaly choose position out of trial positions

    rel_weight  = data.RosenbluthWeight/data.id### weight per Residue
    if data.RosenbluthWeight< perm.MinWeight[data.id]
        if rand() < 1/perm.K
            prune(data, param, perm)
        else
            data.RosenbluthWeight *= perm.K
        end
    else
        if data.RosenbluthWeight > perm.MaxWeight[data.id]
            enrich(data, param, perm)
        end
    end
    perm.Z[data.id] += data.RosenbluthWeight
    perm.N_meas[data.id] += 1

    nothing
end


@inline function UpdateAcceptanceWeights(data::SimData{T,I}, perm::PermData) where {T<:Real,I<:Integer}
    if perm.Measure 
        NMeasured = data.batch_id * data.BatchSize + data.id_in_batch -1 ### the main loop updated too early for this use case
        #perm.Z += data.RosenbluthWeight
        if NMeasured< perm.PreAverage && NMeasured>1
            val = log10(data.RosenbluthWeight)
            old = deepcopy(perm.μ)
            perm.μ = perm.μ + (val - perm.μ)/NMeasured
            perm.σ_sqr +=val^2

            if data.RosenbluthWeight==0
                println("Bla: NMeasured: $(NMeasured)")
            end
            #perm.Z_Min = min(perm.Z_Min, data.RosenbluthWeight)
        end
        if NMeasured == perm.PreAverage
            ### assume that we have a log gaussian distribution

            σ = √(perm.σ_sqr/perm.PreAverage - perm.μ^2) 
            println("μ $(perm.μ) ± $(σ)")
            perm.cMaxWeight = 10^(2.0*σ  )#10^(perm.μ+5*σ)
            perm.cMinWeight = 10^(-1.5*σ ) #perm.μ-0.5*σ)
  
            range= data.NBeads./collect(1:data.NBeads)

            tmp = perm.Z./perm.N_meas
            #tmp ./= tmp[end]
            perm.MaxWeight = @. tmp.* 10^(2*σ/range) #perm.cMaxWeight
            perm.MinWeight = @. tmp.*  10^(-2*σ/range) # perm.cMinWeight./range#/perm.Z[1]


            println("min/max = $(perm.cMinWeight) $(perm.cMaxWeight) , ratio: $(perm.cMaxWeight/perm.cMinWeight)")
        end

        if NMeasured > perm.PreAverage ### start updating after one chain is grown
            #perm.Z += data.RosenbluthWeight# (data.RosenbluthWeight+ T(NMeasured-1)*perm.Z)/T(NMeasured) ### stream computation of mean of rosenbluth weights

            #=
            range= data.NBeads./collect(1:data.NBeads)
            tmp = perm.Z./perm.N_meas
            tmp ./= tmp[end]
            perm.MaxWeight = @. tmp.* perm.cMaxWeight./range#/perm.Z[1]
            perm.MinWeight = @. tmp.* perm.cMinWeight./range#/perm.Z
            [1]
            =#

            #perm.MaxWeight = @. perm.Z*perm.cMaxWeight/perm.N_meas
            #perm.MinWeight = @. perm.Z*perm.cMinWeight/perm.N_meas
        end
    end
    nothing
end