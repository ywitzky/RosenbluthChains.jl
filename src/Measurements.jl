#import Base.close

#### measurement for radius of gyration and end-to-end distance length
mutable struct RG_Measurement{T<:Real} <: AbstractMeasurement 
    Path::String
    RGs::Vector{T}
    REEs::Vector{T}
    Weights::Vector{T}

    io_RG::IOStream
    io_REE::IOStream
    io_weights::IOStream
    com::Vector3{T}
    tmp::Vector3{T}
    cnt::Int64
    RG_Measurement(FolderPath::String, NProteins::Int; T=Float64) = begin 
        StoragePath = "$(FolderPath)tmp/"
        mkpath(StoragePath)

        io_rg  = open("$(StoragePath)rg.bin", "w+")
        io_ree     = open("$(StoragePath)ree.bin", "w+")
        io_weights = open("$(StoragePath)weights.bin", "w+")
    
        rg      = Mmap.mmap(io_rg      , Vector{T}, (NProteins))
        ree     = Mmap.mmap(io_ree     , Vector{T}, (NProteins))
        weights = Mmap.mmap(io_weights , Vector{T}, (NProteins))

        new{T}(FolderPath, rg, ree, weights, io_rg, io_ree, io_weights, zeros(T, 3),zeros(T, 3), 1) 
    end 
end


### initialise measurments with empty vectors
RosenbluthChains.InitMeasurement(data::SimData, param::SimulationParameters, Tmp::RG_Measurement) =  RG_Measurement(data.FolderPath, data.BatchSize*data.NumberOfBatches)

function RosenbluthChains.MeasureAfterChainGrowth(data::SimData, param::SimulationParameters,Measurement::RG_Measurement) 
    for pos in data.xyz
        Measurement.com .+= pos
    end
    Measurement.com ./= data.NBeads

    for pos in data.xyz
        Measurement.tmp .= pos-Measurement.com
        Measurement.RGs[Measurement.cnt] += sqr_norm(Measurement.tmp)
    end
    Measurement.RGs[Measurement.cnt] /= data.NBeads

    Measurement.REEs[Measurement.cnt] = sqr_norm(data.xyz[end]-data.xyz[1])
    Measurement.Weights[Measurement.cnt] =  getRosenbluthWeigth(data, param)

    Measurement.cnt+= 1
end

RosenbluthChains.MeasureAfterBatch(data::SimData, param::SimulationParameters, Measurement::RG_Measurement) = nothing


function close(Measurement::RG_Measurement{T}) where {T<:Real}
    StoragePath = "$(Measurement.Path)tmp/"

    Base.close(Measurement.io_RG)
    Base.close(Measurement.io_REE)
    Base.close(Measurement.io_weights)

    rm("$StoragePath/rg.bin")
    rm("$StoragePath/ree.bin")
    rm("$StoragePath/weights.bin")
end



#### measurement for radius of gyration and end-to-end distance length
mutable struct RG_Mass_Measurement{T<:Real} <: AbstractMeasurement 
    Path::String
    RGs::Vector{T}
    REEs::Vector{T}
    Weights::Vector{T}
    Masses::Vector{T}
    Total_Mass::T

    io_RG::IOStream
    io_REE::IOStream
    io_weights::IOStream
    com::Vector3{T}
    tmp::Vector3{T}
    cnt::Int64
    RG_Mass_Measurement(FolderPath::String, NProteins::I, Masses::Vector{T}) where {T<:Real, I<:Integer}= begin 
        StoragePath = "$(FolderPath)tmp/"
        mkpath(StoragePath)

        io_rg  = open("$(StoragePath)rg.bin", "w+")
        io_ree     = open("$(StoragePath)ree.bin", "w+")
        io_weights = open("$(StoragePath)weights.bin", "w+")
    
        rg      = Mmap.mmap(io_rg      , Vector{T}, (NProteins))
        ree     = Mmap.mmap(io_ree     , Vector{T}, (NProteins))
        weights = Mmap.mmap(io_weights , Vector{T}, (NProteins))

        new{T}(FolderPath, rg, ree, weights, Masses, sum(Masses), io_rg, io_ree, io_weights, zeros(T, 3),zeros(T, 3), 1) 
    end 
end


### initialise measurments with empty vectors
RosenbluthChains.InitMeasurement(data::SimData, param::SimulationParameters, Tmp::RG_Mass_Measurement) =  RG_Mass_Measurement(data.FolderPath, data.BatchSize*data.NumberOfBatches, Tmp.Masses)

function RosenbluthChains.MeasureAfterChainGrowth(data::SimData, param::SimulationParameters,Measurement::RG_Mass_Measurement) 
    for (i,pos) in enumerate(data.xyz)
        Measurement.com .+= pos*Measurement.Masses[i]
    end
    Measurement.com ./= Measurement.Total_Mass

    for (i,pos) in enumerate(data.xyz)
        Measurement.tmp .= pos-Measurement.com
        Measurement.RGs[Measurement.cnt] += sqr_norm(Measurement.tmp) *Measurement.Masses[i]
    end
    Measurement.RGs[Measurement.cnt] /= Measurement.Total_Mass

    Measurement.REEs[Measurement.cnt] = sqr_norm(data.xyz[end]-data.xyz[1])
    Measurement.Weights[Measurement.cnt] =  getRosenbluthWeigth(data, param)

    Measurement.cnt+= 1
end

RosenbluthChains.MeasureAfterBatch(data::SimData, param::SimulationParameters, Measurement::RG_Mass_Measurement) = nothing


function close(Measurement::RG_Mass_Measurement{T}) where {T<:Real}
    StoragePath = "$(Measurement.Path)tmp/"

    Base.close(Measurement.io_RG)
    Base.close(Measurement.io_REE)
    Base.close(Measurement.io_weights)

    rm("$StoragePath/rg.bin")
    rm("$StoragePath/ree.bin")
    rm("$StoragePath/weights.bin")
end

