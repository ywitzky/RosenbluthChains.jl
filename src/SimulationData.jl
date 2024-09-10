export SimData

mutable struct SimData{T<:Real, I<:Integer}
    FolderPath::String
    NBeads::I
    NTrials::I
    BatchSize::I
    NumberOfBatches::I
    TType::T
    IType::I
    id::I ### current bead 
    tid::I ### trial id
    id_in_batch::I
    batch_id::I
    xyz::Vector{Vector3{T}}
    max_interaction_length::T ### used for getPotentialNeighbors to hash non bonded interactions

    ### Rosebluth Weight stuff
    LogRosenbluthWeight::T   ### Log of weight for the whole chain at current step;
    RosenbluthWeight::T   ### weight for the whole chain only updated at the end;
    BoltzmannFaktor::Vector{T}
    LogBoltzmannFaktor::Vector{T}

    ### Data for Trial positions
    trial_radius::Vector{T} ### radius for the trial positions
    trial_angle::Vector{T} ### trial angle for bonds
    sin_trial_angle::Vector{T} ### sine of trial angle for bonds
    cos_trial_angle::Vector{T} ### cosine of trial angle for bonds
    trial_torsion_angle::Vector{T} ### trial torsion angles for bonds
    sin_trial_torsion_angle::Vector{T} ### sine of trial torsion angles for bonds
    cos_trial_torsion_angle::Vector{T} ### cosine of trial torsion angles for bonds
    trial_positions::Vector{Vector3{T}} ### trial positions for Rosenbluth
    
    
    ### Allocations for computations, that will be reused
    RotationMatrix::Matrix3{T} ### RotationMatrix
    crossproduct::Vector3{T}
    current::Vector3{T}
    new_vec::Vector3{T}
    tmp1::Vector3{T}
    tmp2::Vector3{T}
    tmp3::Vector3{T}
    tmp4::Vector{T}
    tmp5::Vector{T}
    btmp::Vector{T}
    btmp2::Vector{T}
    btmp3::Vector{T}
    id_arr::Vector{I}
    id_arr2::Vector{I}
    id_arr3::Vector{I}

    val_arr::Vector{T}
    x_arr::Vector{T}
    y_arr::Vector{T}
    z_arr::Vector{T}


    ### Constructor
    SimData(FolderPath::String, type::T, NBeads::I, NTrials::I, BatchSize::I, NumberOfBatches::I) where {T<: Real, I<: Integer} =  
    new{T,I}(FolderPath, NBeads, NTrials, BatchSize, NumberOfBatches, 1,1, 1,-1,0,0, Vector{Vector3{T}}(NBeads) ,0.0,
    0.0,0.0,zeros(NTrials), zeros(NTrials),
    zeros(T, NTrials), zeros(T, NTrials),  zeros(T, NTrials),zeros(T, NTrials),zeros(T, NTrials),zeros(T, NTrials),zeros(T, NTrials), Vector{Vector3{T}}(NTrials), 
    Matrix3{T}(), Vector3{T}(0,0,0),Vector3{T}(0,0,0), Vector3{T}(0,0,0) , Vector3{T}(0,0,0), Vector3{T}(0,0,0), Vector3{T}(0,0,0) , zeros(T, NTrials), zeros(T, NTrials), zeros(T, NTrials),  zeros(T, NTrials),  zeros(T, NTrials) ,  zeros(I, NBeads) ,  zeros(I, NBeads) ,  zeros(I, NBeads) ,  zeros(T, NBeads),  zeros(T, NTrials),  zeros(T, NTrials),  zeros(T, NTrials)) #MVector{NTrials, T}(zeros(T, NTrials)), MVector{NTrials, T}(zeros(T, NTrials)), MVector{NTrials, T}(zeros(T, NTrials)))

end