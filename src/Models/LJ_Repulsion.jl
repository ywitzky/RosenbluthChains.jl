export LJ_Repulsion, clear
using Printf

using LoopVectorization

import Base.@kwdef, Base.show, Base.iterate

#=
CELL LIST seems to have too much overhead for short polymers; Long polymers dont work well with rosenbluth-rosenbluth method anyway.

export LinkedList, add!

abstract type AbstractSortedLinkedList end

@kwdef mutable struct Link{I <: Integer, T<:Any}
    prev::Union{Nothing, Link{I,T}, AbstractSortedLinkedList, Link{I,Nothing}}
    next::Union{Nothing, Link{I,T}, AbstractSortedLinkedList, Link{I,Nothing}}
    key::I
    value::T
end

mutable struct LinkedList{I <: Union{Integer, Nothing, AbstractSortedLinkedList}, T<:Any}  <: AbstractSortedLinkedList
    start::Union{Link{I,T}, Nothing, Link{I,Nothing}}
    stop::Union{Link{I,T}, Nothing, Link{I,Nothing}}
    current::Union{Link{I,T}, Nothing, Link{I,Nothing}}
end

function LinkedList(key::I, X::Any) where {I<:Integer}  
    NewLink = Link{I,typeof(X)}(nothing, nothing, 0, X)
    LinkedList{I, typeof(X)}(NewLink, NewLink, NewLink)
end

function LinkedList(_::I) where {I<:Integer}
    LinkedList{I, I}(nothing, nothing, nothing)
end

function add!(list::LinkedList, value::Any)
    NewLink = Link{Int64,typeof(value)}(list.stop, nothing, 0, value)
    if isnothing(list.start) ### checks for the initial state
        setfield!(list, :start, NewLink)
        setfield!(list, :stop , NewLink)
    else
        setfield!(list.stop,:next, NewLink)
        setfield!(list, :stop, NewLink)
    end
    nothing
end

function add!(list::LinkedList{I, LinkedList{I,I}}, value::Any) where {I<:Integer}
    NewLink = Link{Int64,typeof(value)}(list.stop, nothing, 0, value)
    if isnothing(list.start.value.start) ### checks for the initial state
        setfield!(list, :start, NewLink)
        setfield!(list, :stop , NewLink)
    else
        setfield!(list.stop,:next, NewLink)
        setfield!(list, :stop, NewLink)
    end
    nothing
end

function append!(A::LinkedList, B::LinkedList) ### outdated technically probably not needed
    NewList = deepcopy(B)
    setfield!(A.stop,:next, NewList.start)
    setfield!(A, :stop, NewList.stop)
    nothing
end

function insert!(A::LinkedList{I, LinkedList{I,T}}, B::LinkedList) where {I<: Integer , T<:Any}
    add!(A, B)
    nothing
end

function Base.iterate(A::LinkedList{I,T}) where {I<: Integer , T<:Any}
    A.current= A.start
    if isnothing(A.current)
        return nothing
    else
        val = A.current.value
        setfield!(A, :current, A.start.next)
        return (val, A )
    end
end

function Base.iterate(_::LinkedList, A::Union{LinkedList, Nothing})
    if isnothing(A.current)
        return nothing
    else
        val = A.current.value
        setfield!(A, :current, A.start.next)
        return (val, A )
    end
end

function Base.iterate(A::LinkedList{I, LinkedList{I,T}}) where {I<: Integer , T<:Any}
    A.current= A.start
    if isnothing(A.current) || isnothing(A.current.value.start)
        return nothing
    else
        A.current.value.current = A.current.value.start ### start searching the sublist at the start
        val = A.current.value.current.value
        A.current.value.current = A.current.value.current.next
        return (val, val) #Base.iterate(A.current.value)
    end
end

function Base.iterate(A::LinkedList{I, LinkedList{I,T}}, _::T) where {I<: Integer , T<:Any}
    if isnothing(A.current) ### check if outer list is done
        return nothing
    end
    if isnothing(A.current.value.current) ### check if inner list is empty and choose the next one
        A.current = A.current.next
        if isnothing(A.current) 
            return nothing
        else
            A.current.value.current = A.current.value.start ### start at the start of the inner list
        end
    end
    val = A.current.value.current.value
    setfield!(A.current.value, :current, A.current.value.current.next)
    return (val, val )
end


mutable struct SortedList{I <: Union{Integer, Nothing, AbstractSortedLinkedList}, T<:Any}  <: AbstractSortedLinkedList
    start::Link{I,T}
    stop::Link{I,T}
    last::Link{I,T}
end

function SortedList(key::I, X::Any) where {I<:Integer}  
    NewLink = Link{I,typeof(X)}(nothing, nothing, key, X)
    SortedList{I, typeof(X)}(NewLink, NewLink, NewLink)
end

function SortedList3D(keys::Vector{I},value::Any) where {I<:Integer}
    SortedList(keys[1], SortedList(keys[2], SortedList(keys[3], LinkedList(0, value))))
end

function GenerationFunc(value::Any, nextInd::Vector{I}) where {I<:Integer}
    len = length(nextInd)
    if len==0
        return LinkedList(0,value)
    elseif len==1
        return SortedList(nextInd[1], LinkedList(0,value))
    elseif len==2
        return SortedList(nextInd[1], SortedList(nextInd[2], LinkedList(0,value)))
    end
end


contains(list::AbstractSortedLinkedList, ind::Integer) = list.stop.key<ind || list.start.key>ind ? false : true

function add!(list::AbstractSortedLinkedList, inds::Union{Vector{I}, Vector3{I}}, value::Any) where {I<:Integer} add!(list, inds, value, GenerationFunc) end

function add!(list::AbstractSortedLinkedList, inds::Union{Vector{I}, Vector3{I}}, value::Any, GenFunc::Function)  where {I<:Integer}
    if(inds[1]> list.stop.key)
        NewLink = Link(list.stop, nothing, inds[1], GenFunc(value, inds[2:end]))
        #list.stop.next = NewLink
        #list.stop = NewLink
        setfield!(list.stop, :next, NewLink)
        setfield!(list, :stop, NewLink)
    elseif ( inds[1]<list.start.key)
        NewLink = Link(nothing, list.start, inds[1], GenFunc(value,  inds[2:end]))
        #list.start.prev = NewLink
        #list.start = NewLink
        setfield!(list.start, :prev, NewLink)
        setfield!(list, :start, NewLink)
    else
        current =  list.start
        while(current.key < inds[1])
            current = current.next
        end
        if inds[1] == current.key
            if length(inds)>= 2
                add!(current.value, inds[2:end], value, GenFunc)
            else ### assuming that the lowest list is of type LinkedList
                add!(current.value, value)
            end
        else
            if length(inds)>= 2
                NewLink = Link{I,typeof(current.value)}(current.prev, current,inds[1], GenFunc(value,  inds[2:end]))
                setfield!(current.prev, :next, NewLink)
                setfield!( current, :prev, NewLink)
            elseif length(inds)== 1
                NewLink = Link{I,typeof(current.value)}(current.prev, current,inds[1], LinkedList(0,value))
                setfield!(current.prev, :next, NewLink)
                setfield!( current, :prev, NewLink)
            else ### assuming that the lowest list is of type LinkedList
                add!(current.prev.value, value)
            end

        end
    end
end



function Base.show(io::IO,  l::Link)  @printf(io, "%s, ",l.value); if !isnothing(l.next) show(io,  l.next) end end

function Base.show(io::IO, X::MIME"text/plain", l::Link)  @printf(io, "%s, ",l.value); if !isnothing(l.next) show(io, X , l.next) end end

function Base.show(io::IO,  l::LinkedList) 
    @printf(io, "[" ) 
    if !isnothing(l.start)
        show(io, l.start)
    end
    @printf(io, "]" ) 
end

function Base.show(io::IO,  l::SortedList) 
    @printf(io, "[" ) 
    if !isnothing(l.start)
        show(io, l.start)
    end
    @printf(io, "]" ) 
end

#=
### version with deepcopy
function neighbours(list::AbstractSortedLinkedList, inds::Union{Vector{I}, Vector3{I}},resList::LinkedList) where {I <: Integer}
    if typeof(list)===LinkedList{I,I}
        resList = append!(resList, list)
        return nothing
    end
    current =  list.start
    while(current!= nothing && current.key < inds[1]-1 && current.key < inds[1]+2 )
        current = current.next
    end
    if current== nothing || current.key >= inds[1]+2 ### no solution is in sortedlinkedlist
        return nothing
    else
        neighbours(current.value, inds[2:end], resList)
        current = current.next
        if current!= nothing && current.key < inds[1]+2
            neighbours(current.value, inds[2:end], resList)
            current = current.next
            if current!= nothing && current.key < inds[1]+2
                neighbours(current.value, inds[2:end], resList)
            end
        end
    end
end
=#

function neighbours(list::AbstractSortedLinkedList, inds::Union{Vector{I}, Vector3{I}},resList::LinkedList) where {I <: Integer}
    if typeof(list)===LinkedList{I,I} ### check if final list is reached -> append list
        resList = add!(resList, list)
        return nothing
    end
    current = list.last ### search from the last hit we had since, the searches will be correlated to the last hit
    if current.key>=inds[1]-1
        while(~isnothing(current)  &&  current.key >= inds[1]-1  && ~isnothing(current.prev) )
            current = current.prev
        end
        if ~isnothing(current) && ~isnothing(current.next)  && ~isnothing(current.prev)
            current = current.next
        end
    elseif current.key< inds[1]-1
        while(~isnothing(current)  && current.key < inds[1]-1 && current.key < inds[1]+2  && ~isnothing(current.next))
            current = current.next
        end
    end
    if isnothing(current) || current.key >= inds[1]+2 ### no solution is in sortedlinkedlist
        return nothing
    else
        neighbours(current.value, inds[2:end], resList)
        list.last= current
        current = current.next
        if ~isnothing(current) && current.key < inds[1]+2
            neighbours(current.value, inds[2:end], resList)
            current = current.next
            if ~isnothing(current) && current.key < inds[1]+2
                neighbours(current.value, inds[2:end], resList)
            end
        end
    end
end

function neighbours(list::AbstractSortedLinkedList, inds::Union{Vector{I}, Vector3{I}}) where {I <: Integer}
    resList = LinkedList(0,LinkedList(zero(I)))
    neighbours(list, inds, resList)
    resList
end


###############################
###############################
###############################



### naiv implementation, probably possibilities to speed up
mutable struct CellList{I <: Integer}
    dict::Dict{Vector3{I}, Vector{I}}
    HashFunc::Function
    offsets::Vector{Vector3{I}}
    neighbourcells::SVector{27,Vector3{I}}
    CellList(GridSpacing::T ) where{T<:Real} = begin
        #d = collect(floor(-c/s):1:ceil(c/s))
        #bla = [[x,y,z] for x in d for y in d for z in d if x^2+y^2+z^2 < ceil((c/s+0.001)^2)] #### c/s is slightly smaller because of numerics

        new{Int64}(Dict{Vector3{Int64}, Vector{Int64}}(),     x-> Get3DHash(x, GridSpacing), [ [x,y,z] for x in [-1,0,1] for y in [-1,0,1] for z in [-1,0,1]], [[0,0,0] for i in 1:27])
    end
end



function add(CL::CellList{I}, pos::Vector3{T}, id::I) where{I<:Integer, T<:Real}
    ind = CL.HashFunc(pos)
    if ind in keys(CL.dict)
        append!(CL.dict[ind], id)
    else
        CL.dict[ind]=[id]
    end
end

#=
function neighbours(CL::CellList{I}, pos::Vector3{T}) where{I<:Integer, T<:Real} 
    ind = CL.HashFunc(pos) 
    for id in 1:27
        CL.neighbourcells[id] .= ind + CL.offsets[id]
    end
    key = copy(pos)
    return [getkey(CL.dict, i, empty) end for i in CL.neighbourcells if i in keys(CL.dict)]
    #return [getkey(CL.dict, i, empty) for i in CL.neighbourcells if i in keys(CL.dict)]
end

function neighbours(CL::CellList{I}, pos::Vector3{T}) where{I<:Integer, T<:Real} 
    ind = CL.HashFunc(pos) 
    for id in 1:27
        for i in 1:3
            CL.neighbourcells[id][i] = ind[i] + CL.offsets[id][i]
        end
    end
   # key = copy(pos)
   # return [getkey(CL.dict, [ind[1]+], empty) end for i in CL.neighbourcells if i in keys(CL.dict)]
    return [getkey(CL.dict, i, empty) for i in CL.neighbourcells if i in keys(CL.dict)]
end
=#


mutable struct CellLinkedList{I <: Integer}
    HashFunc::Function
    data::SortedList{Int64, SortedList{Int64, SortedList{Int64, LinkedList{Int64, Int64}}}}
    CellLinkedList(GridSpacing::T ) where{T<:Real} = new{Int64}(x-> Get3DHash(x, GridSpacing), SortedList3D([0,0,0], 1))
end

@inline function clear(CL::CellLinkedList{I}) where {I<:Integer} 
    CL.data = SortedList3D([0,0,0], 1) ### add first bead
    nothing
end

=#

struct LJ_Repulsion{T<:Real, I<:Integer} <: AbstractSelfAvoidanceParameters
    ϵ::T #### saved as 4*epsilon
    type::Vector{I}
    sqr_cutoffs::Matrix{T}
    σ_sqr_mean::Matrix{T}
    #CL::CellLinkedList{I}
    LJ_Repulsion(BeadType::Vector{I}, TypeToSigma::Dict{I, T}, ϵ::T) where {T<:Real, I<:Integer} =begin
        BeadType_ = BeadType
        BeadType = collect(Set(BeadType))

        σ = zeros(T, length(BeadType))
        σ_mean = zeros(T, (length(BeadType),length(BeadType)))
        for (id, type) in enumerate(BeadType)
            σ[id] = TypeToSigma[type]
        end

        for (id, type) in enumerate(BeadType)
            for (id2, type2) in enumerate(BeadType)
                σ_mean[id, id2] = (σ[id]+σ[id2])/2.0
            end
        end
        #SMatrix{length(BeadType),length(BeadType)}((
        sqr_cutoffs =(2.0^(1.0/6.0).*σ_mean).^2

        return new{T,I}(4*ϵ,BeadType_, sqr_cutoffs, σ_mean.^2)#, CellLinkedList(maximum(sqr_cutoffs)) )
    end
end

function InitSimParam(data::SimData,param::LJ_Repulsion ) 
    nothing
end

#=
@inline function clear(CL::CellList{I}) where{I<:Integer} 
    CL.dict = Dict{Vector3{I}, Vector{I}}()
    nothing
end
=#

@inline function clear(data::SimData,param::LJ_Repulsion)
    #clear(param.CL)
    nothing
end

function get_non_bonded_cutoff(data::SimData,NonBonded::LJ_Repulsion)
    return sqrt(maximum(NonBonded.sqr_cutoffs))
end


###naiv implementation
@inline function GetTrialBoltzmannWeight(data::SimData,param::LJ_Repulsion)
    fill!(data.tmp4, 0.0) ### saves energy for every test position  

    @inbounds data.x_arr .= getindex.(data.trial_positions, 1) 
    @inbounds data.y_arr .= getindex.(data.trial_positions, 2) 
    @inbounds data.z_arr .= getindex.(data.trial_positions, 3) 
    ### preselect the indices according to locality sensitive hashing
    #off = (data.id-1) * size(param.sqr_cutoffs)[1]
    #println("id $(data.id), size( $(size(param.sqr_cutoffs)[1])), off = $(off)")

    ind2=param.type[data.id]

    for id in getPotentialNeighbors(data)
 
        LoopVectorization.@avx data.btmp  .= (data.x_arr.-data.xyz[id][1]).^2 
        LoopVectorization.@avx data.btmp .+= (data.y_arr.-data.xyz[id][2]).^2 
        LoopVectorization.@avx data.btmp .+= (data.z_arr.-data.xyz[id][3]).^2


        ind1= param.type[id]
        #data.rand_val =  @inbounds param.sqr_cutoffs[id , data.id]
        data.rand_val =  @inbounds param.sqr_cutoffs[ind1,ind2]

        @inbounds fill!(data.btmp3, data.rand_val)

        LoopVectorization.@avx data.tmp5 .= data.btmp3.>data.btmp ### mask for cutoff, 0 if larger than cutoff

        #LoopVectorization.@avx data.tmp5 .= param.sqr_cutoffs[id, data.id].>data.btmp ### mask for cutoff, 0 if larger than cutoff
        #if ~all(data.tmp5 .==0)

        data.rand_val =  @inbounds param.sqr_cutoffs[ind1,ind2]

        @inbounds  fill!(data.btmp3, param.σ_sqr_mean[ind1,ind2])

        LoopVectorization.@avx data.btmp .= data.btmp3./data.btmp ### ^2

        LoopVectorization.@avx data.btmp .= data.btmp.*data.btmp.*data.btmp ### ^6

        LoopVectorization.@avx data.tmp4 .+= data.tmp5.* param.ϵ .*(data.btmp.*data.btmp .-data.btmp.+0.25)### +0.25 is offset at cutoff, ϵ is actually 4*ϵ
        #end
    end
    LoopVectorization.@avx data.LogBoltzmannFaktor .-=  data.tmp4
    nothing
end


#=
###CellList implementation
@inline function GetTrialBoltzmannWeight(data::SimData,param::LJ_Repulsion)
    fill!(data.tmp4, 0.0) ### saves energy for every test position  

    ### add last added bead that is not direct neighbour
    if  data.id-2>0
        add(param.CL, data.xyz[data.id-2], data.id-2)
    end

    #for tid in 1:data.NTrials
        for neigh in neighbours(param.CL, data.xyz[data.id]) 
            for nid in param.CL.dict[neigh]
                data.btmp[tid] = sqr_norm(data.trial_positions[tid]-data.xyz[nid])### square distance
                if param.sqr_cutoffs[nid, data.id]>data.btmp[tid]
                    data.btmp[tid] = param.σ_sqr_mean[nid, data.id]/data.btmp[tid] ### ^2
                    data.btmp[tid] = data.btmp[tid]*data.btmp[tid]*data.btmp[tid] ### ^6
                    data.tmp4[tid] += param.ϵ*(data.btmp[tid]*data.btmp[tid]-data.btmp[tid]+0.25)### +0.25 is offset at cutoff, ϵ is actually 4*ϵ
                end
            end
        end
    #end
    data.LogBoltzmannFaktor .-=  data.tmp4
    nothing
end
=#

#=
###CellList implementation
@inline function GetTrialBoltzmannWeight(data::SimData,param::LJ_Repulsion)
    fill!(data.tmp4, 0.0) ### saves energy for every test position  

    ### add last added bead that is not direct neighbour
    if  data.id-2>0
        add(param.CL, data.xyz[data.id-2], data.id-2)
    end

    if data.id>1
    for neigh in neighbours(param.CL, data.xyz[data.id-1]) 
        for nid in param.CL.dict[neigh]
            for tid in 1:data.NTrials
            data.btmp[tid] = sqr_norm(data.trial_positions[tid]-data.xyz[nid]) ### square distance
            end
            data.tmp5 .= param.sqr_cutoffs[nid, data.id].>data.btmp ### mask for cutoff, 0 if larger than cutoff
            data.btmp .= param.σ_sqr_mean[nid, data.id]./data.btmp ### ^2
            data.btmp .= data.btmp.*data.btmp.*data.btmp ### ^6
            data.tmp4 .+= data.tmp5.* param.ϵ .*(data.btmp.*data.btmp .-data.btmp.+0.25)### +0.25 is offset at cutoff, ϵ is actually 4*ϵ
        end
    end
    data.LogBoltzmannFaktor .-=  data.tmp4
    end
    nothing
end
=#



@inline function Get3DHash(Position::Vector3{T}, Spacing::T) where {T<: Real}
    return floor.(Int32, Position./Spacing)
end

function sqr_norm(Vec::Vector{Vector3{T}}, A::Vector3{T}) where {T<:Number}
    return [sqr_norm(x.-A) for x in Vec]
end


#=
###SortedList implementation
@inline function GetTrialBoltzmannWeight(data::SimData,param::LJ_Repulsion)
    fill!(data.tmp4, 0.0) ### saves energy for every test position  

    ### add last added bead that is not direct neighbour
    if  data.id-3>0 ### first bead is already added
        #add(param.CL, data.xyz[data.id-2], data.id-2)
        add!(param.CL.data, param.CL.HashFunc(data.xyz[data.id-2]), data.id-2)
        #println("data.id: $(data.id) CL:")
        #println(param.CL.data)
    end


    if data.id>2
        

        #=
        check = [i for i in 1:data.id-2 if any(sqr_norm(data.trial_positions,data.xyz[i]).<(3.8*2.0^(1.0/6.0))^2)]
        for c in check
            res =false
            for i in neighbours(param.CL.data, param.CL.HashFunc(data.xyz[data.id-1]))
                if i ≈ c 
                    res=true
                end
            end
            if res==false
                println()

                printstyled("$c  was not found \n "; color=:red)
                println("data.id: $(data.id) , chain: $(data.id_in_batch)  , $(data.xyz[data.id-1]), $(param.CL.HashFunc(data.xyz[data.id-1]))")
                #println("rel. hash: $([x.-param.CL.HashFunc(data.xyz[data.id-1]) for x in param.CL.HashFunc.(data.xyz[1:data.id-2])])")
                println("neightbours: ", neighbours(param.CL.data, param.CL.HashFunc(data.xyz[data.id-1])) )
                #for i in neighbours(param.CL.data, param.CL.HashFunc(data.xyz[data.id-1])) println(" i $(i), $(data.xyz[i])"); end
                #println("pos: $(data.xyz[1:data.id-1])")
                println("dist neigh: ", [i for i in 1:data.id-2 if sqr_norm(data.xyz[data.id-1]-data.xyz[i])<(2.0*3.8*2.0^(1.0/6.0))^2])
                #println("trial: $(data.trial_positions)")

                println("dist neigh2: ", check)
                println("pos:  ", data.xyz[check])
                println("inds: ", param.CL.HashFunc.(data.xyz[check])) 
                
            end
        end
        =#
        #=
        for nid in neighbours(param.CL.data, param.CL.HashFunc(data.xyz[data.id-1])) ### take all neighbours based on the last bead since all trial positions are around that point, cutoff is increased
            for tid in 1:data.NTrials
            data.btmp[tid] = sqr_norm(data.trial_positions[tid]-data.xyz[nid]) ### square distance
            end
            data.tmp5 .= param.sqr_cutoffs[nid, data.id].>data.btmp ### mask for cutoff, 0 if larger than cutoff
            if ~all(data.tmp5 .==0)
                data.btmp .= param.σ_sqr_mean[nid, data.id]./data.btmp ### ^2
                data.btmp .= data.btmp.*data.btmp.*data.btmp ### ^6
                data.tmp4 .+= data.tmp5.* param.ϵ .*(data.btmp.*data.btmp .-data.btmp.+0.25)### +0.25 is offset at cutoff, ϵ is actually 4*ϵ
            end
        end
        =#
        #println("energy: $(data.tmp4)")

               for tid in 1:data.NTrials
                for nid in neighbours(param.CL.data, param.CL.HashFunc(data.trial_positions[tid]))
                    data.btmp[tid] = sqr_norm(data.trial_positions[tid]-data.xyz[nid])### square distance
                    if param.sqr_cutoffs[nid, data.id]>data.btmp[tid]
                        data.btmp[tid] = param.σ_sqr_mean[nid, data.id]/data.btmp[tid] ### ^2
                        data.btmp[tid] = data.btmp[tid]*data.btmp[tid]*data.btmp[tid] ### ^6
                        data.tmp4[tid] += param.ϵ*(data.btmp[tid]*data.btmp[tid]-data.btmp[tid]+0.25)### +0.25 is offset at cutoff, ϵ is actually 4*ϵ
                    end
                end
            end
            
    #end 
    data.LogBoltzmannFaktor .-=  data.tmp4
    end
    nothing
end
=#




