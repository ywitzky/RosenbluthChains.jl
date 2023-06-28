using Mmap
import Base.finalizer

mutable struct Output{T<: Real}
    xFileName::String
    yFileName::String
    zFileName::String
    xio::IOStream
    yio::IOStream
    zio::IOStream
    x::Matrix{T}
    y::Matrix{T}
    z::Matrix{T}
    Output(Path::String, NBeads::I,NSteps::I, _::T)  where {T<:Real, I<:Integer} = begin 
        xName = Path*"x.bin"
        yName = Path*"y.bin"
        zName = Path*"z.bin"
        xio = open(xName,"w+")
        yio = open(yName,"w+")
        zio = open(zName,"w+")
        xout =  Mmap.mmap(xio, Matrix{eltype(T)}, (NBeads,NSteps))
        yout =  Mmap.mmap(yio, Matrix{eltype(T)}, (NBeads,NSteps))
        zout =  Mmap.mmap(zio, Matrix{eltype(T)}, (NBeads,NSteps))
        new{T}(xName,yName,zName, xio,yio,zio, xout, yout, zout)
    end
end 

### closes Syncs and Closes IOStream, as soon as Output struct is not reachable, similar to deconstructor
finalizer(x::Output{T}) where {T<:Real} =  CloseOutput(x)

function CloseOutputOut(Out::Output)
    Mmap.sync!(Out.x)
    Mmap.sync!(Out.y)
    Mmap.sync!(Out.z)
    close(Out.xio)
    close(Out.yio)
    close(Out.zio)
end

#Out= Output(Path::String, data.NBeads::I,Iterations, one(T))