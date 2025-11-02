struct Labels
    x::Vector{Symbol}
    u::Vector{Symbol}
    y::Vector{Symbol}
    d::Vector{Symbol}
end

function Labels(nx::Int,nu::Int,ny::Int,nd::Int)
    xlabel = [Symbol("x"*string(i)) for i in 1:nx]
    ulabel = [Symbol("u"*string(i)) for i in 1:nu]
    ylabel = [Symbol("y"*string(i)) for i in 1:ny]
    dlabel = [Symbol("d"*string(i)) for i in 1:nd]
    return Labels(xlabel,ulabel,ylabel,dlabel)
end

struct Model
    F::Matrix{Float64}
    G::Matrix{Float64}
    Gd::Matrix{Float64}
    h::Vector{Float64}

    wmin::Vector{Float64}
    wmax::Vector{Float64}

    C::Matrix{Float64}
    Dd::Matrix{Float64}
    
    nx::Int
    nu::Int
    ny::Int
    nd::Int

    Ts::Float64

    labels::Labels
end

function Model(F,G,Gd,C,Dd;Ts=-1.0, h=zeros(0))
    Model(F,G;Gd,h,C,Dd,Ts)
end

function Model(F,G;Ts=-1.0, C = zeros(0,0), Gd = zeros(0,0), h=zeros(0), Dd = zeros(0,0), wmin=zeros(0), wmax=zeros(0))
    G = reshape(G,size(G,1),:) 
    nx,nu = size(G)
    C = isempty(C) ? Matrix{Float64}(I,nx,nx) : float(C)
    ny = size(C,1);
    (size(C,2)==size(F,1)==nx) || throw(ArgumentError("Dimensions of ss-model incompatible"))
    # disturbance
    Gd = isempty(Gd) ? zeros(nx,0) : Gd
    Dd = isempty(Dd) ? zeros(ny,0) : Dd 
    h = isempty(h) ? zeros(nx) : h
    wmin = isempty(wmin) ? zeros(nx) : wmin
    wmax = isempty(wmax) ? zeros(nx) : wmax
    nd = max(size(Gd,2),size(Dd,2))
    Gd = [Gd zeros(nx,nd-size(Gd,2))]
    Dd = [Dd zeros(ny,nd-size(Dd,2))]
    Model(float(F),float(G),float(Gd), float(h),float(wmin), float(wmax), float(C),float(Dd),nx,nu,ny,nd,Ts,Labels(nx,nu,ny,nd))
end

function Model(A,B,Ts; Bd = zeros(0,0), C = zeros(0,0), Dd = zeros(0,0), h=zeros(0))
    dims = size(B);
    nx,nu = length(dims)==1 ? (dims[1],1) : dims
    (size(A,1) == nx) || throw(ArgumentError("Dimensions of ss-model incompatible"))
    Bd = isempty(Bd) ? zeros(nx,0) : Bd
    h  = isempty(h) ? zeros(nx) : h 
    F,Gext =zoh(A,[B Bd h],Ts)
    G,Gd,h = Gext[:,1:nu], Gext[:,nu+1:nu+size(Bd,2)], Gext[:,end]
    return Model(F,G;Ts,Gd,C,Dd,h)
end

function Model(A,B,Bd,C,Dd,Ts::AbstractFloat;h=zeros(0))
    Model(A,B,Ts;Bd,C,Dd,h)
end

using ForwardDiff

function linearize(f,h,x,u;d=zeros(0))

    nx,nu,nd = length(x),length(u),length(d)
    fz = z->f(z[1:nx],z[nx+1:nx+nu],z[nx+nu+1:nx+nu+nd])
    F = ForwardDiff.jacobian(fz,[x;u;d])
    A,B,Bd = F[:,1:nx],F[:,nx+1:nx+nu],F[:,nx+nu+1:end]
    offset = f(x,u,d)-A*x-B*u-Bd*d

    hz = z->h(z[1:nx],z[nx+1:nx+nu],z[nx+nu+1:nx+nu+nd])
    H = ForwardDiff.jacobian(hz,[x;u;d])
    C,D,Dd = H[:,1:nx],H[:,nx+1:nx+nu],H[:,nx+nu+1:end]
    return A,B,Bd,C,D,Dd,offset
end

function Model(f,h,x::AbstractVector,u::AbstractVector,Ts;d=zeros(0))
    A,B,Bd,C,D,Dd,offset = linearize(f,h,x,u;d)
    iszero(D) || throw(ArgumentError("Non-proper system"))
    return Model(A,B,Ts;Bd,C,Dd,h=offset)
end

function Model(f,h,x::AbstractVector,u::AbstractVector;d=zeros(0),Ts=-1.0)
    F,G,Gd,C,D,Dd,offset = linearize(f,h,x,u;d)
    iszero(D) || throw(ArgumentError("Non-proper system"))
    return Model(F,G;Gd,C,Dd,Ts,h=offset)
end
