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
    f_offset::Vector{Float64}

    xo::Vector{Float64}
    uo::Vector{Float64}

    wmin::Vector{Float64}
    wmax::Vector{Float64}

    C::Matrix{Float64}
    Dd::Matrix{Float64}

    true_dynamics::Function
    
    nx::Int
    nu::Int
    ny::Int
    nd::Int

    Ts::Float64

    labels::Labels
end

function Model(F,G,Gd,C,Dd;Ts=-1.0, f_offset=zeros(0),xo=zeros(0),uo=zeros(0),true_dynamics=nothing)
    Model(F,G;Gd,f_offset,xo,uo,C,Dd,Ts,true_dynamics)
end

function Model(F,G;Ts=-1.0, C = zeros(0,0), Gd = zeros(0,0), f_offset=zeros(0), xo=zeros(0),uo=zeros(0),
        Dd = zeros(0,0), wmin=zeros(0), wmax=zeros(0), true_dynamics=nothing)
    G = reshape(G,size(G,1),:) 
    nx,nu = size(G)
    C = isempty(C) ? Matrix{Float64}(I,nx,nx) : float(C)
    ny = size(C,1);
    (size(C,2)==size(F,1)==nx) || throw(ArgumentError("Dimensions of ss-model incompatible"))
    # disturbance
    Gd = isempty(Gd) ? zeros(nx,0) : Gd
    Dd = isempty(Dd) ? zeros(ny,0) : Dd 
    f_offset = isempty(f_offset) ? zeros(nx) : f_offset
    xo = isempty(xo) ? zeros(nx) : xo
    uo = isempty(uo) ? zeros(nu) : uo
    wmin = isempty(wmin) ? zeros(nx) : wmin
    wmax = isempty(wmax) ? zeros(nx) : wmax
    nd = max(size(Gd,2),size(Dd,2))
    Gd = [Gd zeros(nx,nd-size(Gd,2))]
    Dd = [Dd zeros(ny,nd-size(Dd,2))]
    true_dynamics = isnothing(true_dynamics) ? (x,u,d)->F*x+G*u+Gd*d+f_offset : true_dynamics
    Model(float(F),float(G),float(Gd), float(f_offset), float(xo), float(uo),
          float(wmin), float(wmax), float(C),float(Dd), true_dynamics,
          nx,nu,ny,nd,Ts,Labels(nx,nu,ny,nd))
end

function Model(A,B,Ts; Bd = zeros(0,0), C = zeros(0,0), Dd = zeros(0,0), f_offset=zeros(0),xo=zeros(0),uo=zeros(0),
        true_dynamics=nothing)
    dims = size(B);
    nx,nu = length(dims)==1 ? (dims[1],1) : dims
    (size(A,1) == nx) || throw(ArgumentError("Dimensions of ss-model incompatible"))
    Bd = isempty(Bd) ? zeros(nx,0) : Bd
    f_offset  = isempty(f_offset) ? zeros(nx) : f_offset 
    F,Gext =zoh(A,[B Bd f_offset],Ts)
    G,Gd,f_offset = Gext[:,1:nu], Gext[:,nu+1:nu+size(Bd,2)], Gext[:,end]
    f  = isnothing(true_dynamics) ?  nothing : (x,u,d)->x+Ts*true_dynamics(x,u,d)
    return Model(F,G;Ts,Gd,C,Dd,f_offset,xo,uo,true_dynamics=f)
end

function Model(A,B,Bd,C,Dd,Ts::AbstractFloat;f_offset=zeros(0),xo=zeros(0),uo=zeros(0), true_dynamics=nothing)
    Model(A,B,Ts;Bd,C,Dd,f_offset,xo,uo,true_dynamics)
end

using ForwardDiff

function linearize(f,h,x,u;d=zeros(0))

    nx,nu,nd = length(x),length(u),length(d)
    fz = z->f(z[1:nx],z[nx+1:nx+nu],z[nx+nu+1:nx+nu+nd])
    F = ForwardDiff.jacobian(fz,[x;u;d])
    A,B,Bd = F[:,1:nx],F[:,nx+1:nx+nu],F[:,nx+nu+1:end]
    f_offset = f(x,u,d)-A*x-B*u-Bd*d

    hz = z->h(z[1:nx],z[nx+1:nx+nu],z[nx+nu+1:nx+nu+nd])
    H = ForwardDiff.jacobian(hz,[x;u;d])
    C,D,Dd = H[:,1:nx],H[:,nx+1:nx+nu],H[:,nx+nu+1:end]
    return A,B,Bd,C,D,Dd,f_offset
end

function Model(f,h,x::AbstractVector,u::AbstractVector,Ts;d=zeros(0))
    A,B,Bd,C,D,Dd,f_offset = linearize(f,h,x,u;d)
    iszero(D) || throw(ArgumentError("Non-proper system"))
    return Model(A,B,Ts;Bd,C,Dd,f_offset,true_dynamics=f,xo=x,uo=u)
end

function Model(f,h,x::AbstractVector,u::AbstractVector;d=zeros(0),Ts=-1.0)
    F,G,Gd,C,D,Dd,f_offset = linearize(f,h,x,u;d)
    iszero(D) || throw(ArgumentError("Non-proper system"))
    return Model(F,G;Gd,C,Dd,Ts,f_offset,true_dynamics=f,xo=x,uo=u)
end
