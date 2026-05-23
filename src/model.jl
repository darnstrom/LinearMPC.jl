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
    D::Matrix{Float64}
    Dd::Matrix{Float64}
    h_offset::Vector{Float64}

    true_dynamics::Function
    true_h::Function
    
    nx::Int
    nu::Int
    ny::Int
    nd::Int

    Ts::Float64

    labels::Labels
end

function Model(F,G,Gd,C,Dd;Ts=-1.0, f_offset=zeros(0), h_offset=zeros(0),
        xo=zeros(0),uo=zeros(0),true_dynamics=nothing,true_h=nothing)
    Model(F,G;Gd,f_offset,h_offset,xo,uo,C,Dd,Ts,true_dynamics,true_h)
end

function Model(F,G,Gd,C,D,Dd;Ts=-1.0, f_offset=zeros(0), h_offset=zeros(0),
        xo=zeros(0),uo=zeros(0),true_dynamics=nothing,true_h=nothing)
    Model(F,G;Gd,f_offset,h_offset,xo,uo,C,D,Dd,Ts,true_dynamics,true_h)
end

function Model(F,G;Ts=-1.0, C = zeros(0,0), D = zeros(0,0), Gd = zeros(0,0), f_offset=zeros(0), h_offset=zeros(0), 
        xo=zeros(0),uo=zeros(0), Dd = zeros(0,0), wmin=zeros(0), wmax=zeros(0), 
        true_dynamics=nothing, true_h=nothing)
    G = reshape(G,size(G,1),:) 
    nx,nu = size(G)
    C = isempty(C) ? Matrix{Float64}(I,nx,nx) : float(C)
    ny = size(C,1);
    (size(C,2)==size(F,1)==nx) || throw(ArgumentError("Dimensions of ss-model incompatible"))
    D = isempty(D) ? zeros(ny,nu) : float(D)
    size(D) == (ny,nu) || throw(ArgumentError("D must have size ($ny, $nu)"))
    # disturbance
    Gd = isempty(Gd) ? zeros(nx,0) : Gd
    Dd = isempty(Dd) ? zeros(ny,0) : Dd 
    f_offset = isempty(f_offset) ? zeros(nx) : f_offset
    h_offset = isempty(h_offset) ? zeros(ny) : h_offset
    xo = isempty(xo) ? zeros(nx) : xo
    uo = isempty(uo) ? zeros(nu) : uo
    wmin = isempty(wmin) ? zeros(nx) : wmin
    wmax = isempty(wmax) ? zeros(nx) : wmax
    nd = max(size(Gd,2),size(Dd,2))
    Gd = [Gd zeros(nx,nd-size(Gd,2))]
    Dd = [Dd zeros(ny,nd-size(Dd,2))]
    true_dynamics = isnothing(true_dynamics) ? (x,u,d)->F*x+G*u+Gd*d+f_offset : true_dynamics
    true_h = isnothing(true_h) ? (x,u,d)->C*x+D*u+Dd*d+h_offset : true_h
    Model(float(F),float(G),float(Gd), float(f_offset), float(xo), float(uo),
          float(wmin), float(wmax), float(C),float(D),float(Dd), float(h_offset),
          true_dynamics,true_h,
          nx,nu,ny,nd,Ts,Labels(nx,nu,ny,nd))
end

function Model(A,B,Ts; Bd=zeros(0,0), C=zeros(0,0), D=zeros(0,0), Dd=zeros(0,0), f_offset=zeros(0), h_offset=zeros(0),
        xo=zeros(0),uo=zeros(0),true_dynamics=nothing, true_h = nothing)
    dims = size(B);
    nx,nu = length(dims)==1 ? (dims[1],1) : dims
    (size(A,1) == nx) || throw(ArgumentError("Dimensions of ss-model incompatible"))
    Bd = isempty(Bd) ? zeros(nx,0) : Bd
    f_offset  = isempty(f_offset) ? zeros(nx) : f_offset 
    F,Gext =zoh(A,[B Bd f_offset],Ts)
    G,Gd,f_offset = Gext[:,1:nu], Gext[:,nu+1:nu+size(Bd,2)], Gext[:,end]
    f  = isnothing(true_dynamics) ?  nothing : (x,u,d)->x+Ts*true_dynamics(x,u,d)
    true_h  = isnothing(true_h) ?  nothing : true_h 
    return Model(F,G;Ts,Gd,C,D,Dd,f_offset,h_offset,xo,uo,true_dynamics=f,true_h)
end

function Model(A,B,Bd,C,Dd,Ts::AbstractFloat;f_offset=zeros(0),h_offset=zeros(0),
        xo=zeros(0),uo=zeros(0), true_dynamics=nothing, true_h=nothing)
    Model(A,B,Ts;Bd,C,Dd,f_offset,h_offset,xo,uo,true_dynamics,true_h)
end

function Model(A,B,Bd,C,D,Dd,Ts::AbstractFloat;f_offset=zeros(0),h_offset=zeros(0),
        xo=zeros(0),uo=zeros(0), true_dynamics=nothing, true_h=nothing)
    Model(A,B,Ts;Bd,C,D,Dd,f_offset,h_offset,xo,uo,true_dynamics,true_h)
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
    h_offset = h(x,u,d)-C*x-D*u-Dd*d
    return A,B,Bd,C,D,Dd,f_offset,h_offset
end

function Model(f,h,x::AbstractVector,u::AbstractVector,Ts;d=zeros(0))
    A,B,Bd,C,D,Dd,f_offset,h_offset = linearize(f,h,x,u;d)
    return Model(A,B,Ts;Bd,C,D,Dd,f_offset,h_offset,true_dynamics=f,true_h=h,xo=x,uo=u)
end

function Model(f,h,x::AbstractVector,u::AbstractVector;d=zeros(0),Ts=-1.0)
    F,G,Gd,C,D,Dd,f_offset,h_offset = linearize(f,h,x,u;d)
    return Model(F,G;Gd,C,D,Dd,Ts,f_offset,h_offset,true_dynamics=f,true_h=h,xo=x,uo=u)
end

struct MLDModel
    model::Model
    ncontrols::Int
    ndelta::Int
    nz::Int
    Eu::Matrix{Float64}
    Edelta::Matrix{Float64}
    Ez::Matrix{Float64}
    Ex::Matrix{Float64}
    be::Vector{Float64}
    zmin::Vector{Float64}
    zmax::Vector{Float64}
end

function MLDModel(F,Bu,Bdelta,Bz;
        C=zeros(0,0), Du=zeros(0,0), Ddelta=zeros(0,0), Dz=zeros(0,0),
        Eu=zeros(0,0), Edelta=zeros(0,0), Ez=zeros(0,0), Ex=zeros(0,0), be=zeros(0),
        bx=zeros(0), by=zeros(0), Ts=-1.0, xo=zeros(0), uo=zeros(0),
        zmin=zeros(0), zmax=zeros(0), delta_labels=nothing, z_labels=nothing)

    Bu = reshape(Bu, size(Bu,1), :)
    Bdelta = reshape(Bdelta, size(Bdelta,1), :)
    Bz = reshape(Bz, size(Bz,1), :)
    nx = size(F, 1)
    size(Bu, 1) == nx || throw(ArgumentError("Bu must have $nx rows"))
    size(Bdelta, 1) == nx || throw(ArgumentError("Bdelta must have $nx rows"))
    size(Bz, 1) == nx || throw(ArgumentError("Bz must have $nx rows"))

    ncontrols = size(Bu, 2)
    ndelta = size(Bdelta, 2)
    nz = size(Bz, 2)
    ny = isempty(C) ? nx : size(C, 1)

    Du = isempty(Du) ? zeros(ny, ncontrols) : float(Du)
    Ddelta = isempty(Ddelta) ? zeros(ny, ndelta) : float(Ddelta)
    Dz = isempty(Dz) ? zeros(ny, nz) : float(Dz)
    size(Du) == (ny, ncontrols) || throw(ArgumentError("Du must have size ($ny, $ncontrols)"))
    size(Ddelta) == (ny, ndelta) || throw(ArgumentError("Ddelta must have size ($ny, $ndelta)"))
    size(Dz) == (ny, nz) || throw(ArgumentError("Dz must have size ($ny, $nz)"))

    Eu = isempty(Eu) ? zeros(0, ncontrols) : float(Eu)
    Edelta = isempty(Edelta) ? zeros(size(Eu,1), ndelta) : float(Edelta)
    Ez = isempty(Ez) ? zeros(size(Eu,1), nz) : float(Ez)
    Ex = isempty(Ex) ? zeros(size(Eu,1), nx) : float(Ex)
    be = isempty(be) ? zeros(size(Eu,1)) : float(be)
    size(Edelta) == (size(Eu,1), ndelta) || throw(ArgumentError("Edelta must have size ($(size(Eu,1)), $ndelta)"))
    size(Ez) == (size(Eu,1), nz) || throw(ArgumentError("Ez must have size ($(size(Eu,1)), $nz)"))
    size(Ex) == (size(Eu,1), nx) || throw(ArgumentError("Ex must have size ($(size(Eu,1)), $nx)"))
    length(be) == size(Eu,1) || throw(ArgumentError("be must have length $(size(Eu,1))"))

    bx = isempty(bx) ? zeros(nx) : float(bx)
    by = isempty(by) ? zeros(ny) : float(by)

    zmin = isempty(zmin) ? fill(-1e30, nz) : float(zmin)
    zmax = isempty(zmax) ? fill(1e30, nz) : float(zmax)
    length(zmin) == nz || throw(ArgumentError("zmin must have length $nz"))
    length(zmax) == nz || throw(ArgumentError("zmax must have length $nz"))

    G = [Bu Bdelta Bz]
    D = [Du Ddelta Dz]
    uo_full = isempty(uo) ? zeros(size(G,2)) : float(uo)
    length(uo_full) == size(G,2) || throw(ArgumentError("uo must have length $(size(G,2))"))
    model = Ts < 0 ?
        Model(F, G; C, D, f_offset=bx, h_offset=by, xo, uo=uo_full) :
        Model(F, G, Ts; C, D, f_offset=bx, h_offset=by, xo, uo=uo_full)

    if !isnothing(delta_labels)
        length(delta_labels) == ndelta || throw(ArgumentError("Need $ndelta delta labels"))
        model.labels.u[ncontrols+1:ncontrols+ndelta] .= Symbol.(delta_labels)
    end
    if !isnothing(z_labels)
        length(z_labels) == nz || throw(ArgumentError("Need $nz auxiliary continuous labels"))
        model.labels.u[ncontrols+ndelta+1:end] .= Symbol.(z_labels)
    end

    return MLDModel(model, ncontrols, ndelta, nz, Eu, Edelta, Ez, Ex, be, zmin, zmax)
end
