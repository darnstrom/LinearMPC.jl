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

    C::Matrix{Float64}
    Dd::Matrix{Float64}
    
    nx::Int
    nu::Int
    ny::Int
    nd::Int

    Ts::Float64

    labels::Labels
end

function Model(F,G,Gd,C,Dd;Ts=1.0)
    Model(F,G;Gd,C,Dd,Ts)
end

function Model(F,G;Ts=1.0, C = zeros(0,0), Gd = zeros(0,0), Dd = zeros(0,0))
    G = reshape(G,size(G,1),:) 
    nx,nu = size(G)
    C = isempty(C) ? Matrix{Float64}(I,nx,nx) : float(C)
    ny = size(C,1);
    (size(C,2)==size(F,1)==nx) || throw(ArgumentError("Dimensions of ss-model incompatible"))
    # disturbance
    Gd = isempty(Gd) ? zeros(nx,0) : Gd
    Dd = isempty(Dd) ? zeros(ny,0) : Dd 
    nd = max(size(Gd,2),size(Dd,2))
    Gd = [Gd zeros(nx,nd-size(Gd,2))]
    Dd = [Dd zeros(ny,nd-size(Dd,2))]
    Model(float(F),float(G),float(Gd),float(C),float(Dd),nx,nu,ny,nd,Ts,Labels(nx,nu,ny,nd))
end

function Model(A,B,Ts; Bd = zeros(0,0), C = zeros(0,0), Dd = zeros(0,0))
    (size(A,1)==size(B,1)) || throw(ArgumentError("Dimensions of ss-model incompatible"))
    F,G,Gd=zoh(A,B,Ts;Bd)
    return Model(F,G;Ts,Gd,C,Dd)
end

function Model(A,B,Bd,C,Dd,Ts::AbstractFloat)
    Model(A,B,Ts;Bd,C,Dd)
end
