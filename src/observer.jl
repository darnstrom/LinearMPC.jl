struct KalmanFilter
    F::Matrix{Float64}
    G::Matrix{Float64}
    C::Matrix{Float64}
    K::Matrix{Float64}
    x::Vector{Float64} 
end

function KalmanFilter(model::Model;x0=nothing, u0=nothing, Q=nothing,R=nothing)
    # Solve equation 
    F = model.F
    G = model.G
    C = model.C
    ny,nx = size(C) 
    nu = size(G,2)

    x0 = isnothing(x0) ? zeros(nx) : x0;
    Q = isnothing(Q) ? Matrix{Float64}(I,nx,nx) : matrixify(Q,nx);
    R = isnothing(R) ? Matrix{Float64}(I,ny,ny) : matrixify(R,nu);

    P,_ = ared(F,C',R,Q);
    K = P*C'*inv(C*P*C'+R) 
    return KalmanFilter(F,G,C,K,x0)

end

function KalmanFilter(mpc;x0=nothing,Q=I,R=I)
    return Kalmanfilter(mpc.model;x0,Q,R)
end

function set_state!(kf::KalmanFilter,x)
    kf.x[:] = x
end
function predict!(kf::KalmanFilter,u)
    kf.x[:] = kf.F*kf.x + kf.G*u  
end

function correct!(kf::KalmanFilter,y)
    inov = y - kf.C*kf.x
    kf.x[:] += kf.K*inov
end
