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

function codegen(kf::KalmanFilter,fh,fsrc)
    ny,nx = size(kf.C)
    nu = size(kf.G,2)
    @printf(fh, "#define N_MEASUREMENT %d\n",ny);
    @printf(fh, "extern F_OBSERVER[%d]\n",nx*nx);
    @printf(fh, "extern G_OBSERVER[%d]\n",nx*nu);
    @printf(fh, "extern C_OBSERVER[%d]\n",ny*nx);
    @printf(fh, "extern K_TRANSPOSE_OBSERVER[%d]\n",ny*nx);
    fmpc_h = open(joinpath(dirname(pathof(LinearMPC)),"../codegen/mpc_observer.h"), "r");
    write(fh, read(fmpc_h))
    close(fmpc_h)

    write_float_array(fsrc,kf.F[:],"F_OBSERVER");
    write_float_array(fsrc,kf.G[:],"G_OBSERVER");
    write_float_array(fsrc,kf.C[:],"C_OBSERVER");
    write_float_array(fsrc,kf.K'[:],"K_TRANSPOSE_OBSERVER");
    fmpc_src = open(joinpath(dirname(pathof(LinearMPC)),"../codegen/mpc_observer.c"), "r");
    write(fsrc, read(fmpc_src))
    close(fmpc_src)
end
