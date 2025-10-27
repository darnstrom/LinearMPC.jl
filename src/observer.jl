struct KalmanFilter
    F::Matrix{Float64}
    G::Matrix{Float64}
    C::Matrix{Float64}
    K::Matrix{Float64}
    x::Vector{Float64} 
end

function KalmanFilter(F,G,C;x0=nothing, Q=nothing,R=nothing)
    # Solve equation 
    ny,nx = size(C) 
    nu = size(G,2)

    x0 = isnothing(x0) ? zeros(nx) : x0;
    Q = isnothing(Q) ? Matrix{Float64}(I,nx,nx) : matrixify(Q,nx);
    R = isnothing(R) ? Matrix{Float64}(I,ny,ny) : matrixify(R,nu);

    P,_ = ared(F',C',R,Q);
    K = P*C'/(C*P*C'+R) 
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
    @printf(fh, "extern c_float F_OBSERVER[%d];\n",nx*nx);
    @printf(fh, "extern c_float G_OBSERVER[%d];\n",nx*nu);
    @printf(fh, "extern c_float C_OBSERVER[%d];\n",ny*nx);
    @printf(fh, "extern c_float K_TRANSPOSE_OBSERVER[%d];\n",ny*nx);
    fmpc_h = open(joinpath(dirname(pathof(LinearMPC)),"../codegen/mpc_observer.h"), "r");
    write(fh, read(fmpc_h))
    close(fmpc_h)

    write_float_array(fsrc,kf.F'[:],"F_OBSERVER");
    write_float_array(fsrc,kf.G[:],"G_OBSERVER");
    write_float_array(fsrc,kf.C'[:],"C_OBSERVER");
    write_float_array(fsrc,kf.K[:],"K_TRANSPOSE_OBSERVER");
    fmpc_src = open(joinpath(dirname(pathof(LinearMPC)),"../codegen/mpc_observer.c"), "r");
    write(fsrc, read(fmpc_src))
    close(fmpc_src)
end

function update_state!(mpc,u,y)
    isnothing(u) || predict!(mpc.state_observer,u)
    isnothing(y) || correct!(mpc.state_observer,y)
    return mpc.state_observer.x
end

function set_state_estimate!(mpc,x)
    set_state!(mpc.state_observer,x)
end
