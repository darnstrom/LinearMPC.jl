struct KalmanFilter
    F::Matrix{Float64}
    G::Matrix{Float64}
    Gd::Matrix{Float64}
    f_offset::Vector{Float64}
    C::Matrix{Float64}
    Dd::Matrix{Float64}
    h_offset::Vector{Float64}
    K::Matrix{Float64}
    x::Vector{Float64} 
end

struct OffsetFreeObserver
    estimator::KalmanFilter
    C::Matrix{Float64}
    Dd::Matrix{Float64}
    h_offset::Vector{Float64}
    nx::Int
    nd_measured::Int
    nd_offsetfree::Int
    nd_control::Int
    period::Int
    formulation::Symbol
end

function Base.propertynames(observer::OffsetFreeObserver, private::Bool=false)
    public_names = (:estimator, :C, :Dd, :h_offset, :nx, :nd_measured,
                    :nd_offsetfree, :nd_control, :period, :formulation, :x, :d, :K)
    return private ? (public_names..., fieldnames(KalmanFilter)...) : public_names
end

function Base.getproperty(observer::OffsetFreeObserver, name::Symbol)
    if name === :x
        nx = getfield(observer, :nx)
        return @view getfield(observer, :estimator).x[1:nx]
    elseif name === :d
        nx = getfield(observer, :nx)
        ndo = getfield(observer, :nd_offsetfree)
        return @view getfield(observer, :estimator).x[nx+1:nx+ndo]
    elseif name === :K
        return getfield(observer, :estimator).K
    elseif name in fieldnames(OffsetFreeObserver)
        return getfield(observer, name)
    elseif name in fieldnames(KalmanFilter)
        return getproperty(getfield(observer, :estimator), name)
    else
        return getfield(observer, name)
    end
end

get_estimated_disturbance(::KalmanFilter) = zeros(0)
get_estimated_disturbance(observer::OffsetFreeObserver) = collect(observer.d)

function get_current_offset_free_disturbance(observer::OffsetFreeObserver)
    ndc = observer.nd_control
    ndc == 0 && return zeros(0)
    return collect(@view observer.d[1:ndc])
end

function get_offset_free_disturbance_preview(observer::OffsetFreeObserver, Np::Int)
    ndc = observer.nd_control
    ndc == 0 && return zeros(0, Np)
    profile = collect(observer.d)
    preview = zeros(ndc, Np)
    for k in 1:Np
        block = mod(k - 1, observer.period) + 1
        rows = (block - 1) * ndc + 1:block * ndc
        preview[:, k] .= profile[rows]
    end
    return preview
end

function KalmanFilter(F,G,C;Gd=nothing,Dd=nothing,f_offset=nothing, h_offset=nothing, x0=nothing, Q=nothing,R=nothing)
    # Solve equation 
    ny,nx = size(C) 
    nu = size(G,2)

    Gd = isnothing(Gd) ? zeros(nx,0) : Gd
    Dd = isnothing(Dd) ? zeros(ny,0) : Dd
    f_offset = isnothing(f_offset) ? zeros(nx) : f_offset;
    h_offset = isnothing(h_offset) ? zeros(ny) : h_offset;
    x0 = isnothing(x0) ? zeros(nx) : x0;
    Q = isnothing(Q) ? Matrix{Float64}(I,nx,nx) : matrixify(Q,nx);
    R = isnothing(R) ? Matrix{Float64}(I,ny,ny) : matrixify(R,ny);

    P,_ = ared(F',C',R,Q);
    K = P*C'/(C*P*C'+R) 
    return KalmanFilter(F,G,Gd,f_offset,C,Dd,h_offset,K,x0)

end

function set_state!(kf::KalmanFilter,x)
    kf.x .= x
end
function set_state!(observer::OffsetFreeObserver, x, d0=nothing)
    xaug = observer.estimator.x
    if length(x) == length(xaug)
        xaug .= x
    elseif length(x) == observer.nx
        xaug[1:observer.nx] .= x
        if isnothing(d0)
            xaug[observer.nx+1:end] .= 0
        else
            length(d0) == observer.nd_offsetfree || throw(ArgumentError("Offset-free disturbance estimate must have length $(observer.nd_offsetfree)"))
            xaug[observer.nx+1:end] .= d0
        end
    else
        throw(ArgumentError("Observer state must have length $(observer.nx) or $(length(xaug))"))
    end
    return observer.x
end

function get_measured_disturbance(observer::OffsetFreeObserver, d)
    ndm = observer.nd_measured
    isnothing(d) && return ndm == 0 ? nothing : zeros(ndm)
    if length(d) == ndm
        return d
    elseif length(d) == ndm + observer.nd_offsetfree
        return d[1:ndm]
    else
        throw(ArgumentError("Disturbance vector must have length $ndm or $(ndm + observer.nd_offsetfree)"))
    end
end

function predict!(kf::KalmanFilter,u,d=nothing)
    kf.x .= kf.F*kf.x + kf.G*u +kf.f_offset
    isnothing(d) || (kf.x .+= kf.Gd*d)
    return kf.x
end
function predict!(observer::OffsetFreeObserver,u,d=nothing)
    predict!(observer.estimator,u,get_measured_disturbance(observer,d))
    return observer.x
end

function correct!(kf::KalmanFilter,y,d=nothing)
    inov = y - kf.C*kf.x - kf.h_offset
    isnothing(d) || (inov .-= kf.Dd*d)
    kf.x .+= kf.K*inov
end
function correct!(observer::OffsetFreeObserver,y,d=nothing)
    correct!(observer.estimator,y,get_measured_disturbance(observer,d))
    return observer.x
end

function render_observer_codegen(kf::KalmanFilter,fh,fsrc)
    ny,nx = size(kf.C)
    nu = size(kf.G,2)
    nd = size(kf.Gd,2)
    @printf(fh, "#define N_MEASUREMENT %d\n",ny);
    @printf(fh, "#define N_OBSERVER_STATE %d\n",nx);
    @printf(fh, "#define N_OBSERVER_CONTROL %d\n",nu);
    @printf(fh, "#define N_OBSERVER_DISTURBANCE %d\n",nd);
    @printf(fh, "extern c_float MPC_PLANT_DYNAMICS[%d];\n",nx*(1+nx+nu+nd));
    @printf(fh, "extern c_float MPC_MEASUREMENT_FUNCTION[%d];\n",ny*(1+nx+nd));
    @printf(fh, "extern c_float K_TRANSPOSE_OBSERVER[%d];\n",ny*nx);
    fmpc_h = open(joinpath(dirname(pathof(LinearMPC)),"../codegen/mpc_observer.h"), "r");
    write(fh, read(fmpc_h))
    close(fmpc_h)

    write_float_array(fsrc,[kf.f_offset kf.F kf.G kf.Gd]'[:],"MPC_PLANT_DYNAMICS");
    write_float_array(fsrc,[kf.h_offset kf.C kf.Dd]'[:],"MPC_MEASUREMENT_FUNCTION");
    write_float_array(fsrc,kf.K[:],"K_TRANSPOSE_OBSERVER");
    fmpc_src = open(joinpath(dirname(pathof(LinearMPC)),"../codegen/mpc_observer.c"), "r");
    write(fsrc, read(fmpc_src))
    close(fmpc_src)
end

codegen(kf::KalmanFilter,fh,fsrc) = render_observer_codegen(kf,fh,fsrc)
codegen(kf::KalmanFilter, mpc::Union{MPC,ExplicitMPC}, fh, fsrc) = render_observer_codegen(kf,fh,fsrc)

function codegen(observer::OffsetFreeObserver, mpc::Union{MPC,ExplicitMPC}, fh, fsrc)
    render_observer_codegen(observer.estimator, fh, fsrc)

    @printf(fh, "#define N_MEASURED_DISTURBANCE %d\n", observer.nd_measured)
    @printf(fh, "#define N_OFFSET_FREE_DISTURBANCE %d\n", observer.nd_control)
    @printf(fh, "#define N_OFFSET_FREE_PROFILE %d\n", observer.nd_offsetfree)
    @printf(fh, "#define N_OFFSET_FREE_PERIOD %d\n", observer.period)
    @printf(fh, "void mpc_get_estimated_state(c_float* state, c_float* observer_state);\n")
    @printf(fh, "void mpc_get_estimated_disturbance(c_float* disturbance, c_float* observer_state, c_float* measured_disturbance);\n")
    mpc.settings.disturbance_preview &&
        @printf(fh, "void mpc_get_estimated_disturbance_preview(c_float* disturbance, c_float* observer_state, c_float* measured_disturbance);\n")
    if mpc.np > 0
        @printf(fh, "int mpc_compute_control_observer(c_float* control, c_float* observer_state, c_float* reference, c_float* measured_disturbance, c_float* affine_parameter);\n")
    else
        @printf(fh, "int mpc_compute_control_observer(c_float* control, c_float* observer_state, c_float* reference, c_float* measured_disturbance);\n")
    end

    write(fsrc, """
void mpc_get_estimated_state(c_float* state, c_float* observer_state){
    int i;
    for(i=0;i<N_STATE;i++) state[i] = observer_state[i];
}

void mpc_get_estimated_disturbance(c_float* disturbance, c_float* observer_state, c_float* measured_disturbance){
    int i;
    for(i=0;i<N_MEASURED_DISTURBANCE;i++) disturbance[i] = measured_disturbance ? measured_disturbance[i] : 0;
    for(i=0;i<N_OFFSET_FREE_DISTURBANCE;i++) disturbance[N_MEASURED_DISTURBANCE+i] = observer_state[N_STATE+i];
}
""")

    if mpc.settings.disturbance_preview
        write(fsrc, """
void mpc_get_estimated_disturbance_preview(c_float* disturbance, c_float* observer_state, c_float* measured_disturbance){
    int k, i, block, row, offset;
    for(k=0;k<N_DISTURBANCE_PREVIEW_HORIZON;k++){
        offset = k*N_DISTURBANCE_BASE;
        for(i=0;i<N_MEASURED_DISTURBANCE;i++) disturbance[offset+i] = measured_disturbance ? measured_disturbance[i] : 0;
        block = k % N_OFFSET_FREE_PERIOD;
        for(i=0;i<N_OFFSET_FREE_DISTURBANCE;i++){
            row = block*N_OFFSET_FREE_DISTURBANCE + i;
            disturbance[offset+N_MEASURED_DISTURBANCE+i] = observer_state[N_STATE+row];
        }
    }
}
""")
    end

    if mpc.np > 0
        get_disturbance_call = mpc.settings.disturbance_preview ?
            "mpc_get_estimated_disturbance_preview(disturbance, observer_state, measured_disturbance);" :
            "mpc_get_estimated_disturbance(disturbance, observer_state, measured_disturbance);"
        write(fsrc, """
int mpc_compute_control_observer(c_float* control, c_float* observer_state, c_float* reference, c_float* measured_disturbance, c_float* affine_parameter){
    c_float state[N_STATE];
    c_float disturbance[N_DISTURBANCE];
    mpc_get_estimated_state(state, observer_state);
    $get_disturbance_call
    return mpc_compute_control(control, state, reference, disturbance, affine_parameter);
}
""")
    else
        get_disturbance_call = mpc.settings.disturbance_preview ?
            "mpc_get_estimated_disturbance_preview(disturbance, observer_state, measured_disturbance);" :
            "mpc_get_estimated_disturbance(disturbance, observer_state, measured_disturbance);"
        write(fsrc, """
int mpc_compute_control_observer(c_float* control, c_float* observer_state, c_float* reference, c_float* measured_disturbance){
    c_float state[N_STATE];
    c_float disturbance[N_DISTURBANCE];
    mpc_get_estimated_state(state, observer_state);
    $get_disturbance_call
    return mpc_compute_control(control, state, reference, disturbance);
}
""")
    end
end

function codegen(observer::OffsetFreeObserver, fh, fsrc)
    throw(ArgumentError("Need the MPC to generate code for OffsetFreeObserver"))
end

function get_control_disturbance(mpc::Union{MPC,ExplicitMPC}, d=nothing)
    observer = mpc.state_observer
    !(observer isa OffsetFreeObserver) && return d
    d_est = mpc.settings.disturbance_preview ?
            get_offset_free_disturbance_preview(observer, mpc.Np) :
            get_current_offset_free_disturbance(observer)
    ndm = observer.nd_measured

    if isnothing(d)
        d_measured = mpc.settings.disturbance_preview ? zeros(ndm, mpc.Np) : zeros(ndm)
        return [d_measured; d_est]
    elseif d isa AbstractMatrix && size(d, 1) == ndm
        d_measured = format_measured_disturbance_preview(d, ndm, size(d_est, 2))
        return [d_measured; d_est]
    elseif d isa AbstractMatrix && size(d, 1) == mpc.model.nd
        return d
    elseif d isa AbstractVector && length(d) == ndm
        if mpc.settings.disturbance_preview
            return [repeat(float(d), 1, mpc.Np); d_est]
        else
            return [d; d_est]
        end
    elseif length(d) == mpc.model.nd
        return d
    else
        throw(ArgumentError("Disturbance must have $(observer.nd_measured) measured rows or $(mpc.model.nd) controller rows"))
    end
end

function format_measured_disturbance_preview(d::AbstractMatrix, nd::Int, Np::Int)
    nd == 0 && return zeros(0, Np)
    if size(d, 2) >= Np
        return float(d[:, 1:Np])
    end
    d_extended = zeros(nd, Np)
    d_extended[:, 1:size(d, 2)] .= d
    d_extended[:, size(d, 2)+1:end] .= repeat(d[:, end], 1, Np - size(d, 2))
    return d_extended
end

function periodic_offset_free_target_preview(observer::OffsetFreeObserver, F, G, Cz, r;
        Bd=nothing, Cd=nothing, Np=observer.period, atol=1e-8)
    observer.formulation == :periodic ||
        throw(ArgumentError("periodic_offset_free_target_preview requires a periodic offset-free observer"))
    nd, period = observer.nd_control, observer.period
    nx, nu = size(G)
    nz = size(Cz, 1)
    Bd = isnothing(Bd) ? observer.estimator.F[1:nx, nx+1:nx+nd] : float(Bd)
    Cd = isnothing(Cd) ? zeros(size(observer.C, 1), nd) : float(Cd)
    rmat = r isa AbstractVector ? reshape(float(r), nz, :) : float(r)
    size(rmat, 1) == nz || throw(ArgumentError("Reference must have $nz rows"))
    if size(rmat, 2) < period
        r_full = zeros(nz, period)
        r_full[:, 1:size(rmat, 2)] .= rmat
        r_full[:, size(rmat, 2)+1:end] .= repeat(rmat[:, end], 1, period - size(rmat, 2))
    else
        r_full = rmat[:, 1:period]
    end

    S = cyclic_shift_matrix(period)
    Sx = kron(S, Matrix{Float64}(I, nx, nx))
    FN = kron(Matrix{Float64}(I, period, period), F)
    GN = kron(Matrix{Float64}(I, period, period), G)
    BN = kron(Matrix{Float64}(I, period, period), Bd)
    CN = kron(Matrix{Float64}(I, period, period), observer.C)
    CDN = kron(Matrix{Float64}(I, period, period), Cd)
    HN = kron(Matrix{Float64}(I, period, period), Cz)

    A = [FN - Sx GN; HN * CN zeros(nz * period, nu * period)]
    b = [-BN * collect(observer.d); vec(r_full) - HN * CDN * collect(observer.d)]
    z = pinv(A; atol) * b
    xbar = reshape(z[1:nx*period], nx, period)
    ubar = reshape(z[nx*period+1:end], nu, period)
    return get_periodic_preview(xbar, Np), get_periodic_preview(ubar, Np)
end

function get_periodic_preview(profile::AbstractMatrix, Np::Int)
    n, period = size(profile)
    out = zeros(n, Np)
    for k in 1:Np
        out[:, k] .= profile[:, mod(k - 1, period) + 1]
    end
    return out
end

simulation_disturbance_dim(mpc::Union{MPC,ExplicitMPC}) = mpc.state_observer isa OffsetFreeObserver ? mpc.state_observer.nd_measured : mpc.model.nd
get_estimated_disturbance(mpc::Union{MPC,ExplicitMPC}) = isnothing(mpc.state_observer) ? zeros(0) : get_estimated_disturbance(mpc.state_observer)

"""
    predict_state!(mpc,u,d=nothing)
Predict the state at the next time step if the control `u` is applied.
This updates the state of `state_observer`
"""
function predict_state!(mpc::Union{MPC,ExplicitMPC},u,d=nothing)
    predict!(mpc.state_observer,u,d)
end

"""
    correct_state!(mpc,y,d=nothing)
Correct the state estimated based on measurement `u`.
This updates the state of `state_observer`
"""
function correct_state!(mpc::Union{MPC,ExplicitMPC},y,d=nothing)
    correct!(mpc.state_observer,y,d)
end

"""
    set_state!(mpc,x)
Set the state of `state_observer` to `x`  
"""
function set_state!(mpc::Union{MPC,ExplicitMPC},x)
    set_state!(mpc.state_observer,x)
end

"""
    get_state!(mpc)
Get the current state of the observer
"""
function get_state(mpc::Union{MPC,ExplicitMPC})
    return mpc.state_observer.x
end

function update_state!(mpc::Union{MPC,ExplicitMPC},u,y)
    isnothing(u) || predict!(mpc.state_observer,u)
    isnothing(y) || correct!(mpc.state_observer,y)
    return mpc.state_observer.x
end
