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
    formulation::Symbol
end

function Base.propertynames(observer::OffsetFreeObserver, private::Bool=false)
    public_names = (:estimator, :C, :Dd, :h_offset, :nx, :nd_measured,
                    :nd_offsetfree, :formulation, :x, :d, :K)
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

function codegen(kf::KalmanFilter,fh,fsrc)
    ny,nx = size(kf.C)
    nu = size(kf.G,2)
    nd = size(kf.Gd,2)
    @printf(fh, "#define N_MEASUREMENT %d\n",ny);
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

function codegen(::OffsetFreeObserver, fh, fsrc)
    throw(ArgumentError("Code generation for OffsetFreeObserver is not yet supported"))
end

function get_control_disturbance(mpc::Union{MPC,ExplicitMPC}, d=nothing)
    observer = mpc.state_observer
    !(observer isa OffsetFreeObserver) && return d

    if isnothing(d)
        d_measured = observer.nd_measured == 0 ? zeros(0) : zeros(observer.nd_measured)
        return [d_measured; get_estimated_disturbance(observer)]
    elseif length(d) == observer.nd_measured
        return [d; get_estimated_disturbance(observer)]
    elseif length(d) == mpc.model.nd
        return d
    else
        throw(ArgumentError("Disturbance vector must have length $(observer.nd_measured) or $(mpc.model.nd)"))
    end
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
