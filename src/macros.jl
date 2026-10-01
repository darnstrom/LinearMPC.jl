# ============================================================================
# JuMP-like macro interface for LinearMPC
# ============================================================================
#
# Two usage styles are supported:
#
# --- Style 1: full-matrix (known dimensions upfront) ---
#
#   mpc = MPC(F, G; C=C, Np=10)
#   @control mpc u umin=-ones(2) umax=ones(2)
#   @state   mpc x
#   @output  mpc y
#   @disturbance mpc d
#   @dynamics mpc x_next = F_new*x + G_new*u   # full F and G matrices
#   @objective mpc Q=I R=0.1*I
#   @constraint mpc -1 <= u <= 1
#   @constraint mpc 0 <= y <= 5
#   setup!(mpc)
#
# --- Style 2: incremental (dimensions built row by row) ---
#
#   mpc = MPC(Np=10)                  # zero-dimensional MPC
#   @state   mpc x1                   # x1 = StateRef(mpc,1), x1_next also defined
#   @state   mpc v[1:4]               # v  = StateVec(mpc,[2,3,4,5]), v_next also defined
#   @control mpc u[1:2]               # u  = ControlVec(mpc,[1,2])
#   @dynamics mpc x1 = 0.9*x1 + 0.5*v[1]           # set row 1 (can omit _next)
#   @dynamics mpc v[1] = 0.8*v[1] + 0.2*v[2] + u[1] # set row for v[1]
#   @dynamics mpc [v[1]; v[2]] = A*v[1:2] + B*u      # set multiple rows at once
#   @objective mpc Q=I R=0.1
#   @constraint mpc -1 <= u <= 1
#   setup!(mpc)

# ---- Module-level registry for incremental-mode MPCs ----------------------
const _incremental_registry = WeakKeyDict{MPC, NTuple{3,Int}}()

"""
Return true if `mpc` was created in incremental mode (via `MPC(; Np=...)`).
"""
_is_incremental(mpc::MPC) = haskey(_incremental_registry, mpc)

# ---- Signal reference types ------------------------------------------------

"""
    ControlRef(mpc, idx)
    ControlRef(mpc)

Reference to one control input (`idx > 0`) or the full control vector (`idx == 0`).
Created by [`@control`](@ref).
"""
struct ControlRef
    mpc::MPC
    idx::Int
end
ControlRef(mpc::MPC) = ControlRef(mpc, 0)

"""
    StateRef(mpc, idx)
    StateRef(mpc)

Reference to one state (`idx > 0`) or the full state vector (`idx == 0`).
Created by [`@state`](@ref).
"""
struct StateRef
    mpc::MPC
    idx::Int
end
StateRef(mpc::MPC) = StateRef(mpc, 0)

"""
    OutputRef(mpc)

Reference to the output signal. Created by [`@output`](@ref).
"""
struct OutputRef
    mpc::MPC
end

"""
    DisturbanceRef(mpc, idx)
    DisturbanceRef(mpc)

Reference to one disturbance (`idx > 0`) or the full disturbance vector (`idx == 0`).
Created by [`@disturbance`](@ref).
"""
struct DisturbanceRef
    mpc::MPC
    idx::Int
end
DisturbanceRef(mpc::MPC) = DisturbanceRef(mpc, 0)

# ---- Vector signal types ---------------------------------------------------

"""
    StateVec(mpc, indices)

A named vector of state variables covering the given indices.
Created by `@state mpc v[1:4]`.

Supports scalar and range indexing:
- `v[3]`   → `StateRef(mpc, indices[3])`
- `v[2:4]` → `StateVec(mpc, indices[2:4])`
"""
struct StateVec
    mpc::MPC
    indices::Vector{Int}
end
Base.length(sv::StateVec) = length(sv.indices)
Base.getindex(sv::StateVec, i::Int)           = StateRef(sv.mpc, sv.indices[i])
Base.getindex(sv::StateVec, r::AbstractRange) = StateVec(sv.mpc, sv.indices[r])
Base.iterate(sv::StateVec, s=1) = s > length(sv) ? nothing : (sv[s], s+1)

"""
    ControlVec(mpc, indices)

A named vector of control variables. Created by `@control mpc u[1:2]`.
"""
struct ControlVec
    mpc::MPC
    indices::Vector{Int}
end
Base.length(cv::ControlVec) = length(cv.indices)
Base.getindex(cv::ControlVec, i::Int)           = ControlRef(cv.mpc, cv.indices[i])
Base.getindex(cv::ControlVec, r::AbstractRange) = ControlVec(cv.mpc, cv.indices[r])
Base.iterate(cv::ControlVec, s=1) = s > length(cv) ? nothing : (cv[s], s+1)

"""
    DisturbanceVec(mpc, indices)

A named vector of disturbance variables. Created by `@disturbance mpc d[1:2]`.
"""
struct DisturbanceVec
    mpc::MPC
    indices::Vector{Int}
end
Base.length(dv::DisturbanceVec) = length(dv.indices)
Base.getindex(dv::DisturbanceVec, i::Int)           = DisturbanceRef(dv.mpc, dv.indices[i])
Base.getindex(dv::DisturbanceVec, r::AbstractRange) = DisturbanceVec(dv.mpc, dv.indices[r])
Base.iterate(dv::DisturbanceVec, s=1) = s > length(dv) ? nothing : (dv[s], s+1)

# ---- SignalTerm and SignalExpr ---------------------------------------------

"""
Internal: A * [signal_kind] term.

`col_indices`:
- `Int[]`   (empty)  -- dense term spanning all columns.
- `[i]`     (single) -- scalar at column i.
- `[i,j,.]` (multi)  -- sparse block at those columns.
"""
struct SignalTerm
    A::Matrix{Float64}
    kind::Symbol       # :u, :x, :r, :d, :uprev
    col_indices::Vector{Int}
end
SignalTerm(A, kind)           = SignalTerm(A, kind, Int[])
SignalTerm(A, kind, col::Int) = SignalTerm(A, kind, col == 0 ? Int[] : [col])

"""
    SignalExpr(mpc, terms)

A linear combination of weighted MPC signals. Built automatically by the
arithmetic overloads on signal references.
"""
struct SignalExpr
    mpc::MPC
    terms::Vector{SignalTerm}
end

# ---- DynamicsExpr ----------------------------------------------------------

"""
    DynamicsExpr(signals, offset)

RHS of a dynamics equation: weighted signals + optional constant offset.
"""
struct DynamicsExpr
    signals::SignalExpr
    offset::Vector{Float64}
end
DynamicsExpr(expr::SignalExpr) = DynamicsExpr(expr, zeros(0))

# ---- Arithmetic overloads --------------------------------------------------

# Dense matrix * full-vector ref (idx == 0)
Base.:*(A::AbstractMatrix, x::StateRef)       = SignalExpr(x.mpc, [SignalTerm(float(A), :x)])
Base.:*(A::AbstractMatrix, u::ControlRef)     = SignalExpr(u.mpc, [SignalTerm(float(A), :u)])
Base.:*(A::AbstractMatrix, d::DisturbanceRef) = SignalExpr(d.mpc, [SignalTerm(float(A), :d)])

# Dense matrix * StateVec/ControlVec/DisturbanceVec (sparse block at col_indices)
Base.:*(A::AbstractMatrix, sv::StateVec)       = SignalExpr(sv.mpc, [SignalTerm(float(A), :x, sv.indices)])
Base.:*(A::AbstractMatrix, cv::ControlVec)     = SignalExpr(cv.mpc, [SignalTerm(float(A), :u, cv.indices)])
Base.:*(A::AbstractMatrix, dv::DisturbanceVec) = SignalExpr(dv.mpc, [SignalTerm(float(A), :d, dv.indices)])

# Vector coefficient * full-vector ref
Base.:*(A::AbstractVector, x::StateRef)       = SignalExpr(x.mpc, [SignalTerm(reshape(float(A),:,1), :x)])
Base.:*(A::AbstractVector, u::ControlRef)     = SignalExpr(u.mpc, [SignalTerm(reshape(float(A),:,1), :u)])
Base.:*(A::AbstractVector, d::DisturbanceRef) = SignalExpr(d.mpc, [SignalTerm(reshape(float(A),:,1), :d)])

# Scalar * scalar ref (idx > 0) or full-vector ref (idx == 0)
function Base.:*(a::Number, x::StateRef)
    x.idx == 0 ?
        SignalExpr(x.mpc, [SignalTerm(float(a)*Matrix{Float64}(I, x.mpc.model.nx, x.mpc.model.nx), :x)]) :
        SignalExpr(x.mpc, [SignalTerm(fill(Float64(a), 1, 1), :x, x.idx)])
end
function Base.:*(a::Number, u::ControlRef)
    u.idx == 0 ?
        SignalExpr(u.mpc, [SignalTerm(float(a)*Matrix{Float64}(I, u.mpc.model.nu, u.mpc.model.nu), :u)]) :
        SignalExpr(u.mpc, [SignalTerm(fill(Float64(a), 1, 1), :u, u.idx)])
end
function Base.:*(a::Number, d::DisturbanceRef)
    d.idx == 0 ?
        SignalExpr(d.mpc, [SignalTerm(float(a)*Matrix{Float64}(I, d.mpc.model.nd, d.mpc.model.nd), :d)]) :
        SignalExpr(d.mpc, [SignalTerm(fill(Float64(a), 1, 1), :d, d.idx)])
end

# Scalar * StateVec/ControlVec/DisturbanceVec: one sparse term per element
function Base.:*(a::Number, sv::StateVec)
    SignalExpr(sv.mpc, [SignalTerm(fill(Float64(a), 1, 1), :x, [j]) for j in sv.indices])
end
function Base.:*(a::Number, cv::ControlVec)
    SignalExpr(cv.mpc, [SignalTerm(fill(Float64(a), 1, 1), :u, [j]) for j in cv.indices])
end
function Base.:*(a::Number, dv::DisturbanceVec)
    SignalExpr(dv.mpc, [SignalTerm(fill(Float64(a), 1, 1), :d, [j]) for j in dv.indices])
end

# Matrix * Vector{StateRef} / Vector{ControlRef} / Vector{DisturbanceRef}
# (e.g., A * [x1; x2] where [x1;x2] is a vcat of indexed refs)
function Base.:*(A::AbstractMatrix, refs::AbstractVector{<:StateRef})
    isempty(refs) && error("Empty StateRef vector")
    mpc = refs[1].mpc
    all(r -> r.mpc === mpc && r.idx > 0, refs) ||
        error("@dynamics: mixed MPC or unindexed StateRefs in vector expression")
    SignalExpr(mpc, [SignalTerm(float(A), :x, [r.idx for r in refs])])
end
function Base.:*(A::AbstractMatrix, refs::AbstractVector{<:ControlRef})
    isempty(refs) && error("Empty ControlRef vector")
    mpc = refs[1].mpc
    all(r -> r.mpc === mpc && r.idx > 0, refs) ||
        error("@dynamics: mixed MPC or unindexed ControlRefs in vector expression")
    SignalExpr(mpc, [SignalTerm(float(A), :u, [r.idx for r in refs])])
end
function Base.:*(a::Number, refs::AbstractVector{<:StateRef})
    isempty(refs) && error("Empty StateRef vector")
    SignalExpr(refs[1].mpc, [SignalTerm(fill(Float64(a),1,1), :x, [r.idx for r in refs])])
end
function Base.:*(a::Number, refs::AbstractVector{<:ControlRef})
    isempty(refs) && error("Empty ControlRef vector")
    SignalExpr(refs[1].mpc, [SignalTerm(fill(Float64(a),1,1), :u, [r.idx for r in refs])])
end

# Bare signal -> identity SignalExpr
_bare(x::StateRef) = x.idx == 0 ?
    SignalExpr(x.mpc, [SignalTerm(Matrix{Float64}(I, x.mpc.model.nx, x.mpc.model.nx), :x)]) :
    SignalExpr(x.mpc, [SignalTerm(fill(1.0, 1, 1), :x, x.idx)])
_bare(u::ControlRef) = u.idx == 0 ?
    SignalExpr(u.mpc, [SignalTerm(Matrix{Float64}(I, u.mpc.model.nu, u.mpc.model.nu), :u)]) :
    SignalExpr(u.mpc, [SignalTerm(fill(1.0, 1, 1), :u, u.idx)])
_bare(d::DisturbanceRef) = d.idx == 0 ?
    SignalExpr(d.mpc, [SignalTerm(Matrix{Float64}(I, d.mpc.model.nd, d.mpc.model.nd), :d)]) :
    SignalExpr(d.mpc, [SignalTerm(fill(1.0, 1, 1), :d, d.idx)])
_bare(sv::StateVec)       = SignalExpr(sv.mpc, [SignalTerm(fill(1.0,1,1), :x, [j]) for j in sv.indices])
_bare(cv::ControlVec)     = SignalExpr(cv.mpc, [SignalTerm(fill(1.0,1,1), :u, [j]) for j in cv.indices])
_bare(dv::DisturbanceVec) = SignalExpr(dv.mpc, [SignalTerm(fill(1.0,1,1), :d, [j]) for j in dv.indices])
function _bare(refs::AbstractVector{<:StateRef})
    isempty(refs) && error("Empty StateRef vector")
    SignalExpr(refs[1].mpc, [SignalTerm(fill(1.0,1,1), :x, [r.idx for r in refs])])
end
function _bare(refs::AbstractVector{<:ControlRef})
    isempty(refs) && error("Empty ControlRef vector")
    SignalExpr(refs[1].mpc, [SignalTerm(fill(1.0,1,1), :u, [r.idx for r in refs])])
end

# Combine SignalExprs
function Base.:+(e1::SignalExpr, e2::SignalExpr)
    @assert e1.mpc === e2.mpc "Signal expressions must refer to the same MPC"
    SignalExpr(e1.mpc, [e1.terms; e2.terms])
end

# Negation and subtraction
Base.:-(e::SignalExpr) =
    SignalExpr(e.mpc, [SignalTerm(-t.A, t.kind, t.col_indices) for t in e.terms])
Base.:-(e1::SignalExpr, e2::SignalExpr)         = e1 + (-e2)
Base.:-(e::SignalExpr,  x::StateRef)            = e + (-_bare(x))
Base.:-(e::SignalExpr,  u::ControlRef)          = e + (-_bare(u))
Base.:-(e::SignalExpr,  d::DisturbanceRef)      = e + (-_bare(d))
Base.:-(e::SignalExpr,  sv::StateVec)           = e + (-_bare(sv))
Base.:-(e::SignalExpr,  cv::ControlVec)         = e + (-_bare(cv))
Base.:-(x::StateRef,    e::SignalExpr)          = _bare(x) + (-e)
Base.:-(u::ControlRef,  e::SignalExpr)          = _bare(u) + (-e)
Base.:-(sv::StateVec,   e::SignalExpr)          = _bare(sv) + (-e)

# Bare ref + something
for (T, bare_fn) in [(:StateRef, :_bare), (:ControlRef, :_bare),
                      (:DisturbanceRef, :_bare), (:StateVec, :_bare),
                      (:ControlVec, :_bare), (:DisturbanceVec, :_bare)]
    @eval Base.:+(x::$T, e::SignalExpr)  = $bare_fn(x) + e
    @eval Base.:+(e::SignalExpr, x::$T)  = e + $bare_fn(x)
end
Base.:+(x::StateRef,   u::ControlRef)       = _bare(x)  + _bare(u)
Base.:+(u::ControlRef, x::StateRef)         = _bare(u)  + _bare(x)
Base.:+(x::StateRef,   d::DisturbanceRef)   = _bare(x)  + _bare(d)
Base.:+(d::DisturbanceRef, x::StateRef)     = _bare(d)  + _bare(x)
Base.:+(sv::StateVec,  cv::ControlVec)      = _bare(sv) + _bare(cv)
Base.:+(cv::ControlVec, sv::StateVec)       = _bare(cv) + _bare(sv)

# Support + between Vector{StateRef}/Vector{ControlRef} and SignalExpr
Base.:+(e::SignalExpr, refs::AbstractVector{<:StateRef})   = e + _bare(refs)
Base.:+(refs::AbstractVector{<:StateRef}, e::SignalExpr)   = _bare(refs) + e
Base.:+(e::SignalExpr, refs::AbstractVector{<:ControlRef}) = e + _bare(refs)
Base.:+(refs::AbstractVector{<:ControlRef}, e::SignalExpr) = _bare(refs) + e

# SignalExpr + constant offset vector -> DynamicsExpr
Base.:+(e::SignalExpr, ofs::AbstractVector)   = DynamicsExpr(e, float(ofs))
Base.:+(ofs::AbstractVector, e::SignalExpr)   = DynamicsExpr(e, float(ofs))
Base.:+(e::DynamicsExpr, ofs::AbstractVector) = DynamicsExpr(e.signals, e.offset + ofs)
Base.:+(ofs::AbstractVector, e::DynamicsExpr) = DynamicsExpr(e.signals, e.offset + ofs)
Base.:+(e1::DynamicsExpr, e2::SignalExpr)     = DynamicsExpr(e1.signals + e2, e1.offset)
Base.:+(e1::SignalExpr, e2::DynamicsExpr)     = DynamicsExpr(e1 + e2.signals, e2.offset)

# ---- Model expansion helpers (incremental mode) ---------------------------

function _expand_state!(mpc::MPC)
    m = mpc.model
    nx, nu, nd = m.nx, m.nu, m.nd
    new_nx = nx + 1
    F_new  = [m.F          zeros(nx, 1);   zeros(1, nx) zeros(1, 1)]
    G_new  = [m.G;  zeros(1, nu)]
    Gd_new = [m.Gd; zeros(1, nd)]
    C_new  = Matrix{Float64}(I, new_nx, new_nx)   # full-state output
    mpc.model = Model(F_new, G_new;
        Gd=Gd_new, C=C_new, Dd=m.Dd,
        f_offset=[m.f_offset; 0.0], h_offset=zeros(new_nx),
        xo=[m.xo; 0.0], uo=m.uo,
        wmin=[m.wmin; 0.0], wmax=[m.wmax; 0.0], Ts=m.Ts)
    mpc.mpqp_issetup = false
    nx_s, nu_s, nd_s = _incremental_registry[mpc]
    _incremental_registry[mpc] = (nx_s + 1, nu_s, nd_s)
    return new_nx
end

function _expand_control!(mpc::MPC)
    m = mpc.model
    mpc.model = Model(m.F, [m.G zeros(m.nx, 1)];
        Gd=m.Gd, C=m.C, Dd=m.Dd,
        f_offset=m.f_offset, h_offset=m.h_offset,
        xo=m.xo, uo=[m.uo; 0.0], wmin=m.wmin, wmax=m.wmax, Ts=m.Ts)
    mpc.mpqp_issetup = false
    nx_s, nu_s, nd_s = _incremental_registry[mpc]
    _incremental_registry[mpc] = (nx_s, nu_s + 1, nd_s)
    return m.nu + 1
end

function _expand_disturbance!(mpc::MPC)
    m = mpc.model
    mpc.model = Model(m.F, m.G;
        Gd=[m.Gd zeros(m.nx,1)], C=m.C, Dd=[m.Dd zeros(m.ny,1)],
        f_offset=m.f_offset, h_offset=m.h_offset,
        xo=m.xo, uo=m.uo, wmin=m.wmin, wmax=m.wmax, Ts=m.Ts)
    mpc.mpqp_issetup = false
    nx_s, nu_s, nd_s = _incremental_registry[mpc]
    _incremental_registry[mpc] = (nx_s, nu_s, nd_s + 1)
    return m.nd + 1
end

# ---- Internal helpers ------------------------------------------------------

_bound(::Nothing)           = zeros(0)
_bound(x)                   = x
_to_vec_bound(::Nothing, n) = zeros(0)
_to_vec_bound(x::AbstractVector, n) = x
_to_vec_bound(x::Number, n) = fill(Float64(x), n)

_is_signal(::ControlRef)     = true
_is_signal(::OutputRef)      = true
_is_signal(::StateRef)       = true
_is_signal(::DisturbanceRef) = true
_is_signal(::SignalExpr)     = true
_is_signal(::ControlVec)     = true
_is_signal(::StateVec)       = true
_is_signal(::DisturbanceVec) = true
_is_signal(::Any)            = false

# ---- Core constraint dispatch ----------------------------------------------

function _add_constraint!(mpc::MPC, lb, ::ControlRef, ub; kw...)
    set_input_bounds!(mpc; umin=_bound(lb), umax=_bound(ub))
end

function _add_constraint!(mpc::MPC, lb, ::ControlVec, ub; kw...)
    set_input_bounds!(mpc; umin=_bound(lb), umax=_bound(ub))
end

function _add_constraint!(mpc::MPC, lb, ::OutputRef, ub;
        ks=2:mpc.Np, soft=true, binary=false, prio=0, kw...)
    set_output_bounds!(mpc; ymin=_to_vec_bound(lb, mpc.model.ny),
                            ymax=_to_vec_bound(ub, mpc.model.ny),
                            ks=ks, soft=soft, binary=binary, prio=prio)
end

function _add_constraint!(mpc::MPC, lb, ::StateRef, ub;
        ks=2:mpc.Np, soft=false, binary=false, prio=0, kw...)
    add_constraint!(mpc;
        Ax=Matrix{Float64}(I, mpc.model.nx, mpc.model.nx),
        lb=_to_vec_bound(lb, mpc.model.nx),
        ub=_to_vec_bound(ub, mpc.model.nx),
        ks=ks, soft=soft, binary=binary, prio=prio)
end

function _add_constraint!(mpc::MPC, lb, ::DisturbanceRef, ub; kw...)
    set_disturbance!(mpc, _bound(lb), _bound(ub))
end

function _add_constraint!(mpc::MPC, lb, expr::SignalExpr, ub;
        ks=2:mpc.Np, soft=false, binary=false, prio=0, kw...)
    Au=zeros(0,0); Ax=zeros(0,0); Ar=zeros(0,0); Ad=zeros(0,0); Aup=zeros(0,0)
    for t in expr.terms
        A = t.A
        if     t.kind == :u;     Au  = isempty(Au)  ? A : Au  + A
        elseif t.kind == :x;     Ax  = isempty(Ax)  ? A : Ax  + A
        elseif t.kind == :r;     Ar  = isempty(Ar)  ? A : Ar  + A
        elseif t.kind == :d;     Ad  = isempty(Ad)  ? A : Ad  + A
        elseif t.kind == :uprev; Aup = isempty(Aup) ? A : Aup + A
        end
    end
    add_constraint!(mpc;
        Au=isempty(Au) ? nothing : Au, Ax=isempty(Ax) ? nothing : Ax,
        Ar=Ar, Ad=Ad, Aup=Aup,
        lb=_bound(lb), ub=_bound(ub),
        ks=ks, soft=soft, binary=binary, prio=prio)
end

function _add_constraint_onesided!(mpc::MPC, lhs, op::Symbol, rhs; kw...)
    if op == :(<=)
        _is_signal(lhs) ? _add_constraint!(mpc, nothing, lhs, rhs; kw...) :
        _is_signal(rhs) ? _add_constraint!(mpc, lhs, rhs, nothing; kw...) :
            error("@constraint: no signal reference found")
    elseif op == :(>=)
        _is_signal(lhs) ? _add_constraint!(mpc, rhs, lhs, nothing; kw...) :
        _is_signal(rhs) ? _add_constraint!(mpc, nothing, rhs, lhs; kw...) :
            error("@constraint: no signal reference found")
    else
        error("@constraint: unsupported operator $op")
    end
end

# ---- Macro helpers ---------------------------------------------------------

function _parse_macro_kwargs(args)
    kw_pairs = Expr[]
    for arg in args
        (arg isa Expr && arg.head == :(=)) || continue
        push!(kw_pairs, Expr(:kw, arg.args[1], esc(arg.args[2])))
    end
    return kw_pairs
end

function _make_call(fn, kw_pairs, pos_args...)
    escaped = [esc(a) for a in pos_args]
    isempty(kw_pairs) ?
        Expr(:call, fn, escaped...) :
        Expr(:call, fn, Expr(:parameters, kw_pairs...), escaped...)
end

# ---- Dynamics helpers ------------------------------------------------------

"""
Apply one `SignalTerm` to a row `r` of matrices `F_new`, `G_new`, `Gd_new`.
Handles dense (`col_indices` empty), single-column, and multi-column terms.
"""
function _apply_term_to_row!(F_new, G_new, Gd_new, term::SignalTerm, r::Int)
    A, cols = term.A, term.col_indices
    mat = term.kind == :x ? F_new :
          term.kind == :u ? G_new :
          term.kind == :d ? Gd_new : nothing
    mat === nothing && return
    if isempty(cols)
        # Dense: use row r of A (or row 1 if A has only 1 row)
        row_vec = size(A, 1) == 1 ? vec(A) : A[r, :]
        mat[r, :] .+= row_vec
    elseif length(cols) == 1
        mat[r, cols[1]] += A[1, 1]
    else
        # Sparse block: A is (1 x ncols) or (nrows x ncols)
        row_vec = size(A, 1) == 1 ? vec(A) : A[r, :]
        mat[r, cols] .+= row_vec
    end
end

"""
Apply one `SignalTerm` to a block of rows `rows` of F, G, Gd.
"""
function _apply_term_to_rows!(F_new, G_new, Gd_new, term::SignalTerm, rows::Vector{Int})
    A, cols = term.A, term.col_indices
    mat = term.kind == :x ? F_new :
          term.kind == :u ? G_new :
          term.kind == :d ? Gd_new : nothing
    mat === nothing && return
    n = length(rows)
    if isempty(cols)
        # Dense: A is (n x ncols) or (1 x ncols)
        for (i, r) in enumerate(rows)
            row_vec = size(A, 1) == 1 ? vec(A) : A[i, :]
            mat[r, :] .+= row_vec
        end
    elseif length(cols) == 1
        for (i, r) in enumerate(rows)
            mat[r, cols[1]] += size(A,1) == 1 ? A[1,1] : A[i,1]
        end
    else
        for (i, r) in enumerate(rows)
            row_vec = size(A, 1) == 1 ? vec(A) : A[i, :]
            mat[r, cols] .+= row_vec
        end
    end
end

"""
    _set_dynamics!(mpc, expr)

Replace full F, G, Gd, f_offset from a `SignalExpr` or `DynamicsExpr`.
Preserves C, Dd, h_offset, Ts, operating point.
"""
function _set_dynamics!(mpc::MPC, expr::SignalExpr)
    _set_dynamics!(mpc, DynamicsExpr(expr))
end

function _set_dynamics!(mpc::MPC, dexpr::DynamicsExpr)
    m = mpc.model
    nx, nu, nd = m.nx, m.nu, m.nd
    F = zeros(nx, nx); G = zeros(nx, nu); Gd = zeros(nx, nd)
    f_offset = isempty(dexpr.offset) ? zeros(nx) : Vector{Float64}(dexpr.offset)
    for t in dexpr.signals.terms
        if     t.kind == :x; F  .+= t.A
        elseif t.kind == :u; G  .+= t.A
        elseif t.kind == :d; Gd .+= t.A
        end
    end
    mpc.model = Model(F, G;
        Gd=Gd, C=m.C, Dd=m.Dd, f_offset=f_offset, h_offset=m.h_offset,
        Ts=m.Ts, xo=m.xo, uo=m.uo, wmin=m.wmin, wmax=m.wmax)
    mpc.mpqp_issetup = false
    return mpc
end

"""
    _set_dynamics_row!(mpc, row_ref, expr)

Set a single row of F, G, Gd from `expr`. `row_ref` must have `idx > 0`.
"""
function _set_dynamics_row!(mpc::MPC, row_ref::StateRef, expr)
    r = row_ref.idx
    r > 0 || error("_set_dynamics_row!: expected indexed StateRef (idx > 0)")
    m = mpc.model
    F_new = copy(m.F); G_new = copy(m.G); Gd_new = copy(m.Gd); fo_new = copy(m.f_offset)
    signals = expr isa DynamicsExpr ? expr.signals : expr
    if expr isa DynamicsExpr && !isempty(expr.offset)
        ofs = expr.offset
        fo_new[r] += length(ofs) >= r ? ofs[r] : ofs[end]
    end
    for t in signals.terms
        _apply_term_to_row!(F_new, G_new, Gd_new, t, r)
    end
    mpc.model = Model(F_new, G_new;
        Gd=Gd_new, C=m.C, Dd=m.Dd, f_offset=fo_new, h_offset=m.h_offset,
        Ts=m.Ts, xo=m.xo, uo=m.uo, wmin=m.wmin, wmax=m.wmax)
    mpc.mpqp_issetup = false
    return mpc
end

"""
    _set_dynamics_rows!(mpc, rows, expr)

Set multiple rows of F, G, Gd simultaneously from `expr`.
`rows` is a `Vector{Int}` of row indices.
"""
function _set_dynamics_rows!(mpc::MPC, rows::Vector{Int}, expr)
    m = mpc.model
    F_new = copy(m.F); G_new = copy(m.G); Gd_new = copy(m.Gd); fo_new = copy(m.f_offset)
    signals = expr isa DynamicsExpr ? expr.signals : expr
    if expr isa DynamicsExpr && !isempty(expr.offset)
        ofs = expr.offset
        for (i, r) in enumerate(rows)
            fo_new[r] += length(ofs) >= i ? ofs[i] : 0.0
        end
    end
    for t in signals.terms
        _apply_term_to_rows!(F_new, G_new, Gd_new, t, rows)
    end
    mpc.model = Model(F_new, G_new;
        Gd=Gd_new, C=m.C, Dd=m.Dd, f_offset=fo_new, h_offset=m.h_offset,
        Ts=m.Ts, xo=m.xo, uo=m.uo, wmin=m.wmin, wmax=m.wmax)
    mpc.mpqp_issetup = false
    return mpc
end

# Helper: extract row indices from various LHS types
_dynamics_rows(sv::StateVec)  = sv.indices
function _dynamics_rows(v::AbstractVector)
    all(x -> x isa StateRef && x.idx > 0, v) ||
        error("@dynamics: vector LHS must contain only indexed StateRefs")
    [x.idx for x in v]
end

# ---- Public macros ---------------------------------------------------------

"""
    @control(mpc, name)
    @control(mpc, name[1:n])
    @control(mpc, name, umin=val, umax=val)
    @control(mpc, name[1:n], umin=val, umax=val)

Declare control signal reference(s).

In **incremental mode** (`MPC(Np=...)`):
- `@control mpc u`      — adds 1 control; `u` is a `ControlRef`.
- `@control mpc u[1:3]` — adds 3 controls; `u` is a `ControlVec`.

In **style-1 mode** (`MPC(F,G;...)`):
- `@control mpc u`      — `u` is a full-vector `ControlRef`.
- `@control mpc u[1:3]` — `u` is a `ControlVec` over indices `[1,2,3]`.
"""
macro control(mpc_ex, name_ex, args...)
    umin_ex = nothing; umax_ex = nothing
    for arg in args
        (arg isa Expr && arg.head == :(=)) || continue
        k = Symbol(arg.args[1])
        k == :umin && (umin_ex = arg.args[2])
        k == :umax && (umax_ex = arg.args[2])
    end
    um = isnothing(umin_ex) ? :(zeros(0)) : esc(umin_ex)
    uM = isnothing(umax_ex) ? :(zeros(0)) : esc(umax_ex)
    apply_bounds = !isnothing(umin_ex) || !isnothing(umax_ex)

    if name_ex isa Expr && name_ex.head == :ref
        # Vector form: u[1:n]
        var_sym   = name_ex.args[1]
        range_ex  = name_ex.args[2]
        result = Expr(:block)
        push!(result.args, quote
            if LinearMPC._is_incremental($(esc(mpc_ex)))
                _n = length($(esc(range_ex)))
                _si = $(esc(mpc_ex)).model.nu + 1
                for _ in 1:_n; LinearMPC._expand_control!($(esc(mpc_ex))); end
                $(esc(var_sym)) = LinearMPC.ControlVec($(esc(mpc_ex)), collect(_si:_si+_n-1))
            else
                $(esc(var_sym)) = LinearMPC.ControlVec($(esc(mpc_ex)), collect($(esc(range_ex))))
            end
        end)
        apply_bounds && push!(result.args,
            :(LinearMPC.set_input_bounds!($(esc(mpc_ex)); umin=$um, umax=$uM)))
        push!(result.args, :($(esc(var_sym))))
        return result
    end

    # Scalar form: u
    var_sym = name_ex isa Symbol ? name_ex : error("@control: expected symbol or symbol[range]")
    result = Expr(:block)
    push!(result.args, quote
        if LinearMPC._is_incremental($(esc(mpc_ex)))
            _idx = LinearMPC._expand_control!($(esc(mpc_ex)))
            $(esc(var_sym)) = LinearMPC.ControlRef($(esc(mpc_ex)), _idx)
        else
            $(esc(var_sym)) = LinearMPC.ControlRef($(esc(mpc_ex)))
        end
    end)
    apply_bounds && push!(result.args,
        :(LinearMPC.set_input_bounds!($(esc(mpc_ex)); umin=$um, umax=$uM)))
    push!(result.args, :($(esc(var_sym))))
    return result
end

"""
    @state(mpc, name)
    @state(mpc, name[1:n])

Declare state signal reference(s).

In **incremental mode**:
- `@state mpc x1`      — adds 1 state; defines `x1` and `x1_next` (both `StateRef`).
- `@state mpc v[1:4]`  — adds 4 states; defines `v` and `v_next` (both `StateVec`).

In **style-1 mode**:
- `@state mpc x`       — `x` and `x_next` are full-vector `StateRef`s.
- `@state mpc v[1:4]`  — `v` and `v_next` are `StateVec`s over indices `[1,2,3,4]`.

The `_next` variant is provided for readability in `@dynamics` LHS, but the plain
name also works there.
"""
macro state(mpc_ex, name_ex)
    if name_ex isa Expr && name_ex.head == :ref
        var_sym  = name_ex.args[1]
        range_ex = name_ex.args[2]
        next_sym = Symbol(string(var_sym) * "_next")
        return quote
            if LinearMPC._is_incremental($(esc(mpc_ex)))
                _n = length($(esc(range_ex)))
                _si = $(esc(mpc_ex)).model.nx + 1
                for _ in 1:_n; LinearMPC._expand_state!($(esc(mpc_ex))); end
                $(esc(var_sym))  = LinearMPC.StateVec($(esc(mpc_ex)), collect(_si:_si+_n-1))
                $(esc(next_sym)) = LinearMPC.StateVec($(esc(mpc_ex)), collect(_si:_si+_n-1))
            else
                $(esc(var_sym))  = LinearMPC.StateVec($(esc(mpc_ex)), collect($(esc(range_ex))))
                $(esc(next_sym)) = LinearMPC.StateVec($(esc(mpc_ex)), collect($(esc(range_ex))))
            end
        end
    end

    name_ex isa Symbol || error("@state: expected symbol or symbol[range]")
    next_sym = Symbol(string(name_ex) * "_next")
    quote
        if LinearMPC._is_incremental($(esc(mpc_ex)))
            _idx = LinearMPC._expand_state!($(esc(mpc_ex)))
            $(esc(name_ex))  = LinearMPC.StateRef($(esc(mpc_ex)), _idx)
            $(esc(next_sym)) = LinearMPC.StateRef($(esc(mpc_ex)), _idx)
        else
            $(esc(name_ex))  = LinearMPC.StateRef($(esc(mpc_ex)))
            $(esc(next_sym)) = LinearMPC.StateRef($(esc(mpc_ex)))
        end
    end
end

"""
    @output(mpc, name)

Declare an output signal reference.
"""
macro output(mpc_ex, name_ex)
    quote; $(esc(name_ex)) = LinearMPC.OutputRef($(esc(mpc_ex))); end
end

"""
    @disturbance(mpc, name)
    @disturbance(mpc, name[1:n])
    @disturbance(mpc, name, wmin=val, wmax=val)

Declare disturbance signal reference(s). Like `@state` and `@control`, supports
both scalar (`d`) and vector (`d[1:2]`) forms in both style-1 and incremental mode.
"""
macro disturbance(mpc_ex, name_ex, args...)
    wmin_ex = nothing; wmax_ex = nothing
    for arg in args
        (arg isa Expr && arg.head == :(=)) || continue
        k = Symbol(arg.args[1])
        k == :wmin && (wmin_ex = arg.args[2])
        k == :wmax && (wmax_ex = arg.args[2])
    end
    wm = isnothing(wmin_ex) ? :(zeros(0)) : esc(wmin_ex)
    wM = isnothing(wmax_ex) ? :(zeros(0)) : esc(wmax_ex)
    apply_bounds = !isnothing(wmin_ex) || !isnothing(wmax_ex)

    if name_ex isa Expr && name_ex.head == :ref
        var_sym  = name_ex.args[1]
        range_ex = name_ex.args[2]
        next_sym = Symbol(string(var_sym) * "_next")
        result = Expr(:block)
        push!(result.args, quote
            if LinearMPC._is_incremental($(esc(mpc_ex)))
                _n = length($(esc(range_ex)))
                _si = $(esc(mpc_ex)).model.nd + 1
                for _ in 1:_n; LinearMPC._expand_disturbance!($(esc(mpc_ex))); end
                $(esc(var_sym))  = LinearMPC.DisturbanceVec($(esc(mpc_ex)), collect(_si:_si+_n-1))
                $(esc(next_sym)) = LinearMPC.DisturbanceVec($(esc(mpc_ex)), collect(_si:_si+_n-1))
            else
                $(esc(var_sym))  = LinearMPC.DisturbanceVec($(esc(mpc_ex)), collect($(esc(range_ex))))
                $(esc(next_sym)) = LinearMPC.DisturbanceVec($(esc(mpc_ex)), collect($(esc(range_ex))))
            end
        end)
        apply_bounds && push!(result.args,
            :(LinearMPC.set_disturbance!($(esc(mpc_ex)), $wm, $wM)))
        push!(result.args, :($(esc(var_sym))))
        return result
    end

    name_ex isa Symbol || error("@disturbance: expected symbol or symbol[range]")
    next_sym = Symbol(string(name_ex) * "_next")
    result = Expr(:block)
    push!(result.args, quote
        if LinearMPC._is_incremental($(esc(mpc_ex)))
            _idx = LinearMPC._expand_disturbance!($(esc(mpc_ex)))
            $(esc(name_ex))  = LinearMPC.DisturbanceRef($(esc(mpc_ex)), _idx)
            $(esc(next_sym)) = LinearMPC.DisturbanceRef($(esc(mpc_ex)), _idx)
        else
            $(esc(name_ex)) = LinearMPC.DisturbanceRef($(esc(mpc_ex)))
        end
    end)
    apply_bounds && push!(result.args,
        :(LinearMPC.set_disturbance!($(esc(mpc_ex)), $wm, $wM)))
    push!(result.args, :($(esc(name_ex))))
    return result
end

"""
    @constraint(mpc, lb <= signal <= ub)
    @constraint(mpc, signal <= ub)
    @constraint(mpc, signal >= lb)
    @constraint(mpc, expr, soft=true, ks=2:5, binary=false, prio=0)

Add a constraint. Works with `ControlRef`, `ControlVec`, `OutputRef`, `StateRef`,
`DisturbanceRef`, and `SignalExpr`.
"""
macro constraint(mpc_ex, expr, args...)
    kw_pairs = _parse_macro_kwargs(args)
    if expr isa Expr && expr.head == :comparison && length(expr.args) == 5
        lb_ex, op1, mid_ex, op2, ub_ex = expr.args
        if op1 == :(<=) && op2 == :(<=)
            return _make_call(:(LinearMPC._add_constraint!), kw_pairs,
                              mpc_ex, lb_ex, mid_ex, ub_ex)
        elseif op1 == :(>=) && op2 == :(>=)
            return _make_call(:(LinearMPC._add_constraint!), kw_pairs,
                              mpc_ex, ub_ex, mid_ex, lb_ex)
        end
    end
    if expr isa Expr && expr.head == :call && length(expr.args) == 3
        op, lhs_ex, rhs_ex = expr.args
        op in (:(<=), :(>=)) &&
            return _make_call(:(LinearMPC._add_constraint_onesided!), kw_pairs,
                              mpc_ex, lhs_ex, QuoteNode(op), rhs_ex)
    end
    return :(error("@constraint: unsupported expression `" * string($(QuoteNode(expr))) * "`"))
end

"""
    @dynamics(mpc, lhs = rhs)

Set (or update) system dynamics.

**LHS forms:**
- `x_next` or `x`       — single `StateRef` with `idx > 0`; sets one row of F, G.
- `v[3]`                — index into `StateVec`; sets the row for `v.indices[3]`.
- `v[1:2]` or `[x; z]` — `StateVec` or vector of `StateRef`s; sets multiple rows.
- `x_next` with `idx==0` — full-matrix assignment (style 1).

**RHS forms:**
- `F*x + G*u`          — style 1 full-matrix.
- `0.9*x1 - 0.2*x2 + 0.5*u` — scalar coefficient terms (style 2 single row).
- `A*v[1:3] + B*u`    — matrix times `StateVec`/`ControlVec` (style 2 multi-row).

Undefined rows of F, G remain zero.

# Examples
```julia
# Style 1
mpc = MPC(2, 1; Np=10);  @state mpc x;  @control mpc u
A = [1.0 0.1; 0.0 1.0];  B = [0.0; 1.0]
@dynamics mpc x_next = A*x + B*u

# Style 2 -- scalar
mpc = MPC(Np=10)
@state mpc x1;  @state mpc x2;  @control mpc u
@dynamics mpc x1 = 0.9*x1 - 0.2*x2 + 0.5*u   # row 1 (can use x1 or x1_next)
@dynamics mpc x2 = 0.1*x1 + 0.8*x2            # row 2

# Style 2 -- vector
@state mpc v[1:3];  @control mpc u
A3 = rand(3,3)
@dynamics mpc v = A3*v + 0.5*u[1]   # sets rows for all elements of v

# Style 2 -- vcat LHS
@state mpc x1;  @state mpc x2
@dynamics mpc [x1; x2] = A*[x1; x2]   # alternative multi-row form
```
"""
macro dynamics(mpc_ex, eq_ex)
    (eq_ex isa Expr && eq_ex.head == :(=)) ||
        error("@dynamics: expected assignment, e.g. `x1 = 0.9*x1 + 0.5*u`")
    lhs = eq_ex.args[1]
    rhs = eq_ex.args[2]
    return quote
        let _lhs_val = $(esc(lhs)), _rhs_val = $(esc(rhs))
            if _lhs_val isa LinearMPC.StateRef && _lhs_val.idx > 0
                # Single row (scalar StateRef, includes x1 and x1_next)
                LinearMPC._set_dynamics_row!($(esc(mpc_ex)), _lhs_val, _rhs_val)
            elseif _lhs_val isa LinearMPC.StateVec
                # Multiple rows via StateVec (e.g. v, v[1:2])
                LinearMPC._set_dynamics_rows!($(esc(mpc_ex)),
                    LinearMPC._dynamics_rows(_lhs_val), _rhs_val)
            elseif _lhs_val isa AbstractVector
                # Multiple rows via [x1; x2] vcat
                LinearMPC._set_dynamics_rows!($(esc(mpc_ex)),
                    LinearMPC._dynamics_rows(_lhs_val), _rhs_val)
            else
                # Full-matrix (style 1)
                LinearMPC._set_dynamics!($(esc(mpc_ex)), _rhs_val)
            end
        end
    end
end

"""
    @objective(mpc, Q=Q_val, R=R_val, Rr=Rr_val, S=S_val, Qf=Qf_val, Qfx=Qfx_val)

Set the objective function weights.
"""
macro objective(mpc_ex, args...)
    _make_call(:(LinearMPC.set_objective!), _parse_macro_kwargs(args), mpc_ex)
end

# ---- MPC constructor for incremental mode ----------------------------------

"""
    MPC(; Np=10, Nc=Np)

Create a zero-dimensional MPC in **incremental mode**. States, controls, and
disturbances are added one at a time (or in batches) with [`@state`](@ref),
[`@control`](@ref), and [`@disturbance`](@ref). Dynamics are set row-by-row
with [`@dynamics`](@ref); undefined rows default to zero.

# Example
```julia
mpc = MPC(Np=10)
@state mpc x1;  @state mpc x2;  @control mpc u

@dynamics mpc x1 = 0.9*x1 - 0.2*x2 + 0.5*u
@dynamics mpc x2 = 0.1*x1 + 0.8*x2

@objective mpc Q=I R=0.1
@constraint mpc -1 <= u <= 1
setup!(mpc)
```
"""
function MPC(; Np::Int=10, Nc::Int=Np)
    mpc = MPC(Model(zeros(0,0), zeros(0,0)); Np=Np, Nc=Nc)
    _incremental_registry[mpc] = (0, 0, 0)
    return mpc
end
