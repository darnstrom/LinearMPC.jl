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
#   mpc = MPC(Np=10)          # zero-dimensional MPC
#   @state   mpc x1           # x1 = StateRef(mpc,1), x1_next also defined
#   @state   mpc x2           # x2 = StateRef(mpc,2), x2_next also defined
#   @control mpc u1           # u1 = ControlRef(mpc,1)
#   @dynamics mpc x1_next = 0.9*x1 - 0.2*x2 + 0.5*u1   # sets row 1 of F,G
#   @dynamics mpc x2_next = 0.1*x1 + 0.8*x2             # row 2; G[2,:] stays 0
#   @objective mpc Q=I R=0.1
#   @constraint mpc -1 <= u1 <= 1
#   setup!(mpc)

# ---- Module-level registry for incremental-mode MPCs ----------------------
# Maps each incremental MPC (by object identity) to its (nx_count, nu_count, nd_count).
const _incremental_registry = WeakKeyDict{MPC, NTuple{3,Int}}()

"""
Return true if `mpc` was created in incremental mode (via `MPC(; Np=...)`).
"""
_is_incremental(mpc::MPC) = haskey(_incremental_registry, mpc)

# ---- Signal reference types ------------------------------------------------

"""
    ControlRef(mpc, idx)
    ControlRef(mpc)

Reference to the control input signal of an MPC controller.

- `idx = 0` -- refers to the full control vector (style 1).
- `idx > 0` -- refers to a single control variable at column `idx` (style 2).

Created by the [`@control`](@ref) macro.
"""
struct ControlRef
    mpc::MPC
    idx::Int
end
ControlRef(mpc::MPC) = ControlRef(mpc, 0)

"""
    StateRef(mpc, idx)
    StateRef(mpc)

Reference to the state signal of an MPC controller.

- `idx = 0` -- refers to the full state vector (style 1).
- `idx > 0` -- refers to a single state variable at row/column `idx` (style 2).

Created by the [`@state`](@ref) macro.
"""
struct StateRef
    mpc::MPC
    idx::Int
end
StateRef(mpc::MPC) = StateRef(mpc, 0)

"""
    OutputRef(mpc)

Reference to the output signal of an MPC controller.
Created by the [`@output`](@ref) macro for use in [`@constraint`](@ref) expressions.
"""
struct OutputRef
    mpc::MPC
end

"""
    DisturbanceRef(mpc, idx)
    DisturbanceRef(mpc)

Reference to the disturbance signal of an MPC controller.

- `idx = 0` -- refers to the full disturbance vector (style 1).
- `idx > 0` -- refers to a single disturbance variable at column `idx` (style 2).

Created by the [`@disturbance`](@ref) macro.
"""
struct DisturbanceRef
    mpc::MPC
    idx::Int
end
DisturbanceRef(mpc::MPC) = DisturbanceRef(mpc, 0)

# ---- SignalTerm and SignalExpr ---------------------------------------------

"""
Internal: A * [signal_kind] term, possibly at a specific column.

- `col_idx == 0`: dense matrix coefficient (full vector signal).
- `col_idx > 0`: scalar coefficient at column `col_idx` (single signal variable).
"""
struct SignalTerm
    A::Matrix{Float64}
    kind::Symbol       # :u, :x, :r, :d, :uprev
    col_idx::Int       # 0 = dense full matrix; i > 0 = scalar at column i
end
SignalTerm(A, kind) = SignalTerm(A, kind, 0)

"""
    SignalExpr(mpc, terms)

A linear combination of weighted MPC signals (e.g. `Au*u + Ax*x`).
Built automatically when multiplying a matrix by a signal reference.
"""
struct SignalExpr
    mpc::MPC
    terms::Vector{SignalTerm}
end

# ---- DynamicsExpr: SignalExpr + optional constant offset -------------------

"""
    DynamicsExpr(signals, offset)

Right-hand side of a dynamics equation: F*x + G*u + Gd*d + offset.
"""
struct DynamicsExpr
    signals::SignalExpr
    offset::Vector{Float64}
end
DynamicsExpr(expr::SignalExpr) = DynamicsExpr(expr, zeros(0))

# ---- Arithmetic overloads --------------------------------------------------

# Dense matrix * full-state/control/disturbance ref (style 1)
Base.:*(A::AbstractMatrix, x::StateRef) =
    SignalExpr(x.mpc, [SignalTerm(float(A), :x, 0)])
Base.:*(A::AbstractMatrix, u::ControlRef) =
    SignalExpr(u.mpc, [SignalTerm(float(A), :u, 0)])
Base.:*(A::AbstractMatrix, d::DisturbanceRef) =
    SignalExpr(d.mpc, [SignalTerm(float(A), :d, 0)])

# Column-vector * full-state/control/disturbance ref (style 1, row-vector coefficient)
Base.:*(A::AbstractVector, x::StateRef) =
    SignalExpr(x.mpc, [SignalTerm(reshape(float(A), :, 1), :x, 0)])
Base.:*(A::AbstractVector, u::ControlRef) =
    SignalExpr(u.mpc, [SignalTerm(reshape(float(A), :, 1), :u, 0)])
Base.:*(A::AbstractVector, d::DisturbanceRef) =
    SignalExpr(d.mpc, [SignalTerm(reshape(float(A), :, 1), :d, 0)])

# Scalar * ref: for indexed refs (style 2) -> sparse column term;
#               for full-vector refs (style 1) -> scaled identity
function Base.:*(a::Number, x::StateRef)
    if x.idx == 0
        nx = x.mpc.model.nx
        return SignalExpr(x.mpc, [SignalTerm(float(a) * Matrix{Float64}(I, nx, nx), :x, 0)])
    else
        return SignalExpr(x.mpc, [SignalTerm(fill(Float64(a), 1, 1), :x, x.idx)])
    end
end
function Base.:*(a::Number, u::ControlRef)
    if u.idx == 0
        nu = u.mpc.model.nu
        return SignalExpr(u.mpc, [SignalTerm(float(a) * Matrix{Float64}(I, nu, nu), :u, 0)])
    else
        return SignalExpr(u.mpc, [SignalTerm(fill(Float64(a), 1, 1), :u, u.idx)])
    end
end
function Base.:*(a::Number, d::DisturbanceRef)
    if d.idx == 0
        nd = d.mpc.model.nd
        return SignalExpr(d.mpc, [SignalTerm(float(a) * Matrix{Float64}(I, nd, nd), :d, 0)])
    else
        return SignalExpr(d.mpc, [SignalTerm(fill(Float64(a), 1, 1), :d, d.idx)])
    end
end

# Bare signal as identity (for use in `x + G*u` etc.)
_bare(x::StateRef) = x.idx == 0 ?
    SignalExpr(x.mpc, [SignalTerm(Matrix{Float64}(I, x.mpc.model.nx, x.mpc.model.nx), :x, 0)]) :
    SignalExpr(x.mpc, [SignalTerm(fill(1.0, 1, 1), :x, x.idx)])
_bare(u::ControlRef) = u.idx == 0 ?
    SignalExpr(u.mpc, [SignalTerm(Matrix{Float64}(I, u.mpc.model.nu, u.mpc.model.nu), :u, 0)]) :
    SignalExpr(u.mpc, [SignalTerm(fill(1.0, 1, 1), :u, u.idx)])
_bare(d::DisturbanceRef) = d.idx == 0 ?
    SignalExpr(d.mpc, [SignalTerm(Matrix{Float64}(I, d.mpc.model.nd, d.mpc.model.nd), :d, 0)]) :
    SignalExpr(d.mpc, [SignalTerm(fill(1.0, 1, 1), :d, d.idx)])

# Combine two SignalExprs
function Base.:+(e1::SignalExpr, e2::SignalExpr)
    @assert e1.mpc === e2.mpc "Signal expressions must refer to the same MPC controller"
    SignalExpr(e1.mpc, [e1.terms; e2.terms])
end

# Unary negation and subtraction for SignalExpr
Base.:-(e::SignalExpr) =
    SignalExpr(e.mpc, [SignalTerm(-term.A, term.kind, term.col_idx) for term in e.terms])
Base.:-(e1::SignalExpr, e2::SignalExpr) = e1 + (-e2)
Base.:-(e::SignalExpr, x::StateRef)     = e + (-_bare(x))
Base.:-(e::SignalExpr, u::ControlRef)   = e + (-_bare(u))
Base.:-(e::SignalExpr, d::DisturbanceRef) = e + (-_bare(d))
Base.:-(x::StateRef,   e::SignalExpr)   = _bare(x) + (-e)
Base.:-(u::ControlRef, e::SignalExpr)   = _bare(u) + (-e)

# Allow bare signal + something: x + G*u, u + F*x, etc.
Base.:+(x::StateRef, e::SignalExpr) = _bare(x) + e
Base.:+(e::SignalExpr, x::StateRef) = e + _bare(x)
Base.:+(u::ControlRef, e::SignalExpr) = _bare(u) + e
Base.:+(e::SignalExpr, u::ControlRef) = e + _bare(u)
Base.:+(d::DisturbanceRef, e::SignalExpr) = _bare(d) + e
Base.:+(e::SignalExpr, d::DisturbanceRef) = e + _bare(d)
Base.:+(x::StateRef, u::ControlRef) = _bare(x) + _bare(u)
Base.:+(u::ControlRef, x::StateRef) = _bare(u) + _bare(x)
Base.:+(x::StateRef, d::DisturbanceRef) = _bare(x) + _bare(d)
Base.:+(d::DisturbanceRef, x::StateRef) = _bare(d) + _bare(x)

# DynamicsExpr: SignalExpr + constant offset vector
Base.:+(e::SignalExpr, offset::AbstractVector) = DynamicsExpr(e, float(offset))
Base.:+(offset::AbstractVector, e::SignalExpr) = DynamicsExpr(e, float(offset))
Base.:+(e::DynamicsExpr, offset::AbstractVector) = DynamicsExpr(e.signals, e.offset + offset)
Base.:+(offset::AbstractVector, e::DynamicsExpr) = DynamicsExpr(e.signals, e.offset + offset)
Base.:+(e1::DynamicsExpr, e2::SignalExpr) = DynamicsExpr(e1.signals + e2, e1.offset)
Base.:+(e1::SignalExpr, e2::DynamicsExpr) = DynamicsExpr(e1 + e2.signals, e2.offset)

# ---- Model expansion helpers (incremental mode) ---------------------------

"""
Expand mpc.model to add one more state dimension. Returns the new state index.
The new state has zero rows/columns in F, G, Gd and is unobserved in C by default.
"""
function _expand_state!(mpc::MPC)
    m = mpc.model
    nx, nu, nd = m.nx, m.nu, m.nd
    new_nx = nx + 1
    F_new  = [m.F          zeros(nx, 1);
              zeros(1, nx) zeros(1, 1)]
    G_new  = [m.G;  zeros(1, nu)]
    Gd_new = [m.Gd; zeros(1, nd)]
    # Default to full-state output (C = I) in incremental mode
    C_new  = Matrix{Float64}(I, new_nx, new_nx)
    mpc.model = Model(F_new, G_new;
        Gd       = Gd_new,
        C        = C_new,
        Dd       = m.Dd,
        f_offset = [m.f_offset; 0.0],
        h_offset = zeros(new_nx),   # h_offset grows with ny
        xo       = [m.xo; 0.0],
        uo       = m.uo,
        wmin     = [m.wmin; 0.0],
        wmax     = [m.wmax; 0.0],
        Ts       = m.Ts)
    mpc.mpqp_issetup = false
    nx_s, nu_s, nd_s = _incremental_registry[mpc]
    _incremental_registry[mpc] = (nx_s + 1, nu_s, nd_s)
    return new_nx   # new state index
end

"""
Expand mpc.model to add one more control dimension. Returns the new control index.
"""
function _expand_control!(mpc::MPC)
    m = mpc.model
    G_new = [m.G zeros(m.nx, 1)]
    mpc.model = Model(m.F, G_new;
        Gd       = m.Gd,
        C        = m.C,
        Dd       = m.Dd,
        f_offset = m.f_offset,
        h_offset = m.h_offset,
        xo       = m.xo,
        uo       = [m.uo; 0.0],
        wmin     = m.wmin,
        wmax     = m.wmax,
        Ts       = m.Ts)
    mpc.mpqp_issetup = false
    nx_s, nu_s, nd_s = _incremental_registry[mpc]
    _incremental_registry[mpc] = (nx_s, nu_s + 1, nd_s)
    return m.nu + 1   # new control index
end

"""
Expand mpc.model to add one more disturbance dimension. Returns the new disturbance index.
"""
function _expand_disturbance!(mpc::MPC)
    m = mpc.model
    Gd_new = [m.Gd zeros(m.nx, 1)]
    Dd_new = [m.Dd zeros(m.ny, 1)]
    mpc.model = Model(m.F, m.G;
        Gd       = Gd_new,
        C        = m.C,
        Dd       = Dd_new,
        f_offset = m.f_offset,
        h_offset = m.h_offset,
        xo       = m.xo,
        uo       = m.uo,
        wmin     = m.wmin,
        wmax     = m.wmax,
        Ts       = m.Ts)
    mpc.mpqp_issetup = false
    nx_s, nu_s, nd_s = _incremental_registry[mpc]
    _incremental_registry[mpc] = (nx_s, nu_s, nd_s + 1)
    return m.nd + 1   # new disturbance index
end

# ---- Internal helpers ------------------------------------------------------

_bound(::Nothing) = zeros(0)
_bound(x) = x

_to_vec_bound(::Nothing, n)           = zeros(0)
_to_vec_bound(x::AbstractVector, n)   = x
_to_vec_bound(x::Number, n)           = fill(Float64(x), n)

_is_signal(::ControlRef)     = true
_is_signal(::OutputRef)      = true
_is_signal(::StateRef)       = true
_is_signal(::DisturbanceRef) = true
_is_signal(::SignalExpr)     = true
_is_signal(::Any)            = false

# ---- Core constraint dispatch ----------------------------------------------

function _add_constraint!(mpc::MPC, lb, ::ControlRef, ub; kw...)
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
        Ax = Matrix{Float64}(I, mpc.model.nx, mpc.model.nx),
        lb = _to_vec_bound(lb, mpc.model.nx),
        ub = _to_vec_bound(ub, mpc.model.nx),
        ks=ks, soft=soft, binary=binary, prio=prio)
end

function _add_constraint!(mpc::MPC, lb, ::DisturbanceRef, ub; kw...)
    set_disturbance!(mpc, _bound(lb), _bound(ub))
end

function _add_constraint!(mpc::MPC, lb, expr::SignalExpr, ub;
        ks=2:mpc.Np, soft=false, binary=false, prio=0, kw...)
    Au = zeros(0,0); Ax = zeros(0,0)
    Ar = zeros(0,0); Ad = zeros(0,0); Aup = zeros(0,0)
    for term in expr.terms
        A = term.A
        if     term.kind == :u;     Au  = isempty(Au)  ? A : Au  + A
        elseif term.kind == :x;     Ax  = isempty(Ax)  ? A : Ax  + A
        elseif term.kind == :r;     Ar  = isempty(Ar)  ? A : Ar  + A
        elseif term.kind == :d;     Ad  = isempty(Ad)  ? A : Ad  + A
        elseif term.kind == :uprev; Aup = isempty(Aup) ? A : Aup + A
        end
    end
    add_constraint!(mpc;
        Au  = isempty(Au)  ? nothing : Au,
        Ax  = isempty(Ax)  ? nothing : Ax,
        Ar  = Ar, Ad = Ad, Aup = Aup,
        lb  = _bound(lb), ub = _bound(ub),
        ks=ks, soft=soft, binary=binary, prio=prio)
end

function _add_constraint_onesided!(mpc::MPC, lhs, op::Symbol, rhs; kw...)
    if op == :(<=)
        if _is_signal(lhs)
            _add_constraint!(mpc, nothing, lhs, rhs; kw...)
        elseif _is_signal(rhs)
            _add_constraint!(mpc, lhs, rhs, nothing; kw...)
        else
            error("@constraint: no signal reference found in expression")
        end
    elseif op == :(>=)
        if _is_signal(lhs)
            _add_constraint!(mpc, rhs, lhs, nothing; kw...)
        elseif _is_signal(rhs)
            _add_constraint!(mpc, nothing, rhs, lhs; kw...)
        else
            error("@constraint: no signal reference found in expression")
        end
    else
        error("@constraint: unsupported operator $op")
    end
end

# ---- Macro helper ----------------------------------------------------------

function _parse_macro_kwargs(args)
    kw_pairs = Expr[]
    for arg in args
        if arg isa Expr && arg.head == :(=)
            push!(kw_pairs, Expr(:kw, arg.args[1], esc(arg.args[2])))
        end
    end
    return kw_pairs
end

function _make_call(fn, kw_pairs, pos_args...)
    escaped = [esc(a) for a in pos_args]
    if isempty(kw_pairs)
        return Expr(:call, fn, escaped...)
    else
        return Expr(:call, fn, Expr(:parameters, kw_pairs...), escaped...)
    end
end

# ---- Dynamics: full-matrix and row-by-row ----------------------------------

"""
    _set_dynamics!(mpc, expr)

Set the full F, G, Gd, f_offset from a `SignalExpr` or `DynamicsExpr`.
"""
function _set_dynamics!(mpc::MPC, expr::SignalExpr)
    _set_dynamics!(mpc, DynamicsExpr(expr))
end

function _set_dynamics!(mpc::MPC, dexpr::DynamicsExpr)
    model = mpc.model
    nx, nu, nd = model.nx, model.nu, model.nd

    F  = zeros(nx, nx)
    G  = zeros(nx, nu)
    Gd = zeros(nx, nd)
    f_offset = isempty(dexpr.offset) ? zeros(nx) : Vector{Float64}(dexpr.offset)

    for term in dexpr.signals.terms
        A = term.A
        if     term.kind == :x;  F  .+= A
        elseif term.kind == :u;  G  .+= A
        elseif term.kind == :d;  Gd .+= A
        end
    end

    mpc.model = Model(F, G;
        Gd       = Gd,
        C        = model.C,
        Dd       = model.Dd,
        f_offset = f_offset,
        h_offset = model.h_offset,
        Ts       = model.Ts,
        xo       = model.xo,
        uo       = model.uo,
        wmin     = model.wmin,
        wmax     = model.wmax)
    mpc.mpqp_issetup = false
    return mpc
end

"""
    _set_dynamics_row!(mpc, row_ref, expr)

Set a single row of F, G, and Gd.  `row_ref` must be a `StateRef` with
`idx > 0` (created by `@state` in incremental mode).

Each `SignalTerm` in `expr` contributes:
- `col_idx == 0` (dense A): if A has one row, that row is used; otherwise row `r` of A.
- `col_idx > 0` (sparse): places scalar `A[1,1]` at the given column.

Rows not specified via this function remain zero.
"""
function _set_dynamics_row!(mpc::MPC, row_ref::StateRef, expr)
    r = row_ref.idx
    r > 0 || error("_set_dynamics_row!: expected an indexed StateRef (idx > 0), got idx=0")
    m = mpc.model

    F_new  = copy(m.F)
    G_new  = copy(m.G)
    Gd_new = copy(m.Gd)
    fo_new = copy(m.f_offset)

    signals = expr isa DynamicsExpr ? expr.signals : expr
    if expr isa DynamicsExpr && !isempty(expr.offset)
        ofs = expr.offset
        fo_new[r] += length(ofs) >= r ? ofs[r] : ofs[end]
    end

    for term in signals.terms
        A, col = term.A, term.col_idx
        if term.kind == :x
            if col == 0
                row_vec = size(A, 1) == 1 ? vec(A) : A[r, :]
                F_new[r, :] .+= row_vec
            else
                F_new[r, col] += A[1, 1]
            end
        elseif term.kind == :u
            if col == 0
                row_vec = size(A, 1) == 1 ? vec(A) : A[r, :]
                G_new[r, :] .+= row_vec
            else
                G_new[r, col] += A[1, 1]
            end
        elseif term.kind == :d
            if col == 0
                row_vec = size(A, 1) == 1 ? vec(A) : A[r, :]
                Gd_new[r, :] .+= row_vec
            else
                Gd_new[r, col] += A[1, 1]
            end
        end
    end

    mpc.model = Model(F_new, G_new;
        Gd       = Gd_new,
        C        = m.C,
        Dd       = m.Dd,
        f_offset = fo_new,
        h_offset = m.h_offset,
        Ts       = m.Ts,
        xo       = m.xo,
        uo       = m.uo,
        wmin     = m.wmin,
        wmax     = m.wmax)
    mpc.mpqp_issetup = false
    return mpc
end

# ---- Public macros ---------------------------------------------------------

"""
    @control(mpc, name)
    @control(mpc, name, umin=umin_val, umax=umax_val)

Declare a control input signal reference.

**Style 1** (`MPC(F, G; ...)`): assigns a full-vector `ControlRef` with no model change.

**Style 2** (`MPC(Np=...)`): increments `nu` by 1 and assigns an indexed `ControlRef`.

Optionally set input bounds via `umin`/`umax`.

# Examples
```julia
# Style 1
mpc = MPC(F, G; Np=10);  @control mpc u umin=-1 umax=1
# Style 2
mpc = MPC(Np=10);  @state mpc x1;  @control mpc u
@dynamics mpc x1_next = 0.9*x1 + 0.5*u
```
"""
macro control(mpc_ex, name_ex, args...)
    umin_ex = nothing; umax_ex = nothing
    for arg in args
        (arg isa Expr && arg.head == :(=)) || continue
        k = Symbol(arg.args[1])
        k == :umin && (umin_ex = arg.args[2])
        k == :umax && (umax_ex = arg.args[2])
    end

    result = Expr(:block)
    push!(result.args, quote
        if LinearMPC._is_incremental($(esc(mpc_ex)))
            _idx = LinearMPC._expand_control!($(esc(mpc_ex)))
            $(esc(name_ex)) = LinearMPC.ControlRef($(esc(mpc_ex)), _idx)
        else
            $(esc(name_ex)) = LinearMPC.ControlRef($(esc(mpc_ex)))
        end
    end)
    if !isnothing(umin_ex) || !isnothing(umax_ex)
        um = isnothing(umin_ex) ? :(zeros(0)) : esc(umin_ex)
        uM = isnothing(umax_ex) ? :(zeros(0)) : esc(umax_ex)
        push!(result.args,
              :(LinearMPC.set_input_bounds!($(esc(mpc_ex)); umin=$um, umax=$uM)))
    end
    push!(result.args, :($(esc(name_ex))))
    return result
end

"""
    @state(mpc, name)

Declare a state signal reference.

**Style 1** (`MPC(F, G; ...)`): assigns a full-vector `StateRef` with no model change.

**Style 2** (`MPC(Np=...)`): increments `nx` by 1 and assigns an indexed `StateRef`.
Also defines `name_next` (e.g., `x1_next`) for use as the LHS of [`@dynamics`](@ref).
Unspecified rows of F, G remain zero.

# Examples
```julia
# Style 1
mpc = MPC(F, G; Np=10);  @state mpc x;  @constraint mpc -5 <= x <= 5
# Style 2
mpc = MPC(Np=10)
@state mpc x1   # defines x1 (StateRef index 1) and x1_next
@state mpc x2   # defines x2 (StateRef index 2) and x2_next
@control mpc u
@dynamics mpc x1_next = 0.9*x1 - 0.2*x2 + 0.5*u
@dynamics mpc x2_next = 0.1*x1 + 0.8*x2   # row 2 of G stays zero
```
"""
macro state(mpc_ex, name_ex)
    name_ex isa Symbol || error("@state: expected a symbol for the state name")
    next_sym = Symbol(string(name_ex) * "_next")
    quote
        if LinearMPC._is_incremental($(esc(mpc_ex)))
            _idx = LinearMPC._expand_state!($(esc(mpc_ex)))
            $(esc(name_ex))  = LinearMPC.StateRef($(esc(mpc_ex)), _idx)
            $(esc(next_sym)) = LinearMPC.StateRef($(esc(mpc_ex)), _idx)
        else
            # Style 1: no model change; define both name and name_next as full-state refs
            $(esc(name_ex))  = LinearMPC.StateRef($(esc(mpc_ex)))
            $(esc(next_sym)) = LinearMPC.StateRef($(esc(mpc_ex)))
        end
    end
end

"""
    @output(mpc, name)

Declare an output signal reference for use in [`@constraint`](@ref) expressions.

# Example
```julia
@output mpc y
@constraint mpc 0 <= y <= 5
```
"""
macro output(mpc_ex, name_ex)
    quote
        $(esc(name_ex)) = LinearMPC.OutputRef($(esc(mpc_ex)))
    end
end

"""
    @disturbance(mpc, name)
    @disturbance(mpc, name, wmin=wmin_val, wmax=wmax_val)

Declare a disturbance signal reference.

**Style 2** (`MPC(Np=...)`): increments `nd` by 1 and assigns an indexed
`DisturbanceRef`. Also defines `name_next` for use in [`@dynamics`](@ref).

# Example
```julia
@disturbance mpc d wmin=-0.1*ones(1) wmax=0.1*ones(1)
```
"""
macro disturbance(mpc_ex, name_ex, args...)
    name_ex isa Symbol || error("@disturbance: expected a symbol for the disturbance name")
    next_sym = Symbol(string(name_ex) * "_next")
    wmin_ex = nothing; wmax_ex = nothing
    for arg in args
        (arg isa Expr && arg.head == :(=)) || continue
        k = Symbol(arg.args[1])
        k == :wmin && (wmin_ex = arg.args[2])
        k == :wmax && (wmax_ex = arg.args[2])
    end

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
    if !isnothing(wmin_ex) || !isnothing(wmax_ex)
        wm = isnothing(wmin_ex) ? :(zeros(0)) : esc(wmin_ex)
        wM = isnothing(wmax_ex) ? :(zeros(0)) : esc(wmax_ex)
        push!(result.args,
              :(LinearMPC.set_disturbance!($(esc(mpc_ex)), $wm, $wM)))
    end
    push!(result.args, :($(esc(name_ex))))
    return result
end

"""
    @constraint(mpc, lb <= signal <= ub)
    @constraint(mpc, signal <= ub)
    @constraint(mpc, signal >= lb)
    @constraint(mpc, lb <= A*signal1 + B*signal2 <= ub)
    @constraint(mpc, expr, soft=true, ks=2:5, binary=false, prio=0)

Add a constraint to the MPC controller `mpc`.

| Signal type  | Calls                      | Default `soft` |
|:-------------|:---------------------------|:---------------|
| `ControlRef` | `set_input_bounds!`        | N/A            |
| `OutputRef`  | `set_output_bounds!`       | `true`         |
| `StateRef`   | `add_constraint!` (Ax = I) | `false`        |
| `SignalExpr` | `add_constraint!`          | `false`        |

# Examples
```julia
@constraint mpc -1 <= u <= 1
@constraint mpc 0 <= y <= 5
@constraint mpc u <= 0.5
@constraint mpc lb <= Au*u + Ax*x <= ub
@constraint mpc -1 <= u <= 1 soft=true ks=3:8
```
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
        if op in (:(<=), :(>=))
            return _make_call(:(LinearMPC._add_constraint_onesided!), kw_pairs,
                              mpc_ex, lhs_ex, QuoteNode(op), rhs_ex)
        end
    end

    return :(error("@constraint: unsupported expression `" * string($(QuoteNode(expr))) * "`\n" *
                   "Expected: lb <= signal <= ub, signal <= ub, or signal >= lb"))
end

"""
    @dynamics(mpc, x_next = F*x + G*u)
    @dynamics(mpc, x_next = F*x + G*u + Gd*d + f_offset)
    @dynamics(mpc, x1_next = 0.9*x1 - 0.2*x2 + 0.5*u1)

Set (or update) system dynamics.

**Style 1** (full-matrix): the LHS `x_next` is a `StateRef` with `idx == 0`
(created by `@state mpc x` on a dimensioned MPC).  The RHS provides F, G, Gd,
and an optional constant offset vector all at once.

**Style 2** (row-by-row): the LHS `x1_next` is an indexed `StateRef` (created
by `@state mpc x1` on an incremental `MPC(Np=...)`).  The RHS expression
defines *one row* of F and G using scalar or row-vector coefficients.
Calling `@dynamics` for each state sets the corresponding row; rows not
explicitly set remain zero.

The C, Dd, h_offset, Ts, and operating-point settings of the existing model
are always preserved.

# Examples
```julia
# Style 1 -- full matrices
mpc = MPC(2, 1; Np=10)
@state mpc x;  @control mpc u
A = [1.0 0.1; 0.0 1.0];  B = [0.0; 1.0]
@dynamics mpc x_next = A*x + B*u

# Style 2 -- row by row (incremental)
mpc = MPC(Np=10)
@state mpc x1;  @state mpc x2;  @control mpc u
@dynamics mpc x1_next = 0.9*x1 - 0.2*x2 + 0.5*u   # row 1 of F and G
@dynamics mpc x2_next = 0.1*x1 + 0.8*x2            # row 2; G[2,:] stays 0

# With disturbance and affine offset (style 2)
@disturbance mpc d
@dynamics mpc x1_next = 0.9*x1 + 0.5*u + 0.1*d
```
"""
macro dynamics(mpc_ex, eq_ex)
    (eq_ex isa Expr && eq_ex.head == :(=)) ||
        error("@dynamics: expected an assignment expression, e.g. `x1_next = 0.9*x1 + 0.5*u`")
    lhs = eq_ex.args[1]
    rhs = eq_ex.args[2]
    return quote
        let _lhs_val = $(esc(lhs)), _rhs_val = $(esc(rhs))
            if _lhs_val isa LinearMPC.StateRef && _lhs_val.idx > 0
                LinearMPC._set_dynamics_row!($(esc(mpc_ex)), _lhs_val, _rhs_val)
            else
                LinearMPC._set_dynamics!($(esc(mpc_ex)), _rhs_val)
            end
        end
    end
end

"""
    @objective(mpc, Q=Q_val, R=R_val, Rr=Rr_val, S=S_val, Qf=Qf_val, Qfx=Qfx_val)

Set the objective function weights for the MPC controller.
Equivalent to [`set_objective!`](@ref)`(mpc; Q=Q_val, R=R_val, ...)`.

# Example
```julia
using LinearAlgebra
@objective mpc Q=I R=0.1*I Rr=0.01*I
```
"""
macro objective(mpc_ex, args...)
    kw_pairs = _parse_macro_kwargs(args)
    return _make_call(:(LinearMPC.set_objective!), kw_pairs, mpc_ex)
end

# ---- MPC constructor for incremental mode ----------------------------------

"""
    MPC(; Np=10, Nc=Np)

Create a zero-dimensional MPC controller in **incremental mode**.

States, controls, and disturbances are added one at a time using the
[`@state`](@ref), [`@control`](@ref), and [`@disturbance`](@ref) macros.
Each `@state mpc x_i` call increments the state dimension by 1 and defines
both `x_i` (for RHS use) and `x_i_next` (for the LHS of [`@dynamics`](@ref)).

Dynamics are specified row-by-row via [`@dynamics`](@ref); undefined rows
default to zero.

# Example
```julia
mpc = MPC(Np=10)
@state mpc x1;  @state mpc x2
@control mpc u

@dynamics mpc x1_next = 0.9*x1 - 0.2*x2 + 0.5*u
@dynamics mpc x2_next = 0.1*x1 + 0.8*x2   # G[2,:] stays zero

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
