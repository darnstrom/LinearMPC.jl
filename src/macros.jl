# ============================================================================
# JuMP-like macro interface for LinearMPC
# ============================================================================
#
# Usage example:
#
#   mpc = MPC(2, 1; Np=10)       # create with just state/control dimensions
#   x = @state   mpc x           # declare signal references
#   u = @control mpc u umin=-ones(1) umax=ones(1)
#   d = @disturbance mpc d
#
#   A = [1.0 0.1; 0.0 1.0]; B = [0.0; 1.0]
#   @dynamics mpc x_next = A*x + B*u        # set dynamics
#
#   y = @output mpc y
#   @objective mpc Q=I R=0.1*I Rr=0.01*I
#   @constraint mpc -1 <= u <= 1
#   @constraint mpc 0 <= y <= 5
#   @constraint mpc lb <= Au*u + Ax*x <= ub   # general linear constraint
#   @constraint mpc u <= 0.5
#   @constraint mpc -1 <= u <= 1 soft=true
#   setup!(mpc)

# ---- Signal reference types ------------------------------------------------

"""
    ControlRef(mpc)

Reference to the control input signal of an MPC controller.
Created by the [`@control`](@ref) macro for use in [`@constraint`](@ref) expressions.
"""
struct ControlRef
    mpc::MPC
end

"""
    StateRef(mpc)

Reference to the state signal of an MPC controller.
Created by the [`@state`](@ref) macro for use in [`@constraint`](@ref) expressions.
"""
struct StateRef
    mpc::MPC
end

"""
    OutputRef(mpc)

Reference to the output signal of an MPC controller.
Created by the [`@output`](@ref) macro for use in [`@constraint`](@ref) expressions.
"""
struct OutputRef
    mpc::MPC
end

"""
    DisturbanceRef(mpc)

Reference to the disturbance signal of an MPC controller.
Created by the [`@disturbance`](@ref) macro for use in [`@constraint`](@ref) expressions.
"""
struct DisturbanceRef
    mpc::MPC
end

# ---- SignalExpr: weighted combination of signals ---------------------------

# Internal: A * [signal_kind] term
struct SignalTerm
    A::Matrix{Float64}
    kind::Symbol  # :u, :x, :r, :d, :uprev
end

"""
    SignalExpr(mpc, terms)

A linear combination of weighted MPC signals (e.g. `Au*u + Ax*x`),
used inside [`@constraint`](@ref) expressions.

Built automatically when multiplying a matrix by a signal reference:
```julia
Au*u + Ax*x   # where u::ControlRef, x::StateRef
```
"""
struct SignalExpr
    mpc::MPC
    terms::Vector{SignalTerm}
end

# ---- Arithmetic overloads --------------------------------------------------

"""
    A * u  where u::ControlRef

Produce a weighted control signal expression for use in `@constraint` and `@dynamics`.
"""
Base.:*(A::AbstractMatrix, u::ControlRef) =
    SignalExpr(u.mpc, [SignalTerm(float(A), :u)])
Base.:*(A::AbstractVector, u::ControlRef) =
    SignalExpr(u.mpc, [SignalTerm(reshape(float(A), :, 1), :u)])

"""
    A * x  where x::StateRef

Produce a weighted state signal expression for use in `@constraint` and `@dynamics`.
"""
Base.:*(A::AbstractMatrix, x::StateRef) =
    SignalExpr(x.mpc, [SignalTerm(float(A), :x)])
Base.:*(A::AbstractVector, x::StateRef) =
    SignalExpr(x.mpc, [SignalTerm(reshape(float(A), :, 1), :x)])

"""
    A * d  where d::DisturbanceRef

Produce a weighted disturbance signal expression for use in `@constraint` and `@dynamics`.
"""
Base.:*(A::AbstractMatrix, d::DisturbanceRef) =
    SignalExpr(d.mpc, [SignalTerm(float(A), :d)])
Base.:*(A::AbstractVector, d::DisturbanceRef) =
    SignalExpr(d.mpc, [SignalTerm(reshape(float(A), :, 1), :d)])

"""Combine two `SignalExpr` values (sum)."""
function Base.:+(e1::SignalExpr, e2::SignalExpr)
    @assert e1.mpc === e2.mpc "Signal expressions must refer to the same MPC controller"
    SignalExpr(e1.mpc, [e1.terms; e2.terms])
end

# ---- DynamicsExpr: SignalExpr + optional constant offset -------------------

"""
    DynamicsExpr(signals, offset)

Internal type representing the right-hand side of a dynamics equation:
`x_{k+1} = F*x_k + G*u_k + Gd*d_k + f_offset`

Built automatically by expressions like `F*x + G*u + f_offset` where `x`,
`u` are signal references.  Passed to [`_set_dynamics!`](@ref) by the
[`@dynamics`](@ref) macro.
"""
struct DynamicsExpr
    signals::SignalExpr
    offset::Vector{Float64}
end

DynamicsExpr(expr::SignalExpr) = DynamicsExpr(expr, zeros(0))

Base.:+(e::SignalExpr, offset::AbstractVector) = DynamicsExpr(e, float(offset))
Base.:+(offset::AbstractVector, e::SignalExpr) = DynamicsExpr(e, float(offset))
Base.:+(e::DynamicsExpr, offset::AbstractVector) = DynamicsExpr(e.signals, e.offset + offset)
Base.:+(offset::AbstractVector, e::DynamicsExpr) = DynamicsExpr(e.signals, e.offset + offset)
Base.:+(e1::DynamicsExpr, e2::SignalExpr) = DynamicsExpr(e1.signals + e2, e1.offset)
Base.:+(e1::SignalExpr, e2::DynamicsExpr) = DynamicsExpr(e1 + e2.signals, e2.offset)

# Scalar multiplication: allow e.g. `0.5*x` (treats signal as identity-weighted)
Base.:*(a::Number, u::ControlRef) =
    SignalExpr(u.mpc, [SignalTerm(float(a) * Matrix{Float64}(I, u.mpc.model.nu, u.mpc.model.nu), :u)])
Base.:*(a::Number, x::StateRef) =
    SignalExpr(x.mpc, [SignalTerm(float(a) * Matrix{Float64}(I, x.mpc.model.nx, x.mpc.model.nx), :x)])
Base.:*(a::Number, d::DisturbanceRef) =
    SignalExpr(d.mpc, [SignalTerm(float(a) * Matrix{Float64}(I, d.mpc.model.nd, d.mpc.model.nd), :d)])

# Bare signal with no matrix weight → identity
_bare(u::ControlRef) = SignalExpr(u.mpc, [SignalTerm(Matrix{Float64}(I, u.mpc.model.nu, u.mpc.model.nu), :u)])
_bare(x::StateRef)   = SignalExpr(x.mpc, [SignalTerm(Matrix{Float64}(I, x.mpc.model.nx, x.mpc.model.nx), :x)])
_bare(d::DisturbanceRef) = SignalExpr(d.mpc, [SignalTerm(Matrix{Float64}(I, d.mpc.model.nd, d.mpc.model.nd), :d)])

# Allow bare signal + something:  x + G*u, etc.
Base.:+(x::StateRef, e::SignalExpr) = _bare(x) + e
Base.:+(e::SignalExpr, x::StateRef) = e + _bare(x)
Base.:+(u::ControlRef, e::SignalExpr) = _bare(u) + e
Base.:+(e::SignalExpr, u::ControlRef) = e + _bare(u)
Base.:+(d::DisturbanceRef, e::SignalExpr) = _bare(d) + e
Base.:+(e::SignalExpr, d::DisturbanceRef) = e + _bare(d)

# ---- _set_dynamics!: update mpc.model from a DynamicsExpr -----------------

"""
    _set_dynamics!(mpc, expr)

Internal function used by [`@dynamics`](@ref).  Extracts F, G, Gd, and
f_offset from `expr` (a `SignalExpr` or `DynamicsExpr`) and rebuilds
`mpc.model` accordingly, preserving the output matrix C, offset h_offset,
and sample time Ts.
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

# ---- Internal helpers ------------------------------------------------------

# Convert `nothing` → empty Float64 vector (for optional bound arguments)
_bound(::Nothing) = zeros(0)
_bound(x) = x

# Expand a bound to a vector of length `n` when a scalar is provided.
# Useful so that `@constraint mpc 0 <= y <= 5` works for any ny.
_to_vec_bound(::Nothing, n)           = zeros(0)
_to_vec_bound(x::AbstractVector, n)   = x
_to_vec_bound(x::Number, n)           = fill(Float64(x), n)

# Check whether a value is a signal reference (for one-sided constraint dispatch)
_is_signal(::ControlRef)     = true
_is_signal(::OutputRef)      = true
_is_signal(::StateRef)       = true
_is_signal(::DisturbanceRef) = true
_is_signal(::SignalExpr)     = true
_is_signal(::Any)            = false

# ---- Core constraint dispatch ----------------------------------------------

# Control: lb ≤ u ≤ ub  →  set_input_bounds!
function _add_constraint!(mpc::MPC, lb, ::ControlRef, ub; kw...)
    set_input_bounds!(mpc; umin=_bound(lb), umax=_bound(ub))
end

# Output: lb ≤ y ≤ ub  →  set_output_bounds!
function _add_constraint!(mpc::MPC, lb, ::OutputRef, ub;
        ks=2:mpc.Np, soft=true, binary=false, prio=0, kw...)
    set_output_bounds!(mpc; ymin=_to_vec_bound(lb, mpc.model.ny),
                            ymax=_to_vec_bound(ub, mpc.model.ny),
                            ks=ks, soft=soft, binary=binary, prio=prio)
end

# State: lb ≤ x ≤ ub  →  add_constraint! with Ax = I
function _add_constraint!(mpc::MPC, lb, ::StateRef, ub;
        ks=2:mpc.Np, soft=false, binary=false, prio=0, kw...)
    add_constraint!(mpc;
        Ax = Matrix{Float64}(I, mpc.model.nx, mpc.model.nx),
        lb = _to_vec_bound(lb, mpc.model.nx),
        ub = _to_vec_bound(ub, mpc.model.nx),
        ks=ks, soft=soft, binary=binary, prio=prio)
end

# Disturbance: lb ≤ d ≤ ub  →  set_disturbance!
function _add_constraint!(mpc::MPC, lb, ::DisturbanceRef, ub; kw...)
    set_disturbance!(mpc, _bound(lb), _bound(ub))
end

# General SignalExpr: lb ≤ Au*u + Ax*x + … ≤ ub  →  add_constraint!
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

# One-sided: dispatch based on which side carries the signal
function _add_constraint_onesided!(mpc::MPC, lhs, op::Symbol, rhs; kw...)
    if op == :(<=)
        if _is_signal(lhs)
            # signal <= bound  →  upper bound
            _add_constraint!(mpc, nothing, lhs, rhs; kw...)
        elseif _is_signal(rhs)
            # bound <= signal  →  lower bound
            _add_constraint!(mpc, lhs, rhs, nothing; kw...)
        else
            error("@constraint: no signal reference found in expression")
        end
    elseif op == :(>=)
        if _is_signal(lhs)
            # signal >= bound  →  lower bound
            _add_constraint!(mpc, rhs, lhs, nothing; kw...)
        elseif _is_signal(rhs)
            # bound >= signal  →  upper bound
            _add_constraint!(mpc, nothing, rhs, lhs; kw...)
        else
            error("@constraint: no signal reference found in expression")
        end
    else
        error("@constraint: unsupported operator $op")
    end
end

# ---- Macro helper (parse keyword args from raw macro arg list) -------------

function _parse_macro_kwargs(args)
    kw_pairs = Expr[]
    for arg in args
        if arg isa Expr && arg.head == :(=)
            push!(kw_pairs, Expr(:kw, arg.args[1], esc(arg.args[2])))
        end
    end
    return kw_pairs
end

# Build a function-call Expr with optional keyword arguments
function _make_call(fn, kw_pairs, pos_args...)
    escaped = [esc(a) for a in pos_args]
    if isempty(kw_pairs)
        return Expr(:call, fn, escaped...)
    else
        return Expr(:call, fn, Expr(:parameters, kw_pairs...), escaped...)
    end
end

# ---- Public macros ---------------------------------------------------------

"""
    @control(mpc, name)
    @control(mpc, name, umin=umin_val, umax=umax_val)

Declare a control input signal reference named `name` for the MPC controller
`mpc`.  The variable `name` is assigned a [`ControlRef`](@ref) that can be
used in [`@constraint`](@ref) expressions.

Optionally set the input bounds `umin ≤ u ≤ umax` (equivalent to calling
[`set_input_bounds!`](@ref)).

# Examples
```julia
@control mpc u
@control mpc u umin=-ones(2) umax=ones(2)
@constraint mpc -1 <= u <= 1
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
    push!(result.args, :($(esc(name_ex)) = LinearMPC.ControlRef($(esc(mpc_ex)))))
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

Declare a state signal reference named `name` for the MPC controller `mpc`.
The variable `name` is assigned a [`StateRef`](@ref) that can be used in
[`@constraint`](@ref) expressions.

# Example
```julia
@state mpc x
@constraint mpc -5 <= x <= 5
```
"""
macro state(mpc_ex, name_ex)
    quote
        $(esc(name_ex)) = LinearMPC.StateRef($(esc(mpc_ex)))
    end
end

"""
    @output(mpc, name)

Declare an output signal reference named `name` for the MPC controller `mpc`.
The variable `name` is assigned an [`OutputRef`](@ref) that can be used in
[`@constraint`](@ref) expressions.

Output constraints are soft by default (consistent with [`set_output_bounds!`](@ref)).

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

Declare a disturbance signal reference named `name` for the MPC controller
`mpc`.  The variable `name` is assigned a [`DisturbanceRef`](@ref) that can
be used in [`@constraint`](@ref) expressions.

Optionally set disturbance bounds (equivalent to calling
[`set_disturbance!`](@ref)).

# Example
```julia
@disturbance mpc d wmin=-0.1*ones(2) wmax=0.1*ones(2)
```
"""
macro disturbance(mpc_ex, name_ex, args...)
    wmin_ex = nothing; wmax_ex = nothing
    for arg in args
        (arg isa Expr && arg.head == :(=)) || continue
        k = Symbol(arg.args[1])
        k == :wmin && (wmin_ex = arg.args[2])
        k == :wmax && (wmax_ex = arg.args[2])
    end

    result = Expr(:block)
    push!(result.args, :($(esc(name_ex)) = LinearMPC.DisturbanceRef($(esc(mpc_ex)))))
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

The `signal` must be a reference created by [`@control`](@ref),
[`@state`](@ref), [`@output`](@ref), or [`@disturbance`](@ref).

| Signal type  | Calls                      | Default `soft` |
|:-------------|:---------------------------|:---------------|
| `ControlRef` | `set_input_bounds!`        | N/A            |
| `OutputRef`  | `set_output_bounds!`       | `true`         |
| `StateRef`   | `add_constraint!` (Ax = I) | `false`        |
| `SignalExpr` | `add_constraint!`          | `false`        |

Optional keyword arguments (for `StateRef` and `SignalExpr` constraints):
- `soft=false`: Penalise violations instead of enforcing them hard
- `binary=false`: Enforce with equality
- `prio=0`: Priority level for hierarchical optimisation
- `ks=2:mpc.Np`: Time steps at which the constraint is active

# Examples
```julia
@control mpc u
@output  mpc y
@state   mpc x

@constraint mpc -1 <= u <= 1           # input bounds (hard)
@constraint mpc 0 <= y <= 5            # output bounds (soft by default)
@constraint mpc u <= 0.5               # one-sided upper bound
@constraint mpc -1 <= u <= 1 soft=true # explicit soft constraint
@constraint mpc lb <= Au*u + Ax*x <= ub             # general linear
@constraint mpc lb <= Au*u + Ax*x <= ub soft=true ks=3:8
```
"""
macro constraint(mpc_ex, expr, args...)
    kw_pairs = _parse_macro_kwargs(args)

    # Double-sided: lb <= mid <= ub  or  lb >= mid >= ub
    if expr isa Expr && expr.head == :comparison && length(expr.args) == 5
        lb_ex, op1, mid_ex, op2, ub_ex = expr.args
        if op1 == :(<=) && op2 == :(<=)
            return _make_call(:(LinearMPC._add_constraint!), kw_pairs,
                              mpc_ex, lb_ex, mid_ex, ub_ex)
        elseif op1 == :(>=) && op2 == :(>=)
            # a >= mid >= b  ≡  b <= mid <= a
            return _make_call(:(LinearMPC._add_constraint!), kw_pairs,
                              mpc_ex, ub_ex, mid_ex, lb_ex)
        end
    end

    # Single-sided: a <= b  or  a >= b
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
    @dynamics(mpc, x_next = F*x + G*u + Gd*d)
    @dynamics(mpc, x_next = F*x + G*u + f_offset)
    @dynamics(mpc, x_next = F*x + G*u + Gd*d + f_offset)

Set (or update) the system dynamics on the MPC controller `mpc`.

The right-hand side must be an expression built from signal references
(`x::StateRef`, `u::ControlRef`, `d::DisturbanceRef`) created by the
[`@state`](@ref), [`@control`](@ref), and [`@disturbance`](@ref) macros, plus
optional constant vectors for the affine offset.

The state, control, and disturbance matrices are extracted automatically from
the expression:

| Expression term | Interpretation      |
|:----------------|:--------------------|
| `F*x`           | State matrix F       |
| `G*u`           | Input matrix G       |
| `Gd*d`          | Disturbance matrix Gd|
| `f_offset`      | Affine offset vector |

The left-hand side name is ignored (any symbol may be used).
The C / Dd / h_offset / Ts / operating-point settings on the existing model
are preserved.

This macro is particularly useful when combined with [`MPC(nx, nu)`](@ref),
which creates an MPC controller from state and control dimensions alone
without requiring the matrices upfront.

# Examples
```julia
# Basic usage with an existing MPC
mpc = LinearMPC.MPC(F, G; Np=10)
@state mpc x;  @control mpc u
@dynamics mpc x_next = F_new*x + G_new*u     # update dynamics in place

# Build an MPC purely through macros (no matrices needed upfront)
mpc = LinearMPC.MPC(2, 1; Np=10)             # 2 states, 1 control
@state mpc x;  @control mpc u

A = [1.0 0.1; 0.0 1.0];  B = [0.0; 1.0]
@dynamics mpc x_next = A*x + B*u

@objective mpc Q=I R=0.1
@constraint mpc -1 <= u <= 1
setup!(mpc)

# With disturbance and affine offset
mpc = LinearMPC.MPC(2, 1; nd=1, Np=10)
@state mpc x;  @control mpc u;  @disturbance mpc d
Gd = [0.1; 0.0]
fo = [0.01; 0.0]
@dynamics mpc x_next = A*x + B*u + Gd*d + fo
```
"""
macro dynamics(mpc_ex, eq_ex)
    (eq_ex isa Expr && eq_ex.head == :(=)) ||
        error("@dynamics: expected an assignment expression, e.g. `x_next = F*x + G*u`")
    rhs = eq_ex.args[2]
    return :(LinearMPC._set_dynamics!($(esc(mpc_ex)), $(esc(rhs))))
end

"""
    @objective(mpc, Q=Q_val, R=R_val, Rr=Rr_val, S=S_val, Qf=Qf_val, Qfx=Qfx_val)
Set the objective function weights for the MPC controller.
Equivalent to [`set_objective!`](@ref)`(mpc; Q=Q_val, R=R_val, …)`.

A vector is interpreted as a diagonal weight matrix.

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
