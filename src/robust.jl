"""
    constraint_tightening(Ax,F,ks,wmin,wmax,x0_uncertainty;Aup=zeros(0,0),Bw=I,F0=F)

Tightening of the constraint rows `lb ≤ Ax eₜ + Aup eₜ₋₁ ≤ ub` (in terms of the deviation `eₜ` of
the state from its nominal prediction) at the time steps `t = k-1`, `k ∈ ks`, for the error dynamics

    e₁ = F0 e₀ + Bw w₀,   eₜ₊₁ = F eₜ + Bw wₜ (t ≥ 1),   -δ ≤ e₀ ≤ δ,   wmin ≤ wₜ ≤ wmax,

where `F` is the prestabilized dynamics `F0 - G*K`, `F0` the open-loop dynamics and
`δ = x0_uncertainty`. The control `u₀` is computed from the estimate `x̂₀` and therefore carries no
error, which is why the first step is open loop. `Ax` holds the coefficients of the error of the state
and of the control of the same time step (`Ax - Au*K` for the prestabilizing feedback `u = -K x + v`),
`Aup` those of the previous control (`-Aup*K`); the `Aup` term is present from `t = 2` on. Rows at
`k = 1` (`t = 0`) are not tightened, since `x₀` cannot be influenced.

Returns the amounts `(tight_upper, tight_lower)` by which `ub` is to be decreased and `lb` increased,
with the rows ordered as `ks`. They are the support functions of the reachable error set,

    tight_upper(t) = |Mₜ₋₁ F0| δ + Σ_{s=0}^{t-1} max_{w} Mₛ Bw w,
    tight_lower(t) = |Mₜ₋₁ F0| δ - Σ_{s=0}^{t-1} min_{w} Mₛ Bw w,

with `Mₛ = Ax Fˢ + Aup Fˢ⁻¹`, the `Aup` term present for `s ≥ 1`.
"""
function constraint_tightening(Ax,F,ks,wmin,wmax,x0_uncertainty;Aup=zeros(0,0),Bw=I,F0=F)
    m,nx = size(Ax)
    Bw = Bw isa UniformScaling ? Matrix{Float64}(Bw,nx,nx) : Bw
    has_up = !isempty(Aup) && !iszero(Aup)
    kmax = maximum(ks; init=1)
    upper, lower = zeros(m, kmax), zeros(m, kmax)
    accum_upper, accum_lower = zeros(m), zeros(m)

    Fs = Matrix{Float64}(I,nx,nx) # Fˢ
    Fsm1 = zeros(nx,nx)           # Fˢ⁻¹
    for t in 1:kmax-1
        s = t-1
        Ms = Ax*Fs
        has_up && s >= 1 && (Ms += Aup*Fsm1)
        # The disturbance wₜ₋₁₋ₛ, which reaches eₜ through Fˢ Bw
        for (i,ci) in enumerate(eachrow(Ms*Bw))
            accum_upper[i] += sum(max(ci[j]*wmax[j], ci[j]*wmin[j]) for j in eachindex(ci); init=0.0)
            accum_lower[i] -= sum(min(ci[j]*wmax[j], ci[j]*wmin[j]) for j in eachindex(ci); init=0.0)
        end
        # The error of the initial state, which reaches eₜ through Fᵗ⁻¹ F0
        e0 = abs.(Ms*F0)*x0_uncertainty
        upper[:,t+1] = accum_upper + e0
        lower[:,t+1] = accum_lower + e0
        Fsm1, Fs = Fs, Fs*F
    end
    # Rows in the order of ks
    tight_upper, tight_lower = vec(upper[:,ks]), vec(lower[:,ks])
    return tight_upper,tight_lower
end

# Whether the constraints are tightened for a bounded disturbance or an uncertain initial state
is_robust(mpc) = !iszero(mpc.model.wmin) || !iszero(mpc.model.wmax) || !iszero(mpc.Δx0)
