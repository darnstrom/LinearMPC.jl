function constraint_tightening(Ax,F,ks,wmin,wmax,x0_uncertainty) 
    m,nx = size(Ax)
    accum_upper, accum_lower = zeros(m), zeros(m)

    Ck = Ax
    # x0 uncertainty
    for (i,ci) in enumerate(eachrow(Ck))
        accum_upper[i] = sum(abs(ci[j] * x0_uncertainty[j]) for j in 1:nx)
        accum_lower[i] = accum_upper[i]
    end

    # Tightening at each time step up to the last one in ks (none at k = 1)
    kmax = maximum(ks; init=1)
    upper, lower = zeros(m, kmax), zeros(m, kmax)
    for k in 2:kmax
        Ck *= F
        for (i,ci) in enumerate(eachrow(Ck))
            accum_upper[i] += sum(ci[j] * ( ci[j] > 0 ? wmax[j] : wmin[j]) for j in 1:nx)
            accum_lower[i] -= sum(ci[j] * ( ci[j] < 0 ? wmax[j] : wmin[j]) for j in 1:nx)
        end
        upper[:,k] = accum_upper
        lower[:,k] = accum_lower
    end
    # Rows in the order of ks
    tight_upper, tight_lower = vec(upper[:,ks]), vec(lower[:,ks])
    return tight_upper,tight_lower
end
