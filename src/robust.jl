function constraint_tightening(Ax,F,ks,wmin,wmax,x0_uncertainty) 
    m,nx = size(Ax)
    nk = length(ks)
    tight_upper, tight_lower = zeros(m*nk), zeros(m*nk)
    accum_upper, accum_lower = zeros(m), zeros(m)

    Ck = Ax
    # x0 uncertainty
    for (i,ci) in enumerate(eachrow(Ck))
        accum_upper[i] = sum(abs(ci[j] * x0_uncertainty[j]) for j in 1:nx)
        accum_lower[i] = accum_upper[i]
    end

    ki = 1+sum(ks .< 2); # Ignore before k=2 
    for k in 2:maximum(ks) 
        Ck *= F
        for (i,ci) in enumerate(eachrow(Ck))
            accum_upper[i] += sum(ci[j] * ( ci[j] > 0 ? wmax[j] : wmin[j]) for j in 1:nx)
            accum_lower[i] -= sum(ci[j] * ( ci[j] < 0 ? wmax[j] : wmin[j]) for j in 1:nx)
        end
        if k == ks[ki]
            tight_upper[m*(ki-1)+1:m*ki] = accum_upper
            tight_lower[m*(ki-1)+1:m*ki] = accum_lower
            ki += 1
            ki > nk && break;
        end
    end
    return tight_upper,tight_lower
end
