function constraint_tightening(Ax,F,ks,wmin,wmax) 
    m,nx = size(Ax)
    nk = length(ks)
    tight_upper, tight_lower = zeros(m*nk), zeros(m*nk)
    accum_upper, accum_lower = zeros(m), zeros(m)

    # TODO: add special case for uncertainty in x0 (corresponding to w in k=1)

    ki = 1+sum(ks .< 2); # Ignore before k=2 
    Ck = Ax*F
    for k in 2:maximum(ks) 
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
        Ck *=  F # Assume that ks is a UnitRange
    end
    return tight_upper,tight_lower
end
