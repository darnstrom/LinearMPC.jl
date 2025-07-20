function constraint_tightening(Ax,F,ks,wmin,wmax) 
    m,nx = size(Ax)
    nk = length(ks)
    accum_upper, accum_lower = zeros(m), zeros(m)
    tight_upper, tight_lower = zeros(m*nk), zeros(m*nk)
    FK = F 
    Ck = Ax* FK^(first(ks)-1)
    for (ki,k) in enumerate(ks)
        for (i,ci) in enumerate(eachrow(Ck))
            accum_upper[i] += sum(ci[j] * ( ci[j] > 0 ? wmax[j] : wmin[j]) for j in 1:nx)
            accum_lower[i] -= sum(ci[j] * ( ci[j] < 0 ? wmax[j] : wmin[j]) for j in 1:nx)
        end
        tight_upper[m*(ki-1)+1:m*ki] = accum_upper
        tight_lower[m*(ki-1)+1:m*ki] = accum_lower
        if(ki < nk)
            Δk = ks[ki+1] - k
            Ck *= Δk == 1 ? FK : FK^Δk
        end
    end
    return tight_upper,tight_lower
end
