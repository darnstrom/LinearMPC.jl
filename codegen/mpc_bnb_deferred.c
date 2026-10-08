
// Deferred binary controls (codegen of an MPC with groups of deferred binary controls), as in
// LinearMPC.solve with groups declared by defer_binary_controls!:
// 1. With the warm start (DAQP_BNB_WARMSTART), the candidate: the binary decision variables that
//    are not deferred are fixed at the bound that is nearest to the value of their control one
//    step later in the previous solution (bnb_xprev_shift), the deferred ones are relaxed, and they
//    are resolved for the solution of this QP as in 3.
// 2. The relaxed search: daqp_bnb with daqp_work.bnb->bin_ids and nb replaced by
//    bnb_relaxed_bin_ids and N_BNB_RELAXED_BIN, which do not contain the deferred binary decision
//    variables, and with the objective of the candidate as cutoff (see mpc_bnb_warm_start.c). If
//    it finds no solution, the candidate is returned, or, without a candidate, the result of the
//    branch and bound without relaxation.
// 3. The resolution of the relaxed solution: the binary decision variables that are not deferred
//    are fixed at their values, and the candidates of the groups are combined (see
//    mpc_bnb_prepare_group). Each combination is a QP with all binary decision variables fixed,
//    which is solved as the candidate of the warm start: both LDP bounds of a binary decision
//    variable are set to the selected one, it is added to the working set as an immutable
//    constraint, and daqp_node_cleanup_workspace removes it after the solve.
// 4. The best of these and the candidate is returned if its objective exceeds that of the relaxed
//    search by at most bnb_deferred_tol. Otherwise, the branch and bound without relaxation runs
//    with its objective as cutoff, and its result is returned, or the best of these if it finds
//    no better solution.
// The candidates, the limit BNB_MAX_COMBINATIONS and the order of the combinations are those of
// LinearMPC.bnb_resolve_deferred. Objectives are compared in the internal objective of DAQP
// (daqp_work.fval), in which the difference between two solutions for the same parameter is
// twice the difference between their objectives. A cutoff with the objective J = 0.5*fval of a
// solution sets fval_bound to J-abs_subopt-rel_subopt*|J|, the value that daqp_bnb sets when it
// finds that solution itself, lowered by the margin MPC_BNB_CUTOFF_TOL*(1+|J|) (see
// mpc_bnb_warm_start.c). The time limit of DAQP does not apply (see mpc_bnb_warm_start.c).
//
// The positions of binary decision variables refer to bnb_binary_ids. A binary decision variable
// is normalized by its bounds: 0 at bnb_lower and 1 at bnb_upper. Its value in an LDP solution is
// x_j = upper+s_j*(r_j-dupper[j]), with the row value r_j and the scaling s_j of its LDP row (see
// mpc_bnb_warm_start.c).

#define MPC_BNB_ROUND_TOL ((c_float)1e-6) // Tolerance of the rounding (LinearMPC.BNB_ROUND_TOL)
#ifndef MPC_BNB_CUTOFF_TOL
#define MPC_BNB_CUTOFF_TOL ((c_float)1e-9) // Margin of the cutoff relative to 1+|J| (LinearMPC.BNB_CUTOFF_TOL)
#endif

int bnb_source = MPC_BNB_SOURCE_NONE;
static c_float bnb_relaxed_u[NX]; // LDP solution of the relaxed search
static c_float bnb_best_u[NX]; // LDP solution of the best combination
static c_float bnb_best_fval; // Its internal objective
static int bnb_best_flag; // Its exit flag (0 if there is none)
static c_float bnb_z[N_BNB_BINARY]; // Normalized values of the binary decision variables in a relaxed solution
static unsigned char bnb_val[N_BNB_BINARY]; // Values of the binary decision variables (1: upper bound)
static unsigned char bnb_base[N_BNB_DEFERRED]; // Rounded values of the binary decision variables of the groups
static int bnb_code_near[N_BNB_DEFERRED]; // Integer encoding: the code of the nearest integer at each step
static int bnb_code_alt[N_BNB_DEFERRED]; // and of the other integer next to the relaxed value (-1 if there is none)
static int bnb_flip[2*N_BNB_GROUPS]; // Sum-up rounding: the steps that the other candidates change
static int bnb_ncand[N_BNB_GROUPS]; // Number of candidates of each group
static int bnb_choice[N_BNB_GROUPS]; // Candidate of each group in a combination

// Discard the stored solution of the warm start (if it is generated) and the source of the latest solution
void mpc_reset_bnb_deferred(void){
#ifdef DAQP_BNB_WARMSTART
    mpc_reset_bnb_warm_start();
#endif
    bnb_source = MPC_BNB_SOURCE_NONE;
}

// Row value of simple bound id of the LDP for the LDP solution u
static c_float mpc_bnb_deferred_row_value(const int id, const c_float* u){
    int j,disp;
    c_float val = 0;
    if(daqp_work.Rinv == NULL) return u[id]; // Hessian is identity
    for(j=id,disp=id+DAQP_R_OFFSET(id,daqp_work.n);j<daqp_work.n;j++)
        val += daqp_work.Rinv[disp++]*u[j];
    return val;
}

// Normalized values of the binary decision variables in the LDP solution u (with the bounds of
// the LDP for the parameter)
static void mpc_bnb_normalized_values(const c_float* u){
    int i, id;
    c_float x;
    for(i = 0; i < N_BNB_BINARY; i++){
        id = bnb_binary_ids[i];
        x = bnb_upper[i]+bnb_binary_scaling[i]*(mpc_bnb_deferred_row_value(id,u)-daqp_work.dupper[id]);
        x = (x-bnb_lower[i])/(bnb_upper[i]-bnb_lower[i]);
        bnb_z[i] = x < 0 ? 0 : (x > 1 ? 1 : x);
    }
}

// Solve with the binary decision variables at the values in bnb_val, the deferred ones only if
// fix_deferred is nonzero; otherwise they are relaxed. If the exit flag is positive, the LDP
// solution is in daqp_work.u and its internal objective in daqp_work.fval.
static int mpc_bnb_fixed_solve(const int fix_deferred){
    int i, id, exitflag = 0;
    const int n_keep = daqp_work.n_active;
    int* const bin_ids0 = daqp_work.bnb->bin_ids;
    const int nb0 = daqp_work.bnb->nb;
    int sense_save[N_BNB_BINARY];
    c_float dupper_save[N_BNB_BINARY], dlower_save[N_BNB_BINARY];
    for(i = 0; i < N_BNB_BINARY; i++){
        id = bnb_binary_ids[i];
        sense_save[i] = daqp_work.sense[id];
        dupper_save[i] = daqp_work.dupper[id];
        dlower_save[i] = daqp_work.dlower[id];
        if(!fix_deferred && bnb_is_deferred[i]) continue;
        if(daqp_work.sense[id] & DAQP_ACTIVE) continue; // Already fixed (equal bounds)
        if(daqp_work.sing_ind != DAQP_EMPTY_IND) continue; // A previous one was dependent
        if(bnb_val[i]){
            daqp_work.dlower[id] = daqp_work.dupper[id];
            daqp_add_upper_lower(id,&daqp_work);
        }
        else{
            daqp_work.dupper[id] = daqp_work.dlower[id];
            daqp_add_upper_lower(DAQP_ADD_LOWER_FLAG(id),&daqp_work);
        }
        daqp_work.sense[id] |= DAQP_IMMUTABLE;
    }
    // Solve, unless a fixed binary decision variable is linearly dependent on the other
    // constraints of the working set
    if(daqp_work.sing_ind == DAQP_EMPTY_IND){
        if(!fix_deferred){
            daqp_work.bnb->bin_ids = bnb_relaxed_bin_ids;
            daqp_work.bnb->nb = N_BNB_RELAXED_BIN;
        }
        exitflag = daqp_bnb(&daqp_work);
        daqp_work.bnb->bin_ids = bin_ids0;
        daqp_work.bnb->nb = nb0;
    }
    // Remove the fixed binary decision variables from the working set and restore them
    daqp_node_cleanup_workspace(n_keep,&daqp_work);
    for(i = 0; i < N_BNB_BINARY; i++){
        id = bnb_binary_ids[i];
        daqp_work.sense[id] = sense_save[i];
        daqp_work.dupper[id] = dupper_save[i];
        daqp_work.dlower[id] = dlower_save[i];
    }
    daqp_work.reuse_ind = 0;
    return exitflag;
}

// Solve the combination in bnb_val and keep it if it is the best so far. Returns its internal
// objective, or DAQP_INF if it is infeasible.
static c_float mpc_bnb_evaluate(void){
    int i;
    const int flag = mpc_bnb_fixed_solve(1);
    if(flag < 1) return DAQP_INF;
    if(bnb_best_flag < 1 || daqp_work.fval < bnb_best_fval){
        for(i = 0; i < NX; i++) bnb_best_u[i] = daqp_work.u[i];
        bnb_best_fval = daqp_work.fval;
        bnb_best_flag = flag;
    }
    return daqp_work.fval;
}

// Whether the values in bnb_val satisfy the logic constraints of group g
static int mpc_bnb_logic_feasible(const int g){
    int r, j, i;
    c_float act;
    const c_float tol = daqp_work.settings->primal_tol;
    for(r = bnb_group_rstart[g]; r < bnb_group_rstart[g+1]; r++){
        for(j = bnb_lrow_start[r], act = 0; j < bnb_lrow_start[r+1]; j++){
            i = bnb_lrow_pos[j];
            act += bnb_lrow_coef[j]*(bnb_val[i] ? bnb_upper[i] : bnb_lower[i]);
        }
        if(act > bnb_lrow_upper[r]+tol || act < bnb_lrow_lower[r]-tol) return 0;
    }
    return 1;
}

// The candidates of group g from the normalized relaxed values in bnb_z, with the binary decision
// variables that are not deferred at their values in bnb_val. The first candidate rounds the
// relaxed values.
// - Integer encoding (kind 0): at each binary step, the codes of the encodable integers just
//   below and above the relaxed value (the first code in the order of the binary numbers for an
//   integer), the nearest first, combined over the steps (candidate c takes the other integer at
//   the t-th step that has two if bit t of c is set). If there are more than
//   BNB_MAX_COMBINATIONS, only the first.
// - Single control (kind 1): sum-up rounding, and the sequences in which the first step that is
//   off is on and the last step that is on is off, if the logic constraints allow it.
// - Enumeration (kind 2): all combinations (candidate c changes the rounded values at the bits
//   of c), or only the rounded values if there are more than BNB_MAX_COMBINATIONS.
static void mpc_bnb_prepare_group(const int g){
    int i, k, c, t, m, steps, code_below, code_above;
    const int start = bnb_group_start[g];
    const int n = bnb_group_start[g+1]-start;
    const int* pos = bnb_group_pos+start;
    const c_float* w;
    c_float v, val, lo, hi, wmax, tol, below, above, d;
    switch(bnb_group_kind[g]){
    case 0:
        m = bnb_group_nctrl[g];
        steps = n/m;
        w = bnb_group_weights+bnb_group_wstart[g];
        for(i = 0, lo = 0, hi = 0, wmax = 0; i < m; i++){
            if(w[i] < 0) lo += w[i]; else hi += w[i];
            if((w[i] < 0 ? -w[i] : w[i]) > wmax) wmax = w[i] < 0 ? -w[i] : w[i];
        }
        tol = (c_float)1e-6*(1+wmax);
        for(k = 0, t = 0; k < steps; k++){
            for(i = 0, v = 0; i < m; i++) v += w[i]*bnb_z[pos[k*m+i]];
            v = v < lo ? lo : (v > hi ? hi : v);
            code_below = -1; code_above = -1; below = 0; above = 0;
            for(c = 0; c < (1 << m); c++){
                for(i = 0, val = 0; i < m; i++) if((c >> i) & 1) val += w[i];
                if(val <= v+tol && (code_below < 0 || val > below)){ below = val; code_below = c; }
                if(val >= v-tol && (code_above < 0 || val < above)){ above = val; code_above = c; }
            }
            if(v-below <= above-v+tol){
                bnb_code_near[start+k] = code_below;
                bnb_code_alt[start+k] = code_above;
            }
            else{
                bnb_code_near[start+k] = code_above;
                bnb_code_alt[start+k] = code_below;
            }
            if(below == above) bnb_code_alt[start+k] = -1;
            else t++;
        }
        bnb_ncand[g] = (t < 30 && (1L << t) <= BNB_MAX_COMBINATIONS) ? (1 << t) : 1;
        break;
    case 1:
        for(k = 0, d = 0; k < n; k++){
            d += bnb_z[pos[k]]*bnb_group_dt[start+k];
            bnb_base[start+k] = d >= (0.5-MPC_BNB_ROUND_TOL)*bnb_group_dt[start+k];
            if(bnb_base[start+k]) d -= bnb_group_dt[start+k];
            bnb_val[pos[k]] = bnb_base[start+k];
        }
        bnb_ncand[g] = 1;
        for(k = 0; k < n; k++){ // The first step that is off
            if(bnb_base[start+k]) continue;
            bnb_val[pos[k]] = 1;
            c = mpc_bnb_logic_feasible(g);
            bnb_val[pos[k]] = 0;
            if(c){ bnb_flip[2*g+bnb_ncand[g]-1] = k; bnb_ncand[g]++; break; }
        }
        for(k = n-1; k >= 0; k--){ // The last step that is on
            if(!bnb_base[start+k]) continue;
            bnb_val[pos[k]] = 0;
            c = mpc_bnb_logic_feasible(g);
            bnb_val[pos[k]] = 1;
            if(c){ bnb_flip[2*g+bnb_ncand[g]-1] = k; bnb_ncand[g]++; break; }
        }
        if(bnb_ncand[g] > BNB_MAX_COMBINATIONS) bnb_ncand[g] = 1;
        break;
    default:
        for(k = 0; k < n; k++) bnb_base[start+k] = bnb_z[pos[k]] > 0.5+MPC_BNB_ROUND_TOL;
        bnb_ncand[g] = (n < 30 && (1L << n) <= BNB_MAX_COMBINATIONS) ? (1 << n) : 1;
    }
}

// Set the values in bnb_val of the binary decision variables of group g to its candidate c
static void mpc_bnb_set_candidate(const int g, const int c){
    int i, k, t, m, code;
    const int start = bnb_group_start[g];
    const int n = bnb_group_start[g+1]-start;
    const int* pos = bnb_group_pos+start;
    switch(bnb_group_kind[g]){
    case 0:
        m = bnb_group_nctrl[g];
        for(k = 0, t = 0; k < n/m; k++){
            code = bnb_code_near[start+k];
            if(bnb_code_alt[start+k] >= 0){
                if((c >> t) & 1) code = bnb_code_alt[start+k];
                t++;
            }
            for(i = 0; i < m; i++) bnb_val[pos[k*m+i]] = (code >> i) & 1;
        }
        break;
    case 1:
        for(k = 0; k < n; k++) bnb_val[pos[k]] = bnb_base[start+k];
        if(c > 0) bnb_val[pos[bnb_flip[2*g+c-1]]] ^= 1;
        break;
    default:
        for(k = 0; k < n; k++) bnb_val[pos[k]] = bnb_base[start+k] ^ ((c >> k) & 1);
    }
}

// Resolve the deferred binary decision variables of the relaxed LDP solution u with internal
// objective fval and exit flag flag (see 3. above). Returns whether a solution has been found;
// it is then in bnb_best_u, with its internal objective bnb_best_fval and its exit flag
// bnb_best_flag. If the deferred binary decision variables of u are integer feasible, it is u.
static int mpc_bnb_resolve(const c_float* u, const c_float fval, const int flag){
    int i, g, c, best_c, total;
    int integral = 1;
    c_float current, f;
    bnb_best_flag = 0;
    mpc_bnb_normalized_values(u);
    for(i = 0; i < N_BNB_BINARY; i++){
        bnb_val[i] = bnb_z[i] > 0.5+MPC_BNB_ROUND_TOL;
        if(bnb_is_deferred[i] &&
           (bnb_z[i] < 0.5 ? bnb_z[i] : 1-bnb_z[i])*(bnb_upper[i]-bnb_lower[i]) > daqp_work.settings->primal_tol)
            integral = 0;
    }
    if(integral){
        for(i = 0; i < NX; i++) bnb_best_u[i] = u[i];
        bnb_best_fval = fval;
        bnb_best_flag = flag;
        return 1;
    }
    for(g = 0, total = 1; g < N_BNB_GROUPS; g++){
        mpc_bnb_prepare_group(g);
        total *= bnb_ncand[g];
        if(total > BNB_MAX_COMBINATIONS) total = BNB_MAX_COMBINATIONS+1;
    }
    for(g = 0; g < N_BNB_GROUPS; g++){
        bnb_choice[g] = 0;
        mpc_bnb_set_candidate(g,0);
    }
    if(total <= BNB_MAX_COMBINATIONS){
        // All combinations, with the first group changing fastest
        while(1){
            mpc_bnb_evaluate();
            for(g = 0; g < N_BNB_GROUPS; g++){
                if(++bnb_choice[g] < bnb_ncand[g]){
                    mpc_bnb_set_candidate(g,bnb_choice[g]);
                    break;
                }
                bnb_choice[g] = 0;
                mpc_bnb_set_candidate(g,0);
            }
            if(g == N_BNB_GROUPS) break;
        }
    }
    else{
        // The groups in turn: the candidates of a group with the groups before it at their best
        // and the groups after it at their first candidate
        current = mpc_bnb_evaluate();
        for(g = 0; g < N_BNB_GROUPS; g++){
            for(c = 1, best_c = 0; c < bnb_ncand[g]; c++){
                mpc_bnb_set_candidate(g,c);
                f = mpc_bnb_evaluate();
                if(f < current){ current = f; best_c = c; }
            }
            mpc_bnb_set_candidate(g,best_c);
        }
    }
    return bnb_best_flag > 0;
}

// Solve the candidate (with the warm start), the relaxed search, the resolution and, if needed,
// the branch and bound without relaxation. On return, daqp_work.u holds the LDP solution of the
// returned solution.
int mpc_bnb_deferred(void){
    int i, exitflag, flag;
    int cand_flag = 0;
    c_float fval_cand = 0, fval_relaxed, bound;
    int* const bin_ids0 = daqp_work.bnb->bin_ids;
    const int nb0 = daqp_work.bnb->nb;
    const c_float fval_bound0 = daqp_work.settings->fval_bound;
#ifdef DAQP_BNB_WARMSTART
    static c_float cand_u[NX];
    c_float dist_lower, dist_upper;
    int id;
#endif

    bnb_best_flag = 0;
#ifdef DAQP_BNB_WARMSTART
    bnb_candidate_used = 0;
    if(bnb_xprev_valid){
        // The binary decision variables that are not deferred at the bound that is nearest to the
        // value of their control one step later in the previous solution (on a tie, the lower one)
        for(i = 0; i < N_BNB_BINARY; i++){
            dist_lower = bnb_xprev_shift[i]-bnb_lower[i];
            dist_upper = bnb_upper[i]-bnb_xprev_shift[i];
            bnb_val[i] = (dist_lower < 0 ? -dist_lower : dist_lower) > (dist_upper < 0 ? -dist_upper : dist_upper);
        }
        flag = mpc_bnb_fixed_solve(0);
        if(flag > 0){
            for(i = 0; i < NX; i++) bnb_relaxed_u[i] = daqp_work.u[i];
            if(mpc_bnb_resolve(bnb_relaxed_u,daqp_work.fval,flag)){
                cand_flag = bnb_best_flag;
                fval_cand = bnb_best_fval;
                for(i = 0; i < NX; i++) cand_u[i] = bnb_best_u[i];
            }
        }
    }
#endif

    // Relaxed search, which only accepts solutions that improve on the candidate (objective
    // 0.5*fval_cand >= 0) by more than the suboptimality tolerances and the margin
    if(cand_flag > 0){
        bound = 0.5*fval_cand*(1-daqp_work.settings->rel_subopt-MPC_BNB_CUTOFF_TOL)
            -daqp_work.settings->abs_subopt-MPC_BNB_CUTOFF_TOL;
        if(bound < fval_bound0) daqp_work.settings->fval_bound = bound;
    }
    daqp_work.bnb->bin_ids = bnb_relaxed_bin_ids;
    daqp_work.bnb->nb = N_BNB_RELAXED_BIN;
    exitflag = daqp_bnb(&daqp_work);
    daqp_work.bnb->bin_ids = bin_ids0;
    daqp_work.bnb->nb = nb0;
    daqp_work.settings->fval_bound = fval_bound0;
    bnb_source = MPC_BNB_SOURCE_SEARCH;

    if(exitflag < 1){
        if(cand_flag > 0){
            // No solution improves on the candidate by more than the suboptimality tolerances
            bnb_best_flag = cand_flag;
            bnb_source = MPC_BNB_SOURCE_CANDIDATE;
        }
        else
            exitflag = daqp_bnb(&daqp_work);
    }
    else{
        // Resolution, and the best of it and the candidate
        flag = exitflag;
        fval_relaxed = daqp_work.fval;
        for(i = 0; i < NX; i++) bnb_relaxed_u[i] = daqp_work.u[i];
        if(mpc_bnb_resolve(bnb_relaxed_u,fval_relaxed,flag)) bnb_source = MPC_BNB_SOURCE_DEFERRED;
        if(cand_flag > 0 && (bnb_best_flag < 1 || fval_cand <= bnb_best_fval)){
            bnb_best_flag = cand_flag;
            bnb_best_fval = fval_cand;
            bnb_source = MPC_BNB_SOURCE_CANDIDATE;
        }
        if(bnb_best_flag < 1 || 0.5*(bnb_best_fval-fval_relaxed) > bnb_deferred_tol){
            // Branch and bound without relaxation, which only accepts solutions that improve on
            // the best one by more than the suboptimality tolerances and the margin
            if(bnb_best_flag > 0){
                bound = 0.5*bnb_best_fval*(1-daqp_work.settings->rel_subopt-MPC_BNB_CUTOFF_TOL)
                    -daqp_work.settings->abs_subopt-MPC_BNB_CUTOFF_TOL;
                if(bound < fval_bound0) daqp_work.settings->fval_bound = bound;
            }
            exitflag = daqp_bnb(&daqp_work);
            daqp_work.settings->fval_bound = fval_bound0;
            if(exitflag > 0 || bnb_best_flag < 1){
                bnb_best_flag = 0;
                bnb_source = MPC_BNB_SOURCE_SEARCH;
            }
        }
        else exitflag = bnb_best_flag; // Accepted, returned below
    }

    // Return the best combination or the candidate unless the branch and bound has found a solution
    if(bnb_best_flag > 0 && bnb_source != MPC_BNB_SOURCE_SEARCH){
#ifdef DAQP_BNB_WARMSTART
        if(bnb_source == MPC_BNB_SOURCE_CANDIDATE){
            for(i = 0; i < NX; i++) daqp_work.u[i] = cand_u[i];
            bnb_candidate_used = 1;
        }
        else
#endif
        for(i = 0; i < NX; i++) daqp_work.u[i] = bnb_best_u[i];
        exitflag = bnb_best_flag;
    }

#ifdef DAQP_BNB_WARMSTART
    // Store the solution for the next call (the controls one step later, in QP variables)
    bnb_xprev_valid = exitflag > 0;
    if(bnb_xprev_valid){
        for(i = 0; i < N_BNB_BINARY; i++){
            id = bnb_shift_ids[i];
            bnb_xprev_shift[i] = bnb_shift_upper[i]
                + bnb_shift_scaling[i]*(mpc_bnb_deferred_row_value(id,daqp_work.u)-daqp_work.dupper[id]);
        }
    }
#endif
    return exitflag;
}
