
// Warm start of the branch and bound (codegen with bnb_warm_start = true), as in LinearMPC.solve
// with the setting bnb_warm_start:
// 1. Each binary decision variable takes the value of its control one step later in the previous
//    solution (bnb_xprev_shift), rounded to the nearest of its bounds.
// 2. The QP with all binary decision variables fixed at these values gives the candidate.
// 3. The branch and bound only accepts solutions with a lower objective than the candidate.
// 4. If it finds none, the candidate is returned with the exit flag of its QP, which is
//    DAQP_EXIT_OPTIMAL (bnb_candidate_used = 1).
// 5. The returned solution is stored for the next call if its exit flag is positive.
//
// Simple bound j of the QP (lower <= x_j <= upper, with bounds that do not depend on the parameter)
// is row j of the LDP, with the row value r_j = (x_j+c_j)/s_j. Here s_j > 0 is the norm of row j of
// [I; A]*inv(R), with H = R'*R, by which the row is normalized, and c_j = (H\f(theta))_j. Hence
// dupper[j] = (upper+c_j)/s_j and dlower[j] = (lower+c_j)/s_j: x_j is at its upper bound if
// r_j = dupper[j] and at its lower bound if r_j = dlower[j], and x_j = upper+s_j*(r_j-dupper[j]).
//
// A binary decision variable is fixed for the candidate by setting both of its bounds to the
// selected one and by adding it to the working set as an immutable constraint, as the branch and
// bound of DAQP fixes binary constraints. On entry, the working set consists of the immutable
// constraints that the previous call of daqp_bnb kept (daqp_node_cleanup_workspace), and the
// fixed binary decision variables are added after them. After the solve of the candidate,
// daqp_node_cleanup_workspace removes them again, so that the branch and bound starts from the
// same working set as without the warm start.
//
// The cutoff: settings->fval_bound refers to the objective J (darnstrom/daqp#214), which for the
// LDP of the generated code (without linear term) is half the internal objective work->fval that
// daqp_bnb returns. When daqp_bnb finds an integer-feasible solution with objective J, it sets
// fval_bound to J-abs_subopt-rel_subopt*|J| and discards the nodes that cannot improve on it.
// Setting fval_bound to this value for the candidate therefore discards the same nodes as if the
// branch and bound had found the candidate itself. Since daqp_bnb discards only the nodes whose
// objective exceeds fval_bound, the bound is lowered by the margin MPC_BNB_CUTOFF_TOL*(1+|J|)
// (LinearMPC.BNB_CUTOFF_TOL), so that with abs_subopt = rel_subopt = 0 the branch and bound does
// not find the candidate again.
//
// DAQP enforces its time limit only if it is compiled with PROFILING and the solve is started
// by daqp_solve, which sets work->timer. The generated code calls daqp_bnb directly, so the time
// limit does not apply to it.

#ifndef MPC_BNB_CUTOFF_TOL
#define MPC_BNB_CUTOFF_TOL ((c_float)1e-9) // Margin of the cutoff relative to 1+|J|
#endif

int bnb_xprev_valid = 0;
c_float bnb_xprev_shift[N_BNB_BINARY];
c_float bnb_candidate_u[NX];
int bnb_candidate_used = 0;

// Discard the stored solution, so that the next call does not form a candidate from it
void mpc_reset_bnb_warm_start(void){
    bnb_xprev_valid = 0;
}

// Row value of simple bound id of the LDP for the solution in daqp_work.u (see daqp_binary_diff)
static c_float mpc_bnb_row_value(const int id){
    int j,disp;
    c_float val = 0;
    if(daqp_work.Rinv == NULL) return daqp_work.u[id]; // Hessian is identity
    for(j=id,disp=id+DAQP_R_OFFSET(id,daqp_work.n);j<daqp_work.n;j++)
        val += daqp_work.Rinv[disp++]*daqp_work.u[j];
    return val;
}

// Solve the candidate (if there is a previous solution) and the branch and bound. On return,
// daqp_work.u holds the LDP solution of the returned solution.
int mpc_bnb_warm_start(void){
    int i, id, exitflag, cand_flag = 0;
    const int n_keep = daqp_work.n_active;
    int sense_save[N_BNB_BINARY];
    c_float dupper_save[N_BNB_BINARY], dlower_save[N_BNB_BINARY];
    c_float dist_lower, dist_upper, fval_cand = 0, bound;
    const c_float fval_bound0 = daqp_work.settings->fval_bound;

    bnb_candidate_used = 0;
    if(bnb_xprev_valid){
        // Fix each binary decision variable at the bound that is nearest to the value of its
        // control one step later in the previous solution (on a tie, at the lower bound)
        for(i = 0; i < N_BNB_BINARY; i++){
            id = bnb_binary_ids[i];
            sense_save[i] = daqp_work.sense[id];
            dupper_save[i] = daqp_work.dupper[id];
            dlower_save[i] = daqp_work.dlower[id];
            if(daqp_work.sense[id] & DAQP_ACTIVE) continue; // Already fixed (equal bounds)
            if(daqp_work.sing_ind != DAQP_EMPTY_IND) continue; // A previous one was dependent
            dist_lower = bnb_xprev_shift[i]-bnb_lower[i];
            dist_upper = bnb_upper[i]-bnb_xprev_shift[i];
            if((dist_lower < 0 ? -dist_lower : dist_lower) <= (dist_upper < 0 ? -dist_upper : dist_upper)){
                daqp_work.dupper[id] = daqp_work.dlower[id];
                daqp_add_upper_lower(DAQP_ADD_LOWER_FLAG(id),&daqp_work);
            }
            else{
                daqp_work.dlower[id] = daqp_work.dupper[id];
                daqp_add_upper_lower(id,&daqp_work);
            }
            daqp_work.sense[id] |= DAQP_IMMUTABLE;
        }
        // Solve the QP of the candidate, unless a fixed binary decision variable is linearly
        // dependent on the other constraints of the working set
        if(daqp_work.sing_ind == DAQP_EMPTY_IND){
            cand_flag = daqp_bnb(&daqp_work);
            if(cand_flag > 0){
                fval_cand = daqp_work.fval;
                for(i = 0; i < NX; i++) bnb_candidate_u[i] = daqp_work.u[i];
            }
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
    }

    // Branch and bound, which only accepts solutions that improve on the candidate (objective
    // 0.5*fval_cand >= 0) by more than the suboptimality tolerances and the margin
    if(cand_flag > 0){
        bound = 0.5*fval_cand*(1-daqp_work.settings->rel_subopt-MPC_BNB_CUTOFF_TOL)
            -daqp_work.settings->abs_subopt-MPC_BNB_CUTOFF_TOL;
        if(bound < fval_bound0) daqp_work.settings->fval_bound = bound;
    }
    exitflag = daqp_bnb(&daqp_work);
    daqp_work.settings->fval_bound = fval_bound0;

    if(exitflag < 1 && cand_flag > 0){
        // No solution with a lower objective than the candidate has been found
        for(i = 0; i < NX; i++) daqp_work.u[i] = bnb_candidate_u[i];
        exitflag = cand_flag;
        bnb_candidate_used = 1;
    }

    // Store the solution for the next call (the controls one step later, in QP variables)
    bnb_xprev_valid = exitflag > 0;
    if(bnb_xprev_valid){
        for(i = 0; i < N_BNB_BINARY; i++){
            id = bnb_shift_ids[i];
            bnb_xprev_shift[i] = bnb_shift_upper[i]
                + bnb_shift_scaling[i]*(mpc_bnb_row_value(id)-daqp_work.dupper[id]);
        }
    }
    return exitflag;
}
