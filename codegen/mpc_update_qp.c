void mpc_update_qp(c_float* th, c_float* dupper, c_float* dlower){
    int i,j,disp;
    c_float b_shift_th;
    for(i =0,disp=0; i < N_CONSTR; i++){
        b_shift_th = 0;
        for(j = 0; j < N_THETA; j++) b_shift_th += Dth[disp++]*th[j];
        dupper[i] = du[i] + b_shift_th;
        dlower[i] = dl[i] + b_shift_th;
    }
}

// Assumes that x is stacked such that the first
// N_CONTROL elements are the controls at the first time step
void mpc_get_solution(c_float* th, c_float* control, c_float* xstar){
    int i,j,disp;
    c_float ctr_shift_th;
    for(i = 0, disp=0; i < N_CONTROL; i++){
        ctr_shift_th = u_offset[i];
        for(j = 0; j < N_THETA; j++) ctr_shift_th += Uth_offset[disp++]*th[j];
        control[i] = uscaling[i]*xstar[i]+ctr_shift_th;
    }
}

#include"daqp.h"
#ifdef DAQP_BNB
#include "bnb.h"
#endif

#if N_AFFINE_PARAMETER > 0
int mpc_compute_control(c_float* control, c_float* state, c_float* reference, c_float* disturbance, c_float* affine_parameter){
    mpc_update_parameter(mpc_parameter, control, state, reference, disturbance, affine_parameter);
#else
int mpc_compute_control(c_float* control, c_float* state, c_float* reference, c_float* disturbance){
    mpc_update_parameter(mpc_parameter, control, state, reference, disturbance);
#endif
    // update problem
    mpc_update_qp(mpc_parameter,daqp_work.dupper,daqp_work.dlower);
    daqp_work.reuse_ind=0; // clear workspace cache

#ifdef DAQP_BNB
    daqp_node_cleanup_workspace(0, &daqp_work);
#if N_EQUALITY > 0
    // The cleanup removes the equality constraints from the working set as well, and daqp_bnb takes
    // the constraints active on entry as the equality constraints of the problem: reactivate them
    for(int i = 0; i < N_EQUALITY; i++) daqp_work.sense[equality_ids[i]] |= DAQP_ACTIVE + DAQP_IMMUTABLE;
    reset_daqp_workspace(&daqp_work);
    daqp_activate_constraints(&daqp_work);
#endif
    int exitflag = daqp_bnb(&daqp_work);
#else
#ifndef DAQP_WARMSTART
    daqp_deactivate_constraints(&daqp_work);
    reset_daqp_workspace(&daqp_work);
#endif
    int exitflag = daqp_ldp(&daqp_work);
#endif

    ldp2qp_solution(&daqp_work);
    mpc_get_solution(mpc_parameter,control,daqp_work.x);
    return exitflag;
}
