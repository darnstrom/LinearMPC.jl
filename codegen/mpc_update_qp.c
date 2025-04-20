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
        ctr_shift_th = 0;
        for(j = 0; j < N_THETA; j++) ctr_shift_th += Xth[disp++]*th[j];
        control[i] = uscaling[i]*xstar[i]+ctr_shift_th;
    }
}

#include"daqp.h"
#ifdef DAQP_BNB
#include "bnb.h"
#endif

int mpc_compute_control(c_float* theta, c_float* control, DAQPWorkspace* work){
    mpc_update_qp(theta,work->dupper,work->dlower);

#ifdef DAQP_BNB
    node_cleanup_workspace(0, work);
    int exitflag = daqp_bnb(work);
#else
#ifndef DAQP_WARMSTART
    deactivate_constraints(work);
    reset_daqp_workspace(work);
#endif
    int exitflag = daqp_ldp(work);
#endif

    ldp2qp_solution(work);
    mpc_get_solution(theta,control,work->x);
    return exitflag;
}
