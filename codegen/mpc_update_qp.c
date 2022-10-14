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
        control[i] = xstar[i]+ctr_shift_th;
    }
}
