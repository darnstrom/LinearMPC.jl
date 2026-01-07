void mpc_predict_state(c_float* state, c_float* control){
    int i,j,disp=0;
    c_float state_old[N_STATE];
    for(i=0;i<N_STATE;i++) state_old[i] = state[i];
    for(i=0,disp=0;i<N_STATE;i++){
        state[i] = MPC_PLANT_DYNAMICS[disp++];
        for(j=0;j<N_STATE;j++) state[i] += MPC_PLANT_DYNAMICS[disp++]*state_old[j];
        for(j=0;j<N_CONTROL;j++) state[i] += MPC_PLANT_DYNAMICS[disp++]*control[j];
    }
}

void mpc_correct_state(c_float* state, c_float* measurement){
    int i,j,disp_C=0,disp_K=0;
    c_float innovation, state_old[N_STATE];
    for(i=0;i<N_STATE;i++) state_old[i] = state[i];
    for(j=0;j<N_MEASUREMENT;j++){
        innovation = measurement[j]-C_OBSERVER[disp_C++];
        for(i=0;i<N_STATE;i++) innovation -= C_OBSERVER[disp_C++]*state_old[i];
        for(i=0;i<N_STATE;i++) state[i] += K_TRANSPOSE_OBSERVER[disp_K++]*innovation;
    }
}
