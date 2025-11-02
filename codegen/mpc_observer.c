void mpc_predict_state(c_float* state, c_float* control){
    int i,j,disp_F,disp_G=0;
    c_float state_old[N_STATE];
    for(i=0;i<N_STATE;i++) state_old[i] = state[i];
    for(i=0,disp_F=0,disp_G=0;i<N_STATE;i++){
        state[i] = OFFSET_OBSERVER[i];
        for(j=0;j<N_STATE;j++) state[i] += F_OBSERVER[disp_F++]*state_old[j];
        for(j=0;j<N_CONTROL;j++) state[i] += G_OBSERVER[disp_G++]*control[j];
    }
}

void mpc_correct_state(c_float* state, c_float* measurement){
    int i,j,disp_C,disp_K=0;
    c_float innovation, state_old[N_STATE];
    for(i=0;i<N_STATE;i++) state_old[i] = state[i];
    for(j=0;j<N_MEASUREMENT;j++){
        innovation = measurement[j];
        for(i=0;i<N_STATE;i++) innovation -= C_OBSERVER[disp_C++]*state_old[i];
        for(i=0;i<N_STATE;i++) state[i] += K_TRANSPOSE_OBSERVER[disp_K++]*innovation;
    }
}
