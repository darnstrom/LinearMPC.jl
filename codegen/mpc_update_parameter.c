void  mpc_update_parameter(c_float* parameter, c_float* control, c_float* state, c_float* reference, c_float* disturbance){
    int i,j;
    // update parameter
    for(i=0,j=0;j<N_STATE;i++, j++) parameter[i] = state[j];
    for(j=0;j<N_REFERENCE;i++, j++) parameter[i] = reference[j];
    for(j=0;j<N_DISTURBANCE;i++, j++) parameter[i] = disturbance[j];
    for(j=0;j<N_CONTROL_PREV;i++, j++) parameter[i] = control[j];
}
