#if N_AFFINE_PARAMETER > 0
void mpc_update_parameter(c_float* parameter, c_float* control, c_float* state, c_float* reference, c_float* disturbance, c_float* affine_parameter){
#else
void mpc_update_parameter(c_float* parameter, c_float* control, c_float* state, c_float* reference, c_float* disturbance){
#endif
    int i,j;
    // update parameter
    for(i=0,j=0;j<N_STATE;i++, j++) parameter[i] = state[j];
#ifdef N_PREVIEW_HORIZON
    int disp = 0;
    for(;i<N_REFERENCE;i++) parameter[i] = 0.0; // reset parameter
    for(i=0;i<N_REFERENCE*N_PREVIEW_HORIZON;i++)
        for(j=0;j<N_REFERENCE;j++)
            parameter[N_STATE+j] += reference[i]*traj2setpoint[disp++];
    i = N_STATE + N_REFERENCE; // Setup i for remaining parameters
#else
    for(j=0;j<N_REFERENCE;i++, j++) parameter[i] = reference[j];
#endif
#ifdef N_DISTURBANCE_PREVIEW_HORIZON
    for(j=0;j<N_DISTURBANCE_BASE*N_DISTURBANCE_PREVIEW_HORIZON;i++, j++) parameter[i] = disturbance[j];
#else
    // If a disturbance trajectory is passed in column-major form, use the first column only.
    for(j=0;j<N_DISTURBANCE_BASE;i++, j++) parameter[i] = disturbance[j];
#endif
    for(j=0;j<N_CONTROL_PREV;i++, j++) parameter[i] = control[j];
#if N_AFFINE_PARAMETER > 0
    for(j=0;j<N_AFFINE_PARAMETER;i++, j++) parameter[i] = affine_parameter[j];
#endif
}
