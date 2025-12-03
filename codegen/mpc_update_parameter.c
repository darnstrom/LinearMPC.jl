#if N_LINEAR_COST > 0
void mpc_update_parameter(c_float* parameter, c_float* control, c_float* state, c_float* reference, c_float* disturbance, c_float* linear_cost){
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
    for(j=0;j<N_DISTURBANCE;i++, j++) parameter[i] = disturbance[j];
    for(j=0;j<N_CONTROL_PREV;i++, j++) parameter[i] = control[j];
#if N_LINEAR_COST > 0
#ifdef N_MOVE_BLOCKS
    // Average linear cost over move blocks
    // linear_cost is (N_CONTROL x N_PREDICTION_HORIZON), column-major
    int block_offset = 0;
    for(int b = 0; b < N_MOVE_BLOCKS; b++) {
        int block_size = move_blocks[b];
        for(int u = 0; u < N_CONTROL; u++) {
            c_float sum = 0.0;
            for(int k = 0; k < block_size; k++) {
                sum += linear_cost[u + (block_offset + k) * N_CONTROL];
            }
            parameter[i++] = sum / block_size;
        }
        block_offset += block_size;
    }
#else
    for(j=0;j<N_LINEAR_COST;i++, j++) parameter[i] = linear_cost[j];
#endif
#endif
}
