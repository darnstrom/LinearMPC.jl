void mpc_update_qp(c_float* th, c_float* dupper, c_float* dlower);
void mpc_get_solution(c_float* th, c_float* control, c_float* xstar);
#if N_LINEAR_COST > 0 && N_AFFINE_PARAMETER > 0
int mpc_compute_control(c_float* control, c_float* state, c_float* reference, c_float* disturbance, c_float* linear_cost, c_float* affine_parameter);
#elif N_LINEAR_COST > 0
int mpc_compute_control(c_float* control, c_float* state, c_float* reference, c_float* disturbance, c_float* linear_cost);
#elif N_AFFINE_PARAMETER > 0
int mpc_compute_control(c_float* control, c_float* state, c_float* reference, c_float* disturbance, c_float* affine_parameter);
#else
int mpc_compute_control(c_float* control, c_float* state, c_float* reference, c_float* disturbance);
#endif
