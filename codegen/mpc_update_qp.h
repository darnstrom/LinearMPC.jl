void mpc_update_qp(c_float* th, c_float* dupper, c_float* dlower);
void mpc_get_solution(c_float* th, c_float* control, c_float* xstar);
#if N_LINEAR_COST > 0
int mpc_compute_control(c_float* control, c_float* state, c_float* reference, c_float* disturbance, c_float* linear_cost);
#else
int mpc_compute_control(c_float* control, c_float* state, c_float* reference, c_float* disturbance);
#endif
