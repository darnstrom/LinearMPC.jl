void mpc_update_qp(c_float* th, c_float* dupper, c_float* dlower);
void mpc_get_solution(c_float* th, c_float* control, c_float* xstar);
int mpc_compute_control(c_float* control, c_float* state, c_float* reference, c_float* in_dist, c_float* out_dist);
