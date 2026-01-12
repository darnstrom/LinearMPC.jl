#ifndef MPC_WORKSPACE_H
#define MPC_WORKSPACE_H

#include "types.h"
#include "constants.h"
// Settings prototype
extern DAQPSettings settings;

#undef NX
#define NX 10
#undef N_CONSTR
#define N_CONSTR 0
#undef N_SIMPLE
#define N_SIMPLE 0 
// Workspace prototypes
extern c_float M[0];
extern c_float dupper[0];
extern c_float dlower[0];
extern c_float Rinv[55];
extern int sense[0];

extern c_float x[11];
extern c_float xold[11];

extern c_float lam[11];
extern c_float lam_star[11];
extern c_float u[11];

extern c_float L[66];
extern c_float D[11];
extern c_float xldl[11];
extern c_float zldl[11];

extern int WS[11];

extern DAQPWorkspace daqp_work;

#endif // ifndef MPC_WORKSPACE_H
#ifndef MPC_WORKSPACE_MPC_H
#define MPC_WORKSPACE_MPC_H

#define N_THETA 5
#define N_STATE 2
#define N_REFERENCE 1
#define N_DISTURBANCE 2
#define N_CONTROL_PREV 0
#define N_LINEAR_COST 0
#define N_CONTROL 1

extern c_float mpc_parameter[5];
extern c_float Dth[0];
extern c_float du[0];
extern c_float dl[0];

extern c_float Uth_offset[5];

extern c_float u_offset[1];

extern c_float uscaling[1];

void mpc_update_qp(c_float* th, c_float* dupper, c_float* dlower);
void mpc_get_solution(c_float* th, c_float* control, c_float* xstar);
#if N_LINEAR_COST > 0
int mpc_compute_control(c_float* control, c_float* state, c_float* reference, c_float* disturbance, c_float* linear_cost);
#else
int mpc_compute_control(c_float* control, c_float* state, c_float* reference, c_float* disturbance);
#endif
#define N_MEASUREMENT 1
extern c_float MPC_PLANT_DYNAMICS[12];
extern c_float MPC_MEASUREMENT_FUNCTION[5];
extern c_float K_TRANSPOSE_OBSERVER[2];
void mpc_predict_state(c_float* state, c_float* control);
void mpc_correct_state(c_float* state, c_float* measurement);
#endif // ifndef MPC_WORKSPACE_MPC_H
