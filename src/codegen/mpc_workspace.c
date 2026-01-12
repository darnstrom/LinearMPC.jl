#include "types.h"
#include "constants.h"
// Settings
DAQPSettings settings = {(c_float)0.00000100000000000000, (c_float)0.00000000000100000000, (c_float)0.00000000001000000000, (c_float)0.00000100000000000000, (c_float)0.00000000000001000000, 10, 1000, (c_float)1000000000000000019884624838656.00000000000000000000, (c_float)0.00000000000000000000, (c_float)0.00000100000000000000, (c_float)0.00000100000000000000,(c_float)0.00000000000000000000,(c_float)0.00000000000000000000,(c_float)0.00000000003700000000,(c_float)0.00000000100000000000};

// Workspace
c_float M[0] = {
};
c_float dupper[0];
c_float dlower[0];
c_float Rinv[55] = {
(c_float)0.05913123959890825843,
(c_float)-0.44219070853771852425,
(c_float)0.07563943930642773317,
(c_float)0.13450308198127863824,
(c_float)0.08296285908315684243,
(c_float)0.03351503650556399150,
(c_float)0.00820342234489719255,
(c_float)0.00006112801796133156,
(c_float)-0.00082011661591903980,
(c_float)0.00000000000000000000,
(c_float)0.52694392767411457612,
(c_float)-0.54460396300627211819,
(c_float)-0.12435190598269314777,
(c_float)0.01319074746928446638,
(c_float)0.03322104495727436907,
(c_float)0.02057807639059361371,
(c_float)0.00803069335971536014,
(c_float)0.00187900136052345586,
(c_float)0.00000000000000000000,
(c_float)0.55649016061156619806,
(c_float)-0.51770997592793199793,
(c_float)-0.13954422322772958021,
(c_float)-0.00058798309659074217,
(c_float)0.02474930809138660079,
(c_float)0.01593913068350648332,
(c_float)0.00539823595288601698,
(c_float)0.00000000000000000000,
(c_float)0.60290734591607353376,
(c_float)-0.47139565956036033612,
(c_float)-0.13464812911885135072,
(c_float)-0.00806438128820495710,
(c_float)0.01569461861166062192,
(c_float)0.00867870241656273217,
(c_float)0.00000000000000000000,
(c_float)0.63523862812625531138,
(c_float)-0.43510749147576188722,
(c_float)-0.13139379857506258698,
(c_float)-0.01673379491700923219,
(c_float)0.00526328005406391259,
(c_float)0.00000000000000000000,
(c_float)0.66662583575762968113,
(c_float)-0.39070536930785826346,
(c_float)-0.12107168057534141437,
(c_float)-0.01998385424866426507,
(c_float)0.00000000000000000000,
(c_float)0.70945699177699161897,
(c_float)-0.32031081411917405877,
(c_float)-0.08746180365973840742,
(c_float)0.00000000000000000000,
(c_float)0.77460660260809366395,
(c_float)-0.19758581709861106068,
(c_float)0.00000000000000000000,
(c_float)0.87697457421527746924,
(c_float)-0.00000000000000000000,
(c_float)1.00000000000000000000,
};
int sense[0] = {
};
c_float x[11];
c_float xold[11];

c_float lam[11];
c_float lam_star[11];
c_float u[11];

c_float L[66];
c_float D[11];
c_float xldl[11];
c_float zldl[11];

int WS[11];

DAQPWorkspace daqp_work= {
NULL,
10, 0, 0,
M, dupper, dlower, Rinv, NULL, sense,
NULL,
NULL,
x, xold,
lam, lam_star, u, -1,
L, D, xldl,zldl,0,
WS, 0,
0,-1,
0.000000,
&settings, 
NULL,
0, NULL};

c_float mpc_parameter[5] = {
(c_float)0.00000000000000000000,
(c_float)0.00000000000000000000,
(c_float)0.00000000000000000000,
(c_float)0.00000000000000000000,
(c_float)0.00000000000000000000,
};
c_float Dth[0] = {
};
c_float du[0] = {
};
c_float dl[0] = {
};
c_float Uth_offset[5] = {
(c_float)-0.48053094928123435414,
(c_float)-1.24961535310130322785,
(c_float)0.48053094928123435414,
(c_float)-1.24961535310130322785,
(c_float)-0.48053094928123435414,
};
c_float u_offset[1] = {
(c_float)-0.00000000000000000000,
};
c_float uscaling[1] = {
(c_float)1.00000000000000000000,
};
#include "mpc_workspace.h"
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
        ctr_shift_th = u_offset[i];
        for(j = 0; j < N_THETA; j++) ctr_shift_th += Uth_offset[disp++]*th[j];
        control[i] = uscaling[i]*xstar[i]+ctr_shift_th;
    }
}

#include"daqp.h"
#ifdef DAQP_BNB
#include "bnb.h"
#endif

#if N_LINEAR_COST > 0
int mpc_compute_control(c_float* control, c_float* state, c_float* reference, c_float* disturbance, c_float* linear_cost){
    mpc_update_parameter(mpc_parameter, control, state, reference, disturbance, linear_cost);
#else
int mpc_compute_control(c_float* control, c_float* state, c_float* reference, c_float* disturbance){
    mpc_update_parameter(mpc_parameter, control, state, reference, disturbance);
#endif
    // update problem
    mpc_update_qp(mpc_parameter,daqp_work.dupper,daqp_work.dlower);
    daqp_work.reuse_ind=0; // clear workspace cache

#ifdef DAQP_BNB
    node_cleanup_workspace(0, &daqp_work);
    int exitflag = daqp_bnb(&daqp_work);
#else
#ifndef DAQP_WARMSTART
    deactivate_constraints(&daqp_work);
    reset_daqp_workspace(&daqp_work);
#endif
    int exitflag = daqp_ldp(&daqp_work);
#endif

    ldp2qp_solution(&daqp_work);
    mpc_get_solution(mpc_parameter,control,daqp_work.x);
    return exitflag;
}
c_float MPC_PLANT_DYNAMICS[12] = {
(c_float)0.00000000000000000000,
(c_float)1.00000000000000000000,
(c_float)1.00000000000000000000,
(c_float)0.00000000000000000000,
(c_float)1.00000000000000000000,
(c_float)0.00000000000000000000,
(c_float)0.00000000000000000000,
(c_float)0.00000000000000000000,
(c_float)1.00000000000000000000,
(c_float)1.00000000000000000000,
(c_float)0.00000000000000000000,
(c_float)0.00000000000000000000,
};
c_float MPC_MEASUREMENT_FUNCTION[5] = {
(c_float)0.00000000000000000000,
(c_float)1.00000000000000000000,
(c_float)0.00000000000000000000,
(c_float)0.00000000000000000000,
(c_float)1.00000000000000000000,
};
c_float K_TRANSPOSE_OBSERVER[2] = {
(c_float)0.99623457684784844091,
(c_float)0.61363043863156097046,
};
void mpc_predict_state(c_float* state, c_float* control, c_float* disturbance){
    int i,j,disp=0;
    c_float state_old[N_STATE];
    for(i=0;i<N_STATE;i++) state_old[i] = state[i];
    for(i=0,disp=0;i<N_STATE;i++){
        state[i] = MPC_PLANT_DYNAMICS[disp++];
        for(j=0;j<N_STATE;j++) state[i] += MPC_PLANT_DYNAMICS[disp++]*state_old[j];
        for(j=0;j<N_CONTROL;j++) state[i] += MPC_PLANT_DYNAMICS[disp++]*control[j];
        for(j=0;j<N_DISTURBANCE;j++) state[i] += MPC_PLANT_DYNAMICS[disp++]*disturbance[j];
    }
}

void mpc_correct_state(c_float* state, c_float* measurement, c_float* disturbance){
    int i,j,disp_C=0,disp_K=0;
    c_float innovation, state_old[N_STATE];
    for(i=0;i<N_STATE;i++) state_old[i] = state[i];
    for(j=0;j<N_MEASUREMENT;j++){
        innovation = measurement[j]-MPC_MEASUREMENT_FUNCTION[disp_C++];
        for(i=0;i<N_STATE;i++) innovation -= MPC_MEASUREMENT_FUNCTION[disp_C++]*state_old[i];
        for(i=0;i<N_DISTURBANCE;i++) innovation -= MPC_MEASUREMENT_FUNCTION[disp_C++]*disturbance[i];
        for(i=0;i<N_STATE;i++) state[i] += K_TRANSPOSE_OBSERVER[disp_K++]*innovation;
    }
}
