extern int bnb_xprev_valid; // Whether bnb_xprev_shift holds a solution of the previous call
extern c_float bnb_xprev_shift[N_BNB_BINARY]; // Control one step later of each binary decision variable in that solution
extern c_float bnb_candidate_u[NX]; // LDP solution of the candidate
extern int bnb_candidate_used; // Whether the latest call returned the candidate
void mpc_reset_bnb_warm_start(void);
int mpc_bnb_warm_start(void);
